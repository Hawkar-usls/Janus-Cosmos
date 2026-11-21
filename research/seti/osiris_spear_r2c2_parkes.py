#!/usr/bin/env python3
"""OSIRIS R2C2: Parkes-only, bounded ATOA TAP metadata footprint diagnostics.
No frequency-time spectra read. Neither s_fov nor s_region proves source pointing.
"""
import concurrent.futures,hashlib,json,math,pathlib,requests
import osiris_spear_r2c_atoa as older
P=pathlib.Path("research/seti/osiris_spear_r2c2_parkes_protocol.json")
Q=json.loads(P.read_text())
OUT=pathlib.Path("runtime/osiris_spear_r2c2");OUT.mkdir(parents=True,exist_ok=True)
HEAD={"User-Agent":"JANUS-COSMOS-OSIRIS-R2C2-PARKES-CATALOG-READONLY/1.0",
      "Accept":"text/csv, application/x-votable+xml"}
def sha(x):return hashlib.sha256(x).hexdigest()
def save(name,v):(OUT/name).write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def num(x):
    try:
        y=float(x)
        return y if math.isfinite(y) else None
    except (ValueError,TypeError):return None
def query(t):
    ra,dec=t["ra_deg"],t["dec_deg"]
    r=Q["radius_deg"];dra=r/max(.12,math.cos(math.radians(dec)))
    if ra-dra<0 or ra+dra>360:raise ValueError("RA_WRAP_NOT_HANDLED")
    select=",".join(Q["fields"])
    adql=(f"SELECT TOP {Q['query_limit_rows']} {select} FROM {Q['table']} "
        f"WHERE instrument_name = 'Parkes' AND "
        f"s_ra BETWEEN {ra-dra:.9f} AND {ra+dra:.9f} AND "
        f"s_dec BETWEEN {dec-r:.9f} AND {dec+r:.9f}")
    receipt={"id":t["id"],"parent":t["parent"],"kind":t["kind"],"ra_deg":ra,"dec_deg":dec,
       "status":"UNAVAILABLE","query_sha256":sha(adql.encode()),
       "radii_deg":[.05,.1,.5,1.5],"errors":[]}
    rows=[]
    for method in ("GET","POST"):
        try:
            data={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":adql,
                   "MAXREC":Q["query_limit_rows"]}
            if method=="GET":r0=requests.get(Q["tap"],params=data,timeout=(8,Q["timeout_seconds"]),headers=HEAD)
            else:r0=requests.post(Q["tap"],data=data,timeout=(8,Q["timeout_seconds"]),headers=HEAD)
            r0.raise_for_status()
            if len(r0.content)>12_000_000:raise ValueError("MAX_RESULT_SIZE_EXCEEDED")
            fieldnames,rows,fmt=older.parse_rows(r0.content)
            if not {"obs_id","s_ra","s_dec","instrument_name"}.issubset(fieldnames):raise ValueError("MISSING_REQUIRED_COLUMNS")
            receipt.update({"status":"OK","method":method,"response_sha256":sha(r0.content),
               "response_columns":fieldnames,"raw_rows":len(rows),"format":fmt})
            break
        except Exception as e:receipt["errors"].append(f"{method}:{type(e).__name__}:{str(e)[:260]}")
    selected=[]
    for x in rows:
        try:
            if x["instrument_name"].strip().lower()!="parkes":continue
            bra=float(x["s_ra"]);bdec=float(x["s_dec"])
            sep=older.angle(ra,dec,bra,bdec)
            if sep>Q["radius_deg"]+1e-8:continue
        except (ValueError,TypeError,KeyError):continue
        s_fov=num(x.get("s_fov"))
        # s_fov in IVOA ObsCore is nominal full field-of-view diameter (degrees).
        # Its association with time-variable observation sky footprint is not proven here.
        fov_possible=bool(s_fov is not None and s_fov>=2*sep)
        mjd_start=num(x.get("t_min"));mjd_end=num(x.get("t_max"))
        state="MISSING_TIME" if mjd_start is None or mjd_end is None else (
            "REVERSED_TIME_FIELDS" if mjd_start>mjd_end else "START_END_ORDER_VALID")
        entry={k:(str(v)[:700] if v is not None else None) for k,v in x.items()}
        entry.update({"angular_separation_deg":sep,"nominal_s_fov_deg":s_fov,
           "nominal_fov_contains_center_by_diameter_rule":fov_possible,
           "time_metadata_quality":state,"catalog_only":True})
        selected.append(entry)
    selected.sort(key=lambda x:(x["angular_separation_deg"],x.get("obs_id") or ""))
    receipt.update({"parkes_rows_in_1p5":len(selected),
       "parkes_rows_under_0p1":sum(x["angular_separation_deg"]<=.1 for x in selected),
       "parkes_rows_under_0p5":sum(x["angular_separation_deg"]<=.5 for x in selected),
       "closest_center_deg":selected[0]["angular_separation_deg"] if selected else None,
       "fov_cover_candidate_rows":sum(x["nominal_fov_contains_center_by_diameter_rule"] for x in selected),
       "missing_fov":sum(x["nominal_s_fov_deg"] is None for x in selected),
       "reversed_timing_rows":sum(x["time_metadata_quality"]=="REVERSED_TIME_FIELDS" for x in selected),
       "possibly_truncated":bool(receipt.get("raw_rows",0)>=Q["query_limit_rows"]),
       "rows":selected[:140]})
    if receipt["possibly_truncated"]:receipt["status"]="TRUNCATED_UNKNOWN"
    print("OSIRIS_PARKES_R2C2",json.dumps({k:v for k,v in receipt.items() if k not in ("rows","response_columns","errors")},ensure_ascii=False),flush=True)
    return receipt
def main():
    frozen=older.targets()
    assert len(frozen)==15
    with concurrent.futures.ThreadPoolExecutor(max_workers=Q["max_parallel_requests"]) as pool:
        r=list(pool.map(query,frozen))
    save("parkes_queries.json",{"data":r,"protocol_sha256":sha(P.read_bytes()),
       "original_spectra_read":0,"actual_pointing_verified":False})
    original=[x for x in r if x["kind"]=="original"]
    summaries=[{"target":x["id"],"status":x["status"],
       "parkes_rows_in_1p5":x["parkes_rows_in_1p5"],
       "parkes_rows_with_nominal_s_fov_cover":x["fov_cover_candidate_rows"],
       "closest_archive_center_sep_deg":x["closest_center_deg"],
       "within_0p1":x["parkes_rows_under_0p1"],
       "within_0p5":x["parkes_rows_under_0p5"],
       "truncation_risk":x["possibly_truncated"],"reversed_timing_rows":x["reversed_timing_rows"]} for x in original]
    potential=[]
    for x in original:
        for y in x["rows"]:
            if y["nominal_fov_contains_center_by_diameter_rule"]:
                potential.append({"target":x["id"],"obs_id":y.get("obs_id"),"obs_collection":y.get("obs_collection"),
                 "s_fov_deg":y["nominal_s_fov_deg"],"separation_deg":y["angular_separation_deg"],
                 "s_region":y.get("s_region"),"instrument_name":y.get("instrument_name"),
                 "frequency":y.get("frequency"),"access_format":y.get("access_format"),
                 "observation_source":"public ATOA ObsCore metadata only","actual_pointing":False})
    save("beam_footprint_diagnostics.json",{"potential_metadata_fov_covers":potential[:150],
       "s_region_is_not_encoder_log":True,"multi_beam_feed_offsets_unknown":True,"response_not_valid_footprint_mask":True})
    qualified=[x for x in original if x["status"]=="OK"]
    state=("CATALOG_NOMINAL_FOV_HINTS_EXIST__HOLD_ACTUAL_POINTING" if potential else
           "NO_CATALOG_NOMINAL_FOV_CONTAINMENT__HOLD_ACTUAL_POINTING")
    if any(x["status"]!="OK" for x in original):state="INCOMPLETE_PARKES_ARCHIVE_SEARCH__HOLD_ACTUAL_POINTING"
    report={"stage":"OSIRIS_SPEAR_R2C2","status":state,
      "method":"Parkes-only ATOA ivoa.obscore, frozen ICRS targets and exact native metadata fields, no spectral power.",
      "query_total":len(r),"original_queries_clean":len(qualified),
      "original_position_summary":summaries,
      "nominal_fov_candidate_rows":len(potential),
      "real_on_source_scans_qualified":0,"real_source_spectra_downloaded":0,
      "real_per_feed_antenna_pointing_verified":False,
      "native_JANUS_radiodetector_applied":False,
      "real_technosignature_detected":False,
      "nondetection_claim_permitted":False,"production_changed":False,
      "control_limit":"RA+40/+80 controls are archive footprint checks, not same frequency/time OFF controls",
      "next":"For metadata possible field coverage obtain original telescope scan path/per-feed pointing, receiver beam and raw spectral data through ATNF support ATNFDATSUP-1670. If none qualified, maintain HOLD and seek truly on-target prospective observation."}
    save("r2c2_result.json",report)
    print("OSIRIS_PARKES_R2C2_RESULT",json.dumps(report,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
