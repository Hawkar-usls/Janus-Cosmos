#!/usr/bin/env python3
"""R2D: authentic CASDA + Parkes DAP pulsar archive read-only metadata audit.

A CASDA image/cube containing an OSIRIS direction in its s_region is NOT an
observation of a SETI narrowband signal, and pulsar search files are not
automatically valid ON/OFF technosignature observations. Never load power or
adapt model thresholds. Archived column values are not verified telescope TCS.
"""
import concurrent.futures,hashlib,json,math,pathlib,requests
import osiris_spear_r2c_atoa as old
P=pathlib.Path("research/seti/osiris_spear_r2d_casda_psrda_protocol.json")
S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/osiris_spear_r2d");OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-COSMOS-OSIRIS-R2D-readonly/1.0","Accept":"text/csv, application/x-votable+xml"}
ROUTES={x["id"]:x for x in S["archive_routes"]}
def digest(b):return hashlib.sha256(b).hexdigest()
def write(n,x):(OUT/n).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def safe(x):
    if x is None:return None
    z=str(x)
    return z[:320] if len(z)>320 else z
def decnum(x):
    try:
        z=float(x)
        return z if math.isfinite(z) else None
    except (TypeError,ValueError):return None
def tap(route,adql,maxrec=200):
    endpoint=route["tap_sync"]
    result={"tap_url":endpoint,"adql":adql,"adql_sha256":digest(adql.encode()),"state":"UNAVAILABLE","attempt_errors":[]}
    for method in ("GET","POST"):
        try:
            args={"REQUEST":"doQuery","LANG":"ADQL","QUERY":adql,"FORMAT":"csv","MAXREC":maxrec}
            if method=="GET":r=requests.get(endpoint,params=args,headers=UA,timeout=(9,S["search_limits"]["max_request_seconds"]))
            else:r=requests.post(endpoint,data=args,headers=UA,timeout=(9,S["search_limits"]["max_request_seconds"]))
            r.raise_for_status()
            raw=r.content
            if len(raw)>S["search_limits"]["max_html_or_vo_bytes"]:raise ValueError("TOO_LARGE_RESPONSE")
            cols,rows,fmt=old.parse_rows(raw)
            result.update({"state":"OK","method":method,"headers":cols,"rows":len(rows),"format":fmt,"body_sha256":digest(raw),"maybe_truncated":len(rows)>=maxrec})
            return result,rows
        except Exception as exc:result["attempt_errors"].append(f"{method}:{type(exc).__name__}:{str(exc)[:200]}")
    return result,[]
def schema(service):
    table=ROUTES[service]["table"]
    query=f"SELECT TOP 160 column_name,unit FROM TAP_SCHEMA.columns WHERE table_name='{table}'"
    rec,rows=tap(ROUTES[service],query,160)
    columns=[{"name":z.get("column_name"),"unit":z.get("unit")} for z in rows if z.get("column_name")]
    # Fail closed for incomplete or broken TAP_SCHEMA.
    if not columns:
        q=f"SELECT TOP 1 * FROM {table}"
        rec2,rows2=tap(ROUTES[service],q,1)
        if rec2["state"]=="OK":columns=[{"name":s,"unit":None} for s in rec2["headers"]]
        else:rec["fallback_errors"]=rec2["attempt_errors"]
    out={"service":service,"table":table,"schema_state":"OK" if columns else "UNAVAILABLE",
         "columns":columns,"tap_receipt":rec}
    print("OSIRIS_R2D_SCHEMA",json.dumps({"service":service,"state":out["schema_state"],"column_names":[x["name"] for x in columns]},ensure_ascii=False),flush=True)
    return out
def trimrow(r,fields):
    return {f:safe(r.get(f)) for f in fields if f in r}
def casda_one(t,schem):
    names={x["name"] for x in schem.get("columns",[])}
    target={"id":t["id"],"ra_deg":t["ra_deg"],"dec_deg":t["dec_deg"],
            "service":"CASDA","state":"SCHEMA_UNAVAILABLE","spatial_region_hits":[],
            "small_center_hits":[],"query_receipts":[],"backend_files_downloaded":0}
    if not {"s_ra","s_dec"}.issubset(names):
        target["state"]="MISSING_SKY_COLUMNS";return target
    cols=[x for x in ("obs_id","obs_collection","target_name","facility_name","instrument_name","dataproduct_type","s_ra","s_dec","s_fov","s_region","t_min","t_max","em_min","em_max","em_res_power","access_format","access_url") if x in names]
    maxrec=S["search_limits"]["max_results_per_target"]
    if "s_region" in names:
        # s_region is a sky footprint, but sensitivity inside an ASKAP imaging mosaic
        # cannot be presumed equal to a native Parkes high-res spectrogram.
        q=f"SELECT TOP {maxrec} {','.join(cols)} FROM ivoa.obscore WHERE 1=CONTAINS(POINT('ICRS',{t['ra_deg']:.10f},{t['dec_deg']:.10f}),s_region)"
        receipt,rows=tap(ROUTES["CASDA"],q,maxrec)
        target["query_receipts"].append({"purpose":"region_containment","receipt":receipt})
        if receipt["state"]=="OK":
            target["spatial_region_hits"]=[trimrow(z,cols) for z in rows[:40]]
    # Fixed narrow search of catalog point centers, not a beam footprint test.
    radius=S["search_limits"]["radius_deg"]
    ra,dec=t["ra_deg"],t["dec_deg"]
    half=radius/max(.12,math.cos(math.radians(dec)))
    if ra-half<0 or ra+half>360:
        target["state"]="RA_WRAP_REQUIRES_REVISED_QUERY";return target
    q=(f"SELECT TOP {maxrec} {','.join(cols)} FROM ivoa.obscore "+
       f"WHERE s_ra BETWEEN {ra-half:.9f} AND {ra+half:.9f} "+
       f"AND s_dec BETWEEN {dec-radius:.9f} AND {dec+radius:.9f}")
    receipt,rows=tap(ROUTES["CASDA"],q,maxrec)
    target["query_receipts"].append({"purpose":"catalog_center_0p2deg","receipt":receipt})
    ranked=[]
    for row in rows:
        rra=decnum(row.get("s_ra"));dde=decnum(row.get("s_dec"))
        if rra is None or dde is None or not(0<=rra<360 and -90<=dde<=90):continue
        sep=old.angle(ra,dec,rra,dde)
        if sep<=radius:ranked.append({**trimrow(row,cols),"distance_deg":sep,
             "catalog_center_is_not_actual_beam":True})
    ranked.sort(key=lambda z:z["distance_deg"])
    target["small_center_hits"]=ranked[:60]
    target["spatial_region_count"]=next((i["receipt"]["rows"] for i in target["query_receipts"] if i["purpose"]=="region_containment" and i["receipt"]["state"]=="OK"),None)
    target["narrow_center_count"]=len(ranked)
    target["state"]="CASDA_METADATA_QUERY_OK" if any(q["receipt"]["state"]=="OK" for q in target["query_receipts"]) else "CASDA_METADATA_QUERY_UNAVAILABLE"
    target["region_coverage_not_verified_native_SETI_spectrum"]=True
    target["center_query_maybe_truncated"]=receipt.get("maybe_truncated",False)
    print("OSIRIS_R2D_CASDA",json.dumps({"id":t["id"],"state":target["state"],"region_count":target["spatial_region_count"],
         "center_count":target["narrow_center_count"],"query_states":[a["receipt"]["state"] for a in target["query_receipts"]]},ensure_ascii=False),flush=True)
    return target
def psrda_one(t,schem):
    names={x["name"] for x in schem.get("columns",[])}
    target={"id":t["id"],"ra_deg":t["ra_deg"],"dec_deg":t["dec_deg"],
            "service":"PSRDA_DAP","state":"SCHEMA_UNAVAILABLE",
            "query_receipts":[],"candidate_rows":[],"native_pulsar_spectra_downloaded":0}
    if not {"ra_angle","dec_angle"}.issubset(names):target["state"]="MISSING_PULSAR_SKY_COLUMNS";return target
    cols=[x for x in ("projid","filename","ra_angle","dec_angle","date_obs","obsfreq","obs_mode","obstype","src_name","collection_fedora_PID","nchan","nsub","npol","tsamp") if x in names]
    maxrec=S["search_limits"]["max_results_per_target"]
    radius=S["search_limits"]["radius_deg"]
    ra,dec=t["ra_deg"],t["dec_deg"]
    half=radius/max(.12,math.cos(math.radians(dec)))
    q=(f"SELECT TOP {maxrec} {','.join(cols)} FROM psrda.pulsar_obs "
       f"WHERE ra_angle BETWEEN {ra-half:.9f} AND {ra+half:.9f} "
       f"AND dec_angle BETWEEN {dec-radius:.9f} AND {dec+radius:.9f}")
    receipt,rows=tap(ROUTES["PSRDA_DAP"],q,maxrec)
    target["query_receipts"].append({"purpose":"pulsar_catalog_center_0p2deg","receipt":receipt})
    valid=[]
    if receipt["state"]=="OK":
        for row in rows:
            pra=decnum(row.get("ra_angle"));pdec=decnum(row.get("dec_angle"))
            if pra is None or pdec is None or not(0<=pra<360 and -90<=pdec<=90):continue
            sep=old.angle(ra,dec,pra,pdec)
            if sep<=radius:valid.append({**trimrow(row,cols),"distance_deg":sep})
    valid.sort(key=lambda z:z["distance_deg"])
    target["catalog_rows_in_0p2deg"]=len(valid)
    target["query_maybe_truncated"]=receipt.get("maybe_truncated",False)
    target["candidate_rows"]=valid[:70]
    target["state"]="PSRDA_PULSAR_METADATA_QUERY_OK" if receipt["state"]=="OK" else "PSRDA_PULSAR_QUERY_UNAVAILABLE"
    target["pulsar_observation_is_not_SETI_narrowband_evidence"]=True
    print("OSIRIS_R2D_PSRDA",json.dumps({"id":t["id"],"state":target["state"],"rows":len(valid),"truncated":target["query_maybe_truncated"]},ensure_ascii=False),flush=True)
    return target
def main():
    schemas={s:schema(s) for s in ("CASDA","PSRDA_DAP")}
    write("tap_schema_receipts.json",schemas)
    t=S["frozen_target_list"]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as e:
        cas=list(e.map(lambda x:casda_one(x,schemas["CASDA"]),t))
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as e:
        psr=list(e.map(lambda x:psrda_one(x,schemas["PSRDA_DAP"]),t))
    write("casda_sky_coverage.json",{"sources":cas,"s_region_records_are_image_footprints_not_radio_detection":True})
    write("psrda_sky_coverage.json",{"sources":psr,"Parkes_pulsar_profiles_do_not_prove_SETI_bandwidth":True})
    # No misleading pass: each source is separate and needs native data.
    report={"stage":"OSIRIS_SPEAR_R2D","status":"TWO_OFFICIAL_ARCHIVE_ROUTES_SURVEYED__HOLD_TRUE_ON_TARGET_NATIVE_SPECTRA",
      "schema_status":{k:v["schema_state"] for k,v in schemas.items()},
      "CASDA_per_target":[{"id":z["id"],"state":z["state"],"region":z.get("spatial_region_count"),"center_0p2":z.get("narrow_center_count",0)} for z in cas],
      "PSRDA_per_target":[{"id":z["id"],"state":z["state"],"pulsar_catalog_0p2":z.get("catalog_rows_in_0p2deg",0)} for z in psr],
      "native_original_spectrum_bytes_downloaded":0,
      "verified_measured_TCS_and_per_feed_pointing":False,
      "matched_real_radio_ONOFF_and_comparable_controls":False,
      "unchanged_JANUS_radio_detector_applied":False,
      "new_real_technosignature_claim":False,"nondetection_claim":False,
      "CASDA_ASKAP_image_not_native_SETI_spectrogram":True,
      "DAP_pulsar_archive_not_narrowband_SETI_by_default":True,
      "external_tickets":{"ATNF":"ATNFDATSUP-1670","NRAO":"58509"},
      "production_changed":False,"reference":{"historical_false_alerts":51,"injected_localizations":318,"weak":69},
      "next":"Investigate archive metadata quality and request actual native source scans, TCS/per-feed footprint and adjacent OFF from data custodians. Do not download full wide field image and label it a validated OSIRIS narrowband signal."}
    write("r2d_gate_report.json",report)
    print("OSIRIS_R2D_RESULT",json.dumps(report,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
