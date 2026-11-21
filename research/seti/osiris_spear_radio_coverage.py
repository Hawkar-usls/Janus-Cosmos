#!/usr/bin/env python3
"""OSIRIS SPEAR Radio R1. Real official NRAO GBT ObsCore exposure query, not a
radio detector run. Old optical alignment is not presumed true. No re-fitting.
"""
import concurrent.futures,csv,hashlib,io,json,math,pathlib,requests,urllib.parse
P=pathlib.Path("research/seti/osiris_spear_radio_protocol.json")
Q=json.loads(P.read_text())
OUT=pathlib.Path("runtime/osiris_spear_radio_r1");OUT.mkdir(parents=True,exist_ok=True)
TAP="https://data-query.nrao.edu/tap/sync"
UA={"User-Agent":"JANUS-COSMOS-OSIRIS-RADIO-R1 scientific-readonly/1.0"}
def sha(x):return hashlib.sha256(x).hexdigest()
def write(name,obj):(OUT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def angular(ra1,dec1,ra2,dec2):
    a,b,c,d=map(math.radians,[ra1,dec1,ra2,dec2])
    s=math.sin((b-d)*.5)**2+math.cos(b)*math.cos(d)*math.sin((a-c)*.5)**2
    return math.degrees(2*math.asin(min(1,math.sqrt(max(0,s)))))
def positions():
    out=[]
    for base in Q["coordinates_icrs_j2000"]:
        ra=base["ra_deg"];dec=base["dec_deg"]
        for mirror in Q["mirror_policy"]["operators"]:
            name=mirror["id"]
            if name=="ORIGINAL":mra,mdec=ra,dec
            elif name=="RA_MERIDIAN_REFLECTION":mra,mdec=(360-ra)%360,dec
            elif name=="EQUATORIAL_REFLECTION":mra,mdec=ra,-dec
            elif name=="ANTIPODE":mra,mdec=(ra+180)%360,-dec
            else:raise AssertionError(name)
            out.append({"id":base["id"]+"__"+name,"base_id":base["id"],"role":base["role"],
             "kind":"original" if name=="ORIGINAL" else "mirror","mirror":name,
             "ra_deg":mra,"dec_deg":mdec,"origin":"frozen historical sky vector transformation; not observation"})
        for d in (40.,80.):
            out.append({"id":base["id"]+"__RA_CONTROL_"+str(int(d)),"base_id":base["id"],"role":"survey_coverage_control",
              "kind":"control","mirror":"RA_ROTATED_NULL","ra_deg":(ra+d)%360,"dec_deg":dec,
              "origin":"predeclared RA rotation; only exposure reference"})
    assert len(out)==len(Q["coordinates_icrs_j2000"])*6
    assert len({x["id"] for x in out})==len(out)
    return out

def adql(ra,dec):
    # Rectangular prefilter for spherical separation within 0.5deg; conservative
    # longitude expansion via cos(dec); no control accesses hidden injection labels.
    half=.5/max(.2,math.cos(math.radians(dec)))
    lo=max(0.,ra-half);hi=min(360.,ra+half)
    low=max(-90.,dec-.5);high=min(90.,dec+.5)
    cols="target_name,s_ra,s_dec,t_min,t_max,obs_publisher_did,instrument_name,access_url"
    return ("SELECT TOP 100 "+cols+" FROM tap_schema.obscore WHERE instrument_name='GBT' "+
            f"AND s_ra BETWEEN {lo:.9f} AND {hi:.9f} AND s_dec BETWEEN {low:.9f} AND {high:.9f}")
def one(item):
    query=adql(item["ra_deg"],item["dec_deg"])
    record={"id":item["id"],"kind":item["kind"],"base_id":item["base_id"],
     "ra_deg":item["ra_deg"],"dec_deg":item["dec_deg"],"source":"NRAO ObsCore GBT target_name/s_ra/s_dec metadata only",
     "query_sha256":sha(query.encode()),"result":"UNAVAILABLE",
     "limited_to_first_100_rows":True,"native_waterfall_downloaded":False,
     "independent_telescope_actual_pointing_verified":False}
    for method in ("GET","POST"):
        try:
            params={"REQUEST":"doQuery","LANG":"ADQL","QUERY":query,"FORMAT":"csv","MAXREC":100}
            if method=="GET":r=requests.get(TAP,params=params,headers=UA,timeout=(8,22))
            else:r=requests.post(TAP,data=params,headers=UA,timeout=(8,22))
            r.raise_for_status()
            raw=r.content
            if len(raw)>650_000:raise ValueError("TAP_CSV_RESPONSE_OVER_650KB")
            txt=raw.decode("utf8-sig","replace").strip()
            if txt.startswith(("<","{","ERROR")):raise ValueError("NONCSV "+txt[:110])
            parser=csv.DictReader(io.StringIO(txt));header=parser.fieldnames
            if not header or not {"target_name","s_ra","s_dec","t_min"}.issubset(set(header)):
                raise ValueError("UNEXPECTED_COLUMNS "+repr(header))
            rows=list(parser)
            if len(rows)>100:raise ValueError("OVER_MAXREC")
            matched=[]
            for row in rows:
                try:
                    bra=float(row["s_ra"]);bdec=float(row["s_dec"])
                    if not(-90<=bdec<=90 and 0<=bra<360):continue
                    angle=angular(item["ra_deg"],item["dec_deg"],bra,bdec)
                    if angle>.50001:continue
                    mjd=float(row["t_min"]) if row.get("t_min") else None
                    m={"separation_deg":angle,"target_name":row.get("target_name"),
                      "s_ra":bra,"s_dec":bdec,"mjd_start":mjd,"mjd_end":row.get("t_max"),
                      "observation_id":row.get("obs_publisher_did"),"instrument":row.get("instrument_name"),
                      "access_url":row.get("access_url") or None,
                      "actual_encoder_pointing":False}
                    matched.append(m)
                except (ValueError,TypeError):continue
            matched.sort(key=lambda x:(x["separation_deg"],x["observation_id"] or ""))
            record.update({"result":"OK","method":method,"csv_sha256":sha(raw),
             "archive_rows_returned":len(rows),"truncated_or_at_risk":len(rows)>=100,
             "count_in_0p1_deg":sum(m["separation_deg"]<=.1 for m in matched),
             "count_in_0p5_deg":len(matched),
             "closest_sep_deg":matched[0]["separation_deg"] if matched else None,
             "closest_archive_targets":matched[:12]})
            break
        except Exception as e:
            record.setdefault("query_errors",[]).append(method+":"+str(e)[:220])
    print("OSIRIS_SPEAR_RADIO_COVERAGE",json.dumps({"id":item["id"],"result":record["result"],
      "within_0p1":record.get("count_in_0p1_deg",0),"within_0p5":record.get("count_in_0p5_deg",0),
      "closest":record.get("closest_sep_deg")},ensure_ascii=False),flush=True)
    return record
def main():
    targets=positions()
    write("original_targets.json",{"frozen_protocol_sha256":sha(P.read_bytes()),"positions":targets,
      "semantic_note":"Mirrors and RA controls are new directions, not additional sightings."})
    # Bounded, read-only official metadata requests, no sky signal FFT or detector.
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        results=list(pool.map(one,targets))
    write("radio_coverage_inventory.json",{"archive":TAP,"responses":results,"read_only":True,
      "native_original_HDF5_or_FIL_files_read":0,"native_R42_detector_scores_created":0,
      "algorithm_used":"frozen coordinate to official GBT ObsCore query + angular separation + mirrored/control footprint audit",
      "not_the_JANUS_SETI_detector":True})
    perbase=[]
    for target in Q["coordinates_icrs_j2000"]:
        d=[z for z in results if z["base_id"]==target["id"]]
        assert len(d)==6
        perbase.append({"target":target["id"],"original":d[0],
         "mirrors":[{"id":v["id"],"status":v["result"],"within_0p1":v.get("count_in_0p1_deg"),
             "within_0p5":v.get("count_in_0p5_deg")} for v in d[1:4]],
         "rotated_controls":[{"id":v["id"],"status":v["result"],"within_0p1":v.get("count_in_0p1_deg"),
             "within_0p5":v.get("count_in_0p5_deg")} for v in d[4:]],
         "science":"No SETI score, detection, nondetection, matched OFF, or significance can be inferred from this metadata inventory."})
    write("mirror_and_control_coverage.json",{"per_target":perbase,"exposure_comparisons_not_calibrated":True})
    ok=sum(r["result"]=="OK" for r in results)
    raw_count=0
    n_original_match=sum(x["original"].get("count_in_0p1_deg",0)>0 for x in perbase)
    possible=[{"target":x["target"],"count_0p1":x["original"].get("count_in_0p1_deg",0),
      "count_0p5":x["original"].get("count_in_0p5_deg",0),
      "closest":x["original"].get("closest_sep_deg")}
      for x in perbase]
    report={"stage":"OSIRIS_SPEAR_RADIO_R1",
      "status":"PASS_AUTHENTIC_NRAO_ARCHIVE_METADATA_COVERAGE__HOLD_NATIVE_JANUS_SPECTRAL_TEST" if ok==len(results) else
      "PARTIAL_NRAO_OBSCORE_COVERAGE__HOLD_NATIVE_JANUS_SPECTRAL_TEST",
      "queried_positions":len(targets),"successful_queries":ok,
      "primary_targets":len(perbase),"original_positions_with_0p1deg_NRAO_metadata_match":n_original_match,
      "primary":possible,"raw_spectral_files_verified":raw_count,"JANUS_COSMOS_detector_applied_to_source_spectra":False,
      "qualified_ON_OFF_pointing_groups":0,"true_technosignature_detected":False,
      "can_make_radio_nondetection_claim":False,
      "prior_optical_gate":"FAIL_2.298deg_vs_0.25deg",
      "external_blind_validation":False,"production_model_changed":False,
      "historical_metrics_untouched":{"false_alerts":51,"localized":318,"weak":69},
      "scientific_next":"If matching observation URLs provide original readable radio spectra, verify raw hashes, original TCS pointing, calibrated frequency/time grid and a matched ON/OFF control. Only then run frozen native JANUS R39/R42 lane, with same controls and honest multiple-testing gates. Do not claim geometry selected a detected radio signal."}
    write("gate_report.json",report)
    print("OSIRIS_SPEAR_RADIO_RESULT",json.dumps(report,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
