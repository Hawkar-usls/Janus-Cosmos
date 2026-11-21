#!/usr/bin/env python3
"""OSIRIS SPEAR R2C: independent ATNF ATOA (Parkes) TAP observation-metadata checks.

Not a radio technosignature search. No native spectra opened, no telescope
pointing assumed from a catalog position. All source and control directions were
frozen before querying. Any unavailable or truncated result is UNKNOWN.
"""
import concurrent.futures,csv,hashlib,io,json,math,pathlib,requests,xml.etree.ElementTree as ET
from urllib.parse import urlparse
P=pathlib.Path("research/seti/osiris_spear_r2c_atoa_protocol.json")
C=json.loads(P.read_text())
OUT=pathlib.Path("runtime/osiris_spear_r2c");OUT.mkdir(parents=True,exist_ok=True)
BASE=C["service"]["sync"];ROOT="https://atoavo.atnf.csiro.au"
UA={"User-Agent":"JANUS-COSMOS-OSIRIS-R2C-ATOA-READONLY/1.0","Accept":"application/x-votable+xml, text/csv"}
def sha(b):return hashlib.sha256(b).hexdigest()
def save(f,x):(OUT/f).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def angle(ra1,de1,ra2,de2):
    a,b,c,d=map(math.radians,(ra1,de1,ra2,de2))
    q=math.sin((b-d)/2)**2+math.cos(b)*math.cos(d)*math.sin((a-c)/2)**2
    return 2*math.degrees(math.asin(min(1.,math.sqrt(max(0,q)))))
def parse_rows(raw):
    txt=raw.decode("utf8-sig","replace").strip()
    if txt.startswith("<"):
        root=ET.fromstring(raw)
        for info in root.iter():
            if info.tag.endswith("INFO") and info.attrib.get("name")=="QUERY_STATUS" and info.attrib.get("value")=="ERROR":
                raise ValueError("TAP_QUERY_STATUS_ERROR: "+(info.text or "")[:260])
        fields=[x.attrib.get("name") for x in root.iter() if x.tag.endswith("FIELD")]
        trs=[x for x in root.iter() if x.tag.endswith("TR")]
        rows=[]
        for tr in trs:
            cells=[x.text if x.text is not None else "" for x in tr if x.tag.endswith("TD")]
            if len(cells)==len(fields):rows.append(dict(zip(fields,cells)))
        if not fields:raise ValueError("TAP_NO_FIELD_METADATA")
        if not trs and any(x.tag.endswith("BINARY") or x.tag.endswith("BINARY2") for x in root.iter()):
            raise ValueError("BINARY_VOTABLE_NOT_PARSED")
        return fields,rows,"votable_xml"
    if txt.startswith(("{","ERROR","<!DOCTYPE")):raise ValueError("TAP_UNEXPECTED_RESPONSE "+txt[:150])
    d=csv.DictReader(io.StringIO(txt))
    if not d.fieldnames:raise ValueError("TAP_CSV_WITHOUT_FIELDS")
    return d.fieldnames,list(d),"csv"
def execute(query,timeout=25):
    results={"query_sha256":sha(query.encode()),"method":None,"result":"UNAVAILABLE","errors":[]}
    for method in ("GET","POST"):
        try:
            params={"REQUEST":"doQuery","LANG":"ADQL","QUERY":query,"FORMAT":"csv","MAXREC":C["search"]["max_rows_per_direction"]}
            if method=="GET":resp=requests.get(BASE,params=params,timeout=(8,timeout),headers=UA)
            else:resp=requests.post(BASE,data=params,timeout=(8,timeout),headers=UA)
            resp.raise_for_status()
            if len(resp.content)>C["search"]["max_raw_response_bytes"]:raise ValueError("TOO_MANY_RESULT_BYTES")
            fields,rows,fmt=parse_rows(resp.content)
            results.update({"result":"OK","method":method,"format":fmt,
                "response_sha256":sha(resp.content),"columns":fields,"rows_returned":len(rows)})
            return results,rows
        except Exception as e:results["errors"].append(f"{method}:{type(e).__name__}:{str(e)[:260]}")
    return results,[]
def schema():
    out={"stage":"SCHEMA_PREFLIGHT","tap":BASE,"declared_table":"ivoa.obscore","status":"UNKNOWN"}
    try:
        r=requests.get(ROOT+"/tap/tables",timeout=(8,22),headers=UA)
        r.raise_for_status()
        if len(r.content)<4*1024*1024:
            root=ET.fromstring(r.content)
            tablematches=[]
            for t in root.iter():
                if not t.tag.endswith("table"):continue
                names=[e.text for e in t if e.tag.endswith("name")]
                if any(x and "obscore" in x.lower() for x in names):
                    cols=[]
                    for c in t:
                        if not c.tag.endswith("column"):continue
                        n=next((v.text for v in c if v.tag.endswith("name")),None)
                        if n:cols.append({"name":n,"unit":next((v.text for v in c if v.tag.endswith("unit")),None)})
                    tablematches.append({"name":names[0] if names else None,"columns":cols})
            out.update({"tables_endpoint":"OK","tables_response_sha256":sha(r.content),
                        "tables_matching_obscore":tablematches})
        else:out["tables_endpoint"]="OVERSIZE"
    except Exception as e:out["tables_error"]=str(e)[:280]
    q="SELECT TOP 100 column_name, table_name, unit FROM TAP_SCHEMA.columns WHERE table_name = 'ivoa.obscore'"
    res,rows=execute(q,timeout=22)
    out["schema_tap_response"]=res
    if res["result"]=="OK" and rows:out["columns"]=[{"name":v.get("column_name"),"unit":v.get("unit")} for v in rows]
    elif out.get("tables_matching_obscore"):
        out["columns"]=out["tables_matching_obscore"][0]["columns"]
    if not out.get("columns"):
        # Read-only sample: documented ATOA users guide names the table explicitly.
        p,rr=execute("SELECT TOP 1 * FROM ivoa.obscore",timeout=22)
        out["sample_response"]=p
        if p["result"]=="OK":out["columns"]=[{"name":x,"unit":None} for x in p["columns"]]
    out["status"]="OK" if out.get("columns") else "UNAVAILABLE"
    save("atoa_schema.json",out)
    print("OSIRIS_ATOA_SCHEMA",json.dumps({"status":out["status"],"columns":[x["name"] for x in out.get("columns",[])],
         "error":out.get("tables_error"),"schema_errors":res.get("errors",[])},ensure_ascii=False),flush=True)
    return out
def targets():
    out=[]
    for x in C["frozen_targets"]:
        out.append({**x,"kind":"original","parent":x["id"]})
        for a in (40,80):out.append({**x,"kind":"rotated_coverage_control","parent":x["id"],
          "id":x["id"]+"_RA_PLUS_"+str(a),"ra_deg":(x["ra_deg"]+a)%360})
    assert len(out)==15 and len({x["id"] for x in out})==15
    return out
def is_parkes(r):
    name=" ".join(str(r.get(k,"")) for k in ("facility_name","instrument_name","telescope","obs_collection"))
    return any(x in name.upper() for x in ("PARKES","MURRIYANG","PKS"))
def inspect(target,columns):
    names=set(x["name"] for x in columns)
    if "s_ra" not in names or "s_dec" not in names:
        return {"id":target["id"],"status":"MISSING_S_RA_OR_S_DEC_IN_SCHEMA","rows":[],"parkes_rows":[]}
    fields=[f for f in C["search"]["prefer_metadata_columns"] if f in names]
    sel=",".join(fields)
    radius=1.5
    ra,dec=target["ra_deg"],target["dec_deg"]
    off=radius/max(.12,math.cos(math.radians(dec)))
    # Avoid wrong RA wrap. None of the frozen directions/control RA is near zero;
    # still fail closed if the circular rectangle crosses coordinate boundary.
    low=ra-off;high=ra+off
    if low<0 or high>360:
        return {"id":target["id"],"status":"RA_WRAP_REQUIRES_SEPARATE_QUERY","rows":[]}
    adql=(f"SELECT TOP {C['search']['max_rows_per_direction']} {sel} FROM ivoa.obscore "
          f"WHERE s_ra BETWEEN {low:.9f} AND {high:.9f} "
          f"AND s_dec BETWEEN {dec-radius:.9f} AND {dec+radius:.9f}")
    receipt,rows=execute(adql,timeout=25)
    result={"id":target["id"],"kind":target["kind"],"parent":target["parent"],
      "ra_deg":ra,"dec_deg":dec,"raw_query_receipt":receipt,
      "row_cap":C["search"]["max_rows_per_direction"],
      "status":"UNAVAILABLE" if receipt["result"]!="OK" else "OK",
      "possible_truncation":len(rows)>=C["search"]["max_rows_per_direction"],
      "matching_rows":0,"matched_parkes":0,"min_deg_parkes":None,"rows":[]}
    valid=[]
    for row in rows:
        try:
            rra=float(row["s_ra"]);rdec=float(row["s_dec"])
            if not (0<=rra<360 and -90<=rdec<=90):continue
            d=angle(ra,dec,rra,rdec)
            if d>radius+1e-8:continue
            parkes=is_parkes(row)
            entry={k:(str(v)[:500] if v is not None else None) for k,v in row.items()}
            entry.update({"separation_deg":d,"instrument_parkes_label":parkes})
            valid.append(entry)
        except (ValueError,TypeError):continue
    park=[x for x in valid if x["instrument_parkes_label"]]
    result.update({"matching_rows":len(valid),"matched_parkes":len(park),
      "min_deg_parkes":min([x["separation_deg"] for x in park],default=None),
      "parkes_under_0p1":sum(x["separation_deg"]<=.1 for x in park),
      "parkes_under_0p5":sum(x["separation_deg"]<=.5 for x in park),
      "parkes_under_1p5":len(park),
      "unresolved_instrument_metadata":sum(not x["instrument_parkes_label"] for x in valid),
      "rows":sorted(park,key=lambda q:q["separation_deg"])[:50],
      "sample_unresolved_instrument_rows":sorted([x for x in valid if not x["instrument_parkes_label"]],key=lambda q:q["separation_deg"])[:5]})
    if result["possible_truncation"]:result["status"]="TRUNCATED_NO_COMPLETENESS_CLAIM"
    print("OSIRIS_ATOA_TARGET",json.dumps({k:v for k,v in result.items() if k not in ("rows","sample_unresolved_instrument_rows","raw_query_receipt")},ensure_ascii=False),flush=True)
    return result
def main():
    sk=schema();positions=targets()
    results=[]
    if sk["status"]=="OK":
        with concurrent.futures.ThreadPoolExecutor(max_workers=C["search"]["max_parallel_requests"]) as pool:
            results=list(pool.map(lambda t:inspect(t,sk["columns"]),positions))
    else:
        results=[{"id":t["id"],"kind":t["kind"],"parent":t["parent"],
            "status":"NOT_QUERIED_SCHEMA_UNAVAILABLE","matching_rows":None} for t in positions]
    save("atoa_candidates.json",{"positions":results,"source":"ATNF ATOA VO ivoa.obscore",
        "source_observation_spectra_read":False,"actual_antenna_encoder_pointing_verified":False})
    save("atoa_query_receipts.json",{"protocol_sha256":sha(P.read_bytes()),"table":C["service"]["table"],
        "schema":sk["status"],"queries":[{"id":r["id"],"receipt":r.get("raw_query_receipt")} for r in results],
        "all_query_failures_are_unknown_not_negative":True})
    good=[r for r in results if r["status"]=="OK" and r.get("kind")=="original"]
    original=[r for r in results if r["kind"]=="original"]
    count=sum(r.get("matched_parkes",0) for r in good)
    report={"stage":"OSIRIS_SPEAR_R2C_ATOA",
        "status":"ATNF_ATOA_TAP_AVAILABLE__ARCHIVE_PARKES_CATALOG_CANDIDATES__HOLD_NATIVE_READING" if count else
         "ATNF_ATOA_TAP_QUERY_INCOMPLETE_OR_NO_CATALOG_MATCH__HOLD_NATIVE_READING",
        "schema_status":sk["status"],"directions":len(results),
        "original_directions_fully_queried":len(good),
        "original_directions_with_Parkes_catalog_rows":sum(r.get("matched_parkes",0)>0 for r in good),
        "parkes_candidate_rows_within_1p5_original":count,
        "original_summary":[{"id":x["id"],"status":x["status"],"parkes_rows":x.get("matched_parkes"),
           "nearest_deg":x.get("min_deg_parkes"),"parkes_under_0p1":x.get("parkes_under_0p1"),
           "parkes_under_0p5":x.get("parkes_under_0p5")} for x in original],
        "original_radio_bytes_fetched":0,"original_tcs_pointing_verified":False,
        "native_JANUS_detector_run":False,"actual_NEW_radio_signals_detected":False,
        "cannot_claim_signal_nondetection":True,"external_requests":{"NRAO":58509,"ATNF":"ATNFDATSUP-1670"},
        "external_data_custodian_ready":False,"production_unchanged":True,
        "scientific_next":"For every qualified Parkes archive row verify actual receiver/beam coverage, raw original data hash and native spectral resolution; follow OPAL-access policy; freeze matched ON/OFF before unchanged JANUS and independent baselines. ATOA VO contains metadata, not necessarily BL backend waterfalls."}
    save("r2c_gate_report.json",report)
    print("OSIRIS_ATOA_R2C_RESULT",json.dumps(report,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
