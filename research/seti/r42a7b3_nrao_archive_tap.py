#!/usr/bin/env python3
"""R42A7B3 NRAO official TAP observation-index check.
Indexes archived scans but cannot substitute for antenna encoder FITS or TCS logs.
Independent pointing is always HOLD unless separate original antenna data are inspected.
"""
import csv,hashlib,io,json,pathlib,requests
from datetime import datetime,timezone

OUT=pathlib.Path("runtime/r42a7b3_tap");OUT.mkdir(parents=True,exist_ok=True)
SERVER="https://data-query.nrao.edu/tap/sync"
UA={"User-Agent":"JANUS-COSMOS-R42A7B3-official-TAP-readonly/1.0"}
FROZEN={
 "project":"AGBT16A_999_42",
 "on_target":"HIP53020","off_target":"HIP53020_OFF",
 "on_mjd":57431.36287037037,"off_mjd":57431.36665509259,
 "on_ra_deg":10.847554166666669*15,"off_ra_deg":10.847555833333333*15,
 "on_dec_deg":6.804433611111112,"off_dec_deg":8.804431666666666,
 "on_sha256":"a4d282eddb536f43d8d296ac80207447cc49940ecda3c9187372d9a520790036",
 "off_sha256":"61cc6221938dd5a41a0c3ab3dfda91762ed9b8a89727d1cf67cfef26e0eef6fd"
}
QUERIES={
 "table_schema":"SELECT TOP 15 table_name,description FROM TAP_SCHEMA.tables WHERE table_name LIKE '%obscore%'",
 "by_native_mjd":"SELECT TOP 100 target_name, s_ra, s_dec, t_min, t_max, obs_publisher_did, instrument_name, access_url FROM tap_schema.obscore WHERE instrument_name='GBT' AND t_min BETWEEN 57431.35 AND 57431.39",
 "by_2016_project":"SELECT TOP 100 target_name, s_ra, s_dec, t_min, t_max, obs_publisher_did, instrument_name, access_url FROM tap_schema.obscore WHERE obs_publisher_did LIKE '%AGBT16A_999_42%'",
 "by_source":"SELECT TOP 100 target_name, s_ra, s_dec, t_min, t_max, obs_publisher_did, instrument_name, access_url FROM tap_schema.obscore WHERE instrument_name='GBT' AND target_name LIKE '%HIP53020%'"
}
def sha(x):return hashlib.sha256(x).hexdigest()
def save(n,x):(OUT/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
def query(label,adql):
  r={"label":label,"server":SERVER,"adql_sha256":sha(adql.encode()),"source":"NRAO official archive TAP ObsCore (metadata; no encoder data)","result":"UNAVAILABLE"}
  for how in ("GET","POST"):
    try:
      params={"REQUEST":"doQuery","LANG":"ADQL","QUERY":adql,"FORMAT":"csv","MAXREC":100}
      if how=="GET": response=requests.get(SERVER,params=params,headers=UA,timeout=(10,45))
      else: response=requests.post(SERVER,data=params,headers=UA,timeout=(10,45))
      response.raise_for_status()
      raw=response.content
      if len(raw)>750_000:raise ValueError("response too large")
      s=raw.decode("utf-8-sig","replace").strip()
      if s.startswith(("<","ERROR", "{")):raise ValueError("not CSV: "+s[:300])
      reader=csv.DictReader(io.StringIO(s))
      cols=reader.fieldnames or []
      if not cols or not any("target_name"==x.strip().lower() or "table_name"==x.strip().lower() for x in cols):
        raise ValueError("unexpected CSV headers: "+repr(cols))
      rows=list(reader);assert len(rows)<=100
      r.update({"result":"OK","method":how,"rows":rows,"columns":cols,"row_count":len(rows),"csv_sha256":sha(raw)})
      break
    except Exception as e:
      r.setdefault("errors",[]).append(how+":"+str(e)[:320])
  print("R42A7B3_NRAO_TAP",json.dumps({k:v for k,v in r.items() if k!="rows"},ensure_ascii=False),flush=True)
  return r
def main():
  output={k:query(k,adql) for k,adql in QUERIES.items()}
  save("nrao_tap_metadata.json",{"frozen":FROZEN,"queries":output,
    "security_scope":"Read-only official ObsCore metadata. No original GBT antenna encoder scan/time series retrieved.",
    "independent_pointing_verified":False})
  candidates=[]
  for key,r in output.items():
    for row in r.get("rows",[]):
      try:
        name=row.get("target_name") or ""
        t=float(row.get("t_min") or "nan")
        target="HIP53020" in name.upper()
        time_ok=57431.35 <= t <= 57431.39
        if target or time_ok:
          candidates.append({"query":key,"name":name,"t_min":t,"ra":row.get("s_ra"),"dec":row.get("s_dec"),
             "id":row.get("obs_publisher_did"),"time_in_window":time_ok})
      except (ValueError,TypeError):pass
  report={"stage":"R42A7B3_TAP","status":"NRAO_PUBLIC_OBSCORE_MATCH_NOT_TCS" if candidates else "NRAO_PUBLIC_OBSCORE_NO_MATCH_OR_UNAVAILABLE",
    "project":"AGBT16A_999_42","read_only":True,
    "source_hdf5_full_bytes_previously_verified":2,"tap_query_results":{k:{"status":v["result"],"rows":v.get("row_count",0)} for k,v in output.items()},
    "rows_matching_target_or_epoch":candidates[:50],
    "metadata_can_independently_prove_actual_antenna_pointing":False,
    "original_Antenna_FITS_or_encoder_logs_available":False,
    "GBT_staff_response_received":False,"independent_custodian_confirmed":False,
    "R42A7B4_blind_exam_permitted":False,
    "science":"HOLD_TCS_AND_EXTERNAL_CUSTODIAN","production_unchanged":True,
    "next":"Review GBT Helpdesk response for Antenna FITS/ScanLog and actual encoder time series. Catalog RA/Dec and metadata-only TAP rows are not independent antenna-encoder evidence."}
  save("gate_report.json",report)
  print("R42A7B3_TAP_RESULT",json.dumps(report,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
