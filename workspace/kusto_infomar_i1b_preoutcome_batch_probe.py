#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1A-PREOUTCOME-FIRST-FIVE-BATCH-FREEZE-2026-09-24-v1.0.json").read_text())
I0F=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I0F-SURVEY-SELECTION-RULE-FREEZE-2026-09-24-v1.0.json").read_text())
EXPECTED=PRE["preoutcome_authority"]["first_five_eligible_SURVEY_IDs"]
LAYER=I0F["feature_service"]
S=requests.Session(); S.headers.update({"User-Agent":"JANUS-KUSTO-INFOMAR-I1B/1.0","Accept-Encoding":"identity"})

q=LAYER+"/query"
r=S.get(q,params={
 "where":"1=1",
 "outFields":"OBJECTID,SURVEY_ID,VESSELNAME,SURVEYYEAR,START_DATE,END_DATE,PROJECT,COLLECTOR,PARTNER,GRIDSTATUS,GRID_URL,BS_STATUS,BS_URL,SURVREPURL,SUMREPURL",
 "returnGeometry":"false","resultRecordCount":"2000","f":"json"
},timeout=180)
r.raise_for_status(); j=r.json()
if "error" in j: raise RuntimeError(j["error"])
rows=[f.get("attributes") or {} for f in j.get("features",[]) or []]

def eligible(x):
    sid=str(x.get("SURVEY_ID") or "")
    gu=str(x.get("GRID_URL") or "")
    bu=str(x.get("BS_URL") or "")
    if not sid or gu in ("","NA") or bu in ("","NA"): return False
    if not gu.lower().startswith(("http://","https://")) or not bu.lower().startswith(("http://","https://")): return False
    if "merged" in gu.lower() or "merged" in bu.lower(): return False
    token="/"+sid.lower()+"/"
    if token not in gu.lower() or token not in bu.lower(): return False
    if not gu.lower().split("?")[0].endswith(".zip") or not bu.lower().split("?")[0].endswith(".zip"): return False
    return True

cand=sorted([x for x in rows if eligible(x)],key=lambda x:(str(x["SURVEY_ID"]),int(x.get("OBJECTID") or 0)))
now=[x["SURVEY_ID"] for x in cand[:5]]
if now!=EXPECTED:
    raise RuntimeError(f"preoutcome first-five drift: current={now} frozen={EXPECTED}")

def probe(url):
    rr=S.get(url,headers={"Range":"bytes=0-0"},timeout=180,allow_redirects=True,stream=True)
    cr=rr.headers.get("content-range") or rr.headers.get("Content-Range")
    total=None
    if cr and "/" in cr:
        try: total=int(cr.rsplit("/",1)[1])
        except: total=None
    return {
      "requested":url,"status":int(rr.status_code),"resolved":rr.url,
      "content_type":rr.headers.get("content-type"),"content_range":cr,
      "archive_size_bytes_from_content_range":total,"body_consumed":False
    }

batch=[]
for rank,row in enumerate(cand[:5],start=1):
    item={"rank":rank,"role":"HISTORICAL_RANK_1" if rank==1 else "PROSPECTIVE","row":row}
    if rank>=2:
        item["transport_probes"]={
          "bathymetry":probe(row["GRID_URL"]),
          "backscatter":probe(row["BS_URL"])
        }
        item["both_routes_range_usable"]=all(
          p["status"]==206 and p["archive_size_bytes_from_content_range"] is not None
          for p in item["transport_probes"].values()
        )
    batch.append(item)

out={
 "artifact_id":"JANUS-KUSTO-INFOMAR-I1B-PREOUTCOME-BATCH-METADATA-AND-TRANSPORT-RECEIPT-2026-09-24-v1.0",
 "batch_freeze":PRE["artifact_id"],
 "selection_rule":I0F["artifact_id"],
 "eligible_row_count":len(cand),
 "current_first_five_match_preoutcome_record":True,
 "batch":batch,
 "prospective_member_count":4,
 "all_prospective_routes_range_usable":all(x.get("both_routes_range_usable",False) for x in batch if x["rank"]>=2),
 "archive_bodies_consumed":False,
 "raster_pixels_read":False,
 "groundtruth_read":False,
 "next_gate":"I1C_BATCH_ARCHIVE_MEMBER_INVENTORY_ONLY",
 "claim_ceiling":"INFOMAR_PREOUTCOME_BATCH_METADATA_AND_TRANSPORT_HEADERS_ONLY"
}
p=OUT/"JANUS-KUSTO-INFOMAR-I1B-PREOUTCOME-BATCH-METADATA-AND-TRANSPORT-RECEIPT-2026-09-24-v1.0.json"
p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps(out,indent=2,ensure_ascii=False))
