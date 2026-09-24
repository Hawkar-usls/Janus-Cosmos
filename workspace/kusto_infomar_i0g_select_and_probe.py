#!/usr/bin/env python3
from __future__ import annotations
import json, re
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
FR=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I0F-SURVEY-SELECTION-RULE-FREEZE-2026-09-24-v1.0.json").read_text())
LAYER=FR["feature_service"]
S=requests.Session();S.headers.update({"User-Agent":"JANUS-KUSTO-INFOMAR-I0G/1.0"})

q=LAYER+"/query"
r=S.get(q,params={"where":"1=1","outFields":"OBJECTID,SURVEY_ID,VESSELNAME,SURVEYYEAR,GRID_URL,BS_URL,GRIDSTATUS,BS_STATUS","returnGeometry":"false","resultRecordCount":"2000","f":"json"},timeout=180)
r.raise_for_status();j=r.json()
if "error" in j: raise RuntimeError(j["error"])
rows=[f.get("attributes") or {} for f in j.get("features",[]) or []]

def eligible(x):
    sid=str(x.get("SURVEY_ID") or "")
    gu=str(x.get("GRID_URL") or "")
    bu=str(x.get("BS_URL") or "")
    if not sid or gu in ("","NA") or bu in ("","NA"): return False
    if not gu.lower().startswith(("http://","https://")) or not bu.lower().startswith(("http://","https://")): return False
    if "merged" in gu.lower() or "merged" in bu.lower(): return False
    # Survey-specific path identity, case-insensitive, slash-delimited enough for current portal.
    token="/"+sid.lower()+"/"
    if token not in gu.lower() or token not in bu.lower(): return False
    if not gu.lower().split("?")[0].endswith(".zip") or not bu.lower().split("?")[0].endswith(".zip"): return False
    return True

cand=sorted([x for x in rows if eligible(x)],key=lambda x:(str(x["SURVEY_ID"]),int(x.get("OBJECTID") or 0)))
if not cand: raise RuntimeError("no eligible INFOMAR survey rows")
sel=cand[0]

def probe(url):
    rr=S.get(url,headers={"Range":"bytes=0-0"},timeout=120,allow_redirects=True,stream=True)
    return {"requested":url,"status":rr.status_code,"resolved":rr.url,
            "content_type":rr.headers.get("content-type"),"content_length":rr.headers.get("content-length"),
            "content_range":rr.headers.get("content-range"),"body_consumed":False}

probes={"bathymetry":probe(sel["GRID_URL"]),"backscatter":probe(sel["BS_URL"])}
usable=all(p["status"] in (200,206) and "text/html" not in (p.get("content_type") or "").lower() for p in probes.values())

out={
 "artifact_id":"JANUS-KUSTO-INFOMAR-I0G-METADATA-SELECTION-AND-TRANSPORT-PROBE-2026-09-24-v1.0",
 "selection_rule_freeze":FR["artifact_id"],
 "eligible_row_count":len(cand),
 "selected_row":sel,
 "first_five_eligible_SURVEY_IDs":[x["SURVEY_ID"] for x in cand[:5]],
 "transport_probes":probes,
 "both_routes_usable":usable,
 "raster_pixels_read":False,
 "archive_bodies_consumed":False,
 "next_gate":"I0H_ARCHIVE_MEMBER_INVENTORY_ONLY" if usable else "SOURCE_TRANSPORT_REPAIR_ONLY",
 "claim_ceiling":"INFOMAR_METADATA_SELECTION_AND_TRANSPORT_HEADERS_ONLY"
}
p=OUT/"JANUS-KUSTO-INFOMAR-I0G-METADATA-SELECTION-AND-TRANSPORT-PROBE-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps(out,indent=2,ensure_ascii=False))
