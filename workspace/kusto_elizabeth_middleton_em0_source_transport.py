#!/usr/bin/env python3
from __future__ import annotations
import json,os
from pathlib import Path
from urllib.parse import urlparse
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EM0-UNCHANGED-V2-PREREG-2026-09-24-v1.0.json").read_text())
BASE="https://www.data.gov.au"; API=BASE+"/data/api/3/action/package_search"
S=requests.Session();S.headers.update({"User-Agent":"Mozilla/5.0 (compatible; JANUS-KUSTO-EM0/1.0)","Accept":"*/*"})

r=S.get(API,params={"q":'"Elizabeth" "Middleton" "GA4848"',"rows":100},timeout=120);r.raise_for_status()
pkgs=(r.json().get("result") or {}).get("results",[]) or []
cand={"bathymetry":[],"backscatter":[]}
for p in pkgs:
    title=p.get("title") or ""; tl=title.lower()
    if "elizabeth" not in tl or "middleton" not in tl: continue
    role=None
    if "bathymetry survey" in tl: role="bathymetry"
    elif "backscatter" in tl: role="backscatter"
    if not role: continue
    for res in p.get("resources",[]) or []:
        rn=res.get("name") or ""; t=(rn+" "+str(res.get("format") or "")).lower()
        if "zip" not in t: continue
        if "metadata" in t: continue
        cand[role].append({
          "package_id":p.get("id"),"package_name":p.get("name"),"package_title":title,
          "metadata_modified":p.get("metadata_modified"),"owner_org":p.get("owner_org"),
          "resource_id":res.get("id"),"resource_name":rn,"format":res.get("format"),
          "url":res.get("url"),"created":res.get("created")
        })
for role in cand:
    seen={}
    for x in cand[role]:seen[(x["resource_id"],x["url"])]=x
    cand[role]=sorted(seen.values(),key=lambda x:x.get("metadata_modified") or "",reverse=True)

def probes_for(role,x):
    page=f'{BASE}/data/dataset/{x["package_name"]}/resource/{x["resource_id"]}'
    routes=[]
    if x.get("url"):routes.append(("declared",x["url"]))
    routes.append(("ckan_download",page+"/download"))
    if x.get("url"):
        bn=os.path.basename(urlparse(x["url"]).path)
        if bn:routes.append(("ckan_download_with_basename",page+"/download/"+bn))
    out=[]
    for kind,u in routes:
        try:
            rr=S.get(u,headers={"Referer":page,"Range":"bytes=0-0"},stream=True,timeout=90,allow_redirects=True)
            out.append({"role":role,"route_kind":kind,"requested":u,"status":rr.status_code,
              "resolved":rr.url,"content_type":rr.headers.get("content-type"),
              "content_range":rr.headers.get("content-range"),"content_length":rr.headers.get("content-length"),
              "body_consumed":False,"resource":x})
        except Exception as e:
            out.append({"role":role,"route_kind":kind,"requested":u,"error":type(e).__name__+": "+str(e),
              "body_consumed":False,"resource":x})
    return out

probes=[]
for role,rows in cand.items():
    for x in rows[:12]:probes.extend(probes_for(role,x))
usable=[p for p in probes if int(p.get("status",0) or 0) in (200,206) and "text/html" not in (p.get("content_type") or "").lower()]
dedup={}
for p in usable:dedup[(p["role"],p.get("resolved") or p["requested"])]=p
usable=list(dedup.values())
roles=set(p["role"] for p in usable)
out={
 "artifact_id":"JANUS-KUSTO-ELIZABETH-MIDDLETON-EM0-SOURCE-TRANSPORT-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],"candidate_resources":cand,"transport_probes":probes,
 "usable_routes":usable,"usable_roles":sorted(roles),
 "status":"PASS_AT_LEAST_ONE_ROUTE_PER_ROLE" if roles=={"bathymetry","backscatter"} else "SOURCE_TRANSPORT_INCOMPLETE",
 "raster_values_read":False,"groundtruth_assets_read":False,"response_bodies_consumed":False,
 "next_gate":"EM0B_ARCHIVE_MEMBER_BINDING" if roles=={"bathymetry","backscatter"} else "SOURCE_ROUTE_RESOLUTION_ONLY",
 "claim_ceiling":"SOURCE_TRANSPORT_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-ELIZABETH-MIDDLETON-EM0-SOURCE-TRANSPORT-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({"artifact_id":out["artifact_id"],"status":out["status"],"usable_roles":out["usable_roles"],
 "usable_routes":usable,"raster_values_read":False},indent=2,ensure_ascii=False))
