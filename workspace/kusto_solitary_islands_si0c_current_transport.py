#!/usr/bin/env python3
from __future__ import annotations
import json,os
from pathlib import Path
from urllib.parse import urlparse
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOLITARY-ISLANDS-SI0-UNCHANGED-V2-PREREG-2026-09-24-v1.0.json").read_text())
API="https://www.data.gov.au/data/api/3/action/package_search"
BASE="https://www.data.gov.au"
S=requests.Session()
S.headers.update({"User-Agent":"Mozilla/5.0 (compatible; JANUS-KUSTO-SI0C/1.0)","Accept":"*/*"})

r=S.get(API,params={"q":'"20220017S" "Solitary Islands Gumbaynggirr Yaegl"',"rows":100},timeout=120)
r.raise_for_status()
pkgs=(r.json().get("result") or {}).get("results",[]) or []
rows=[]
for p in pkgs:
    title=p.get("title") or ""
    if "Solitary Islands Gumbaynggirr Yaegl" not in title: continue
    for res in p.get("resources",[]) or []:
        rn=res.get("name") or ""
        t=(title+" "+rn+" "+str(res.get("format") or "")).lower()
        if "bathymetry" not in t or "backscatter" not in t: continue
        if "zip" not in t: continue
        rows.append({
          "package_id":p.get("id"),
          "package_name":p.get("name"),
          "package_title":title,
          "metadata_modified":p.get("metadata_modified"),
          "owner_org":p.get("owner_org"),
          "resource_id":res.get("id"),
          "resource_name":rn,
          "format":res.get("format"),
          "url":res.get("url"),
          "created":res.get("created"),
          "last_modified":res.get("last_modified"),
        })
# dedupe by physical URL/resource
seen={}
for x in rows:
    seen[(x["resource_id"],x["url"])]=x
rows=sorted(seen.values(),key=lambda x:(x.get("metadata_modified") or ""),reverse=True)

def probe(x):
    resource_page=f'{BASE}/data/dataset/{x["package_name"]}/resource/{x["resource_id"]}'
    probes=[]
    urls=[]
    if x.get("url"): urls.append(("declared_url",x["url"]))
    urls.append(("ckan_download",resource_page+"/download"))
    if x.get("url"):
        bn=os.path.basename(urlparse(x["url"]).path)
        if bn: urls.append(("ckan_download_with_basename",resource_page+"/download/"+bn))
    for kind,u in urls:
        try:
            rr=S.get(u,headers={"Referer":resource_page,"Range":"bytes=0-0"},stream=True,timeout=90,allow_redirects=True)
            probes.append({
              "route_kind":kind,"requested":u,"status":rr.status_code,"resolved":rr.url,
              "content_type":rr.headers.get("content-type"),"content_length":rr.headers.get("content-length"),
              "content_range":rr.headers.get("content-range"),"body_consumed":False
            })
        except Exception as e:
            probes.append({"route_kind":kind,"requested":u,"error":type(e).__name__+": "+str(e),"body_consumed":False})
    return probes

allp=[]
for x in rows[:20]:
    for p in probe(x):
        p["resource"]=x
        allp.append(p)
usable=[p for p in allp if int(p.get("status",0) or 0) in (200,206) and "text/html" not in (p.get("content_type") or "").lower()]
dedup={}
for p in usable:
    dedup[p.get("resolved") or p["requested"]]=p
usable=list(dedup.values())

out={
 "artifact_id":"JANUS-KUSTO-SOLITARY-ISLANDS-SI0C-CURRENT-TRANSPORT-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],
 "candidate_resources":rows,
 "transport_probes":allp,
 "usable_routes":usable,
 "status":"PASS_UNIQUE_CURRENT_COMBINED_ZIP_ROUTE" if len(usable)==1 else "BLOCKED_NONUNIQUE_OR_MISSING_CURRENT_ROUTE",
 "raster_values_read":False,
 "response_bodies_consumed":False,
 "truth_layers_read":False,
 "next_gate":"SI0D_ARCHIVE_MEMBER_HASH_INVENTORY_ONLY" if len(usable)==1 else "SOURCE_ROUTE_RESOLUTION_ONLY",
 "claim_ceiling":"CURRENT_TRANSPORT_ROUTE_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-SOLITARY-ISLANDS-SI0C-CURRENT-TRANSPORT-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"status":out["status"],
 "candidate_resources":rows,
 "usable_routes":usable,
 "raster_values_read":False,"response_bodies_consumed":False
},indent=2,ensure_ascii=False))
