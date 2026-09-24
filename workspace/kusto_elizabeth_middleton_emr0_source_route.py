#!/usr/bin/env python3
from __future__ import annotations
import json, os
from pathlib import Path
from urllib.parse import urlparse
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EMR0-UNCHANGED-V2-PREREG-2026-09-24-v1.0.json").read_text())
API="https://www.data.gov.au/data/api/3/action/package_search"
S=requests.Session();S.headers.update({"User-Agent":"JANUS-KUSTO-ElizabethMiddleton-EMR0/1.0"})

queries=['"Elizabeth and Middleton Reef" GA4848 bathymetry','"Elizabeth and Middleton Reef" GA4848 backscatter']
packages={}
for q in queries:
    r=S.get(API,params={"q":q,"rows":100},timeout=120);r.raise_for_status()
    for p in (r.json().get("result") or {}).get("results",[]) or []:
        packages[p.get("id")]=p

def rows_for(role):
    rows=[]
    for p in packages.values():
        title=(p.get("title") or "").lower()
        if "elizabeth" not in title or "middleton" not in title:continue
        if role not in title:continue
        for res in p.get("resources",[]) or []:
            rn=(res.get("name") or "").lower()
            if "metadata" in rn:continue
            if not res.get("url"):continue
            if not any(k in rn for k in ["geotif","geotiff","5m","zip","download"]):continue
            rows.append({
              "package_id":p.get("id"),"package_name":p.get("name"),"package_title":p.get("title"),
              "metadata_modified":p.get("metadata_modified"),"resource_id":res.get("id"),
              "resource_name":res.get("name"),"format":res.get("format"),"url":res.get("url")
            })
    return rows

cand={"bathymetry":rows_for("bathymetry"),"backscatter":rows_for("backscatter")}

probes=[]
for role,rows in cand.items():
    for x in rows[:30]:
        pkg=x["package_name"];rid=x["resource_id"];decl=x["url"]
        routes=[("declared_url",decl)]
        if pkg and rid:
            page=f"https://www.data.gov.au/data/dataset/{pkg}/resource/{rid}"
            routes.append(("ckan_download",page+"/download"))
            bn=os.path.basename(urlparse(decl).path)
            if bn:routes.append(("ckan_download_with_basename",page+"/download/"+bn))
        for kind,u in routes:
            try:
                rr=S.get(u,headers={"Range":"bytes=0-0","Referer":f"https://www.data.gov.au/data/dataset/{pkg}"},timeout=90,allow_redirects=True,stream=True)
                probes.append({
                  "role":role,"route_kind":kind,"requested":u,"status":rr.status_code,"resolved":rr.url,
                  "content_type":rr.headers.get("content-type"),"content_length":rr.headers.get("content-length"),
                  "content_range":rr.headers.get("content-range"),"resource":x,"body_consumed":False
                })
            except Exception as e:
                probes.append({"role":role,"route_kind":kind,"requested":u,"error":type(e).__name__+": "+str(e),"resource":x,"body_consumed":False})

usable=[]
for p in probes:
    st=int(p.get("status",0) or 0);ct=(p.get("content_type") or "").lower()
    if st in (200,206) and "text/html" not in ct:
        usable.append(p)
# collapse catalogue revisions by physical resolved URL and role
ded={}
for p in usable:ded[(p["role"],p.get("resolved") or p["requested"])]=p
usable=list(ded.values())

by_role={}
for role in ("bathymetry","backscatter"):
    rr=[x for x in usable if x["role"]==role]
    by_role[role]=rr

out={
 "artifact_id":"JANUS-KUSTO-ELIZABETH-MIDDLETON-EMR0-SOURCE-ROUTE-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],
 "queries":queries,
 "package_count":len(packages),
 "semantic_candidates":cand,
 "transport_probes":probes,
 "usable_routes_by_role":by_role,
 "status":"PASS_UNIQUE_PHYSICAL_ROUTE_PER_ROLE" if all(len(by_role[r])==1 for r in by_role) else "SOURCE_ROUTE_NONUNIQUE_OR_BLOCKED",
 "raster_values_read":False,
 "response_bodies_consumed":False,
 "groundtruth_assets_read":False,
 "next_gate":"EMR0B_ARCHIVE_MEMBER_INVENTORY_ONLY" if all(len(by_role[r])==1 for r in by_role) else "SOURCE_ROUTE_RESOLUTION_ONLY",
 "claim_ceiling":"SOURCE_ROUTE_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-ELIZABETH-MIDDLETON-EMR0-SOURCE-ROUTE-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"status":out["status"],
 "usable_routes_by_role":{k:[{"route_kind":x["route_kind"],"resolved":x.get("resolved"),"status":x.get("status"),"content_range":x.get("content_range"),"resource":x["resource"]} for x in v] for k,v in by_role.items()},
 "raster_values_read":False,"groundtruth_assets_read":False
},indent=2,ensure_ascii=False))
