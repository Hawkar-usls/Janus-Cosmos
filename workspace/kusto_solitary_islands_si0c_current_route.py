#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOLITARY-ISLANDS-SI0-UNCHANGED-V2-PREREG-2026-09-24-v1.0.json").read_text())
REC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOLITARY-ISLANDS-SI0-SOURCE-BINDING-RECEIPT-2026-09-24-v1.0.json").read_text())

API="https://www.data.gov.au/data/api/3/action/package_search"
S=requests.Session(); S.headers.update({"User-Agent":"JANUS-KUSTO-SolitaryIslands-SI0C/1.0"})
queries=[
 '"Solitary Islands Gumbaynggirr Yaegl" bathymetry',
 '"Solitary Islands Gumbaynggirr Yaegl" backscatter',
 '"20220017S"'
]
packages={}
for q in queries:
    r=S.get(API,params={"q":q,"rows":100},timeout=120); r.raise_for_status()
    for p in (r.json().get("result") or {}).get("results",[]) or []:
        packages[p.get("id")]=p

def resources_for(role):
    rows=[]
    for p in packages.values():
        title=(p.get("title") or "").lower()
        if "solitary" not in title or "islands" not in title: continue
        if role not in title: continue
        for res in p.get("resources",[]) or []:
            rows.append({
              "package_id":p.get("id"),
              "package_name":p.get("name"),
              "package_title":p.get("title"),
              "metadata_modified":p.get("metadata_modified"),
              "resource_id":res.get("id"),
              "resource_name":res.get("name"),
              "format":res.get("format"),
              "url":res.get("url"),
              "url_type":res.get("url_type"),
              "cache_url":res.get("cache_url")
            })
    return rows

def semantic(rows):
    out=[]
    for x in rows:
        t=" ".join(str(x.get(k) or "") for k in ("resource_name","format","url")).lower()
        if "metadata" in t: continue
        if not x.get("url"): continue
        if any(k in t for k in ["geotif","geotiff",".tif",".tiff","cog","download","zip"]):
            out.append(x)
    return out

cand={"bathymetry":semantic(resources_for("bathymetry")),"backscatter":semantic(resources_for("backscatter"))}

probes=[]
for role,rows in cand.items():
    for x in rows[:25]:
        pkg=x["package_name"]; rid=x["resource_id"]; declared=x["url"]
        routes=[("declared_url",declared)]
        if pkg and rid:
            page=f"https://www.data.gov.au/data/dataset/{pkg}/resource/{rid}"
            routes.append(("ckan_download",page+"/download"))
            import os
            bn=os.path.basename(requests.utils.urlparse(declared).path)
            if bn: routes.append(("ckan_download_with_basename",page+"/download/"+bn))
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
    st=int(p.get("status",0) or 0); ct=(p.get("content_type") or "").lower()
    if st in (200,206) and "text/html" not in ct:
        usable.append(p)

# dedupe by physical resolved URL per role
ded={}
for p in usable: ded[(p["role"],p.get("resolved") or p["requested"])]=p
usable=list(ded.values())
roles=set(x["role"] for x in usable)

out={
 "artifact_id":"JANUS-KUSTO-SOLITARY-ISLANDS-SI0C-CURRENT-ROUTE-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],
 "source_receipt":REC["artifact_id"],
 "queries":queries,
 "package_count":len(packages),
 "semantic_candidates":cand,
 "transport_probes":probes,
 "usable_routes":usable,
 "usable_roles":sorted(roles),
 "status":"PASS_CURRENT_ROUTE_PER_ROLE" if roles=={"bathymetry","backscatter"} else "CURRENT_ROUTE_BLOCKED",
 "raster_values_read":False,
 "response_bodies_consumed":False,
 "truth_layers_read":False,
 "claim_ceiling":"CURRENT_SOURCE_ROUTE_RESOLUTION_ONLY"
}
p=OUT/"JANUS-KUSTO-SOLITARY-ISLANDS-SI0C-CURRENT-ROUTE-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False); p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"status":out["status"],"usable_roles":out["usable_roles"],
 "usable_routes":[{"role":x["role"],"route_kind":x["route_kind"],"resolved":x.get("resolved"),"status":x.get("status"),
                    "content_type":x.get("content_type"),"content_range":x.get("content_range"),"resource":x["resource"]} for x in usable],
 "raster_values_read":False,"truth_layers_read":False
},indent=2,ensure_ascii=False))
