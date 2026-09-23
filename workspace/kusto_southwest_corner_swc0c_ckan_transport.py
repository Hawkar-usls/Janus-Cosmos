#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,os
from pathlib import Path
from urllib.parse import urlparse
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOUTHWEST-CORNER-SWC0-UNCHANGED-V2-PREREG-2026-09-24-v1.0.json").read_text())
PAIR=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOUTHWEST-CORNER-SWC0-SOURCE-PAIR-RECEIPT-2026-09-24-v1.0.json").read_text())

BASE="https://www.data.gov.au"
API=BASE+"/data/api/3/action/package_search"
UA="Mozilla/5.0 (compatible; JANUS-KUSTO-SWC0C/1.0; +https://github.com/hawkar-usls/Janus-Cosmos)"
S=requests.Session()
S.headers.update({"User-Agent":UA,"Accept":"*/*"})

q='"South-west Corner Marine Park survey" GA-4858'
r=S.get(API,params={"q":q,"rows":100},timeout=120);r.raise_for_status()
sj=r.json()
pkgs=(sj.get("result") or {}).get("results",[]) or []

wanted={
 "bathymetry":{"title_token":"bathymetry","name_token":"geotif"},
 "backscatter":{"title_token":"backscatter","name_token":"geotif"},
}
candidates={k:[] for k in wanted}
for p in pkgs:
    title=(p.get("title") or "").lower()
    if "south-west corner marine park survey" not in title: continue
    for role,rule in wanted.items():
        if rule["title_token"] not in title: continue
        for res in p.get("resources",[]) or []:
            rn=(res.get("name") or "").lower()
            if rule["name_token"] not in rn: continue
            candidates[role].append({
              "package_id":p.get("id"),
              "package_name":p.get("name"),
              "package_title":p.get("title"),
              "metadata_modified":p.get("metadata_modified"),
              "owner_org":p.get("owner_org"),
              "resource_id":res.get("id"),
              "resource_name":res.get("name"),
              "format":res.get("format"),
              "url":res.get("url"),
              "url_type":res.get("url_type"),
              "cache_url":res.get("cache_url"),
              "last_modified":res.get("last_modified"),
              "created":res.get("created"),
            })

# dedupe exact package/resource ids, sort newest package revision first
for role in candidates:
    seen={}
    for x in candidates[role]:
        seen[(x["package_id"],x["resource_id"])]=x
    candidates[role]=sorted(seen.values(),key=lambda x:(x.get("metadata_modified") or "",x.get("package_id") or ""),reverse=True)

def probe(url, referer):
    headers={"Referer":referer,"Range":"bytes=0-0","Accept":"application/octet-stream,*/*"}
    try:
        rr=S.get(url,headers=headers,timeout=90,allow_redirects=True,stream=True)
        # never consume the body; this is transport/header resolution only
        return {
          "requested":url,
          "status":rr.status_code,
          "resolved":rr.url,
          "content_type":rr.headers.get("content-type"),
          "content_length":rr.headers.get("content-length"),
          "content_range":rr.headers.get("content-range"),
          "location":rr.history[-1].headers.get("location") if rr.history else None,
          "redirect_chain":[{"status":h.status_code,"url":h.url,"location":h.headers.get("location")} for h in rr.history],
          "body_consumed":False
        }
    except Exception as e:
        return {"requested":url,"error":type(e).__name__+": "+str(e),"body_consumed":False}

probes=[]
for role,rows in candidates.items():
    for x in rows[:12]:
        pkg=x["package_name"]
        rid=x["resource_id"]
        resource_page=f"{BASE}/data/dataset/{pkg}/resource/{rid}"
        urls=[]
        if x.get("url"): urls.append(("declared_url",x["url"]))
        urls.append(("ckan_download_no_filename",resource_page+"/download"))
        # CKAN often accepts /download/<basename>; derive only from declared URL basename.
        if x.get("url"):
            basename=os.path.basename(urlparse(x["url"]).path)
            if basename:
                urls.append(("ckan_download_with_basename",resource_page+"/download/"+basename))
        for route_kind,u in urls:
            p=probe(u,resource_page)
            p.update({"role":role,"route_kind":route_kind,"resource":x})
            probes.append(p)

# A route is transport-usable if Range request gets 200/206 with non-HTML content.
usable=[]
for p in probes:
    st=int(p.get("status",0) or 0)
    ct=(p.get("content_type") or "").lower()
    if st in (200,206) and "text/html" not in ct:
        usable.append(p)

# Deduplicate by role + resolved final URL.
dedup={}
for p in usable:
    dedup[(p["role"],p.get("resolved") or p["requested"])]=p
usable=list(dedup.values())

out={
 "artifact_id":"JANUS-KUSTO-SOUTHWEST-CORNER-SWC0C-CKAN-TRANSPORT-ROUTE-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],
 "source_pair_receipt":PAIR["artifact_id"],
 "query":q,
 "package_result_count":len(pkgs),
 "candidates":candidates,
 "transport_probes":probes,
 "usable_routes":usable,
 "usable_roles":sorted(set(p["role"] for p in usable)),
 "raster_values_read":False,
 "response_bodies_consumed":False,
 "status":"PASS_AT_LEAST_ONE_TRANSPORT_ROUTE_PER_ROLE" if set(p["role"] for p in usable)=={"bathymetry","backscatter"} else "TRANSPORT_ROUTE_STILL_BLOCKED",
 "claim_ceiling":"SOURCE_TRANSPORT_RESOLUTION_ONLY"
}
p=OUT/"JANUS-KUSTO-SOUTHWEST-CORNER-SWC0C-CKAN-TRANSPORT-ROUTE-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "status":out["status"],
 "usable_roles":out["usable_roles"],
 "usable_routes":[{"role":x["role"],"route_kind":x["route_kind"],"requested":x["requested"],"status":x.get("status"),"resolved":x.get("resolved"),"content_type":x.get("content_type"),"content_range":x.get("content_range"),"resource":x["resource"]} for x in usable],
 "bathymetry_probe_summary":[{"route_kind":x["route_kind"],"requested":x["requested"],"status":x.get("status"),"resolved":x.get("resolved"),"content_type":x.get("content_type"),"content_range":x.get("content_range"),"error":x.get("error"),"resource":x["resource"]} for x in probes if x["role"]=="bathymetry"],
 "raster_values_read":False,
 "response_bodies_consumed":False
},indent=2,ensure_ascii=False))
