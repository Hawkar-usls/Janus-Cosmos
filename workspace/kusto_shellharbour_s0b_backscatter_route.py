#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,re
from pathlib import Path
from urllib.parse import urljoin,quote_plus
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
REC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHELLHARBOUR-S0-SOURCE-BINDING-RECEIPT-2026-09-23-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-Shellharbour-S0B-route/1.0"}
S=requests.Session();S.headers.update(UA)

META=REC["pair_identity"]["backscatter_product"]["meta_url"]
BATHY=REC["pair_identity"]["bathymetry_product"]["data_url"]
DERIVED=[
 BATHY.replace("/bathymetry/","/backscatter/").replace("Shellharbour-Bathymetry-2m-5m-2017.zip","Shellharbour-Backscatter-5m-2017.zip"),
 BATHY.replace("/bathymetry/","/backscatter/").replace("Shellharbour-Bathymetry-2m-5m-2017.zip","Shellharbour-Backscatter-Wide-5m-2017.zip"),
]

def fetch(u,method="GET",stream=False):
    r=S.request(method,u,timeout=120,allow_redirects=True,stream=stream)
    return r

meta=fetch(META)
meta_record={"requested":META,"status":meta.status_code,"resolved":meta.url,"bytes":len(meta.content),"sha256":hashlib.sha256(meta.content).hexdigest(),"content_type":meta.headers.get("content-type")}
text=meta.text if "text" in (meta.headers.get("content-type") or "").lower() or "html" in (meta.headers.get("content-type") or "").lower() else ""
urls=set()
for h in re.findall(r'href=["\']([^"\']+)["\']',text,re.I):
    urls.add(urljoin(meta.url,h))
for u in re.findall(r'https?://[^\s<>"\']+',text,re.I):
    urls.add(u.rstrip(").,;"))

# Search current data.gov.au catalogue by exact survey/title identifiers; metadata only.
ckan_url="https://www.data.gov.au/data/api/3/action/package_search"
ckan_params={"q":'"Shellharbour Tharawal" backscatter NSWENV20171130'}
try:
    c=S.get(ckan_url,params=ckan_params,timeout=120)
    ckan_status=c.status_code
    cj=c.json() if c.ok else {}
except Exception as e:
    c=None;cj={};ckan_status=None;ckan_error=type(e).__name__+": "+str(e)
else:
    ckan_error=None

ckan_resources=[]
for pkg in (cj.get("result") or {}).get("results",[]) if isinstance(cj,dict) else []:
    for res in pkg.get("resources",[]) or []:
        u=res.get("url")
        if u:
            ckan_resources.append({
              "package_id":pkg.get("id"),"package_title":pkg.get("title"),"resource_name":res.get("name"),
              "format":res.get("format"),"url":u
            })
            urls.add(u)

# Header-only candidate route probes. No response body for ZIP candidates.
probes=[]
candidates=set(DERIVED)
for u in urls:
    ul=u.lower()
    if ("shellharbour" in ul or "nswenv20171130" in ul or "146308" in ul) and ("backscatter" in ul or ul.endswith(".zip")):
        candidates.add(u)

for u in sorted(candidates):
    try:
        r=fetch(u,method="HEAD",stream=True)
        probes.append({"url":u,"status":r.status_code,"resolved":r.url,"content_type":r.headers.get("content-type"),"content_length":r.headers.get("content-length")})
    except Exception as e:
        probes.append({"url":u,"error":type(e).__name__+": "+str(e)})

# Unique direct-download candidate: successful URL ending zip or resolving to zip/media with non-html content.
success=[]
for p in probes:
    if not (200 <= int(p.get("status",0)) < 400):continue
    ru=(p.get("resolved") or p["url"]).lower()
    ct=(p.get("content_type") or "").lower()
    if ".zip" in ru and "text/html" not in ct:
        success.append(p)
# Canonicalize aliases by resolved URL.
by_resolved={}
for p in success:
    by_resolved[p.get("resolved") or p["url"]]=p
success=list(by_resolved.values())

out={
 "artifact_id":"JANUS-KUSTO-SHELLHARBOUR-S0B-BACKSCATTER-ROUTE-RUN-2026-09-23-v1.0",
 "source_receipt":REC["artifact_id"],
 "meta_page":meta_record,
 "ckan":{"status":ckan_status,"error":ckan_error,"resource_count":len(ckan_resources),"resources":ckan_resources},
 "header_only_probes":probes,
 "successful_direct_zip_routes":success,
 "resolved_unique_direct_zip":success[0]["resolved"] if len(success)==1 else None,
 "status":"PASS_UNIQUE_BACKSCATTER_DOWNLOAD_ROUTE_BOUND" if len(success)==1 else "BLOCKED_NONUNIQUE_OR_MISSING_BACKSCATTER_ROUTE",
 "raster_values_read":False,
 "numeric_zip_body_read":False,
 "next_gate":"S1_BATHYMETRY_ONLY_BLIND_CANDIDATE_FREEZE" if len(success)==1 else "SOURCE_ROUTE_RESOLUTION_ONLY",
 "claim_ceiling":"BACKSCATTER_DOWNLOAD_ROUTE_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-SHELLHARBOUR-S0B-BACKSCATTER-ROUTE-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps({
 "artifact_id":out["artifact_id"],"status":out["status"],"meta_page":meta_record,
 "ckan_resources":ckan_resources,"successful_direct_zip_routes":success,
 "raster_values_read":False,"numeric_zip_body_read":False
},indent=2,ensure_ascii=False))
