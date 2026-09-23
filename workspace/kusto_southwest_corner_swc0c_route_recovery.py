#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOUTHWEST-CORNER-SWC0-UNCHANGED-V2-PREREG-2026-09-24-v1.0.json").read_text())
BLOCK=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOUTHWEST-CORNER-SWC0B-ACCESS-BLOCK-2026-09-24-v1.0.json").read_text())
UA={"User-Agent":"Mozilla/5.0 JANUS-KUSTO-SWC0C-route-recovery/1.0"}
S=requests.Session();S.headers.update(UA)
API="https://www.data.gov.au/data/api/3/action/package_search"

queries=[
 '"South-west Corner Marine Park survey" "bathymetry" "GA-4858"',
 '"South-west Corner Marine Park survey" "backscatter" "GA-4858"'
]
pkgs={}
for q in queries:
    r=S.get(API,params={"q":q,"rows":100},timeout=120);r.raise_for_status()
    for p in (r.json().get("result") or {}).get("results",[]) or []:
        title=(p.get("title") or "").lower()
        if "south-west corner marine park survey" not in title or "ga-4858" not in title:
            continue
        pkgs[p["id"]]=p

records=[]
for p in pkgs.values():
    role="bathymetry" if "bathymetry" in (p.get("title") or "").lower() else ("backscatter" if "backscatter" in (p.get("title") or "").lower() else "other")
    for res in p.get("resources",[]) or []:
        records.append({
          "role":role,
          "package_id":p.get("id"),
          "package_title":p.get("title"),
          "package_metadata_modified":p.get("metadata_modified"),
          "organization":(p.get("organization") or {}).get("title"),
          "resource_id":res.get("id"),
          "resource_name":res.get("name"),
          "format":res.get("format"),
          "url":res.get("url"),
          "url_type":res.get("url_type"),
          "resource_created":res.get("created"),
          "resource_last_modified":res.get("last_modified"),
          "resource_state":res.get("state")
        })

# Deduplicate catalogue mirrors by semantic resource URL and keep newest package record as provenance.
physical={}
for x in records:
    u=x.get("url")
    if not u: continue
    key=(x["role"],u)
    prev=physical.get(key)
    if prev is None or str(x.get("package_metadata_modified") or "")>str(prev.get("package_metadata_modified") or ""):
        physical[key]=x
physical=list(physical.values())

# Identify data resources only. Do not use metadata archives as raster routes.
def is_data(x):
    t=" ".join(str(x.get(k) or "") for k in ("resource_name","format","url")).lower()
    if "metadata" in t:return False
    if x["role"]=="bathymetry":
        return "geotif" in t or "geotiff" in t or ".tif" in t or ("zip" in t and "bathym" in (x.get("package_title") or "").lower())
    if x["role"]=="backscatter":
        return "geotif" in t or "geotiff" in t or ".tif" in t or ("zip" in t and "backscatter" in (x.get("package_title") or "").lower())
    return False

candidates=[x for x in physical if is_data(x)]
probes=[]
for x in candidates:
    headers={"Referer":"https://www.data.gov.au/","Accept":"*/*"}
    try:
        h=S.head(x["url"],headers=headers,timeout=90,allow_redirects=True)
        probes.append({
          "resource":x,
          "method":"HEAD",
          "status":h.status_code,
          "resolved":h.url,
          "content_type":h.headers.get("content-type"),
          "content_length":h.headers.get("content-length"),
          "etag":h.headers.get("etag"),
          "last_modified":h.headers.get("last-modified")
        })
    except Exception as e:
        probes.append({"resource":x,"method":"HEAD","error":type(e).__name__+": "+str(e)})

# Also query GA/AODN catalogue landing pages only; response bodies may be metadata HTML/XML, never raster bytes.
metadata_urls=[
 "https://catalogue.aodn.org.au/geonetwork/srv/eng/csw/dataset/south-west-corner-marine-park-survey-bathymetry-ga-4858",
 "https://catalogue.aodn.org.au/geonetwork/srv/eng/csw/dataset/south-west-corner-marine-park-survey-backscatter-ga-4858",
 "https://pid.geoscience.gov.au/dataset/ga/145281",
 "https://pid.geoscience.gov.au/dataset/ga/145282"
]
metadata_probes=[]
for u in metadata_urls:
    try:
        r=S.get(u,timeout=90,allow_redirects=True)
        ct=(r.headers.get("content-type") or "").lower()
        body=r.content if ("html" in ct or "xml" in ct or "json" in ct or "text" in ct) and len(r.content)<5_000_000 else b""
        metadata_probes.append({
          "requested":u,"status":r.status_code,"resolved":r.url,
          "content_type":r.headers.get("content-type"),"bytes":len(r.content),
          "metadata_body_sha256":hashlib.sha256(body).hexdigest() if body else None,
          "body_sample":body.decode("utf-8","replace")[:20000] if body else None
        })
    except Exception as e:
        metadata_probes.append({"requested":u,"error":type(e).__name__+": "+str(e)})

out={
 "artifact_id":"JANUS-KUSTO-SOUTHWEST-CORNER-SWC0C-CURRENT-ROUTE-RECOVERY-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],
 "prior_access_block":BLOCK["artifact_id"],
 "queries":queries,
 "catalog_package_count":len(pkgs),
 "catalog_records":records,
 "deduped_physical_resources":physical,
 "data_resource_candidates":candidates,
 "header_only_data_route_probes":probes,
 "metadata_endpoint_probes":metadata_probes,
 "raster_values_read":False,
 "raster_body_downloaded":False,
 "underwater_imagery_read":False,
 "scientific_method_changed":False,
 "claim_ceiling":"CURRENT_AUTHORITATIVE_ROUTE_RECOVERY_ONLY"
}
p=OUT/"JANUS-KUSTO-SOUTHWEST-CORNER-SWC0C-CURRENT-ROUTE-RECOVERY-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "catalog_package_count":len(pkgs),
 "candidates":candidates,
 "probes":probes,
 "metadata_endpoints":[{k:v for k,v in x.items() if k!="body_sample"} for x in metadata_probes],
 "raster_values_read":False
},indent=2,ensure_ascii=False))
