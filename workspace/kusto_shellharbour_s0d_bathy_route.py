#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
REC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHELLHARBOUR-S0-SOURCE-BINDING-RECEIPT-2026-09-23-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-Shellharbour-S0D-bathy-route/1.0"}
S=requests.Session();S.headers.update(UA)
API="https://www.data.gov.au/data/api/3/action/package_search"
queries=[
 '"Shellharbour Tharawal" bathymetry NSWENV20171130',
 '"Shellharbour" bathymetry 2m 2017'
]
packages=[]
for q in queries:
    r=S.get(API,params={"q":q},timeout=120);r.raise_for_status()
    j=r.json()
    for pkg in (j.get("result") or {}).get("results",[]) or []:
        packages.append(pkg)

seen={}
for p in packages: seen[p.get("id")]=p
packages=list(seen.values())
resources=[]
for pkg in packages:
    for res in pkg.get("resources",[]) or []:
        resources.append({
          "package_id":pkg.get("id"),
          "package_title":pkg.get("title"),
          "resource_name":res.get("name"),
          "format":res.get("format"),
          "url":res.get("url")
        })

def txt(x): return " ".join(str(x.get(k) or "") for k in ["package_title","resource_name","format","url"]).lower()

eligible=[]
for x in resources:
    t=txt(x)
    if "shellharbour" not in t: continue
    if "bathymetry" not in t: continue
    if "2m" not in t and "2 m" not in t: continue
    if (x.get("format") or "").lower() not in ("tif","tiff","geotiff","cog","zip"): continue
    eligible.append(x)

# Prefer direct raster over ZIP/archive.
direct=[x for x in eligible if (x.get("format") or "").lower() in ("tif","tiff","geotiff","cog") and x.get("url")]
chosen=direct if len(direct)==1 else []

probes=[]
for x in (chosen or eligible):
    u=x.get("url")
    if not u: continue
    try:
        h=S.head(u,timeout=90,allow_redirects=True)
        probes.append({"resource":x,"status":h.status_code,"resolved":h.url,"content_type":h.headers.get("content-type"),"content_length":h.headers.get("content-length")})
    except Exception as e:
        probes.append({"resource":x,"error":type(e).__name__+": "+str(e)})

success=[p for p in probes if 200<=int(p.get("status",0))<400]
unique = success[0] if len(success)==1 else None
out={
 "artifact_id":"JANUS-KUSTO-SHELLHARBOUR-S0D-BATHYMETRY-ROUTE-RUN-2026-09-23-v1.0",
 "source_receipt":REC["artifact_id"],
 "queries":queries,
 "package_count":len(packages),
 "resources":resources,
 "eligible_2m_bathymetry_resources":eligible,
 "header_only_probes":probes,
 "resolved_unique_resource":unique,
 "status":"PASS_UNIQUE_2M_BATHYMETRY_ROUTE_BOUND" if unique else "BLOCKED_NONUNIQUE_OR_MISSING_2M_BATHYMETRY_ROUTE",
 "raster_values_read":False,
 "next_gate":"S1_BATHYMETRY_ONLY_BLIND_CANDIDATE_FREEZE" if unique else "SOURCE_ROUTE_RESOLUTION_ONLY",
 "claim_ceiling":"BATHYMETRY_DOWNLOAD_ROUTE_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-SHELLHARBOUR-S0D-BATHYMETRY-ROUTE-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps({
 "artifact_id":out["artifact_id"],"status":out["status"],
 "eligible":eligible,"probes":probes,"raster_values_read":False
},indent=2,ensure_ascii=False))
