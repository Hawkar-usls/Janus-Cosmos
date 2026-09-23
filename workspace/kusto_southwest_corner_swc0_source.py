#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOUTHWEST-CORNER-SWC0-UNCHANGED-V2-PREREG-2026-09-24-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-SouthwestCorner-SWC0-source/1.0"}
S=requests.Session();S.headers.update(UA)
API="https://www.data.gov.au/data/api/3/action/package_search"

queries=['"South-west Corner Marine Park" "GA-4858"','"South-west Corner" bathymetry backscatter']
packages={}
for q in queries:
    r=S.get(API,params={"q":q},timeout=120);r.raise_for_status()
    for p in (r.json().get("result") or {}).get("results",[]) or []:
        packages[p.get("id")]=p

def collect(role):
    out=[]
    for p in packages.values():
        title=(p.get("title") or "").lower()
        if "south-west corner marine park" not in title: continue
        if role not in title: continue
        for res in p.get("resources",[]) or []:
            out.append({
              "package_id":p.get("id"),"package_title":p.get("title"),
              "resource_id":res.get("id"),"resource_name":res.get("name"),
              "format":res.get("format"),"url":res.get("url")
            })
    return out

bathy=collect("bathymetry");back=collect("backscatter")

def semantic_candidates(rows,role):
    good=[]
    for x in rows:
        t=" ".join(str(x.get(k) or "") for k in ("resource_name","format","url")).lower()
        if role not in t and role not in (x.get("package_title") or "").lower(): continue
        if "metadata" in t: continue
        if not x.get("url"): continue
        if any(k in t for k in ["geotif","geotiff",".tif",".tiff","ascii","geotifs","download"]):
            good.append(x)
    return good

bathy_c=semantic_candidates(bathy,"bathymetry")
back_c=semantic_candidates(back,"backscatter")

probes=[]
for role,rows in [("bathymetry",bathy_c),("backscatter",back_c)]:
    for x in rows:
        try:
            h=S.head(x["url"],timeout=90,allow_redirects=True)
            probes.append({"role":role,"resource":x,"status":h.status_code,"resolved":h.url,
              "content_type":h.headers.get("content-type"),"content_length":h.headers.get("content-length")})
        except Exception as e:
            probes.append({"role":role,"resource":x,"error":type(e).__name__+": "+str(e)})

out={
 "artifact_id":"JANUS-KUSTO-SOUTHWEST-CORNER-SWC0-SOURCE-BINDING-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],
 "queries":queries,
 "package_count":len(packages),
 "bathymetry_resources":bathy,
 "backscatter_resources":back,
 "semantic_bathymetry_candidates":bathy_c,
 "semantic_backscatter_candidates":back_c,
 "header_only_probes":probes,
 "raster_values_read":False,
 "underwater_imagery_read":False,
 "next_gate":"SWC0B_EXACT_RESOURCE_AND_ARCHIVE_MEMBER_BINDING",
 "claim_ceiling":"CURRENT_SOURCE_RESOURCE_DISCOVERY_ONLY"
}
p=OUT/"JANUS-KUSTO-SOUTHWEST-CORNER-SWC0-SOURCE-BINDING-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "package_count":len(packages),
 "bathymetry_candidates":bathy_c,
 "backscatter_candidates":back_c,
 "probes":probes,
 "raster_values_read":False,
 "underwater_imagery_read":False
},indent=2,ensure_ascii=False))
