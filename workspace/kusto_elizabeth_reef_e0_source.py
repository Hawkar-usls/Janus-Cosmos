#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,re,time
from pathlib import Path
from urllib.parse import urljoin
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-REEF-UNSEEN-V2-PREREG-2026-09-23-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-Elizabeth-Reef-E0-source/1.0","Accept-Encoding":"identity"}
S=requests.Session();S.headers.update(UA)

BATH=PRE["dataset"]["bathymetry"]["direct_archive"]
BACK_PID=PRE["dataset"]["backscatter"]["metadata_pid"]

def fetch(u,method="get",tries=5):
    last=None
    for a in range(tries):
        try:
            r=S.head(u,timeout=60,allow_redirects=True) if method=="head" else S.get(u,timeout=120,allow_redirects=True)
            return r
        except Exception as e:
            last=e;time.sleep(min(10,2**a))
    raise last

def probe(u,role):
    try:
        h=fetch(u,"head")
        return {"role":role,"url":u,"status":h.status_code,"resolved":h.url,
                "content_type":h.headers.get("content-type"),"content_length":h.headers.get("content-length")}
    except Exception as e:
        return {"role":role,"url":u,"error":type(e).__name__+": "+str(e)}

bath_probe=probe(BATH,"BATHYMETRY")
pages=[]
hrefs=set()
for u in [BACK_PID]:
    try:
        r=fetch(u)
        pages.append({"url":u,"status":r.status_code,"resolved":r.url,"bytes":len(r.content),
                      "sha256":hashlib.sha256(r.content).hexdigest(),"content_type":r.headers.get("content-type")})
        if r.ok:
            for h in re.findall(r'href=["\']([^"\']+)["\']',r.text,re.I):
                hrefs.add(urljoin(r.url,h))
            for x in re.findall(r'https?://[^\s<>"\']+',r.text,re.I):
                hrefs.add(x.rstrip(").,;"))
    except Exception as e:
        pages.append({"url":u,"error":type(e).__name__+": "+str(e)})

back_routes=[]
for u in sorted(hrefs):
    t=u.lower()
    if ("backscatter" in t or "144416" in t or "ga4848" in t) and (".zip" in t or ".tif" in t or ".tiff" in t):
        back_routes.append(u)

# Metadata-preserving fallback: probe only the exact conventional AusSeabed filename
# derived from the officially published bathymetry filename and the same survey identity.
fallback="https://files.ausseabed.gov.au/survey/Elizabeth%20Middleton%20Reef%20Backscatter%202020%205m.zip"
if fallback not in back_routes:
    back_routes.append(fallback)

back_probes=[probe(u,"BACKSCATTER") for u in back_routes]
reachable_back=[x for x in back_probes if 200<=int(x.get("status",0))<400]
bath_ok=200<=int(bath_probe.get("status",0))<400

out={
 "artifact_id":"JANUS-KUSTO-ELIZABETH-REEF-E0-SOURCE-BINDING-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "selected_subarea":PRE["dataset"]["selected_subarea"],
 "bathymetry_probe":bath_probe,
 "backscatter_metadata_pages":pages,
 "backscatter_candidate_routes":back_routes,
 "backscatter_probes":back_probes,
 "status":"PASS_EXACT_PAIRED_ROUTES_BOUND__NO_RASTER_VALUES_READ" if bath_ok and len(reachable_back)==1 else "SOURCE_BINDING_BLOCKED",
 "bathymetry_values_read":False,"backscatter_values_read":False,"numeric_raster_body_read":False,
 "auv_bruv_sediment_truth_read":False,
 "next_gate":"E1_ELIZABETH_REEF_BATHYMETRY_ONLY_BLIND_CANDIDATE_FREEZE" if bath_ok and len(reachable_back)==1 else "SOURCE_ROUTE_REPAIR_ONLY",
 "claim_ceiling":"SOURCE_FILE_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-ELIZABETH-REEF-E0-SOURCE-BINDING-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps({"artifact_id":out["artifact_id"],"status":out["status"],
 "bathymetry":bath_probe,"backscatter":back_probes,"raster_values_read":False},indent=2,ensure_ascii=False))
