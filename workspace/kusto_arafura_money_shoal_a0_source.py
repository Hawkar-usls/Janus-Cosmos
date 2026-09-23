#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,re,time
from pathlib import Path
from urllib.parse import urljoin
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-MONEY-SHOAL-UNSEEN-V2-PREREG-2026-09-23-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-Arafura-A0-source/1.0","Accept-Encoding":"identity"}
S=requests.Session();S.headers.update(UA)

DATASET_ID="47e1ed7e-9182-4f3a-b8a3-779f7bbc774e"
CKAN=f"https://data.gov.au/data/api/3/action/package_show?id={DATASET_ID}"
LANDING=PRE["dataset"]["bathymetry"]["public_dataset"]
BACK=PRE["dataset"]["backscatter"]["direct_archive"]

def fetch(u,method="get",tries=5):
    last=None
    for a in range(tries):
        try:
            if method=="head":r=S.head(u,timeout=60,allow_redirects=True)
            else:r=S.get(u,timeout=120,allow_redirects=True)
            return r
        except Exception as e:
            last=e;time.sleep(min(10,2**a))
    raise last

ckan={"url":CKAN}
resources=[]
try:
    r=fetch(CKAN)
    ckan.update({"status":r.status_code,"resolved":r.url,"bytes":len(r.content),"sha256":hashlib.sha256(r.content).hexdigest()})
    if r.ok:
        j=r.json()
        for x in (j.get("result") or {}).get("resources") or []:
            resources.append({
              "name":x.get("name"),"format":x.get("format"),"url":x.get("url"),
              "description":x.get("description"),"id":x.get("id"),"size":x.get("size")
            })
except Exception as e:
    ckan["error"]=type(e).__name__+": "+str(e)

landing={"url":LANDING}
hrefs=[]
try:
    r=fetch(LANDING)
    landing.update({"status":r.status_code,"resolved":r.url,"bytes":len(r.content),"sha256":hashlib.sha256(r.content).hexdigest()})
    if r.ok:
        hrefs=sorted(set(urljoin(r.url,h) for h in re.findall(r'href=["\']([^"\']+)["\']',r.text,re.I)))
except Exception as e:
    landing["error"]=type(e).__name__+": "+str(e)

# Use metadata text only. A valid bathymetry resource must be an explicit raster resource,
# not a metadata/report ZIP.
bathy=[]
for x in resources:
    t=" ".join(str(x.get(k) or "") for k in ("name","format","description","url")).lower()
    if ("bathym" in t or "arafura marine park bathymetry" in t) and ("tif" in t or "geotiff" in t):
        bathy.append(x)
if not bathy:
    for u in hrefs:
        t=u.lower()
        if ("bathym" in t or "arafura" in t) and (".tif" in t or "resource" in t):
            bathy.append({"name":"LANDING_DISCOVERED","format":None,"url":u,"description":None,"id":None,"size":None})

# Deduplicate by URL.
uniq={}
for x in bathy:
    if x.get("url"):uniq[x["url"]]=x
bathy=list(uniq.values())

probes=[]
for x in bathy:
    try:
        h=fetch(x["url"],"head")
        probes.append({"role":"BATHYMETRY","url":x["url"],"status":h.status_code,"resolved":h.url,
                       "content_type":h.headers.get("content-type"),"content_length":h.headers.get("content-length")})
    except Exception as e:
        probes.append({"role":"BATHYMETRY","url":x["url"],"error":type(e).__name__+": "+str(e)})
try:
    h=fetch(BACK,"head")
    backprobe={"role":"BACKSCATTER","url":BACK,"status":h.status_code,"resolved":h.url,
               "content_type":h.headers.get("content-type"),"content_length":h.headers.get("content-length")}
except Exception as e:
    backprobe={"role":"BACKSCATTER","url":BACK,"error":type(e).__name__+": "+str(e)}
probes.append(backprobe)

reachable_bathy=[p for p in probes if p["role"]=="BATHYMETRY" and 200<=int(p.get("status",0))<400]
back_ok=200<=int(backprobe.get("status",0))<400

out={
 "artifact_id":"JANUS-KUSTO-ARAFURA-MONEY-SHOAL-A0-SOURCE-BINDING-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "selected_subarea":PRE["dataset"]["selected_subarea"],
 "selection_rule":PRE["dataset"]["selection_rule"],
 "ckan":ckan,"ckan_resources":resources,"landing":landing,
 "candidate_bathymetry_resources":bathy,"header_probes":probes,
 "backscatter_archive":BACK,
 "bathymetry_values_read":False,"backscatter_values_read":False,
 "numeric_raster_body_read":False,"post_survey_feature_truth_read":False,
 "status":"PASS_SOURCE_ROUTES_BOUND__NO_RASTER_VALUES_READ" if len(reachable_bathy)==1 and back_ok else "SOURCE_BINDING_BLOCKED",
 "next_gate":"A1_MONEY_SHOAL_BATHYMETRY_ONLY_BLIND_CANDIDATE_FREEZE" if len(reachable_bathy)==1 and back_ok else "SOURCE_ROUTE_REPAIR_ONLY",
 "claim_ceiling":"SOURCE_FILE_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-ARAFURA-MONEY-SHOAL-A0-SOURCE-BINDING-RUN-2026-09-23-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"status":out["status"],
 "bathymetry_candidates":bathy,"probes":probes,
 "raster_values_read":False
},indent=2,ensure_ascii=False))
