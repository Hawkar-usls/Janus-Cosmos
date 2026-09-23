#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,time,xml.etree.ElementTree as ET
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-MONEY-SHOAL-UNSEEN-V2-PREREG-2026-09-23-v1.0.json").read_text())
REP=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-A0-BATHYMETRY-ROUTE-REPAIR-2026-09-23-v1.0.json").read_text())
WFS="https://warehouse.ausseabed.gov.au/geoserver/ows"
UA={"User-Agent":"JANUS-KUSTO-Arafura-A0B-WFS/1.0","Accept-Encoding":"identity"}
S=requests.Session();S.headers.update(UA)

def get(params,timeout=120):
    last=None
    for a in range(5):
        try:
            r=S.get(WFS,params=params,timeout=timeout,allow_redirects=True)
            r.raise_for_status();return r
        except Exception as e:
            last=e;time.sleep(min(12,2**a))
    raise last

def lname(t):return t.split("}",1)[-1]
def childtext(el,name):
    for c in list(el):
        if lname(c.tag)==name:return (c.text or "").strip()
    return ""

cap=get({"SERVICE":"WFS","REQUEST":"GetCapabilities","VERSION":"2.0.0"})
root=ET.fromstring(cap.content)
fts=[]
for ft in root.iter():
    if lname(ft.tag)=="FeatureType":
        n=childtext(ft,"Name");t=childtext(ft,"Title");a=childtext(ft,"Abstract")
        if n:fts.append({"name":n,"title":t,"abstract":a})

tokens=["arafura","bathymetry","2021","6m"]
matches=[x for x in fts if all(t in (" ".join(x.values())).lower() for t in tokens)]
# Prefer L0 coverage because it carries acquisition/product metadata, not a rendered style.
matches=[x for x in matches if "l0 coverage" in x["title"].lower() or "l0_coverage" in x["name"].lower()]
bindings=[]
for ft in matches:
    r=get({"SERVICE":"WFS","REQUEST":"GetFeature","VERSION":"2.0.0","TYPENAMES":ft["name"],
           "COUNT":"1","OUTPUTFORMAT":"application/json"})
    j=r.json();feats=j.get("features") or []
    bindings.append({
      "feature_type":ft,
      "url":r.url,"sha256":hashlib.sha256(r.content).hexdigest(),"bytes":len(r.content),
      "feature_count":len(feats),
      "properties":(feats[0].get("properties") or {}) if feats else {}
    })

numeric_routes=[]
for b in bindings:
    p=b["properties"]
    for k,v in p.items():
        if not isinstance(v,str):continue
        kl=str(k).lower();vl=v.lower()
        if ("url" in kl or "data" in kl) and ("http://" in vl or "https://" in vl):
            if any(x in vl for x in [".tif",".tiff",".zip","download"]):
                numeric_routes.append({"field":k,"url":v,"layer":b["feature_type"]["name"]})
# deterministic de-dup
seen=set();routes=[]
for x in numeric_routes:
    if x["url"] not in seen:
        seen.add(x["url"]);routes.append(x)

probes=[]
for x in routes:
    try:
        h=S.head(x["url"],timeout=60,allow_redirects=True)
        probes.append({**x,"status":h.status_code,"resolved":h.url,
                       "content_type":h.headers.get("content-type"),"content_length":h.headers.get("content-length")})
    except Exception as e:
        probes.append({**x,"error":type(e).__name__+": "+str(e)})
reachable=[x for x in probes if 200<=int(x.get("status",0))<400]

out={
 "artifact_id":"JANUS-KUSTO-ARAFURA-A0B-AUSSEABED-WFS-ROUTE-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],"repair":REP["artifact_id"],
 "wfs_capabilities":{"url":cap.url,"sha256":hashlib.sha256(cap.content).hexdigest(),"bytes":len(cap.content)},
 "matching_layers":matches,"bindings":bindings,"numeric_route_candidates":routes,"header_probes":probes,
 "bathymetry_values_read":False,"backscatter_values_read":False,
 "status":"PASS_NUMERIC_BATHYMETRY_ROUTE_BOUND__NO_RASTER_VALUES_READ" if len(reachable)==1 else "NUMERIC_ROUTE_NOT_UNIQUELY_REACHABLE",
 "next_gate":"A1_MONEY_SHOAL_BATHYMETRY_ONLY_BLIND_CANDIDATE_FREEZE" if len(reachable)==1 else "SOURCE_ROUTE_REPAIR_ONLY",
 "claim_ceiling":"SOURCE_FILE_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-ARAFURA-A0B-AUSSEABED-WFS-ROUTE-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps({"artifact_id":out["artifact_id"],"status":out["status"],
 "matching_layers":matches,"routes":routes,"probes":probes,"raster_values_read":False},indent=2,ensure_ascii=False))
