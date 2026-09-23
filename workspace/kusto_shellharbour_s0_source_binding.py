#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,xml.etree.ElementTree as ET
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHELLHARBOUR-S0-SOURCE-BINDING-PREREG-2026-09-23-v1.0.json").read_text())
WFS=PRE["source_resolution"]["wfs"]
UA={"User-Agent":"JANUS-KUSTO-Shellharbour-S0-source/1.0"}
S=requests.Session();S.headers.update(UA)

def get(params,timeout=180):
    r=S.get(WFS,params=params,timeout=timeout,allow_redirects=True)
    r.raise_for_status()
    return r

def lname(tag):return tag.split("}",1)[-1]
def childtext(el,name):
    for c in list(el):
        if lname(c.tag)==name:return (c.text or "").strip()
    return ""

cap=get({"SERVICE":"WFS","REQUEST":"GetCapabilities","VERSION":"2.0.0"})
cap_sha=hashlib.sha256(cap.content).hexdigest()
root=ET.fromstring(cap.content)
fts=[]
for ft in root.iter():
    if lname(ft.tag)=="FeatureType":
        row={"name":childtext(ft,"Name"),"title":childtext(ft,"Title"),"abstract":childtext(ft,"Abstract")}
        if row["name"]:fts.append(row)

shell=[x for x in fts if "shellharbour" in (" ".join(x.values())).lower()]
required_titles=PRE["source_resolution"]["required_feature_titles"]
bindings={}
for title in required_titles:
    exact=[x for x in shell if x["title"].strip().lower()==title.strip().lower()]
    if len(exact)!=1:
        bindings[title]={"status":"NOT_UNIQUELY_RESOLVED","matches":exact}
        continue
    ft=exact[0]
    desc=get({"SERVICE":"WFS","REQUEST":"DescribeFeatureType","VERSION":"2.0.0","TYPENAMES":ft["name"]})
    feat=get({"SERVICE":"WFS","REQUEST":"GetFeature","VERSION":"2.0.0","TYPENAMES":ft["name"],"COUNT":"1","OUTPUTFORMAT":"application/json"})
    try:
        fj=feat.json()
        props=(fj.get("features") or [{}])[0].get("properties") or {}
    except Exception:
        fj=None;props={}
    bindings[title]={
      "status":"RESOLVED",
      "feature_type":ft,
      "describe_url":desc.url,
      "describe_sha256":hashlib.sha256(desc.content).hexdigest(),
      "describe_bytes":len(desc.content),
      "feature_url":feat.url,
      "feature_sha256":hashlib.sha256(feat.content).hexdigest(),
      "feature_bytes":len(feat.content),
      "first_feature_properties":props,
      "first_feature_geometry_intentionally_omitted":True
    }

all_resolved=all(v.get("status")=="RESOLVED" for v in bindings.values())
out={
 "artifact_id":"JANUS-KUSTO-SHELLHARBOUR-S0-SOURCE-BINDING-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "wfs_capabilities":{"url":cap.url,"sha256":cap_sha,"bytes":len(cap.content)},
 "shellharbour_feature_types":shell,
 "bindings":bindings,
 "all_required_layers_resolved":all_resolved,
 "raster_values_read":False,
 "landform_classification_read":False,
 "next_gate":"S1_BATHYMETRY_ONLY_BLIND_CANDIDATE_FREEZE" if all_resolved else "SOURCE_ROUTE_RESOLUTION_ONLY",
 "claim_ceiling":"SOURCE_METADATA_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-SHELLHARBOUR-S0-SOURCE-BINDING-RUN-2026-09-23-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "shellharbour_feature_types":shell,
 "resolved":{k:v.get("status") for k,v in bindings.items()},
 "properties":{k:v.get("first_feature_properties") for k,v in bindings.items() if v.get("status")=="RESOLVED"},
 "raster_values_read":False
},indent=2,ensure_ascii=False))
