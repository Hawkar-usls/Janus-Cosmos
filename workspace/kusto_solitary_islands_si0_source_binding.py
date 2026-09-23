#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,xml.etree.ElementTree as ET
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOLITARY-ISLANDS-SI0-UNCHANGED-V2-PREREG-2026-09-24-v1.0.json").read_text())
WFS="https://warehouse.ausseabed.gov.au/geoserver/ows"
UA={"User-Agent":"JANUS-KUSTO-SolitaryIslands-SI0-source/1.0"}
S=requests.Session();S.headers.update(UA)

def get(params):
    r=S.get(WFS,params=params,timeout=180,allow_redirects=True);r.raise_for_status();return r
def lname(tag):return tag.split("}",1)[-1]
def childtext(el,name):
    for c in list(el):
        if lname(c.tag)==name:return (c.text or "").strip()
    return ""

cap=get({"SERVICE":"WFS","REQUEST":"GetCapabilities","VERSION":"2.0.0"})
root=ET.fromstring(cap.content)
fts=[]
for ft in root.iter():
    if lname(ft.tag)!="FeatureType":continue
    row={"name":childtext(ft,"Name"),"title":childtext(ft,"Title"),"abstract":childtext(ft,"Abstract")}
    if row["name"]:fts.append(row)

sol=[x for x in fts if "solitary" in (" ".join(x.values())).lower()]
requirements={
 "BATHYMETRY_5M":["solitary","islands","bathymetry","2022","5m"],
 "BACKSCATTER_5M":["solitary","islands","backscatter","2022","5m"]
}
bindings={}
for role,tokens in requirements.items():
    exact=[x for x in sol if all(t in (" ".join([x["name"],x["title"],x["abstract"]])).lower() for t in tokens)]
    if len(exact)!=1:
        bindings[role]={"status":"NOT_UNIQUELY_RESOLVED","tokens":tokens,"matches":exact}
        continue
    ft=exact[0]
    desc=get({"SERVICE":"WFS","REQUEST":"DescribeFeatureType","VERSION":"2.0.0","TYPENAMES":ft["name"]})
    feat=get({"SERVICE":"WFS","REQUEST":"GetFeature","VERSION":"2.0.0","TYPENAMES":ft["name"],"COUNT":"1","OUTPUTFORMAT":"application/json"})
    props={}
    try:
        fj=feat.json()
        props=(fj.get("features") or [{}])[0].get("properties") or {}
    except Exception:
        pass
    bindings[role]={
      "status":"RESOLVED",
      "feature_type":ft,
      "describe_url":desc.url,
      "describe_sha256":hashlib.sha256(desc.content).hexdigest(),
      "first_feature_properties":props,
      "feature_body_sha256":hashlib.sha256(feat.content).hexdigest(),
      "first_feature_geometry_intentionally_omitted":True
    }

all_resolved=all(v.get("status")=="RESOLVED" for v in bindings.values())
out={
 "artifact_id":"JANUS-KUSTO-SOLITARY-ISLANDS-SI0-SOURCE-BINDING-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],
 "wfs_capabilities":{"url":cap.url,"bytes":len(cap.content),"sha256":hashlib.sha256(cap.content).hexdigest()},
 "solitary_feature_types":sol,
 "bindings":bindings,
 "all_required_layers_resolved":all_resolved,
 "raster_values_read":False,
 "towed_video_or_stills_read":False,
 "sediment_labels_read":False,
 "landform_or_substrate_classification_read":False,
 "next_gate":"SI0B_CURRENT_DOWNLOAD_ROUTE_AND_ARCHIVE_MEMBER_BINDING" if all_resolved else "SOURCE_LAYER_RESOLUTION_ONLY",
 "claim_ceiling":"SOURCE_METADATA_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-SOLITARY-ISLANDS-SI0-SOURCE-BINDING-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "resolved":{k:v.get("status") for k,v in bindings.items()},
 "properties":{k:v.get("first_feature_properties") for k,v in bindings.items() if v.get("status")=="RESOLVED"},
 "raster_values_read":False,
 "truth_layers_read":False
},indent=2,ensure_ascii=False))
