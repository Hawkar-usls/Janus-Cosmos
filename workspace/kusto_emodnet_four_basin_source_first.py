#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,re,time,xml.etree.ElementTree as ET
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-EMODNET-FOUR-BASIN-SOURCE-FIRST-PREREG-2026-09-23-v1.0.json").read_text())
WMS=PRE["services"]["wms"];WFS=PRE["services"]["wfs"]
UA={"User-Agent":"JANUS-KUSTO-EMODnet-four-basin-source-first/1.0","Accept-Encoding":"identity"}
S=requests.Session();S.headers.update(UA)

def get(url,params,timeout=180):
    last=None
    for a in range(5):
        try:
            r=S.get(url,params=params,timeout=timeout,allow_redirects=True)
            r.raise_for_status();return r
        except Exception as e:
            last=e;time.sleep(min(15,2**a))
    raise last

def lname(t):return t.split("}",1)[-1]
def childtext(el,name):
    for c in list(el):
        if lname(c.tag)==name:return (c.text or "").strip()
    return ""

cap=get(WMS,{"SERVICE":"WMS","REQUEST":"GetCapabilities","VERSION":"1.1.1"})
root=ET.fromstring(cap.content)
layers=[]
for el in root.iter():
    if lname(el.tag)=="Layer":
        name=childtext(el,"Name");title=childtext(el,"Title");ab=childtext(el,"Abstract")
        if name:layers.append({"name":name,"title":title,"abstract":ab,"queryable":el.attrib.get("queryable")})
def is_source(x):
    t=(" ".join([x["name"],x["title"],x["abstract"]])).lower()
    return any(k in t for k in ["source","cdi","survey","reference","quality"])
src_layers=[x for x in layers if is_source(x) and str(x.get("queryable","0")).lower() in ("1","true")]

wfs=None;wfs_types=[];wfs_error=None
try:
    wfs=get(WFS,{"SERVICE":"WFS","REQUEST":"GetCapabilities","VERSION":"2.0.0"})
    wr=ET.fromstring(wfs.content)
    for ft in wr.iter():
        if lname(ft.tag)=="FeatureType":
            n=childtext(ft,"Name");t=childtext(ft,"Title");a=childtext(ft,"Abstract")
            if n:wfs_types.append({"name":n,"title":t,"abstract":a})
except Exception as e:
    wfs_error=type(e).__name__+": "+str(e)

probes={}
for p in PRE["fixed_probe_points"]:
    rows=[]
    for l in src_layers:
        pad=.02
        params={
          "SERVICE":"WMS","VERSION":"1.1.1","REQUEST":"GetFeatureInfo",
          "LAYERS":l["name"],"QUERY_LAYERS":l["name"],"SRS":"EPSG:4326",
          "BBOX":f'{p["lon"]-pad},{p["lat"]-pad},{p["lon"]+pad},{p["lat"]+pad}',
          "WIDTH":"101","HEIGHT":"101","X":"50","Y":"50","FORMAT":"image/png",
          "INFO_FORMAT":"text/plain","FEATURE_COUNT":"50"
        }
        try:
            r=get(WMS,params,90);txt=r.text[:50000]
            informative=bool(txt.strip()) and not re.search(r"(no features|serviceexception|application/vnd.ogc.se_xml)",txt,re.I)
            rows.append({"layer":l["name"],"title":l["title"],"status":r.status_code,"bytes":len(r.content),
                         "sha256":hashlib.sha256(r.content).hexdigest(),"informative":informative,
                         "response":txt if informative else txt[:2000]})
        except Exception as e:
            rows.append({"layer":l["name"],"title":l["title"],"error":type(e).__name__+": "+str(e),"informative":False})
    probes[p["id"]]={"point":p,"queries":rows}

out={
 "artifact_id":"JANUS-KUSTO-EMODNET-FOUR-BASIN-SOURCE-FIRST-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "wms":{"url":cap.url,"sha256":hashlib.sha256(cap.content).hexdigest(),"bytes":len(cap.content),
        "layer_count":len(layers),"source_candidate_layers":src_layers},
 "wfs":{"url":wfs.url if wfs else None,"sha256":hashlib.sha256(wfs.content).hexdigest() if wfs else None,
        "feature_types":wfs_types,"error":wfs_error},
 "probes":probes,
 "bathymetric_depth_values_intentionally_not_requested":True,
 "composite_dtm_counted_as_independent_survey":False,
 "next_gate":"EXTRACT_EXPLICIT_CDI_OR_SURVEY_IDENTITIES_FROM_INFORMATIVE_RESPONSES",
 "claim_ceiling":"EUROPEAN_SURVEY_LINEAGE_DISCOVERY_ONLY"
}
p=OUT/"JANUS-KUSTO-EMODNET-FOUR-BASIN-SOURCE-FIRST-RUN-2026-09-23-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "source_layers":[{"name":x["name"],"title":x["title"]} for x in src_layers],
 "wfs_feature_types":[x for x in wfs_types if any(k in (" ".join(x.values())).lower() for k in ["source","cdi","survey","quality"])],
 "informative":{
   k:[{"layer":q["layer"],"title":q["title"],"response":q.get("response","")[:3000]} for q in v["queries"] if q.get("informative")]
   for k,v in probes.items()
 },
 "depth_values_requested":False
},indent=2,ensure_ascii=False))
