#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json, re, xml.etree.ElementTree as ET
from pathlib import Path
from urllib.parse import urlencode
import requests

OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
CANDS=[
 {"id":"KN19207_CAND_001","lat":-4.20546708989972,"lon":-14.871748111744372},
 {"id":"KN19207_CAND_002","lat":-3.9727527956056825,"lon":-12.272824298723462},
 {"id":"KN19207_CAND_003","lat":-4.015075679897318,"lon":-12.29915403590432},
 {"id":"KN19207_CAND_004","lat":-4.2399610081365955,"lon":-15.26855667523257},
 {"id":"KN19207_CAND_005","lat":-4.099231435240242,"lon":-12.332032012186781},
 {"id":"KN19207_CAND_006","lat":-4.309010339378846,"lon":-12.270140997573803},
 {"id":"KN19207_CAND_007","lat":-4.119167753526161,"lon":-14.356560209706116},
 {"id":"KN19207_CAND_008","lat":-4.3376789874153,"lon":-12.319260730325983}
]
WMS="https://geo-service.maris.nl/emodnet_bathymetry/wms"
WFS="https://ows.emodnet-bathymetry.eu/wfs"
UA={"User-Agent":"JANUS-KUSTO-EMODnet-CDI-metadata-audit/1.0"}
S=requests.Session();S.headers.update(UA)

def get(url,params,timeout=120):
    r=S.get(url,params=params,timeout=timeout)
    return {"status":r.status_code,"url":r.url,"headers":dict(r.headers),"bytes":len(r.content),
            "sha256":hashlib.sha256(r.content).hexdigest(),"content":r.content}

def lname(tag): return tag.split("}",1)[-1]

def childtext(el,name):
    for c in list(el):
        if lname(c.tag)==name:
            return (c.text or "").strip()
    return ""

def parse_layers(xml):
    root=ET.fromstring(xml)
    rows=[]
    def walk(el,parent_title=""):
        if lname(el.tag)=="Layer":
            name=childtext(el,"Name");title=childtext(el,"Title");ab=childtext(el,"Abstract")
            q=el.attrib.get("queryable")
            crs=[(c.text or "").strip() for c in list(el) if lname(c.tag) in ("CRS","SRS")]
            bboxes=[]
            for c in list(el):
                if lname(c.tag)=="EX_GeographicBoundingBox":
                    vals={}
                    for x in list(c):vals[lname(x.tag)]=(x.text or "").strip()
                    bboxes.append({"type":"geo",**vals})
                elif lname(c.tag)=="LatLonBoundingBox":
                    bboxes.append({"type":"latlon",**c.attrib})
            if name or title:
                rows.append({"name":name,"title":title,"abstract":ab,"queryable":q,"crs":crs,"bbox":bboxes,"parent_title":parent_title})
            p=title or parent_title
        else:p=parent_title
        for c in list(el):walk(c,p)
    walk(root)
    return rows

def source_candidate(l):
    txt=(" ".join([l.get("name",""),l.get("title",""),l.get("abstract","")])).lower()
    kws=["survey","cdi","source","track","polygon","bathymetric data set","dataset"]
    return bool(l.get("name")) and any(k in txt for k in kws)

cap=get(WMS,{"SERVICE":"WMS","REQUEST":"GetCapabilities","VERSION":"1.1.1"},180)
cap_record={k:v for k,v in cap.items() if k!="content"}
layers=[]
cap_error=None
if cap["status"]==200:
    try:layers=parse_layers(cap["content"])
    except Exception as e:cap_error=repr(e)
cand_layers=[x for x in layers if source_candidate(x)]

# WFS inventory only; no GetFeature yet.
wfs_cap=get(WFS,{"SERVICE":"WFS","REQUEST":"GetCapabilities","VERSION":"2.0.0"},180)
wfs_record={k:v for k,v in wfs_cap.items() if k!="content"}
wfs_featuretypes=[]
if wfs_cap["status"]==200:
    try:
        root=ET.fromstring(wfs_cap["content"])
        for ft in root.iter():
            if lname(ft.tag)=="FeatureType":
                row={"name":childtext(ft,"Name"),"title":childtext(ft,"Title"),"abstract":childtext(ft,"Abstract")}
                if row["name"]:wfs_featuretypes.append(row)
    except Exception as e:
        wfs_record["parse_error"]=repr(e)

# WMS point GetFeatureInfo on metadata-looking queryable layers only.
responses={}
for t in CANDS:
    tr=[]
    for l in cand_layers:
        if str(l.get("queryable","0")).lower() not in ("1","true"):
            continue
        pad=.05
        params={
          "SERVICE":"WMS","VERSION":"1.1.1","REQUEST":"GetFeatureInfo",
          "LAYERS":l["name"],"QUERY_LAYERS":l["name"],
          "SRS":"EPSG:4326","BBOX":f'{t["lon"]-pad},{t["lat"]-pad},{t["lon"]+pad},{t["lat"]+pad}',
          "WIDTH":"101","HEIGHT":"101","X":"50","Y":"50",
          "FORMAT":"image/png","INFO_FORMAT":"text/plain","FEATURE_COUNT":"20"
        }
        rr=get(WMS,params,90)
        text=rr["content"].decode("utf-8","replace")[:50000]
        # keep even blank/exception, but flag whether it looks informative
        info=bool(text.strip()) and not re.search(r'(no features|serviceexception)',text,re.I)
        tr.append({
          "layer":l["name"],"title":l["title"],"status":rr["status"],"sha256":rr["sha256"],
          "bytes":rr["bytes"],"informative":info,"response":text,"url":rr["url"]
        })
    responses[t["id"]]={"target":t,"queries":tr}

out={
 "artifact_id":"JANUS-KUSTO-EMODNET-CDI-ALL8-SOURCE-INVENTORY-RUN-2026-09-22-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-EMODNET-CDI-ALL8-SOURCE-INVENTORY-PREREG-2026-09-22-v1.0.json",
 "wms_capabilities":cap_record,
 "wms_parse_error":cap_error,
 "wms_layer_count":len(layers),
 "wms_candidate_source_layers":cand_layers,
 "wfs_capabilities":wfs_record,
 "wfs_featuretypes":wfs_featuretypes,
 "point_feature_info":responses,
 "depth_values_intentionally_not_requested":True,
 "claim_ceiling":"EXTERNAL_SURVEY_METADATA_DISCOVERY_ONLY"
}
p=OUT/"JANUS-KUSTO-EMODNET-CDI-ALL8-SOURCE-INVENTORY-RUN-2026-09-22-v1.0.json"
p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps({
 "wms_capabilities":cap_record,
 "wms_layer_count":len(layers),
 "candidate_layers":[{"name":x["name"],"title":x["title"],"queryable":x["queryable"]} for x in cand_layers],
 "wfs_capabilities":wfs_record,
 "wfs_featuretypes":wfs_featuretypes,
 "informative_point_queries":{
   k:[{"layer":q["layer"],"title":q["title"],"response":q["response"][:4000]} for q in v["queries"] if q["informative"]]
   for k,v in responses.items()
 }
},indent=2,ensure_ascii=False))
