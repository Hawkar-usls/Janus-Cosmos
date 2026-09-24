#!/usr/bin/env python3
from __future__ import annotations
import json,re
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
WEBMAP_ID="2e97ca9e7d5a4ab38605b79b8fff677e"
BASE="https://www.arcgis.com/sharing/rest/content/items"
S=requests.Session();S.headers.update({"User-Agent":"JANUS-KUSTO-INFOMAR-I0D/1.0"})

def jget(url,params=None):
    r=S.get(url,params=params or {"f":"json"},timeout=120,allow_redirects=True);r.raise_for_status()
    try:return r.json()
    except Exception:return {"_raw_prefix":r.text[:2000]}

meta=jget(f"{BASE}/{WEBMAP_ID}")
data=jget(f"{BASE}/{WEBMAP_ID}/data")

layers=[]
def collect_layers(x,path="root"):
    if isinstance(x,dict):
        # ArcGIS webmap layer-like object
        if any(k in x for k in ("url","itemId","title")) and ("layerType" in x or "url" in x or "itemId" in x):
            layers.append({
              "path":path,"title":x.get("title"),"id":x.get("id"),"itemId":x.get("itemId"),
              "url":x.get("url"),"layerType":x.get("layerType"),"visibility":x.get("visibility")
            })
        for k,v in x.items(): collect_layers(v,path+"."+str(k))
    elif isinstance(x,list):
        for i,v in enumerate(x): collect_layers(v,path+f"[{i}]")
collect_layers(data)

# dedupe layer dicts
uniq=[];seen=set()
for x in layers:
    key=(x.get("title"),x.get("itemId"),x.get("url"),x.get("layerType"))
    if key in seen:continue
    seen.add(key);uniq.append(x)
layers=uniq

item_ids=sorted({x["itemId"] for x in layers if x.get("itemId")})
item_meta={}
for iid in item_ids:
    item_meta[iid]=jget(f"{BASE}/{iid}")

service_urls=sorted({x["url"] for x in layers if x.get("url") and any(t in x["url"].lower() for t in ("featureserver","mapserver","imageserver"))})
service_meta={}
for u in service_urls:
    base=u.split("?")[0]
    service_meta[u]=jget(base,{"f":"json"})

out={
 "artifact_id":"JANUS-KUSTO-INFOMAR-I0D-SURVEY-LEG-WEBMAP-LAYER-GRAPH-2026-09-24-v1.0",
 "webmap":{"id":WEBMAP_ID,"title":meta.get("title"),"owner":meta.get("owner")},
 "layers":layers,
 "referenced_item_metadata":item_meta,
 "service_metadata":service_meta,
 "pixel_values_read":False,
 "feature_rows_queried":False,
 "claim_ceiling":"SURVEY_LEG_WEBMAP_LAYER_GRAPH_ONLY"
}
p=OUT/"JANUS-KUSTO-INFOMAR-I0D-SURVEY-LEG-WEBMAP-LAYER-GRAPH-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "webmap":out["webmap"],
 "layers":layers,
 "item_meta":[{"id":k,"title":v.get("title"),"type":v.get("type"),"url":v.get("url")} for k,v in item_meta.items()],
 "services":[{"url":k,"name":v.get("name"),"type":v.get("type"),"layers":v.get("layers"),"tables":v.get("tables")} for k,v in service_meta.items()],
 "pixel_values_read":False,
 "feature_rows_queried":False
},indent=2,ensure_ascii=False))
