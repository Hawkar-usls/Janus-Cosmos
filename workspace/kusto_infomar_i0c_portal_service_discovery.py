#!/usr/bin/env python3
from __future__ import annotations
import json,re
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
APP_ID="3f2815ec89e745d2b65630429d06385c"
BASE="https://www.arcgis.com/sharing/rest/content/items"
S=requests.Session();S.headers.update({"User-Agent":"JANUS-KUSTO-INFOMAR-I0C/1.0"})

def get_json(url,params=None):
    r=S.get(url,params=params or {"f":"json"},timeout=120,allow_redirects=True);r.raise_for_status()
    return r.json()

meta=get_json(f"{BASE}/{APP_ID}")
data=get_json(f"{BASE}/{APP_ID}/data")

urls=set(); ids=set()
def walk(x):
    if isinstance(x,dict):
        for k,v in x.items():
            if isinstance(v,str):
                if v.startswith("http://") or v.startswith("https://"): urls.add(v)
                for m in re.findall(r"\b[0-9a-fA-F]{32}\b",v): ids.add(m.lower())
            walk(v)
    elif isinstance(x,list):
        for v in x: walk(v)
    elif isinstance(x,str):
        if x.startswith("http://") or x.startswith("https://"): urls.add(x)
        for m in re.findall(r"\b[0-9a-fA-F]{32}\b",x): ids.add(m.lower())
walk(data)
ids.discard(APP_ID.lower())

items=[]
for iid in sorted(ids)[:250]:
    try:
        j=get_json(f"{BASE}/{iid}")
    except Exception as e:
        items.append({"id":iid,"error":type(e).__name__+": "+str(e)});continue
    items.append({
      "id":iid,"title":j.get("title"),"type":j.get("type"),"url":j.get("url"),
      "owner":j.get("owner"),"access":j.get("access"),"description":j.get("description"),
      "tags":j.get("tags")
    })

relevant=[]
for x in items:
    text=" ".join(str(x.get(k) or "") for k in ("title","type","url","description","tags")).lower()
    if any(t in text for t in ["infomar","survey","bathym","backscatter","download","marine"]):
        relevant.append(x)

service_urls=sorted(u for u in urls if any(t in u.lower() for t in ["featureserver","imageserver","mapserver","arcgis.com","marine.ie","infomar"]))
out={
 "artifact_id":"JANUS-KUSTO-INFOMAR-I0C-DOWNLOAD-PORTAL-SERVICE-DISCOVERY-2026-09-24-v1.0",
 "experience_app":{"id":APP_ID,"title":meta.get("title"),"owner":meta.get("owner"),"type":meta.get("type")},
 "referenced_item_count":len(ids),
 "referenced_items":items,
 "relevant_items":relevant,
 "service_urls":service_urls,
 "pixel_values_read":False,
 "downloaded_survey_files":False,
 "claim_ceiling":"OFFICIAL_DOWNLOAD_PORTAL_SERVICE_DISCOVERY_ONLY"
}
p=OUT/"JANUS-KUSTO-INFOMAR-I0C-DOWNLOAD-PORTAL-SERVICE-DISCOVERY-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "app":out["experience_app"],
 "referenced_item_count":len(ids),
 "relevant_items":[{"id":x.get("id"),"title":x.get("title"),"type":x.get("type"),"url":x.get("url")} for x in relevant[:100]],
 "service_urls":service_urls[:100],
 "pixel_values_read":False
},indent=2,ensure_ascii=False))
