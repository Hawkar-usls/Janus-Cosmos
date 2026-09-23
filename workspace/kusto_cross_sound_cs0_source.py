#!/usr/bin/env python3
from __future__ import annotations
import json, re, time
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-CROSS-SOUND-ALASKA-UNSEEN-DUALCHANNEL-PREREG-2026-09-23-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-CrossSound-CS0-source/1.0"}
BASE="https://www.sciencebase.gov/catalog"
PARENT=PRE["dataset"]["release_parent_sciencebase_item"]
KNOWN_BATHY=PRE["dataset"]["known_bathymetry_sciencebase_item"]

def get_json(url,params=None):
    last=None
    for a in range(6):
        try:
            r=requests.get(url,params=params,headers=UA,timeout=90)
            if r.status_code==429:
                time.sleep(min(20,2**a));continue
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last=e
            if a==5:raise
            time.sleep(min(20,2**a))
    raise last

def item(i):
    return get_json(f"{BASE}/item/{i}",{"format":"json"})

def children(i):
    probes=[
      (f"{BASE}/items",{"parentId":i,"format":"json","max":100}),
      (f"{BASE}/items",{"q":"","parentId":i,"format":"json","max":100})
    ]
    for u,p in probes:
        try:
            d=get_json(u,p)
            xs=d.get("items") if isinstance(d,dict) else None
            if isinstance(xs,list):return xs
        except Exception:
            pass
    return []

def file_rows(it):
    rows=[]
    iid=str(it.get("id") or "")
    title=str(it.get("title") or "")
    for f in (it.get("files") or []):
        if not isinstance(f,dict):continue
        name=str(f.get("name") or "")
        url=f.get("url") or f.get("downloadUri") or f.get("downloadURL")
        size=f.get("size")
        csum=f.get("checksum")
        if isinstance(csum,dict):
            checksum=csum
        else:
            checksum={"value":csum} if csum else None
        rows.append({
          "item_id":iid,"item_title":title,"name":name,"url":url,
          "size":size,"contentType":f.get("contentType"),"checksum":checksum
        })
    return rows

root=item(PARENT)
level1=children(PARENT)
# Some ScienceBase releases put actual products one level further down.
seen={PARENT}; items=[root]
for x in level1:
    iid=str(x.get("id") or "")
    if iid and iid not in seen:
        try:items.append(item(iid))
        except Exception:items.append(x)
        seen.add(iid)
for x in list(items[1:]):
    iid=str(x.get("id") or "")
    if not iid:continue
    for y in children(iid):
        yid=str(y.get("id") or "")
        if yid and yid not in seen:
            try:items.append(item(yid))
            except Exception:items.append(y)
            seen.add(yid)

# Guarantee known bathymetry child is inspected even if child API is incomplete.
if KNOWN_BATHY not in seen:
    items.append(item(KNOWN_BATHY));seen.add(KNOWN_BATHY)

files=[]
for it in items:files.extend(file_rows(it))

def norm(s):return re.sub(r"[^a-z0-9]+"," ",str(s).lower())
def is_data_file(r):
    n=str(r["name"]).lower()
    return n.endswith((".zip",".tif",".tiff",".grd",".asc"))

def score(r,kind):
    t=norm(r["item_title"]+" "+r["name"])
    n=str(r["name"]).lower()
    s=0
    if kind=="bathy":
        if "bathy" in t or "bathymetry" in t:s+=10
        if "backscatter" in t:s-=20
    else:
        if "backscatter" in t:s+=12
        if "bathy" in t or "bathymetry" in t:s-=5
    if n.endswith(".zip"):s+=5
    elif n.endswith((".tif",".tiff")):s+=4
    if "metadata" in t or n.endswith((".xml",".txt",".jpg",".jpeg",".png")):s-=20
    try:
        if r.get("size") and int(r["size"])>100000:s+=2
    except:pass
    return s

eligible=[r for r in files if is_data_file(r)]
bathy=sorted(eligible,key=lambda r:(score(r,"bathy"),int(r.get("size") or 0)),reverse=True)
back=sorted(eligible,key=lambda r:(score(r,"backscatter"),int(r.get("size") or 0)),reverse=True)
bsel=bathy[0] if bathy and score(bathy[0],"bathy")>=10 else None
ksel=back[0] if back and score(back[0],"backscatter")>=10 else None

# Distinct product requirement.
if bsel and ksel and bsel["url"]==ksel["url"]:ksel=None

status="PASS_BATHY_BACKSCATTER_BOUND__NO_RASTER_VALUES_READ" if bsel and ksel else "BLOCKED_PRODUCT_BINDING_INCOMPLETE"
out={
 "artifact_id":"JANUS-KUSTO-CROSS-SOUND-CS0-SOURCE-BINDING-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],"status":status,
 "sciencebase_parent_item":PARENT,
 "items_inspected":len(items),"item_ids":[str(x.get("id") or "") for x in items],
 "file_metadata_rows":files,
 "selected":{"bathymetry":bsel,"backscatter":ksel},
 "selection_scores":{
   "bathymetry":score(bsel,"bathy") if bsel else None,
   "backscatter":score(ksel,"backscatter") if ksel else None
 },
 "raster_values_read":False,"numeric_payload_downloaded":False,
 "hard_rules_preserved":["METADATA_ONLY","NO_RASTER_VALUE_READ","NO_FEATURE_LABEL_READ","AMBIGUOUS_BINDING_BLOCKS_EXECUTION"],
 "next_gate":"CS1_BATHYMETRY_ONLY_BLIND_CANDIDATE_FREEZE" if status.startswith("PASS") else "SOURCE_BINDING_REPAIR_ONLY",
 "claim_ceiling":"SOURCE_AND_FILE_IDENTITY_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-CROSS-SOUND-CS0-SOURCE-BINDING-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "artifact_id":out["artifact_id"],"status":status,"items_inspected":len(items),
 "files_seen":len(files),
 "file_inventory":[{"item_id":r["item_id"],"item_title":r["item_title"],"name":r["name"],"size":r["size"],"contentType":r["contentType"],"url":r["url"]} for r in files],
 "bathymetry":bsel,"backscatter":ksel,
 "raster_values_read":False
},indent=2))
