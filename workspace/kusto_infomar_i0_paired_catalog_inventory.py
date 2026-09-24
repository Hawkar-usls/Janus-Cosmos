#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json, math, re
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I0-PAIRED-CATALOG-INVENTORY-PREREG-2026-09-24-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-INFOMAR-I0/1.0"}
S=requests.Session();S.headers.update(UA)

FIELDS="OBJECTID,Name,MinPS,MaxPS,Category,Tag,GroupName,ProductName,CenterX,CenterY"

def fetch_catalog(base):
    # First obtain IDs only, then fetch in deterministic chunks. Geometry metadata only.
    q=base+"/query"
    r=S.get(q,params={"where":"1=1","returnIdsOnly":"true","f":"json"},timeout=120);r.raise_for_status()
    j=r.json()
    ids=sorted(int(x) for x in (j.get("objectIds") or []))
    rows=[]
    for i in range(0,len(ids),200):
        chunk=ids[i:i+200]
        rr=S.get(q,params={
          "objectIds":",".join(map(str,chunk)),
          "outFields":FIELDS,
          "returnGeometry":"true",
          "outSR":"3857",
          "f":"json"
        },timeout=180)
        rr.raise_for_status()
        jj=rr.json()
        if "error" in jj: raise RuntimeError(f"ArcGIS query error {jj['error']}")
        for f in jj.get("features",[]) or []:
            a=f.get("attributes") or {}
            g=f.get("geometry")
            if g and "rings" in g:
                xs=[];ys=[]
                for ring in g["rings"]:
                    for x,y in ring:
                        xs.append(float(x));ys.append(float(y))
                bbox=[min(xs),min(ys),max(xs),max(ys)] if xs else None
            else:
                bbox=None
            rows.append({
              "OBJECTID":a.get("OBJECTID"),
              "Name":a.get("Name"),
              "MinPS":a.get("MinPS"),
              "MaxPS":a.get("MaxPS"),
              "Category":a.get("Category"),
              "Tag":a.get("Tag"),
              "GroupName":a.get("GroupName"),
              "ProductName":a.get("ProductName"),
              "CenterX":a.get("CenterX"),
              "CenterY":a.get("CenterY"),
              "bbox_3857":bbox,
              "geometry_present":bool(g)
            })
    return ids,rows

def norm(s):
    return re.sub(r"\s+"," ",(s or "").strip()).casefold()

bath_ids,bath=fetch_catalog(PRE["sources"]["bathymetry_catalog"])
back_ids,back=fetch_catalog(PRE["sources"]["backscatter_catalog"])

bmap={}
for x in bath:
    n=norm(x["Name"])
    if n:bmap.setdefault(n,[]).append(x)
kmap={}
for x in back:
    n=norm(x["Name"])
    if n:kmap.setdefault(n,[]).append(x)

common=[]
for n in sorted(set(bmap)&set(kmap)):
    for b in bmap[n]:
        for k in kmap[n]:
            common.append({
              "normalized_name":n,
              "display_name_bathymetry":b["Name"],
              "display_name_backscatter":k["Name"],
              "bathymetry":b,
              "backscatter":k,
              "both_geometry_present":bool(b["geometry_present"] and k["geometry_present"]),
              "both_primary_category":bool(b.get("Category")==1 and k.get("Category")==1)
            })

# Also report common GroupName pairs as fallback metadata only; no selection based on these yet.
bg={}
kg={}
for x in bath:
    n=norm(x["GroupName"])
    if n:bg.setdefault(n,[]).append(x)
for x in back:
    n=norm(x["GroupName"])
    if n:kg.setdefault(n,[]).append(x)
group_common=[{"normalized_group":n,"bathymetry_count":len(bg[n]),"backscatter_count":len(kg[n])}
              for n in sorted(set(bg)&set(kg))]

primary_exact=[x for x in common if x["both_geometry_present"] and x["both_primary_category"]]
out={
 "artifact_id":"JANUS-KUSTO-INFOMAR-I0-PAIRED-CATALOG-INVENTORY-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],
 "bathymetry_catalog":{"url":PRE["sources"]["bathymetry_catalog"],"object_id_count":len(bath_ids),"rows":bath},
 "backscatter_catalog":{"url":PRE["sources"]["backscatter_catalog"],"object_id_count":len(back_ids),"rows":back},
 "exact_common_name_pairs":common,
 "exact_common_primary_geometry_pairs":primary_exact,
 "common_group_summary":group_common,
 "counts":{
   "bathymetry_items":len(bath),
   "backscatter_items":len(back),
   "exact_common_name_pairs":len(common),
   "exact_common_primary_geometry_pairs":len(primary_exact)
 },
 "pixel_values_read":False,
 "export_image_called":False,
 "get_samples_called":False,
 "next_gate":"I0B_FREEZE_ONE_EUROPEAN_SURVEY_BY_PREREGISTERED_METADATA_RULE" if primary_exact else "SOURCE_SCHEMA_DIAGNOSTIC_ONLY",
 "claim_ceiling":"INFOMAR_PAIRED_CATALOG_METADATA_ONLY"
}
p=OUT/"JANUS-KUSTO-INFOMAR-I0-PAIRED-CATALOG-INVENTORY-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "counts":out["counts"],
 "primary_candidates":[
   {"name":x["display_name_bathymetry"],
    "bathy_oid":x["bathymetry"]["OBJECTID"],"back_oid":x["backscatter"]["OBJECTID"],
    "bathy_minps":x["bathymetry"]["MinPS"],"back_minps":x["backscatter"]["MinPS"],
    "group_bathy":x["bathymetry"]["GroupName"],"group_back":x["backscatter"]["GroupName"],
    "bbox_bathy":x["bathymetry"]["bbox_3857"],"bbox_back":x["backscatter"]["bbox_3857"]}
   for x in primary_exact[:50]
 ],
 "pixel_values_read":False
},indent=2,ensure_ascii=False))
