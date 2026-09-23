#!/usr/bin/env python3
from __future__ import annotations
import json, hashlib
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"
OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-KUSTO-USGS-natural-controls-source-gate/1.0"}
DATASETS=[
 {"id":"USGS_MIAMI_KEY_BISCAYNE_POCKMARKS","sciencebase_id":"5a96f5cee4b06990606c4ff0","doi":"10.5066/F72J6B4Z"},
 {"id":"USGS_LAKE_CRESCENT_2016","sciencebase_id":"586d3165e4b0f5ce109faa51","doi":"10.5066/F7B56GW5"}
]

def getj(url):
 r=requests.get(url,headers=UA,timeout=90)
 r.raise_for_status()
 return r.json()

def file_role(name):
 n=name.lower()
 roles=[]
 if any(x in n for x in ["bathy","bathym","dem","depth"]):roles.append("BATHYMETRY")
 if any(x in n for x in ["backscatter","back_scatter","reflect"]):roles.append("BACKSCATTER")
 if any(x in n for x in [".xml",".txt","metadata","readme","fgdc"]):roles.append("METADATA_OR_TEXT")
 if any(x in n for x in [".tif",".tiff",".asc",".grd",".xyz",".nc"]):roles.append("NUMERIC_OR_RASTER_PRODUCT")
 return roles

outsets=[]
for ds in DATASETS:
 url=f"https://www.sciencebase.gov/catalog/item/{ds['sciencebase_id']}?format=json"
 try:
  j=getj(url)
  files=[]
  for f in j.get("files") or []:
   name=f.get("name") or ""
   files.append({
    "name":name,
    "url":f.get("url"),
    "size":f.get("size"),
    "contentType":f.get("contentType"),
    "checksum":f.get("checksum"),
    "roles":file_role(name)
   })
  outsets.append({
   **ds,
   "sciencebase_api":url,
   "title":j.get("title"),
   "files_total":len(files),
   "bathymetry_candidates":[f for f in files if "BATHYMETRY" in f["roles"]],
   "backscatter_candidates":[f for f in files if "BACKSCATTER" in f["roles"]],
   "all_files":files,
   "status":"BOUND" if files else "BOUND_NO_FILES_LISTED"
  })
 except Exception as e:
  outsets.append({**ds,"sciencebase_api":url,"status":"ERROR","error":type(e).__name__+": "+str(e)})

out={
 "artifact_id":"JANUS-KUSTO-USGS-NATURAL-CONTROLS-SOURCE-GATE-RUN-2026-09-23-v1.0",
 "stage":"G0_SOURCE_AND_PROVENANCE_INVENTORY",
 "depth_values_read":False,
 "detector_run":False,
 "datasets":outsets,
 "hard_rules":[
  "ScienceBase inventory only; no raster values read.",
  "File-name role classification is routing metadata, not scientific classification.",
  "Processed products may calibrate detection but cannot independently validate their own processing lineage."
 ],
 "next_gate":"FREEZE_NUMERIC_PRODUCTS_AND_NATIVE_RESOLUTION_BEFORE_BLIND_DETECTION",
 "claim_ceiling":"SOURCE_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-USGS-NATURAL-CONTROLS-SOURCE-GATE-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "datasets":[{
  "id":d["id"],"status":d["status"],"files_total":d.get("files_total"),
  "bathymetry_candidates":len(d.get("bathymetry_candidates") or []),
  "backscatter_candidates":len(d.get("backscatter_candidates") or [])
 } for d in outsets],
 "depth_values_read":False
},indent=2))
