#!/usr/bin/env python3
from __future__ import annotations
import json, math
from pathlib import Path
import requests
from pyproj import Geod

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I0N-C047-POSTRESULT-INDEPENDENT-TRUTH-CHARACTERIZATION-PREREG-2026-09-24-v1.0.json").read_text())
T=PRE["target"]; lon=float(T["lon"]); lat=float(T["lat"]); radius=float(T["radius_m"])
SEARCH_M=2000.0
LAYER="https://gsi.geodata.gov.ie/server/rest/services/Marine/IE_GSI_MI_Shipwrecks_IE_Waters_WGS84_LAT/FeatureServer/0"
q=LAYER+"/query"
# Geometry envelope is only a transport prefilter. Final classification uses geodesic distance.
dlat=SEARCH_M/110540.0
dlon=SEARCH_M/(111320.0*math.cos(math.radians(lat)))
env={"xmin":lon-dlon,"ymin":lat-dlat,"xmax":lon+dlon,"ymax":lat+dlat,"spatialReference":{"wkid":4326}}
params={
 "where":"1=1","geometry":json.dumps(env),"geometryType":"esriGeometryEnvelope",
 "inSR":"4326","spatialRel":"esriSpatialRelIntersects","outFields":"*",
 "returnGeometry":"true","outSR":"4326","f":"json"
}
S=requests.Session();S.headers.update({"User-Agent":"JANUS-KUSTO-INFOMAR-I0N-C047/1.0"})
r=S.get(q,params=params,timeout=180);r.raise_for_status();j=r.json()
if "error" in j: raise RuntimeError(j["error"])
geod=Geod(ellps="WGS84")
rows=[]
for f in j.get("features",[]) or []:
    g=f.get("geometry") or {}; x=g.get("x"); y=g.get("y")
    if x is None or y is None: continue
    _,_,dist=geod.inv(lon,lat,float(x),float(y))
    dist=float(abs(dist))
    if dist>SEARCH_M: continue
    if dist<=radius: cls="DIRECT"
    elif dist<=500: cls="NEAR"
    else: cls="REGIONAL"
    rows.append({"distance_m":dist,"distance_class":cls,"geometry":{"lon":float(x),"lat":float(y)},"attributes":f.get("attributes") or {}})
rows.sort(key=lambda x:x["distance_m"])
closest=rows[0] if rows else None
out={
 "artifact_id":"JANUS-KUSTO-INFOMAR-I0N-C047-OFFICIAL-GSI-SHIPWRECK-SPATIAL-QUERY-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],"formal_parent":PRE["formal_parent"],
 "target":T,"official_source":{
   "layer":LAYER,"query_url":r.url,"returned_feature_count_before_radius_filter":len(j.get("features",[]) or []),
   "source_role":"OFFICIAL_GSI_MI_INFOMAR_MAPPED_SHIPWRECK_LAYER"
 },
 "search_radius_m":SEARCH_M,"matched_within_2km_count":len(rows),"matches":rows,
 "closest_match":closest,
 "result_class":(closest["distance_class"] if closest else "NONE"),
 "formal_parent_changed":False,"promotion":False,
 "claim_ceiling":"POSTRESULT_CASE_LEVEL_OFFICIAL_GEOREFERENCED_WRECK_PROXIMITY_ONLY"
}
p=OUT/"JANUS-KUSTO-INFOMAR-I0N-C047-OFFICIAL-GSI-SHIPWRECK-SPATIAL-QUERY-2026-09-24-v1.0.json"
p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps({
 "artifact_id":out["artifact_id"],"matched_within_2km_count":len(rows),
 "result_class":out["result_class"],"closest_distance_m":None if closest is None else closest["distance_m"],
 "closest_attributes":None if closest is None else closest["attributes"]
},indent=2,ensure_ascii=False))
