#!/usr/bin/env python3
import json,requests
from pathlib import Path
from shapely.geometry import shape
from shapely.ops import unary_union

OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
ARC="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
EXCLUDED={"PF0501","EX2105","CD169","KN192-07","KNOX15RR"}
UA={"User-Agent":"JANUS-KUSTO-PF0501-geometry-rank/1.0"}

def q(p):
    p=dict(p);p["f"]="geojson"
    r=requests.get(ARC,params=p,headers=UA,timeout=90);r.raise_for_status()
    d=r.json()
    if "error" in d:raise RuntimeError(d["error"])
    return d

pf=q({"where":"SURVEY_ID='PF0501'","outFields":"SURVEY_ID","returnGeometry":"true","outSR":"4326"})
pg=unary_union([shape(f["geometry"]) for f in pf["features"] if f.get("geometry")])
w,s,e,n=pg.bounds
allf=[];off=0
while True:
    d=q({
      "where":"1=1","geometry":f"{w},{s},{e},{n}",
      "geometryType":"esriGeometryEnvelope","inSR":"4326",
      "spatialRel":"esriSpatialRelIntersects",
      "outFields":"SURVEY_ID,PLATFORM,INSTRUMENT,SURVEY_YEAR,SOURCE,DOWNLOAD_URL",
      "returnGeometry":"true","outSR":"4326",
      "resultOffset":off,"resultRecordCount":1000
    })
    fs=d.get("features",[]);allf.extend(fs)
    if len(fs)<1000:break
    off+=len(fs)

by={}
for f in allf:
    p=f.get("properties",{});sid=p.get("SURVEY_ID")
    if sid and f.get("geometry"):
        by.setdefault(sid,{"p":p,"g":[]})["g"].append(shape(f["geometry"]))
rows=[]
for sid,v in by.items():
    if sid in EXCLUDED:continue
    try:y=int(v["p"].get("SURVEY_YEAR"))
    except:continue
    if y<2006:continue
    g=unary_union(v["g"]);inter=pg.intersection(g)
    if inter.is_empty or inter.area<=0:continue
    rows.append({
      "survey_id":sid,"survey_year":y,"platform":v["p"].get("PLATFORM"),
      "instrument":v["p"].get("INSTRUMENT"),"source":v["p"].get("SOURCE"),
      "download_url":v["p"].get("DOWNLOAD_URL"),
      "intersection_area_deg2":float(inter.area),"intersection_bounds":list(inter.bounds)
    })
rows.sort(key=lambda x:(-x["intersection_area_deg2"],x["survey_id"]))
out={"artifact_id":"JANUS-KUSTO-UNSEEN-PF0501-REFERENCE-GEOMETRY-RANK-2026-09-23-v1.0",
     "depth_values_read":False,"eligible_count":len(rows),"ranked":rows,
     "claim_ceiling":"GEOMETRY_ONLY_REFERENCE_RANKING"}
(OUT/"JANUS-KUSTO-UNSEEN-PF0501-REFERENCE-GEOMETRY-RANK-2026-09-23-v1.0.json").write_text(json.dumps(out,indent=2))
print(json.dumps({"eligible_count":len(rows),"top20":rows[:20]},indent=2))
