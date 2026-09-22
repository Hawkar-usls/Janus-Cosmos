#!/usr/bin/env python3
import json, math
from pathlib import Path
from urllib.parse import urlencode
import requests

OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
ENDPOINT="https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/trackline_combined_dynamic/MapServer/1/query"
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
FIELDS="SURVEY_TYPE,SURVEY_ID,INST_SRC,COUNTRY,PLATFORM,PROJECT,CHIEF,START_YR,END_YR,DOWNLOAD_URL,SURVEY_YEAR,LAT_TOP,LAT_BOTTOM,LON_LEFT,LON_RIGHT"
S=requests.Session();S.headers.update({"User-Agent":"JANUS-KUSTO-global-bathy-trackline-discovery/1.0"})

def local_xy(t,lon,lat):
    return ((lon-t["lon"])*111320*math.cos(math.radians(t["lat"])),(lat-t["lat"])*111320)

def segdist(a,b):
    ax,ay=a;bx,by=b
    vx,vy=bx-ax,by-ay
    den=vx*vx+vy*vy
    if den==0:return math.hypot(ax,ay)
    tt=max(0,min(1,-(ax*vx+ay*vy)/den))
    return math.hypot(ax+tt*vx,ay+tt*vy)

def geom_dist(t,geom):
    best=float("inf")
    for path in geom.get("paths",[]):
        if len(path)==1:
            best=min(best,math.hypot(*local_xy(t,path[0][0],path[0][1])))
        for a,b in zip(path,path[1:]):
            best=min(best,segdist(local_xy(t,a[0],a[1]),local_xy(t,b[0],b[1])))
    return None if not math.isfinite(best) else best

def query(t):
    # ~17km box; final filter exact <=15km.
    pad=0.16
    p={
      "geometry":f'{t["lon"]-pad},{t["lat"]-pad},{t["lon"]+pad},{t["lat"]+pad}',
      "geometryType":"esriGeometryEnvelope","inSR":"4326","outSR":"4326",
      "spatialRel":"esriSpatialRelIntersects","outFields":FIELDS,
      "returnGeometry":"true","resultRecordCount":"2000","f":"json"
    }
    r=S.get(ENDPOINT,params=p,timeout=120);r.raise_for_status();d=r.json()
    if "error" in d:raise RuntimeError(d["error"])
    by={}
    for f in d.get("features",[]):
        a=f.get("attributes",{});sid=a.get("SURVEY_ID")
        if not sid:continue
        dist=geom_dist(t,f.get("geometry",{}))
        if dist is None or dist>15000:continue
        rec=by.get(sid)
        item={"survey_id":sid,"min_trackline_distance_m":dist,**a}
        if rec is None or dist<rec["min_trackline_distance_m"]:
            by[sid]=item
    rows=sorted(by.values(),key=lambda x:(x["min_trackline_distance_m"],str(x["survey_id"])))
    return {"query_url":r.url,"survey_count_within_15km":len(rows),"surveys":rows}

res={t["id"]:{"target":t,**query(t)} for t in CANDS}
known={"KN192-07","KNOX15RR","CD169"}
out={
 "artifact_id":"JANUS-KUSTO-NCEI-GLOBAL-BATHY-TRACKLINE-ALL8-DISCOVERY-RUN-2026-09-22-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-NCEI-GLOBAL-BATHY-TRACKLINE-ALL8-DISCOVERY-PREREG-2026-09-22-v1.0.json",
 "source":ENDPOINT,
 "candidates":res,
 "candidate_new_survey_ids":{
   k:[x["survey_id"] for x in v["surveys"] if x["survey_id"] not in known]
   for k,v in res.items()
 },
 "firewall":"TRACKLINE_PROXIMITY_IS_DISCOVERY_CLUE_ONLY__NO_SWATH_COVERAGE_OR_DEPTH_AUTHORITY",
 "claim_ceiling":"GLOBAL_BATHYMETRY_SURVEY_DISCOVERY_CLUES_ONLY"
}
p=OUT/"JANUS-KUSTO-NCEI-GLOBAL-BATHY-TRACKLINE-ALL8-DISCOVERY-RUN-2026-09-22-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "summary":{k:{
   "survey_count":v["survey_count_within_15km"],
   "surveys":[{"survey_id":x["survey_id"],"distance_m":x["min_trackline_distance_m"],"platform":x.get("PLATFORM"),"year":x.get("SURVEY_YEAR"),"project":x.get("PROJECT"),"download":x.get("DOWNLOAD_URL")} for x in v["surveys"]]
 } for k,v in res.items()},
 "new_ids":out["candidate_new_survey_ids"]
},indent=2))
