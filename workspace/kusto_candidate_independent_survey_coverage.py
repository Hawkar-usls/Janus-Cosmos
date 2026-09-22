#!/usr/bin/env python3
import json, math
from pathlib import Path
from urllib.parse import urlencode
import requests

OUT=Path("workspace/kusto_crosssurvey_out")
OUT.mkdir(parents=True,exist_ok=True)
CANDS=[
 ("KN19207_CAND_001",-4.20546708989972,-14.871748111744372),
 ("KN19207_CAND_002",-3.9727527956056825,-12.272824298723462),
 ("KN19207_CAND_003",-4.015075679897318,-12.29915403590432),
 ("KN19207_CAND_004",-4.2399610081365955,-15.26855667523257),
 ("KN19207_CAND_005",-4.099231435240242,-12.332032012186781),
 ("KN19207_CAND_006",-4.309010339378846,-12.270140997573803),
 ("KN19207_CAND_007",-4.119167753526161,-14.356560209706116),
 ("KN19207_CAND_008",-4.3376789874153,-12.319260730325983),
]
FOOT="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
TRACK="https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/multibeam_dynamic/MapServer/0/query"
FIELDS="NCEI_ID,SURVEY_ID,PLATFORM,SOURCE,CHIEF_SCIENTIST,INSTRUMENT,START_TIME,END_TIME,SURVEY_YEAR,DOWNLOAD_URL,SURVEY_AND_VERSION"
TFIELDS="SURVEY_ID,PLATFORM,SURVEY_YEAR,SOURCE,NGDC_ID,CHIEF_SCIENTIST,INSTRUMENT,FILE_COUNT,TRACK_LENGTH,TOTAL_TIME,BATHY_BEAMS,AMP_BEAMS,SIDESCANS,DOWNLOAD_URL,START_TIME,END_TIME"
S=requests.Session()
S.headers.update({"User-Agent":"JANUS-KUSTO-crosssurvey-coverage/1.0"})

def get_json(url,params):
    r=S.get(url,params=params,timeout=90); r.raise_for_status()
    return r.json(),r.url

def hav(lat1,lon1,lat2,lon2):
    R=6371008.8
    p1,p2=map(math.radians,[lat1,lat2]); dp=math.radians(lat2-lat1); dl=math.radians(lon2-lon1)
    a=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R*math.asin(math.sqrt(a))

def segdist(lat,lon,a,b):
    sx=111320*math.cos(math.radians(lat)); sy=111320
    px,py=lon*sx,lat*sy; ax,ay=a[0]*sx,a[1]*sy; bx,by=b[0]*sx,b[1]*sy
    vx,vy=bx-ax,by-ay; wx,wy=px-ax,py-ay
    d=vx*vx+vy*vy
    t=0 if d==0 else max(0,min(1,(wx*vx+wy*vy)/d))
    return math.hypot(px-(ax+t*vx),py-(ay+t*vy))

def polyline_min_dist(lat,lon,geom):
    best=None
    for path in geom.get("paths",[]):
        for a,b in zip(path,path[1:]):
            d=segdist(lat,lon,a,b)
            if best is None or d<best: best=d
    return best

rows=[]
for cid,lat,lon in CANDS:
    ep={
      "geometry":f"{lon},{lat}","geometryType":"esriGeometryPoint","inSR":"4326","outSR":"4326",
      "spatialRel":"esriSpatialRelIntersects","outFields":FIELDS,"returnGeometry":"false","f":"json"
    }
    exact,exact_url=get_json(FOOT,ep)
    exact_feats=[x.get("attributes",{}) for x in exact.get("features",[])]
    independent=[x for x in exact_feats if (x.get("SURVEY_ID") or "").strip()!="KN192-07"]

    pad=0.15
    np={
      "geometry":f"{lon-pad},{lat-pad},{lon+pad},{lat+pad}","geometryType":"esriGeometryEnvelope",
      "inSR":"4326","outSR":"4326","spatialRel":"esriSpatialRelIntersects",
      "outFields":FIELDS,"returnGeometry":"false","f":"json"
    }
    nearp,near_url=get_json(FOOT,np)
    near_feats=[x.get("attributes",{}) for x in nearp.get("features",[])]
    seen=set(); uniq=[]
    for x in near_feats:
        key=(x.get("SURVEY_ID"),x.get("SURVEY_AND_VERSION"),x.get("NCEI_ID"))
        if key not in seen:
            seen.add(key); uniq.append(x)

    nt={
      "geometry":f"{lon-pad},{lat-pad},{lon+pad},{lat+pad}","geometryType":"esriGeometryEnvelope",
      "inSR":"4326","outSR":"4326","spatialRel":"esriSpatialRelIntersects",
      "outFields":TFIELDS,"returnGeometry":"true","f":"json"
    }
    tr,tr_url=get_json(TRACK,nt)
    trs=[]
    for f in tr.get("features",[]):
        a=f.get("attributes",{})
        a["min_trackline_distance_m"]=polyline_min_dist(lat,lon,f.get("geometry",{}))
        trs.append(a)
    trs.sort(key=lambda x:1e99 if x.get("min_trackline_distance_m") is None else x["min_trackline_distance_m"])

    rows.append({
      "candidate_id":cid,"lat":lat,"lon":lon,
      "exact_footprint_intersections":exact_feats,
      "exact_footprint_count":len(exact_feats),
      "independent_exact_footprint_intersections":independent,
      "independent_exact_footprint_count":len(independent),
      "coverage_gate_pass":len(independent)>0,
      "nearby_0p15deg_footprint_candidates":uniq,
      "nearby_trackline_candidates":trs[:30],
      "queries":{"exact":exact_url,"nearby_footprints":near_url,"nearby_tracklines":tr_url}
    })

receipt={
 "artifact_id":"JANUS-KUSTO-CANDIDATE-INDEPENDENT-SURVEY-COVERAGE-RUN-2026-09-22-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CANDIDATE-INDEPENDENT-SURVEY-COVERAGE-PREREG-2026-09-22-v1.0.json",
 "n_candidates":len(rows),
 "n_with_independent_exact_coverage":sum(r["coverage_gate_pass"] for r in rows),
 "results":rows,
 "firewall":"ONLY_EXACT_NCEI_FOOTPRINT_POINT_INTERSECTION_COUNTS_AS_COVERAGE__NEARBY_RESULTS_ARE_DISCOVERY_ONLY",
 "next_if_pass":"FOR_EACH_PASS_FETCH_INDEPENDENT_SURVEY_BATHYMETRY_AND_REPEAT_FIXED_LOCAL_METRICS",
 "next_if_zero":"PRESERVE_ZERO_AND_USE_CD169_TOBI_AS_SEPARATE_NON_NCEI_INDEPENDENT_LINEAGE_WHERE_NAV_SWATH_OVERLAP_CAN_BE_PROVED"
}
p=OUT/"JANUS-KUSTO-CANDIDATE-INDEPENDENT-SURVEY-COVERAGE-RUN-2026-09-22-v1.0.json"
p.write_text(json.dumps(receipt,indent=2),encoding="utf-8")
print(json.dumps({
 "n_candidates":receipt["n_candidates"],
 "n_with_independent_exact_coverage":receipt["n_with_independent_exact_coverage"],
 "summary":[
   {"id":r["candidate_id"],"exact":[x.get("SURVEY_ID") for x in r["exact_footprint_intersections"]],
    "independent":[x.get("SURVEY_ID") for x in r["independent_exact_footprint_intersections"]],
    "nearest_tracks":[(x.get("SURVEY_ID"),x.get("min_trackline_distance_m")) for x in r["nearby_trackline_candidates"][:5]]}
   for r in rows
 ]
},indent=2))
