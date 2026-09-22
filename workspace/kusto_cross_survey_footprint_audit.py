#!/usr/bin/env python3
import json, math
from pathlib import Path
from urllib.parse import urlencode
import requests

CANDS = [
    {"id":"KN19207_CAND_001","lat":-4.20546708989972,"lon":-14.871748111744372},
    {"id":"KN19207_CAND_002","lat":-3.9727527956056825,"lon":-12.272824298723462},
    {"id":"KN19207_CAND_003","lat":-4.015075679897318,"lon":-12.29915403590432},
    {"id":"KN19207_CAND_004","lat":-4.2399610081365955,"lon":-15.26855667523257},
    {"id":"KN19207_CAND_005","lat":-4.099231435240242,"lon":-12.332032012186781},
    {"id":"KN19207_CAND_006","lat":-4.309010339378846,"lon":-12.270140997573803},
    {"id":"KN19207_CAND_007","lat":-4.119167753526161,"lon":-14.356560209706116},
    {"id":"KN19207_CAND_008","lat":-4.3376789874153,"lon":-12.319260730325983},
]
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session()
S.headers.update({"User-Agent":"JANUS-KUSTO-cross-survey-footprint-audit/1.0"})
FP="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
TRACK="https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/multibeam_dynamic/MapServer/0/query"

FIELDS="NCEI_ID,SURVEY_ID,PLATFORM,SOURCE,CHIEF_SCIENTIST,INSTRUMENT,START_TIME,END_TIME,SURVEY_YEAR,DOWNLOAD_URL,SURVEY_AND_VERSION"

def get(endpoint, params):
    r=S.get(endpoint,params=params,timeout=120)
    r.raise_for_status()
    d=r.json()
    if "error" in d:
        raise RuntimeError(d["error"])
    return d

def norm(a):
    if not a:return {}
    # normalize JSON-serializable attribute values
    return {k:v for k,v in a.items()}

def point_fp(c):
    p={
      "geometry":f'{c["lon"]},{c["lat"]}',
      "geometryType":"esriGeometryPoint","inSR":"4326","outSR":"4326",
      "spatialRel":"esriSpatialRelIntersects","outFields":FIELDS,
      "returnGeometry":"false","returnDistinctValues":"false","f":"json"
    }
    d=get(FP,p)
    return [norm(f.get("attributes",{})) for f in d.get("features",[])], FP+"?"+urlencode(p)

def nearby_fp(c,pad=0.05):
    p={
      "geometry":f'{c["lon"]-pad},{c["lat"]-pad},{c["lon"]+pad},{c["lat"]+pad}',
      "geometryType":"esriGeometryEnvelope","inSR":"4326","outSR":"4326",
      "spatialRel":"esriSpatialRelIntersects","outFields":FIELDS,
      "returnGeometry":"false","f":"json"
    }
    d=get(FP,p)
    return [norm(f.get("attributes",{})) for f in d.get("features",[])], FP+"?"+urlencode(p)

def trackline_near(c,pad=0.05):
    p={
      "geometry":f'{c["lon"]-pad},{c["lat"]-pad},{c["lon"]+pad},{c["lat"]+pad}',
      "geometryType":"esriGeometryEnvelope","inSR":"4326","outSR":"4326",
      "spatialRel":"esriSpatialRelIntersects",
      "outFields":"SURVEY_ID,PLATFORM,SURVEY_YEAR,SOURCE,NGDC_ID,CHIEF_SCIENTIST,INSTRUMENT,DOWNLOAD_URL,START_TIME,END_TIME",
      "returnGeometry":"false","f":"json"
    }
    d=get(TRACK,p)
    return [norm(f.get("attributes",{})) for f in d.get("features",[])]

rows=[]
for c in CANDS:
    exact,url=point_fp(c)
    near,nurl=nearby_fp(c)
    tracks=trackline_near(c)
    exact_ids=sorted(set(x.get("SURVEY_ID") for x in exact if x.get("SURVEY_ID")))
    independent=[x for x in exact if x.get("SURVEY_ID") and x.get("SURVEY_ID")!="KN192-07"]
    indep_ids=sorted(set(x.get("SURVEY_ID") for x in independent if x.get("SURVEY_ID")))
    near_ids=sorted(set(x.get("SURVEY_ID") for x in near if x.get("SURVEY_ID")))
    track_ids=sorted(set(x.get("SURVEY_ID") for x in tracks if x.get("SURVEY_ID")))
    row={
      **c,
      "exact_intersections":exact,
      "exact_survey_ids":exact_ids,
      "exact_intersection_count":len(exact),
      "independent_exact_intersections":independent,
      "independent_exact_survey_ids":indep_ids,
      "independent_exact_survey_count":len(indep_ids),
      "stage_B_eligible":bool(indep_ids),
      "nearby_0p05deg_survey_ids":near_ids,
      "nearby_trackline_survey_ids":track_ids,
      "exact_query_url":url,
      "nearby_query_url":nurl,
      "firewall":"ONLY_EXACT_FOOTPRINT_INTERSECTIONS_COUNT_FOR_STAGE_B"
    }
    rows.append(row)
    print(json.dumps({
      "id":c["id"],
      "exact":exact_ids,
      "independent_exact":indep_ids,
      "nearby_footprints":near_ids,
      "nearby_tracklines":track_ids
    },ensure_ascii=False))

summary={
  "artifact_id":"JANUS-KUSTO-CROSS-SURVEY-EXACT-FOOTPRINT-AUDIT-2026-09-22-v1.0",
  "prereg":"data/cousteau/JANUS-KUSTO-CROSS-SURVEY-SEAFLOOR-VALIDATION-PREREG-2026-09-22-v1.0.json",
  "source":"NOAA_NCEI_MULTIBEAM_FOOTPRINT_POLYGON_LAYER",
  "candidates":rows,
  "stage_B_eligible_candidates":[r["id"] for r in rows if r["stage_B_eligible"]],
  "candidate_count":len(rows),
  "stage_B_eligible_count":sum(r["stage_B_eligible"] for r in rows),
  "claim_ceiling":"EXACT_SURVEY_COVERAGE_INVENTORY_ONLY__NO_MORPHOLOGY_REPLICATION_YET"
}
p=OUT/"JANUS-KUSTO-CROSS-SURVEY-EXACT-FOOTPRINT-AUDIT-2026-09-22-v1.0.json"
p.write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding="utf-8")
print("\nSUMMARY")
print(json.dumps({
 "stage_B_eligible_candidates":summary["stage_B_eligible_candidates"],
 "stage_B_eligible_count":summary["stage_B_eligible_count"]
},indent=2))
