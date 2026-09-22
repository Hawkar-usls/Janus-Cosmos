#!/usr/bin/env python3
import json, requests
from urllib.parse import urlencode
lat,lon=-4.03,-12.25
endpoint="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
params={
 "geometry":f"{lon},{lat}",
 "geometryType":"esriGeometryPoint","inSR":"4326","outSR":"4326",
 "spatialRel":"esriSpatialRelIntersects",
 "outFields":"NCEI_ID,SURVEY_ID,PLATFORM,SOURCE,CHIEF_SCIENTIST,INSTRUMENT,START_TIME,END_TIME,SURVEY_YEAR,DOWNLOAD_URL,SURVEY_AND_VERSION",
 "returnGeometry":"false","f":"json"
}
r=requests.get(endpoint,params=params,timeout=60,headers={"User-Agent":"JANUS-KUSTO-4deg02S-footprint/1.0"})
r.raise_for_status()
d=r.json()
feats=[f.get("attributes",{}) for f in d.get("features",[])]
out={
 "artifact_id":"JANUS-KUSTO-MAR-4DEG02S-EXACT-NCEI-FOOTPRINT-2026-09-22-v1.0",
 "target":{"lat":lat,"lon":lon,"role":"external_known_positive_approximate_published_region"},
 "query_url":endpoint+"?"+urlencode(params),
 "intersections":feats,
 "survey_ids":sorted(set(x.get("SURVEY_ID") for x in feats if x.get("SURVEY_ID"))),
 "claim_ceiling":"EXACT_NCEI_POLYGON_COVERAGE_AT_APPROXIMATE_PUBLISHED_CALIBRATION_COORDINATE"
}
print(json.dumps(out,indent=2))
open("workspace/kusto_cross_survey_out/JANUS-KUSTO-MAR-4DEG02S-EXACT-NCEI-FOOTPRINT-2026-09-22-v1.0.json","w").write(json.dumps(out,indent=2))
