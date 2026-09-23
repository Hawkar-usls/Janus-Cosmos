#!/usr/bin/env python3
import json, requests
from pathlib import Path
from shapely.geometry import shape, Polygon, MultiPolygon, mapping
from shapely.ops import unary_union

OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
URL="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
SURVEYS=["PF0501","EX2105"]

def fetch(sid):
    p={
      "where":f"SURVEY_ID='{sid}'",
      "outFields":"SURVEY_ID,PLATFORM,INSTRUMENT,SURVEY_YEAR,START_TIME,END_TIME,SOURCE,DOWNLOAD_URL",
      "returnGeometry":"true",
      "outSR":"4326",
      "f":"geojson"
    }
    r=requests.get(URL,params=p,timeout=120)
    r.raise_for_status()
    d=r.json()
    if "error" in d: raise RuntimeError(d["error"])
    return d

raw={sid:fetch(sid) for sid in SURVEYS}
geoms={}
features={}
for sid,d in raw.items():
    fs=d.get("features",[])
    features[sid]=[{"properties":f.get("properties",{})} for f in fs]
    gs=[shape(f["geometry"]) for f in fs if f.get("geometry")]
    geoms[sid]=unary_union(gs) if gs else None

a,b=geoms["PF0501"],geoms["EX2105"]
inter=None if a is None or b is None else a.intersection(b)
out={
 "artifact_id":"JANUS-KUSTO-UNSEEN-PF0501-X-EX2105-EXACT-FOOTPRINT-OVERLAP-2026-09-23-v1.0",
 "parent_prereg":"data/cousteau/JANUS-KUSTO-UNSEEN-MULTIBEAM-PF0501-X-EX2105-SOURCE-INVENTORY-PREREG-2026-09-23-v1.0.json",
 "depth_values_read":False,
 "sources":{"arcgis_layer":URL},
 "survey_feature_metadata":features,
 "PF0501_geometry_present":a is not None,
 "EX2105_geometry_present":b is not None,
 "PF0501_bounds":None if a is None else list(a.bounds),
 "EX2105_bounds":None if b is None else list(b.bounds),
 "intersection_empty":True if inter is None else inter.is_empty,
 "intersection_geometry_type":None if inter is None else inter.geom_type,
 "intersection_bounds":None if inter is None or inter.is_empty else list(inter.bounds),
 "intersection_area_deg2":None if inter is None else float(inter.area),
 "intersection_geojson":None if inter is None or inter.is_empty else mapping(inter),
 "exact_footprint_overlap_pass":bool(inter is not None and not inter.is_empty and inter.area>0),
 "next_gate":"FREEZE_DETERMINISTIC_OVERLAP_SUBSET_BEFORE_DEPTH" if inter is not None and not inter.is_empty and inter.area>0 else "PAIR_REJECT_NO_EXACT_FOOTPRINT_OVERLAP",
 "claim_ceiling":"UNSEEN_VALIDATION_EXACT_COVERAGE_ONLY"
}
p=OUT/"JANUS-KUSTO-UNSEEN-PF0501-X-EX2105-EXACT-FOOTPRINT-OVERLAP-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!="intersection_geojson"},indent=2))
