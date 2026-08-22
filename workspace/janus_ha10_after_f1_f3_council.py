#!/usr/bin/env python3
import json, hashlib
from datetime import datetime, timezone
from pathlib import Path
OUT=Path('data/cousteau/JANUS-HA10-AFTER-F1-F3-COUNCIL-RUN-001-2026-08-22-v1.0.json')
result={
 'artifact_id':'JANUS-HA10-AFTER-F1-F3-COUNCIL-RUN-001-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'question':'After F1 found no machine-verifiable modern raster intersection and F3 corrected sensor-depth versus seafloor-depth, what physical key should be pursued next?',
 'new_facts':{
  'modern_h10n1_native_intersection_verified':False,
  'ctbto_route_figures_confirm_array_route_context':True,
  'hydrophone_sensor_depth_is_not_seafloor_depth':True,
  'early_southern_hdas_seafloor_depth_estimate_m':2200,
  'p2548_recovered':False,
  'bas_source_grid_bytes_recovered':False
 },
 'council':{
  'HRain':'RECOVER_DEPTH_SURFACE_AT_FROZEN_COORDINATES_BEFORE_SHAPE_IDENTITY',
  'DemiHead':'DO_NOT_INTERPOLATE_INSTALLED_ANCHOR_DEPTH_FROM_EARLY_DESIGN_ESTIMATE',
  'Fast_CAT':'PREFER_SOURCE_RASTER_OR_DIGITAL_SOUNDINGS_OVER_FIGURE_DIGITIZATION',
  'Cousteau':'PRIORITIZE_PRIMARY_BATHYMETRY_ARCHIVE_AND_INSTALLED_ARRAY_SEAFLOOR_GEOMETRY',
  'Fundamentum':'2200M_IS_HISTORICAL_DESIGN_CONTEXT_NOT_MEASUREMENT_AT_FINAL_NODE',
  'AIFC':'REQUIRE_FINAL_AS_BUILT_OR_SURVEY_PROVENANCE',
  'Voice_of_Janus':'LOCATE_FLOOR_THEN_SHAPE_THEN_ACOUSTICS'
 },
 'janus_answer':{
  'NEXT_PRIORITY':'RECOVER_AS_BUILT_OR_PRIMARY_SURVEY_SEAFLOOR_DEPTHS_AROUND_H10S_TRIAD_AND_H10N1',
  'PARALLEL_PRIORITY':'CONTINUE_P2548_AND_BAS_NATIVE_GRID_RECOVERY',
  'APPROVED_SEARCH_TARGETS':['P2548','OCEANS_2003_ORIGINAL_COLOR_FIGURES','HA10_AS_BUILT_ARRAY_ENGINEERING_REPORT','FINAL_H10N_H10S_ROUTE_SURVEY_PRODUCTS','BAS_NATIVE_NETCDF_OR_ASCII_GRID','BGS_UKHO_NATIVE_ASC_GRIDS_IF_COORDINATE_INTERSECTION_CAN_BE_PROVEN'],
  'NOT_APPROVED':['TREAT_2200M_AS_FINAL_H10S_ANCHOR_DEPTH','INFER_WISHBONE_LOCATION_FROM_SEAMOUNT_CENTROID','SCORE_GEOMETRY_FROM_LOW_RESOLUTION_PROXY'],
  'PROMOTION_CONDITION':'AT_LEAST_ONE_PRIMARY_DIGITAL_OR_GEOREFERENCED_SURFACE_AT_FROZEN_COORDINATES',
  'TARGET_IDENTITY':'UNCONFIRMED'
 },
 'hard_rules':['ASK_JANUS_BEFORE_NEXT_NEW_PHYSICAL_STAGE','AS_BUILT_OUTRANKS_DESIGN_ESTIMATE','NO_VISUAL_GUESSING','NO_POSTHOC_RETARGETING','NEGATIVE_RESULTS_REMAIN_NEGATIVE']
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result,indent=2,ensure_ascii=False))
