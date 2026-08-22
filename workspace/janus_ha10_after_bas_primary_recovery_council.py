#!/usr/bin/env python3
from __future__ import annotations
import json, hashlib
from datetime import datetime, timezone
from pathlib import Path

REC=Path('data/cousteau/JANUS-BAS-H10N1-PRIMARY-RECOVERY-RUN-001-2026-08-22-v1.0.json')
OUT=Path('data/cousteau/JANUS-HA10-AFTER-BAS-PRIMARY-RECOVERY-COUNCIL-RUN-003-2026-08-22-v1.0.json')
r=json.loads(REC.read_text(encoding='utf-8'))
ex=r['fixed_cell_extractions']
vals=[e['nearest_cell']['value_topography_m'] for e in ex]

result={
 'artifact_id':'JANUS-HA10-AFTER-BAS-PRIMARY-RECOVERY-COUNCIL-RUN-003-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'question':'The authoritative BAS native grid was recovered and H10N1 was extracted consistently in NetCDF and ASCII. What is the next authorized physical stage before any morphology scoring?',
 'evidence_snapshot':{
   'primary_recovery_success':r['success_gate'],
   'source_dataset':r['source']['dataset_id'],
   'netcdf_sha256':ex[0]['sha256'] if ex else None,
   'ascii_sha256':ex[1]['sha256'] if len(ex)>1 else None,
   'h10n1_frozen':r['frozen_target'],
   'nearest_grid_cell':ex[0]['nearest_cell'] if ex else None,
   'duplicate_format_value_difference_m':abs(vals[0]-vals[1]) if len(vals)>1 else None,
   'morphology_scoring_performed':r['morphology_scoring_performed'],
   'southern_p2548_recovered':False,
   'target_identity':'UNCONFIRMED'
 },
 'council':{
   'HRain':'PRESERVE_THE_RECOVERED_DEPTH_AS_A_MEASUREMENT_NODE__DO_NOT_JUMP_TO_IDENTITY',
   'iNaiHR':'MAY_CONNECT_THIS_NODE_TO_ARCHIVE_SEARCH_GRAPH_ONLY__NO_EVIDENCE_WEIGHT',
   'DemiHead':'REQUIRE_RAW_WINDOW_HASHES_AND_GRID_GEOMETRY_CHECK_BEFORE_DERIVED_METRICS',
   'Fast_CAT':'EXTRACT_PREREGISTERED_WINDOWS_WITHOUT_SCORING__THEN_FREEZE_CONTROL_SELECTION',
   'Aura':'ZERO_EVIDENCE_AUTHORITY',
   'Janus_Cosmos':'FROZEN_COORDINATE_REMAINS_IMMUTABLE',
   'Cousteau':'RECOVER_LOCAL_DEPTH_SURFACE_FIRST__DERIVATIVES_SECOND',
   'Fundamentum':'VOLCANIC_BASELINE_REMAINS_ACTIVE',
   'AIFC':'HASH_EACH_RAW_TILE_AND_RECORD_CELL_BOUNDS_CRS_RESOLUTION_AND_NODATA',
   'Voice_of_Janus':'WE_NOW_KNOW_WHERE_THE_FLOOR_IS__NEXT_RECOVER_ITS_LOCAL_SURFACE__DO_NOT_NAME_THE_SHAPE'
 },
 'janus_answer':{
   'AUTHORIZED_STAGE':'H10N1_RAW_LOCAL_SURFACE_EXTRACTION_001',
   'ACTIONS':[ 
      'USE_RECOVERED_BAS_NETCDF_AS_CANONICAL_NATIVE_GRID_AND_ASCII_AS_DUPLICATE_FORMAT_CROSSCHECK',
      'EXTRACT_RAW_UNSCORED_WINDOWS_CENTERED_ON_FROZEN_H10N1_AT_1_2_5_10_KM_RADII',
      'RECORD_EXACT_GRID_BOUNDS_DIMENSIONS_NODATA_VALID_FRACTION_MIN_MAX_ONLY_AS_DATA_INTEGRITY_FIELDS',
      'SHA256_EACH_EXTRACTED_RAW_WINDOW',
      'DO_NOT_COMPUTE_SLOPE_CURVATURE_RUGOSITY_PLANARITY_RADIAL_SYMMETRY_OR_FACET_HISTOGRAM_YET'
   ],
   'PARALLEL_ALLOWED':[ 
      'CONTINUE_P2548_AND_2002_H10S_PRIMARY_GEOREFERENCE_RECOVERY_WITHOUT_SCORING',
      'CONTINUE_WAITING_FOR_BAS_BGS_EARTHSCOPE_EXTERNAL_RESPONSES_AS_INDEPENDENT_CALIBRATION'
   ],
   'INTEGRITY_CORRECTION':'GRID_VALUE_IS_TOPOGRAPHY_POSITIVE_UP__NEGATIVE_1927_934M_MEANS_SEAFLOOR_ELEVATION_APPROX_1927_934M_BELOW_MEAN_SEA_LEVEL__NOT_HYDROPHONE_SENSOR_DEPTH',
   'SUCCESS_GATE':'ALL_FOUR_PREREGISTERED_RAW_WINDOWS_EXTRACTED_REPRODUCIBLY_WITH_PROVENANCE_AND_HASHES',
   'ON_SUCCESS':'STOP_AND_ASK_JANUS_AGAIN_TO_FREEZE_BLIND_CONTROL_SELECTION_BEFORE_MORPHOLOGY_METRICS',
   'TARGET_IDENTITY':'UNCONFIRMED'
 },
 'hard_rules':[ 
   'ASK_JANUS_BEFORE_EACH_NEW_PHYSICAL_STAGE',
   'NO_MORPHOLOGY_SCORING_IN_THIS_STAGE',
   'NO_VISUAL_GUESSING',
   'NO_POSTHOC_RETARGETING',
   'VOLCANIC_BASELINE_REMAINS_ACTIVE',
   'NEGATIVE_RESULTS_REMAIN_NEGATIVE'
 ]
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result,indent=2,ensure_ascii=False))
