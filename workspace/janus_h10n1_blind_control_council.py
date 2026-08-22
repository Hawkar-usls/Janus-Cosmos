#!/usr/bin/env python3
from __future__ import annotations
import json, hashlib
from datetime import datetime, timezone
from pathlib import Path

RAW=Path('data/cousteau/JANUS-H10N1-RAW-LOCAL-SURFACE-EXTRACTION-RUN-001-2026-08-22-v1.0.json')
OUT=Path('data/cousteau/JANUS-H10N1-BLIND-CONTROL-COUNCIL-RUN-004-2026-08-22-v1.0.json')
r=json.loads(RAW.read_text(encoding='utf-8'))
valid={w['radius_km']:w['valid_fraction'] for w in r['windows']}

result={
 'artifact_id':'JANUS-H10N1-BLIND-CONTROL-COUNCIL-RUN-004-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'question':'Raw H10N1 windows are recovered. How should blind volcanic controls and scoring scales be frozen before any morphology metric is computed?',
 'evidence_snapshot':{
  'target_cell_topography_m':-1927.934326171875,
  'valid_fraction_by_radius_km':valid,
  'morphology_metrics_computed':False,
  'source_grid_sha256':r['source']['netcdf_sha256'],
  'target_identity':'UNCONFIRMED'
 },
 'council':{
  'HRain':'MATCH_CONTROLS_ON_MEASUREMENT_CONTEXT_NOT_ON_SHAPE',
  'iNaiHR':'NO_ASSOCIATIVE_FEATURE_MAY_ENTER_CONTROL_SELECTION',
  'DemiHead':'EXCLUDE_OVERLAPPING_TARGET_NEIGHBORHOOD_AND_FAIL_IF_TOO_FEW_CONTROLS',
  'Fast_CAT':'DETERMINISTIC_PRESELECTION_BY_DEPTH_VALIDITY_AND_DISTANCE_ONLY__FREEZE_IDS_BEFORE_SCORING',
  'Aura':'ZERO_EVIDENCE_AUTHORITY',
  'Janus_Cosmos':'TARGET_COORDINATE_AND_WINDOWS_REMAIN_FROZEN',
  'Cousteau':'PRIMARY_SCORE_2KM__SECONDARY_1KM__5KM_10KM_CONTEXT_ONLY_DUE_MASKING',
  'Fundamentum':'VOLCANIC_BASELINE_REMAINS_PRIMARY',
  'AIFC':'CONTROL_MANIFEST_MUST_INCLUDE_SOURCE_HASH_SELECTION_RULE_AND_SEED',
  'Voice_of_Janus':'COMPARE_LIKE_WITH_LIKE_BEFORE_YOU_NAME_ANY_PATTERN'
 },
 'janus_answer':{
  'AUTHORIZED_STAGE':'BLIND_VOLCANIC_CONTROL_SELECTION_001',
  'SCORING_SCALE_POLICY':{
    'primary_radius_km':2,
    'secondary_radius_km':1,
    'context_only_radius_km':[5,10],
    'reason':'1km and 2km target windows are >=99.9% valid; larger windows are materially masked before any morphology inspection'
  },
  'CONTROL_SELECTION_RULE':{
    'dataset':'same BAS native grid GB/NERC/BAS/PDC/01236',
    'candidate_lattice_step_deg':0.02,
    'target_exclusion_radius_km':4,
    'minimum_2km_valid_fraction':0.99,
    'minimum_1km_valid_fraction':0.99,
    'matching_variable':'CENTER_CELL_TOPOGRAPHY_ONLY',
    'depth_tiers_m':[200,400,600],
    'desired_control_count':20,
    'selection_order':'FIRST_TIER_WITH_AT_LEAST_20_CANDIDATES__THEN_DETERMINISTIC_SHA256_ORDER',
    'deterministic_seed_material':'SOURCE_GRID_SHA256 + FROZEN_TARGET_COORDINATES + STRING JANUS_H10N1_CONTROL_V1',
    'forbidden_selection_inputs':['slope','curvature','rugosity','planarity','symmetry','aspect','visual_appearance','symbolic_association']
  },
  'ON_SUCCESS':'FREEZE_CONTROL_COORDINATES_AS_BLIND_IDS_C01_C20_AND_ASK_JANUS_AGAIN_BEFORE_MORPHOLOGY_METRICS',
  'ON_FAILURE':'FREEZE_FAILURE_AND_REASK_JANUS__DO_NOT_RELAX_RULES_POSTHOC',
  'TARGET_IDENTITY':'UNCONFIRMED'
 },
 'hard_rules':['ASK_JANUS_BEFORE_EACH_NEW_PHYSICAL_STAGE','NO_MORPHOLOGY_METRICS_BEFORE_CONTROL_MANIFEST_FREEZE','NO_VISUAL_CONTROL_SELECTION','NO_POSTHOC_RETARGETING','NEGATIVE_RESULTS_REMAIN_NEGATIVE']
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result,indent=2,ensure_ascii=False))
