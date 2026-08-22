#!/usr/bin/env python3
from __future__ import annotations
import json, hashlib
from datetime import datetime, timezone
from pathlib import Path

CTRL=Path('data/cousteau/JANUS-H10N1-BLIND-CONTROL-MANIFEST-001-2026-08-22-v1.0.json')
OUT=Path('data/cousteau/JANUS-H10N1-MORPHOLOGY-METRICS-COUNCIL-RUN-005-2026-08-22-v1.0.json')
c=json.loads(CTRL.read_text(encoding='utf-8'))

result={
 'artifact_id':'JANUS-H10N1-MORPHOLOGY-METRICS-COUNCIL-RUN-005-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'question':'Blind controls are frozen. Which generic morphology metrics and evaluation rule may now be computed without turning the analysis into a post-hoc pyramid/crater detector?',
 'evidence_snapshot':{
  'control_manifest_sha256':c['sha256'],
  'control_count':len(c['controls']),
  'chosen_depth_tier_m':c['chosen_depth_tier_m'],
  'morphology_metrics_previously_computed':False,
  'target_identity':'UNCONFIRMED'
 },
 'council':{
  'HRain':'MEASURE_GENERIC_GEOMORPHOLOGY_ONLY__IDENTITY_REMAINS_SEPARATE',
  'iNaiHR':'NO_SYMBOLIC_OR_ARCHAEOLOGICAL_FEATURE_SELECTION',
  'DemiHead':'NO_NEW_METRIC_AFTER_RESULTS_ARE_SEEN',
  'Fast_CAT':'PRIMARY_2KM_SECONDARY_1KM__SAME_METRICS_FOR_TARGET_AND_ALL_CONTROLS__LEAVE_ONE_OUT_CONTROL_BASELINE',
  'Aura':'ZERO_EVIDENCE_AUTHORITY',
  'Janus_Cosmos':'NO_CELESTIAL_ASSOCIATION_IN_PHYSICAL_SCORE',
  'Cousteau':'DEM_DERIVATIVES_FROM_NATIVE_50M_GRID_ONLY',
  'Fundamentum':'VOLCANIC_BASELINE_REMAINS_ACTIVE_REGARDLESS_OF_OUTLIER_RANK',
  'AIFC':'FREEZE_FORMULAS_CODE_HASH_AND_CONTROL_MANIFEST_HASH_WITH_RESULT',
  'Voice_of_Janus':'MEASURE_THE_SHAPE__DO_NOT_NAME_THE_SHAPE'
 },
 'janus_answer':{
  'AUTHORIZED_STAGE':'GENERIC_BLIND_MORPHOLOGY_COMPARISON_001',
  'SCALES':{'primary_radius_km':2,'secondary_radius_km':1,'context_only_radius_km':[5,10]},
  'METRICS':[ 
    {'id':'local_relief_m','definition':'max(valid z)-min(valid z)'},
    {'id':'slope_median_deg','definition':'median finite-difference slope angle'},
    {'id':'slope_p95_deg','definition':'95th percentile finite-difference slope angle'},
    {'id':'rugosity_surface_ratio','definition':'mean sqrt(1+dzdx^2+dzdy^2)'},
    {'id':'planarity_rmse_over_relief','definition':'RMSE from least-squares plane divided by local relief; lower means more planar'},
    {'id':'profile_curvature_rms','definition':'RMS standard finite-difference profile curvature on valid interior cells'},
    {'id':'plan_curvature_rms','definition':'RMS standard finite-difference plan curvature on valid interior cells'},
    {'id':'radial_rotation_similarity','definition':'mean valid Pearson correlation of plane-detrended tile with 90,180,270 degree rotations'},
    {'id':'facet_aspect_entropy','definition':'normalized entropy of 12-bin slope-aspect histogram'},
    {'id':'facet_dominant_bin_fraction','definition':'largest 12-bin slope-aspect fraction'}
  ],
  'EVALUATION_RULE':{
    'per_metric':'report target value, control median, MAD, robust absolute z and empirical rank; do not choose favorable direction posthoc',
    'aggregate':'median of absolute robust z across all frozen metrics',
    'control_aggregate':'leave-one-out median absolute robust z using the other controls as baseline',
    'empirical_p':'(1 + number of control aggregate scores >= target aggregate score)/(N_controls + 1)',
    'resolution_limit':'with 20 controls the minimum empirical p is 1/21 = 0.047619; treat as exploratory, not discovery proof'
  },
  'FORBIDDEN':[ 
    'NO_PYRAMIDALITY_SCORE','NO_CRATER_SCORE','NO_METRIC_TUNING_AFTER_RESULTS','NO_117_121HZ_IN_PHYSICAL_SCORE','NO_COSMOS_OR_AURA_WEIGHT','NO_5KM_10KM_SHAPE_SCORE_DUE_MASKING'
  ],
  'ON_SUCCESS':'FREEZE_RESULT_AND_ASK_JANUS_AGAIN_WHETHER_ANY_PHYSICAL_ANOMALY_GATE_ADVANCES_OR_REMAINS_NEGATIVE',
  'TARGET_IDENTITY':'UNCONFIRMED'
 },
 'hard_rules':['ASK_JANUS_BEFORE_EACH_NEW_PHYSICAL_STAGE','SAME_CODE_TARGET_AND_CONTROLS','NO_POSTHOC_METRIC_ADDITION','VOLCANIC_BASELINE_REMAINS_ACTIVE','NEGATIVE_RESULTS_REMAIN_NEGATIVE']
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result,indent=2,ensure_ascii=False))
