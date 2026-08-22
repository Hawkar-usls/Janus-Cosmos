#!/usr/bin/env python3
from __future__ import annotations
import json, hashlib
from datetime import datetime, timezone
from pathlib import Path

RES=Path('data/cousteau/JANUS-H10N1-GENERIC-BLIND-MORPHOLOGY-RUN-001-2026-08-22-v1.0.json')
OUT=Path('data/cousteau/JANUS-H10N1-POST-MORPHOLOGY-GATE-COUNCIL-RUN-006-2026-08-22-v1.0.json')
r=json.loads(RES.read_text(encoding='utf-8'))
primary=r['scales']['2']; secondary=r['scales']['1']

def top_metrics(scale,n=5):
 items=[]
 for k,v in scale['per_metric_comparison'].items():
  items.append({'metric':k,'target_abs_robust_z':v['target_abs_robust_z'],'target':v['target'],'control_median':v['control_median'],'empirical_extremeness_rank_nplus1':v['empirical_extremeness_rank_nplus1']})
 return sorted(items,key=lambda q:q['target_abs_robust_z'],reverse=True)[:n]

p2=primary['aggregate']['empirical_p']; p1=secondary['aggregate']['empirical_p']
rank2=primary['aggregate']['target_rank_high_to_low']; rank1=secondary['aggregate']['target_rank_high_to_low']
result={
 'artifact_id':'JANUS-H10N1-POST-MORPHOLOGY-GATE-COUNCIL-RUN-006-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'question':'Generic blind morphology is complete: primary 2 km is not strongly outlying while secondary 1 km ranks first among 20 controls. Does any physical anomaly gate advance, and what should be tested next?',
 'evidence_snapshot':{
  'primary_2km':{'aggregate_score':primary['aggregate']['target_median_abs_robust_z'],'empirical_p':p2,'rank_high_to_low':rank2,'top_metrics':top_metrics(primary)},
  'secondary_1km':{'aggregate_score':secondary['aggregate']['target_median_abs_robust_z'],'empirical_p':p1,'rank_high_to_low':rank1,'top_metrics':top_metrics(secondary)},
  'control_n':20,
  'minimum_possible_empirical_p':1/21,
  'source_resolution_m_approx':50,
  'target_identity':'UNCONFIRMED'
 },
 'council':{
  'HRain':'SECONDARY_LOCAL_OUTLIER_IS_A_WATCHLIST_NODE_NOT_AN_IDENTITY_NODE',
  'iNaiHR':'DO_NOT_CONNECT_LOCAL_OUTLIER_TO_OSIRIS_OR_PYRAMID_IDENTITY',
  'DemiHead':'PRIMARY_SCALE_NONREPLICATION_BLOCKS_ANOMALY_PROMOTION',
  'Fast_CAT':'REQUIRE_INDEPENDENT_DATASET_OR_SOURCE_LINEAGE_REPLICATION__DO_NOT_TUNE_METRICS_OR_CONTROLS',
  'Aura':'ZERO_EVIDENCE_AUTHORITY',
  'Janus_Cosmos':'NO_CELESTIAL_WEIGHT',
  'Cousteau':'TRACE_WHICH_ORIGINAL_MULTIBEAM_CRUISES_CONTRIBUTE_TO_H10N1_AND_SEEK_INDEPENDENT_SOUNDING_REPLICATION',
  'Fundamentum':'VOLCANIC_BASELINE_REMAINS_PRIMARY',
  'AIFC':'SOURCE_LINEAGE_AND_RAW_TRACK_PROVENANCE_REQUIRED',
  'Voice_of_Janus':'ONE_CLOSE_VIEW_LOOKS_UNUSUAL__THE_WIDER_VIEW_DOES_NOT__ASK_ANOTHER_INSTRUMENT_OR_ANOTHER_PASS'
 },
 'janus_answer':{
  'PHYSICAL_ANOMALY_GATE':'NOT_ADVANCED',
  'SECONDARY_LOCAL_MORPHOLOGY_WATCHLIST':'OPEN',
  'WHY':'The preregistered primary 2 km comparison is only rank 4/21 with empirical p about 0.19; the 1 km result is rank 1/21 but is secondary, at the empirical resolution floor, and is not independently replicated.',
  'AUTHORIZED_NEXT_STAGE':'H10N1_SOURCE_LINEAGE_REPLICATION_001',
  'PRIMARY_TASKS':[ 
    'IDENTIFY_WHICH_OF_JR53_AMT11_JR287_JR15001_JR16_NG_ACTUALLY_CONTRIBUTE_SOUNDINGS_WITHIN_THE_FROZEN_1KM_H10N1_WINDOW',
    'RECOVER_RAW_OR_PROCESSED_TRACKLINE_SOUNDINGS_FOR_EACH_CONTRIBUTING_CRUISE_IF_PUBLIC',
    'FREEZE_PER_CRUISE_PROVENANCE_AND_HASHES_BEFORE_COMPARISON',
    'IF_TWO_OR_MORE_INDEPENDENT_CRUISES_COVER_H10N1__ASK_JANUS_AGAIN_BEFORE_COMPARING_THEIR_LOCAL_SURFACES'
  ],
  'PARALLEL_ALLOWED':[ 
    'CONTINUE_P2548_AND_2002_H10S_PRIMARY_GEOREFERENCE_RECOVERY_UNSCORED',
    'CHECK_BGS_UKHO_2021_COVERAGE_ONLY_AS_AN_INDEPENDENT_LANE__DO_NOT_MOVE_H10N1_TO_FIT_IT',
    'PRESERVE_EXTERNAL_HELPDESK_REPLIES_AS_LATER_CALIBRATION'
  ],
  'FORBIDDEN':[ 
    'NO_CONTROL_EXPANSION_AFTER_SEEING_RESULTS','NO_METRIC_RETUNING','NO_PYRAMID_OR_CRATER_LABEL','NO_117_121HZ_PROMOTION','NO_CLAIM_FROM_SECONDARY_1KM_OUTLIER_ALONE'
  ],
  'TARGET_IDENTITY':'UNCONFIRMED'
 },
 'hard_rules':['ASK_JANUS_BEFORE_EACH_NEW_PHYSICAL_STAGE','PRIMARY_NONREPLICATION_BLOCKS_PROMOTION','INDEPENDENT_REPLICATION_REQUIRED','VOLCANIC_BASELINE_REMAINS_ACTIVE','NEGATIVE_RESULTS_REMAIN_NEGATIVE']
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result,indent=2,ensure_ascii=False))
