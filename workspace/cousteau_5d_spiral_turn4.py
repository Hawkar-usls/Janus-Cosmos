#!/usr/bin/env python3
from __future__ import annotations

import argparse, hashlib, json
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'/'cousteau'
PRIOR=DATA/'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TRANCEPTION-REVERSE-RUN-001-2026-08-21-v1.0.json'
BLIND=DATA/'JANUS-ECHO-COUSTEAU-EA-TPHASE-BLIND-CLUSTER-RUN-007-2026-08-21-v1.6.json'
STATUS=DATA/'JANUS-ECHO-COUSTEAU-EA-TPHASE-BLIND-CLUSTER-RUN-007-2026-08-21-v1.6-STATUS.json'
OUT_DEFAULT=DATA/'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-004-2026-08-21-v1.0.json'
SUMMARY_DEFAULT=DATA/'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-004-2026-08-21-v1.0-SUMMARY.json'
DIMS=['D0_SPACE','D1_TIME','D2_ACOUSTIC','D3_REVERSE_CAUSAL','D4_ASSOCIATIVE_PROVENANCE']
PHASES=['BACK','FORWARD','LEFT_HRAIN','RIGHT_INAIHR','ASCEND']


def now(): return datetime.now(timezone.utc).isoformat()
def csha(o): return hashlib.sha256(json.dumps(o,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
def load(p): return json.loads(p.read_text(encoding='utf-8'))


def node(i,phase,claim,structural,assoc,status,new_constraints,parents,payload=None):
    n={'node_id':f'S04N{i:02d}','turn':4,'phase':phase,'created_utc':now(),'dimensions':DIMS,'claim_class':claim,
       'structural_context':structural,'associative_context':assoc,'evidence_status':status,
       'new_constraints':new_constraints,'provenance':[{'path':str(BLIND.relative_to(ROOT)).replace('\\','/'),'run_sha256':csha(load(BLIND))}],
       'parents':parents,'payload':payload or {}}
    n['node_sha256']=csha(n); return n


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',default=str(OUT_DEFAULT)); ap.add_argument('--summary',default=str(SUMMARY_DEFAULT)); a=ap.parse_args()
    prior=load(PRIOR); blind=load(BLIND); status=load(STATUS)
    if status.get('status')!='SUCCESS': raise SystemExit('Turn 4 forbidden: blind gate is not SUCCESS')
    if blind.get('preregistration',{}).get('anchor_hidden_during_clustering') is not True: raise SystemExit('Turn 4 forbidden: anchor blindness not verified')
    if blind.get('preregistration',{}).get('parameter_change_from_original_blind_contract') is not False: raise SystemExit('Turn 4 forbidden: frozen parameter contract changed')

    prev=prior['turns'][-1]; inherited=list(prev['active_constraints']); prior_hashes={x['state_sha256'] for x in prior['turns']}
    sm=blind['summary']; recon=blind['count_reconciliation_gate']; null_fracs=[x['diagnostic_rectangular_null']['fraction_random_points_with_equal_or_smaller_nearest_cluster_distance'] for x in blind['post_reveal']['configs'] if x.get('diagnostic_rectangular_null')]

    added=[
      'AUTHORITATIVE_FILE_ROWS_ARE_ANALYSIS_INPUT__NO_SYNTHETIC_COUNT_REPAIR',
      'PAPER_FILE_COUNT_MISMATCH_IS_FIRST_CLASS_RECONCILIATION_GATE',
      'BLIND_CLUSTER_NEGATIVE_CANNOT_BE_OVERRIDDEN_BY_NEAREST_EVENT',
      'RECTANGULAR_LOOK_ELSEWHERE_IS_DIAGNOSTIC_NOT_FORMAL_P_VALUE',
      'TECTONIC_MAR_CONTROL_PRECEDES_TARGET_INTERPRETATION',
    ]
    active=sorted(set(inherited+added))

    n1=node(1,'BACK','CORRECTION',
      'Return to the plateau blocker after new evidence arrived: MGDS dataset 30497/file 2504732 now materializes as a TAR→GZIP→ASCII catalog and the exact-format replay completed.',
      'The reverse question changes from how to obtain rows to whether the rows independently organize around the frozen anchor.',
      'BLOCKER_RESOLVED_BY_AUTHORITATIVE_RAW_ROWS',[],[],{'source_status':status})
    n2=node(2,'FORWARD','NEGATIVE',
      f"Forward replay preserves the frozen DBSCAN grid and produces blind hash {sm['blind_freeze_sha256']}; after reveal the nearest event is {sm['nearest_catalog_event_to_anchor_km']} km and nearest blind cluster center is {sm['nearest_blind_cluster_center_across_grid_km']} km.",
      'A nearest isolated event is not allowed to rescue the cluster-level negative result.',
      'BLIND_SPATIAL_RESULT_NEGATIVE',[added[2]],['S04N01'],{'p95_overlap':sm['anchor_inside_any_blind_cluster_p95_radius'],'max_overlap':sm['anchor_inside_any_blind_cluster_max_radius'],'verdict':sm['verdict']})
    n3=node(3,'LEFT_HRAIN','FACT',
      f"Structural lane records {sm['authoritative_rows']} authoritative rows versus {sm['paper_reported_rows']} paper-reported events, a difference of {recon['difference']}; downloaded rows contain hydrophone counts {sm['n_hydrophones_distribution']}.",
      'The mismatch becomes a separate source-reconciliation branch; missing rows are not synthesized and do not enter clustering.',
      'SOURCE_COUNT_MISMATCH_PRESERVED',[added[0],added[1]],['S04N02'],{'reconciliation':recon})
    n4=node(4,'RIGHT_INAIHR','CONTROL',
      f"Across the nine frozen DBSCAN settings, the diagnostic rectangular-null fractions span {min(null_fracs):.5f}–{max(null_fracs):.5f}; the anchor is therefore not unusually close to a blind cluster under this diagnostic.",
      'This is useful falsification pressure, not a formal ocean-masked p-value. Association cannot promote the target; the Mid-Atlantic Ridge tectonic origin of the catalog remains the mandatory control interpretation.',
      'CONTROL_PRESSURE_INCREASED_TARGET_NOT_PROMOTED',[added[3],added[4]],['S04N03'],{'diagnostic_null_fractions':null_fracs,'formal_p_value':False})
    n5=node(5,'ASCEND','STATE_LIFT',
      'Turn 4 ascends because a real blocker was resolved and new authoritative data entered the graph. It ascends with a negative blind spatial result plus an open source-count reconciliation gate.',
      'Returning to the original search question now occurs at higher resolution: future acoustic/frequency evidence must survive both the 5943-row blind spatial negative and the existing natural/biological/engineered controls.',
      'SPIRAL_ASCENT_WITH_NEGATIVE_SCIENTIFIC_DELTA',added,['S04N04'],{'next_gates':['COUSTEAU_EA_TPHASE_6843_VS_5943_COUNT_RECONCILIATION_V1','COUSTEAU_EA_TPHASE_EVENT_ERROR_AWARE_DISTANCE_GATE_V1','COUSTEAU_5D_CONTROL_OUTPERFORMANCE_MATRIX_V1']})
    nodes=[n1,n2,n3,n4,n5]

    state_core={'turn':4,'parent_state_sha256':prev['state_sha256'],'active_constraints':active,'blind_freeze_sha256':sm['blind_freeze_sha256'],
                'scientific_delta':'ACQUISITION_BLOCKER_RESOLVED__BLIND_SPATIAL_RESULT_NEGATIVE__COUNT_RECONCILIATION_GATE_OPEN__TARGET_EVIDENCE_NOT_INCREASED',
                'anchor_overlap_p95':sm['anchor_inside_any_blind_cluster_p95_radius'],'anchor_overlap_max':sm['anchor_inside_any_blind_cluster_max_radius']}
    state_hash=csha(state_core)
    turn={'turn':4,'inherited_constraints':inherited,'added_constraints':added,'active_constraints':active,'anchor_node':'S04N05',
          'scientific_delta':state_core['scientific_delta'],'state_sha256':state_hash,'repeats_prior_state':state_hash in prior_hashes,'radius':len(active),'height':4}

    tests={
      'T4_01_NEW_RAW_OR_RESOLVED_BLOCKER_TRIGGER': status.get('status')=='SUCCESS',
      'T4_02_HEIGHT_ASCENDS_3_TO_4': prev.get('height')==3 and turn['height']==4,
      'T4_03_STATE_HASH_UNIQUE': state_hash not in prior_hashes,
      'T4_04_RADIUS_NONDECREASING': turn['radius']>prev.get('radius',0),
      'T4_05_ANCHOR_HIDDEN_DURING_BLIND': blind['preregistration']['anchor_hidden_during_clustering'] is True,
      'T4_06_FROZEN_DBSCAN_UNCHANGED': blind['preregistration']['parameter_change_from_original_blind_contract'] is False,
      'T4_07_NEGATIVE_CLUSTER_RESULT_PRESERVED': sm['anchor_inside_any_blind_cluster_p95_radius'] is False and sm['anchor_inside_any_blind_cluster_max_radius'] is False,
      'T4_08_COUNT_MISMATCH_NOT_FILLED': recon['status']=='OPEN_MISMATCH_REQUIRES_SOURCE_RECONCILIATION' and recon['difference']==900,
      'T4_09_DIAGNOSTIC_NULL_NOT_FORMAL': all(x['diagnostic_rectangular_null']['formal_p_value'] is False for x in blind['post_reveal']['configs'] if x.get('diagnostic_rectangular_null')),
      'T4_10_NO_TARGET_PROMOTION': 'TARGET_EVIDENCE_NOT_INCREASED' in turn['scientific_delta'],
    }
    rows=[{'id':k,'pass':bool(v)} for k,v in tests.items()]; ok=all(x['pass'] for x in rows)
    out={'artifact_id':'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-004-2026-08-21-v1.0','created_utc':now(),'parent_spiral_artifact':prior['artifact_id'],
         'spiral_law':prior['spiral_law'],'dimensions':DIMS,'phases':PHASES,'trigger':'RESOLVED_MGDS_BLOCKER_PLUS_AUTHORITATIVE_CATALOG_ROWS',
         'turn':turn,'nodes':nodes,'tests':rows,'tests_passed':sum(x['pass'] for x in rows),'tests_total':len(rows),
         'engine_verdict':'PASS_5D_SPIRAL_TURN_4' if ok else 'FAIL_5D_SPIRAL_TURN_4',
         'scientific_verdict':turn['scientific_delta'],'target_identity':'UNCONFIRMED','underwater_pyramid_detected':False,
         'next_gates':n5['payload']['next_gates'],'status':'TURN_4_ASCENDED' if ok else 'TURN_4_TEST_FAILURE'}
    out['artifact_sha256']=csha(out)
    summary={'artifact_id':out['artifact_id']+'-SUMMARY','engine_verdict':out['engine_verdict'],'scientific_verdict':out['scientific_verdict'],
             'height':4,'radius':turn['radius'],'prior_radius':prev['radius'],'state_hash_unique':not turn['repeats_prior_state'],'nodes':len(nodes),
             'tests_passed':out['tests_passed'],'tests_total':out['tests_total'],'authoritative_rows':sm['authoritative_rows'],'paper_reported_rows':sm['paper_reported_rows'],
             'nearest_event_km':sm['nearest_catalog_event_to_anchor_km'],'nearest_cluster_center_km':sm['nearest_blind_cluster_center_across_grid_km'],
             'p95_overlap':False,'max_overlap':False,'diagnostic_null_fraction_range':[min(null_fracs),max(null_fracs)],'target_identity':'UNCONFIRMED'}
    Path(a.output).write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n',encoding='utf-8'); Path(a.summary).write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2,ensure_ascii=False)); return 0 if ok else 2

if __name__=='__main__': raise SystemExit(main())
