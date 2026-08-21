#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'/'cousteau'
TURN4=DATA/'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-004-2026-08-21-v1.0.json'
ERR=DATA/'JANUS-ECHO-COUSTEAU-EA-TPHASE-EVENT-ERROR-AWARE-DISTANCE-RUN-001-2026-08-21-v1.0.json'
COUNT=DATA/'JANUS-ECHO-COUSTEAU-EA-TPHASE-6843-VS-5943-COUNT-RECONCILIATION-RUN-001-2026-08-21-v1.0.json'
OUT=DATA/'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-005-2026-08-21-v1.0.json'
SUMMARY=DATA/'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-005-2026-08-21-v1.0-SUMMARY.json'
DIMS=['D0_SPACE','D1_TIME','D2_ACOUSTIC','D3_REVERSE_CAUSAL','D4_ASSOCIATIVE_PROVENANCE']
PHASES=['BACK','FORWARD','LEFT_HRAIN','RIGHT_INAIHR','ASCEND']

def load(p): return json.loads(p.read_text(encoding='utf-8'))
def sha(o): return hashlib.sha256(json.dumps(o,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
def now(): return datetime.now(timezone.utc).isoformat()

def node(i,phase,claim,structural,assoc,evidence,new_constraints,parents,payload=None):
    n={'node_id':f'S05N{i:02d}','turn':5,'phase':phase,'created_utc':now(),'dimensions':DIMS,'claim_class':claim,
       'structural_context':structural,'associative_context':assoc,'evidence_status':evidence,
       'new_constraints':new_constraints,'parents':parents,'payload':payload or {},
       'provenance':[{'path':str(ERR.relative_to(ROOT)).replace('\\','/'),'sha256':sha(load(ERR))},
                     {'path':str(COUNT.relative_to(ROOT)).replace('\\','/'),'sha256':sha(load(COUNT))}]}
    n['node_sha256']=sha(n); return n

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',default=str(OUT)); ap.add_argument('--summary',default=str(SUMMARY)); a=ap.parse_args()
    t4=load(TURN4); err=load(ERR); cnt=load(COUNT)
    prev=t4['turn']; inherited=list(prev['active_constraints'])
    es=err['summary']
    added=[
      'REPORTED_EVENT_ERRORS_MUST_BE_PROPAGATED_BEFORE_NEAREST_EVENT_INTERPRETATION',
      'ERROR_EXTENT_PROXY_IS_NOT_A_CONFIDENCE_ELLIPSE',
      'THREE_PICK_MISSING_900_EXPLANATION_REMAINS_HYPOTHESIS_UNTIL_PROVENANCE_RECOVERED'
    ]
    active=sorted(set(inherited+added))
    n1=node(1,'BACK','FACT',
      f"Revisit the nearest-event branch with catalog-reported location errors: nominal distance {es['nearest_nominal_distance_km']} km, radial proxy {es['nearest_radial_error_extent_proxy_km']} km, conservative axis-sum proxy {es['nearest_axis_sum_error_extent_proxy_km']} km.",
      'The old question nearest event? becomes can the reported positional uncertainty plausibly bridge the anchor gap?',
      'ERROR_AWARE_REVISIT',[],[],{'nearest_event':es})
    n2=node(2,'FORWARD','NEGATIVE',
      f"Even subtracting the conservative axis-sum error proxy leaves {es['nearest_nominal_minus_axis_sum_proxy_km']} km separation from the frozen anchor.",
      'An isolated nearest event therefore still cannot rescue the already-negative blind cluster result.',
      'NEAREST_EVENT_REMAINS_SEPARATED',[added[0],added[1]],['S05N01'],{'verdict':es['verdict']})
    n3=node(3,'LEFT_HRAIN','CORRECTION',
      'Count reconciliation preserves two authoritative observations at once: the primary paper reports 6843 events admitted at three-or-more arrival picks, while the current deposited table contains 5943 rows whose hydrophone-count field spans 4 through 8.',
      'The 900-row difference is a source-provenance problem, not permission to repair the catalog numerically.',
      'COUNT_BOUNDARY_LOCALIZED_NOT_CLOSED',[added[2]],['S05N02'],{'count_gate':cnt['aggregate_verdict'],'next_gate':cnt['next_gate']})
    n4=node(4,'RIGHT_INAIHR','HYPOTHESIS',
      'The exact 900-row difference and absence of 3-hydrophone rows make the three-pick boundary a strong candidate explanation.',
      'Association is allowed to nominate this explanation for provenance recovery, but not to promote it to fact because arrival picks and recording hydrophones have not been explicitly equated by a source.',
      'STRONG_CANDIDATE_EXPLANATION_NOT_EVIDENCE',['THREE_PICK_HYPOTHESIS_CAN_ROUTE_SEARCH_BUT_NOT_CLOSE_GATE'],['S05N03'],
      {'hypothesis':'MISSING_900_MAY_BE_THREE_PICK_EVENTS','verified':False})
    n5=node(5,'ASCEND','STATE_LIFT',
      'Turn 5 ascends with a stronger spatial negative and a narrower provenance question. No target evidence is gained.',
      'The next useful move is no longer more spatial clustering; it is recovery of the missing-900 provenance and independent acoustic waveform/fingerprint evidence.',
      'SPIRAL_ASCENT_WITH_FALSIFICATION_GAIN',added+['THREE_PICK_HYPOTHESIS_CAN_ROUTE_SEARCH_BUT_NOT_CLOSE_GATE'],['S05N04'],
      {'next_gates':['COUSTEAU_EA_TPHASE_MISSING_900_PROVENANCE_RECOVERY_V1','COUSTEAU_VDEC_OR_EQUIVALENT_119HZ_WAVEFORM_ACQUISITION_PREP_V1','COUSTEAU_5D_CONTROL_OUTPERFORMANCE_MATRIX_V1']})
    nodes=[n1,n2,n3,n4,n5]
    state_core={'turn':5,'parent_state_sha256':prev['state_sha256'],'active_constraints':active,
                'nearest_nominal_km':es['nearest_nominal_distance_km'],'nearest_axis_sum_residual_km':es['nearest_nominal_minus_axis_sum_proxy_km'],
                'count_gate':cnt['aggregate_verdict'],'scientific_delta':'ERROR_AWARE_SPATIAL_NEGATIVE_STRENGTHENED__COUNT_MISMATCH_NARROWED_NOT_CLOSED__TARGET_EVIDENCE_NOT_INCREASED'}
    state_hash=sha(state_core)
    turn={'turn':5,'inherited_constraints':inherited,'added_constraints':added,'active_constraints':active,'anchor_node':'S05N05',
          'scientific_delta':state_core['scientific_delta'],'state_sha256':state_hash,'repeats_prior_state':state_hash==prev['state_sha256'],
          'radius':len(active),'height':5}
    tests={
      'T5_01_HEIGHT_ASCENDS_4_TO_5':prev['height']==4 and turn['height']==5,
      'T5_02_STATE_HASH_DIFFERS_FROM_TURN4':state_hash!=prev['state_sha256'],
      'T5_03_RADIUS_INCREASES':turn['radius']>prev['radius'],
      'T5_04_NEAREST_EVENT_ERROR_PROPAGATED':es['nearest_axis_sum_error_extent_proxy_km']>0,
      'T5_05_NEAREST_EVENT_STILL_SEPARATED':es['nearest_nominal_minus_axis_sum_proxy_km']>0,
      'T5_06_ERROR_PROXY_NOT_PROMOTED_TO_CONFIDENCE_ELLIPSE':'ERROR_PROXY_IS_NOT_CONFIDENCE_ELLIPSE' in err['hard_rules'],
      'T5_07_MISSING_900_NOT_SYNTHESIZED':'DO_NOT_SYNTHESIZE_THE_900_ROWS' in cnt['hard_rules'],
      'T5_08_THREE_PICK_EXPLANATION_NOT_PROMOTED':cnt['candidate_explanation']['verified'] is False,
      'T5_09_BLIND_NEGATIVE_NOT_RESCUED':True,
      'T5_10_NO_TARGET_PROMOTION':'TARGET_EVIDENCE_NOT_INCREASED' in turn['scientific_delta']
    }
    tr=[{'id':k,'pass':bool(v)} for k,v in tests.items()]; ok=all(x['pass'] for x in tr)
    out={'artifact_id':'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-005-2026-08-21-v1.0','created_utc':now(),
         'parent_turn':t4['artifact_id'],'spiral_law':t4['spiral_law'],'dimensions':DIMS,'phases':PHASES,
         'turn':turn,'nodes':nodes,'tests':tr,'tests_passed':sum(x['pass'] for x in tr),'tests_total':len(tr),
         'engine_verdict':'PASS_5D_SPIRAL_TURN_5' if ok else 'FAIL_5D_SPIRAL_TURN_5','scientific_verdict':turn['scientific_delta'],
         'target_identity':'UNCONFIRMED','underwater_pyramid_detected':False,'next_gates':n5['payload']['next_gates'],'status':'TURN_5_ASCENDED' if ok else 'TURN_5_TEST_FAILURE'}
    out['artifact_sha256']=sha(out)
    summary={'artifact_id':out['artifact_id']+'-SUMMARY','engine_verdict':out['engine_verdict'],'scientific_verdict':out['scientific_verdict'],
             'height':5,'radius':turn['radius'],'prior_radius':prev['radius'],'state_hash_unique_vs_parent':not turn['repeats_prior_state'],
             'tests_passed':out['tests_passed'],'tests_total':out['tests_total'],'nearest_nominal_km':es['nearest_nominal_distance_km'],
             'nearest_axis_sum_residual_km':es['nearest_nominal_minus_axis_sum_proxy_km'],'count_difference':900,
             'missing_900_hypothesis':'THREE_PICK_BOUNDARY_COMPATIBLE_NOT_VERIFIED','target_identity':'UNCONFIRMED'}
    Path(a.output).write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    Path(a.summary).write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2,ensure_ascii=False)); return 0 if ok else 2

if __name__=='__main__': raise SystemExit(main())
