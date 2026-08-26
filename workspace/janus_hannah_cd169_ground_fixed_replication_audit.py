#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path

def gate(state,evidence=None,failure=None):return {'state':state,'evidence':evidence or {},'failure':failure}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('receipt',type=Path);ap.add_argument('--output',required=True,type=Path);a=ap.parse_args();d=json.loads(a.receipt.read_text(encoding='utf-8'));g={}
    v1=d['v1_wireout_track_result'];v2=d['veh_nav_directional_repair_v2'];t=d['v2_target_mapping'];em=d['independent_em12_ground_cell'];k=d['killer_gate_result']
    mag=(v1['abs_layback_delta_median_m']<v1['all_anchor_residual_median_m'])
    g['G5A_LAYBACK_MAGNITUDE_INFORMATION']=gate('INFORMATIVE_PASS' if mag else 'FAIL',{'abs_layback_delta_median_m':v1['abs_layback_delta_median_m'],'position_residual_median_m':v1['all_anchor_residual_median_m']})
    dom=(v2['v2_heldout']['median_residual_m']<v2['v1_same_heldout']['median_residual_m'] and v2['v2_heldout']['p90_residual_m']<v2['v1_same_heldout']['p90_residual_m'] and v2['parameters_tuned_after_heldout_results'] is False)
    g['G5B_V2_HELDOUT_DIRECTIONAL_REPAIR']=gate('PASS_DOMINANCE' if dom else 'FAIL',{'v1':v2['v1_same_heldout'],'v2':v2['v2_heldout'],'parameters_tuned_after_results':v2['parameters_tuned_after_heldout_results']})
    native=(d['janus_gate_state']['G5_NATIVE_PER_PING_VEH_NAV']=='RECOVERED')
    g['G5C_NATIVE_VEH_NAV']=gate('PASS' if native else 'BLOCKED',{'receipt_state':d['janus_gate_state']['G5_NATIVE_PER_PING_VEH_NAV']},None if native else 'NATIVE_WIREOUT_DERIVED_VEH_NAV_NOT_RECOVERED')
    mapok=(t['mapped_raw_sample_index']==821 and t['mapping_used_intensity'] is False and t['distance_reconstructed_vehicle_to_frozen_coordinate_m']<250)
    g['G6_V2_TARGET_MAPPING']=gate('APPROX_PASS_SAMPLE_821' if mapok else 'FAIL',{'sample_index':t['mapped_raw_sample_index'],'vehicle_to_frozen_m':t['distance_reconstructed_vehicle_to_frozen_coordinate_m'],'mapping_used_intensity':t['mapping_used_intensity'],'absolute_ground_truth':False})
    winners=em['winning_files_geometry_only'];cov=(em['coverage_pass'] is True and len(winners)>=1 and all(x['nearest_distance_m']<=250 for x in winners) and em['phase_A_depth_values_used_for_file_selection'] is False)
    g['G7B_EM12_LOCAL_COVERAGE']=gate('PASS' if cov else 'FAIL',{'winner_count':len(winners),'nearest_distances_m':[x['nearest_distance_m'] for x in winners],'depth_based_selection':em['phase_A_depth_values_used_for_file_selection']})
    independent=(cov and len(winners)>=2)
    g['G7B_EM12_MULTILINE_BATHYMETRY']=gate('PASS_INDEPENDENT_LINE_GROUP_COVERAGE' if independent else 'PARTIAL',{'line_groups':[x['file'] for x in winners],'surface_consistency_authority':'DESCRIPTIVE_ONLY_NOT_PREREGISTERED_AS_FEATURE_PASS'})
    feature=bool(k['same_structure_replicated_at_same_Earth_location'])
    g['G7_GROUND_FIXED_FEATURE_REPLICATION']=gate('PASS' if feature else 'HOLD_NOT_ESTABLISHED',{'strong_tobi_feature_at_frozen_cell':k['strong_tobi_feature_at_frozen_cell'],'same_structure_replicated':feature},None if feature else 'NO_STRONG_GEOMETRY_SELECTED_TOBI_FEATURE_AT_FROZEN_CELL_TO_REPLICATE')
    promote=(feature and native)
    g['G8_PROMOTION_TO_GEOMETRY_TEST']=gate('PASS' if promote else 'BLOCKED',{'native_veh_nav':native,'feature_replication':feature},None if promote else 'REQUIRES_DEFENSIBLE_GROUND_FIXED_TOBI_FEATURE_REPLICATION')
    fatal=any(x['state']=='FAIL' for x in g.values());status='FAIL' if fatal else 'PASS_EARTH_CELL_ANCHORED__HOLD_GROUND_FIXED_FEATURE_REPLICATION'
    out={'artifact_id':'JANUS-HANNAH-CD169-GROUND-FIXED-REPLICATION-AUDIT-RUN-001-2026-08-26-v1.0','schema':'janus.cosmos.cousteau.hannah_cd169.ground_fixed_replication.audit.v1','source_receipt':str(a.receipt),'audit_status':status,'gates':g,'janus_state':{'EARTH_CELL_IN_INDEPENDENT_EM12':'ANCHORED' if cov else 'NOT_ANCHORED','EM12_INDEPENDENT_LINE_GROUPS':len(winners),'TOBI_HISTORICAL_NAV_V2':'HELDOUT_DOMINANCE_PASS' if dom else 'FAIL','NATIVE_VEH_NAV':'RECOVERED' if native else 'BLOCKED','V2_FROZEN_CELL_SAMPLE':821 if mapok else None,'TOBI_STRONG_FEATURE_AT_FROZEN_CELL':False,'GROUND_FIXED_FEATURE_REPLICATION':'NOT_ESTABLISHED','PROMOTE_TO_GEOMETRY_TEST':False,'UNDERWATER_PYRAMID_DETECTED':False},'posthoc_nav_fit':False,'recentered':False,'next_gate':d['next_gate'],'hard_rules':d['hard_rules']}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8');print(json.dumps({'audit_status':status,'janus_state':out['janus_state']},indent=2));return 0 if not fatal else 2
if __name__=='__main__':raise SystemExit(main())
