#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path

def main():
    ap=argparse.ArgumentParser();ap.add_argument('checkpoint',type=Path);ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
    raw=a.checkpoint.read_bytes();d=json.loads(raw);checks=[]
    def ck(name,cond,observed=None):checks.append({'name':name,'pass':bool(cond),'observed':observed})
    ck('status_partial_blocked',d.get('status')=='PARTIAL_PASS_PROVENANCE_DIAGNOSTICS__NATIVE_OR_EXACT_WIREOUT_BLOCKED',d.get('status'))
    ck('wording_boundary_corrected','unresolved link is frozen Earth-cell -> defensible TOBI feature replication' in d['scientific_boundary_correction']['canonical_wording'])
    ck('earth_cell_anchored',d['independent_earth_fixed_anchor']['state']=='PASS')
    ck('G3A_native_blocked',d['v3_gates']['G3A_NATIVE_VEH_NAV']['state']=='BLOCKED_CURRENT_ARCHIVE_NEGATIVE')
    ck('G3B_cdf_fail',d['v3_gates']['G3B_CALIBRATED_NAV_PRODUCT']['state']=='FAIL_CDF_PRE_MRGNAV')
    ck('cdf_latlon_all_zero',d['v3_gates']['G3B_CALIBRATED_NAV_PRODUCT']['evidence']['latlon_exact_zero_zero_rows']==102040 and d['v3_gates']['G3B_CALIBRATED_NAV_PRODUCT']['evidence']['latlon_valid_nonzero_rows']==0)
    ck('false_positive_rejected',d['v3_gates']['G3B_CALIBRATED_NAV_PRODUCT']['evidence']['formal_v3b_pass_disposition'].startswith('REJECTED_FALSE_POSITIVE'))
    ck('G3C_exact_wireout_blocked',d['v3_gates']['G3C_EXACT_WIREOUT_REPRODUCTION']['state']=='BLOCKED_SOURCE_CONFIG_NOT_RECOVERED')
    ck('G3D_diagnostics_pass',d['v3_gates']['G3D_OUTLIER_DIAGNOSTICS']['state']=='PASS_DIAGNOSTIC')
    ck('0631_raw_exact',d['v3_gates']['G3D_OUTLIER_DIAGNOSTICS']['cases']['06:31']['utc']=='2005-02-28T06:31:00Z' and d['v3_gates']['G3D_OUTLIER_DIAGNOSTICS']['cases']['06:31']['record_index']==129)
    ck('G3E_blocked',d['v3_gates']['G3E_FOOTPRINT_RETEST']['state']=='BLOCKED')
    va=d['v2_support_corrected_audit'];ck('v2_support_corrected_dominance',va['relative_dominance_pass'] and va['v2_support_corrected_median_m']<va['v1_comparator_median_m'] and va['v2_support_corrected_p90_m']<va['v1_comparator_p90_m'],{'v1_median':va['v1_comparator_median_m'],'v2_median':va['v2_support_corrected_median_m'],'v1_p90':va['v1_comparator_p90_m'],'v2_p90':va['v2_support_corrected_p90_m']})
    ck('unchanged_model',va['model_parameters_changed'] is False)
    ck('sample_not_promoted',d['frozen_cell_feature_state']['prior_v2_geometry_selected_raw_sample']==821 and d['frozen_cell_feature_state']['v3_authority_to_remap_sample'] is False)
    ck('feature_replication_not_established',d['frozen_cell_feature_state']['ground_fixed_tobi_feature_replication']=='NOT_ESTABLISHED')
    p=d['process_integrity'];ck('integrity_no_recenter_no_fit_no_drop_no_v4',not p['recentered'] and not p['posthoc_nav_fit'] and not p['outliers_dropped'] and not p['V4_navigation_model_created'])
    ck('no_pyramid_claim',d['janus_state']['UNDERWATER_PYRAMID_DETECTED'] is False)
    passed=all(x['pass'] for x in checks)
    out={'schema':'janus.cosmos.cousteau.hannah_cd169.v3_checkpoint.audit.v1','source_checkpoint':str(a.checkpoint),'source_checkpoint_sha256':hashlib.sha256(raw).hexdigest(),'audit_status':'PASS_V3_CHECKPOINT_INVARIANTS' if passed else 'FAIL_V3_CHECKPOINT_INVARIANTS','checks':checks,'passed':sum(x['pass'] for x in checks),'total':len(checks)}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8');print(json.dumps({'audit_status':out['audit_status'],'passed':out['passed'],'total':out['total'],'failed':[x['name'] for x in checks if not x['pass']]},indent=2));return 0 if passed else 2
if __name__=='__main__':raise SystemExit(main())
