#!/usr/bin/env python3
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"data/cousteau"

def load(name):
    return json.loads((D/name).read_text(encoding="utf-8"))

reg=load("JANUS-KUSTO-EXTERNAL-REAL-DATA-CALIBRATION-REGISTRY-2026-09-21-v1.0.json")
cov=load("JANUS-KUSTO-COVERAGE-AND-PROVENANCE-CONTROL-PACK-2026-09-21-v1.0.json")
mb=load("JANUS-KUSTO-MULTIBEAM-PROCESSING-CALIBRATION-PACK-2026-09-21-v1.0.json")
ha=load("JANUS-KUSTO-HYDROACOUSTIC-INSTRUMENT-AWARENESS-CALIBRATION-PACK-2026-09-21-v1.0.json")
mx=load("JANUS-KUSTO-CALIBRATION-LAYER-MATRIX-2026-09-21-v1.0.json")
pre=load("JANUS-KUSTO-SYNTHETIC-FORWARD-MODEL-V2-PREREG-2026-09-21-v1.1.json")

sources={x["id"]:x for x in reg["sources"]}
cases={x["case_id"]:x for x in cov["frozen_cases"]}

checks={
 "cd169_is_not_v2_unseen_validation": sources["BODC_CD169_TOBI_BODCREQ_9406"]["independent_validation_eligible_for_v2_navigation"] is False,
 "jr15001_not_v2_bathy_unseen_validation": sources["BODC_JR15001_ORIGINATOR_EM122_BODCREQ_9408"]["independent_validation_eligible_for_v2_bathymetry_selection"] is False,
 "ukho_negative_is_exactly_scoped": sources["UKHO_HI1751_AUTHORITATIVE_COVERAGE_NEGATIVE"]["frozen_facts"]["fixed_point_intersects_HI1751"] is False,
 "bgs_negative_is_institution_scoped": sources["BGS_ASCENSION_COVERAGE_REFERRAL_IDA310886"]["important_rule"]=="Institution-specific no-coverage is not a global no-data result.",
 "bas_h10s2_no_soundings_retained": "There are no soundings around the secondary H10S2 point." in sources["BAS_PDC_JR15001_RAW_EM122_LINE_0019"]["custodian_facts"],
 "p2548_not_numeric_calibration": sources["FUGRO_P2548_ARCHIVE_ROUTE"]["usable_numeric_calibration_data"] is False,
 "titanic_no_reply_not_calibration": sources["RMS_TITANIC_BLIND_BENCHMARK_REQUEST"]["usable_numeric_calibration_data"] is False,
 "bbox_not_coverage_control": cases["BOUNDING_BOX_NE_EXACT_SURVEY_COVERAGE"]["expected_class"]=="NO_HI1751_INTERSECTION",
 "coverage_hierarchy_ends_ping_beam": cov["coverage_boolean_hierarchy"][-1]=="PING_BEAM_FOOTPRINT_INTERSECTS_POINT",
 "raw_not_corrected": any(x["state"]=="RAW_KONGSBERG_ALL" and "final bathymetric truth" in x["forbidden_authority_without_correction"] for x in mb["processing_state_controls"]),
 "no_invented_svp_sigma": mb["svp_rule"]["numeric_uncertainty_not_provided"] is True,
 "parser_crossformat_control_tight": mb["parser_invariance_control"]["absolute_cross_format_difference_m"] < 1e-6,
 "ha_negative_preserved": ha["frozen_real_negative_control"]["station_event_passes"]==0,
 "ha_no_sideband_passes": ha["frozen_real_negative_control"]["target_sideband_contrast_pass_fraction"]==0.0,
 "ha_calibration_only": "not unseen validation" in ha["reuse_policy"].lower(),
 "layer_matrix_no_same_layer_revalidation": mx["global_validation_firewall"]["same_dataset_may_validate_a_layer_it_calibrated"] is False,
 "same_unseen_data_v1_v2": mx["global_validation_firewall"]["v1_and_v2_must_be_scored_on_same_unseen_data"] is True,
 "no_single_score": mx["global_validation_firewall"]["no_single_overall_score"] is True,
 "v2_qualitative_uncertainty_guard": "QUALITATIVE_UNCERTAINTY_NE_INVENTED_NUMERIC_PRIOR" in pre["new_hard_rules_from_multisource_audit"],
 "v2_claim_ceiling_unseen": "NO_V2_PREDICTIVE_SKILL_CLAIM_UNTIL_UNSEEN_VALIDATION" in pre["claim_ceiling"]
}
failed=[k for k,v in checks.items() if not v]
print(json.dumps({"passed":len(checks)-len(failed),"total":len(checks),"failed":failed,"checks":checks},indent=2))
raise SystemExit(1 if failed else 0)
