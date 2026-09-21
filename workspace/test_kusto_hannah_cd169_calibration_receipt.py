#!/usr/bin/env python3
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
receipt=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-HANNAH-CD169-REAL-WORLD-CALIBRATION-RECEIPT-2026-09-21-v1.0.json").read_text())
contract=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-HANNAH-CD169-REAL-WORLD-CALIBRATION-CORPUS-CONTRACT-2026-09-21-v1.0.json").read_text())
prereg=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SYNTHETIC-FORWARD-MODEL-V2-PREREG-2026-09-21-v1.0.json").read_text())

checks={
 "receipt_calibration_only": receipt["role"]=="CALIBRATION_ONLY_NOT_INDEPENDENT_VALIDATION",
 "raw_bytes_not_redistributed": receipt["source_integrity"]["raw_source_bytes_redistributed"] is False,
 "v1_negative_retained": receipt["line19_synthetic_v1_vs_native"]["calibration_verdict"]=="FALSIFIED_AS_NATIVE_NAVIGATION_MODEL",
 "v1_native_residual_outside_rounding_envelope": receipt["line19_synthetic_v1_vs_native"]["position_residual_m"]["min"] > receipt["frozen_baseline"]["synthetic_rounding_only_mc_abs_cross_track_p995_m"],
 "old_sample_821_superseded": receipt["frozen_target_native_geometry"]["old_v2_sample_821_target_mapping_disposition"]=="SUPERSEDED_BY_NATIVE_VEH_NAV",
 "no_sonar_intensity": receipt["epistemic_firewall"]["no_sonar_intensity_used_in_this_calibration"] is True,
 "no_morphology_claim": receipt["epistemic_firewall"]["no_morphology_claim"] is True,
 "cd169_cannot_validate_v2": receipt["epistemic_firewall"]["cd169_may_validate_v2_after_calibration"] is False,
 "unseen_validation_required": receipt["epistemic_firewall"]["unseen_real_survey_required_for_v2_validation"] is True,
 "contract_freeze": contract["mandatory_epistemic_firewalls"]["CD169_CANNOT_VALIDATE_SYNTHETIC_V2"] is True,
 "prereg_disjoint_validation": prereg["calibration_dataset"]["may_be_used_for_v2_blind_validation"] is False,
 "no_single_vanity_score": prereg["success_rule"]["overall_v2_winner_score_forbidden"] is True
}
failed=[k for k,v in checks.items() if not v]
print(json.dumps({"checks":checks,"passed":len(checks)-len(failed),"total":len(checks),"failed":failed},indent=2))
raise SystemExit(1 if failed else 0)
