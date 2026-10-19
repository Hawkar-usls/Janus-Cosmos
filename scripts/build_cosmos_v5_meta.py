#!/usr/bin/env python3
"""Create detailed V5 meta-registry JSON after validated Cosmos PR merge."""

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = {
    "caldb_gbm_manifest_builder": "scripts/build_cosmos_v5_dataset_manifest.py",
    "caldb_gbm_manifest": "data/cosmos/COSMOS-ALI-V5-FERMI-LAT-GBM-CALDB-MANIFEST-2026-10-19.json",
    "instrument_calibration_replay_engine": "workspace/janus_cosmos_fermi_calibration_v5.py",
    "real_instrument_pulse_receipt": "data/cosmos/COSMOS-ALI-V5-CALDB-GBM-LAT-INSTRUMENT-PULSE-2026-10-19.json",
    "joint_source_lag_toy_engine": "workspace/janus_cosmos_fermi_joint_v5.py",
    "source_lag_toy_receipt": "data/cosmos/COSMOS-ALI-V5-SOURCE-LAG-TOY-SENSITIVITY-2026-10-19.json",
    "instrument_science_guide": "docs/JANUS-COSMOS-FERMI-LAT-GBM-V5-2026-10-19.md",
    "v5_tests": "tests/test_cosmos_fermi_calibration_v5.py",
    "v4_manifest": "data/cosmos/COSMOS-ALI-V4-FERMI-RAW-FITS-MANIFEST-2026-10-18.json",
    "v4_event_receipt": "data/cosmos/COSMOS-ALI-V4-FERMI-RAW-EVENT-DIAGNOSTICS-2026-10-18.json",
    "v4_toy": "data/cosmos/COSMOS-ALI-V4-JOINT-TOY-LIKELIHOOD-2026-10-18.json",
    "frozen_v3_result": "data/cosmos/COSMOS-ALI-V3-FERMI-PUBLISHED-TOF-METHOD-RESULT-2026-10-18.json",
    "frozen_v2_result": "data/cosmos/COSMOS-ALI-V2-LATTICE-TOF-REVEAL-2026-10-18.json",
    "frozen_v1_result": "data/cosmos/COSMOS-ALI-SIMULATION-AUDIT-RESULT-2026-10-18.json",
}


def sha(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def digest(x: dict) -> str:
    return sha(json.dumps(x, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8"))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cosmos-merge-sha", required=True)
    p.add_argument("--cosmos-pr", required=True, type=int)
    p.add_argument("--ci-passed", required=True, type=int)
    p.add_argument("--tests-passed", required=True, type=int)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.ci_passed < 1 or args.tests_passed < 87 or len(args.cosmos_merge_sha) != 40:
        raise ValueError("V5 requires verified Cosmos merge and >=87 passing unit tests")
    files = {k: {"path": rel, "sha256": sha((ROOT / rel).read_bytes()),
                 "github_url": "https://github.com/Hawkar-usls/Janus-Cosmos/blob/main/" + rel}
             for k, rel in FILES.items()}
    manifest = json.loads((ROOT / FILES["caldb_gbm_manifest"]).read_text())
    obs = json.loads((ROOT / FILES["real_instrument_pulse_receipt"]).read_text())
    toy = json.loads((ROOT / FILES["source_lag_toy_receipt"]).read_text())
    for data in (obs, toy):
        value = dict(data)
        asserted = value.pop("result_sha256")
        if digest(value) != asserted:
            raise ValueError("V5 frozen digest mismatch")
    if manifest["v4_frozen_manifest_sha256"] != files["v4_manifest"]["sha256"] or manifest["v4_frozen_observed_json_sha256"] != files["v4_event_receipt"]["sha256"]:
        raise ValueError("frozen V4 lineage mismatch")
    if len(manifest["lat_caldb"]["caldb_files"]) != 9 or len(manifest["gbm_bursts"]) != 3:
        raise ValueError("missing P8R3_SOURCE_V3/GBM input families")
    if obs["counts"]["gbm_tte_events_all_selected_detectors"] != 6962621 or obs["counts"]["lat_v4_selected_roi_candidates"] != 175:
        raise ValueError("observed event counts changed")
    if obs["scientific_gates"]["joint_observational_liv_hypothesis_test"] != "BLOCKED_NO_VALIDATED_LIKELIHOOD" or toy["observational_liv_p_value"] is not None:
        raise ValueError("hypothesis promoted without observed physical model")
    if toy["unconstrained_source_lag_true_identifiability"]["eta_mle"] is not None:
        raise ValueError("exact intrinsic lag degeneracy broken")

    registry = {
        "artifact_uuid": "JANUS-COSMOS-ALI-FERMI-LAT-GBM-CALDB-MULTIBURST-V5-2026-10-19",
        "registry_schema_id": "janus.registry.physics.fermi_lat_gbm_caldb_source_lag.v5",
        "date_utc": "2026-10-19", "version": "1.0.0", "owner": "Hawkar-usls",
        "title": "JANUS Cosmos V5 — source-verified LAT P8R3_SOURCE_V3 CALDB grids, 6.96M real GBM TTE events, CSPEC/RSP2 metadata, intrinsic-lag identifiability control",
        "experiment_id": manifest["experiment_id"],
        "classification": "OBSERVATIONAL_INSTRUMENT_AND_GBM_LIGHTCURVE_DATA_WITH_SYNTHETIC_ONLY_LIV_PRIOR_SENSITIVITY",
        "scientific_status": "REAL_INSTRUMENT_INGEST_PASS_TRUE_LIV_FORWARD_MODEL_NOT_YET_BUILT_NO_DETECTION",
        "research_lineage": {
            "original_video": "https://www.youtube.com/watch?v=QBpuorUxuSM",
            "source_status": "CONCEPTUAL_MOTIVATION_NOT_EMPIRICAL_EVIDENCE",
            "predecessor_cosmos_v4_pr": "https://github.com/Hawkar-usls/Janus-Cosmos/pull/40",
            "predecessor_registry_v4_pr": "https://github.com/Hawkar-usls/janus-meta-registry/pull/259",
            "predecessor_cosmos_v4_commit": "034c36488b85abbd4c67d75a426af9e3f8af16a2",
            "v4_observed_sha256": files["v4_event_receipt"]["sha256"],
            "v4_manifest_sha256": files["v4_manifest"]["sha256"],
            "historical_receipt_hashes": {k: files[k]["sha256"] for k in ("frozen_v1_result", "frozen_v2_result", "frozen_v3_result", "v4_toy")},
            "prior_frozen_results_unchanged": True,
            "osiris_authority": "UNMODIFIED"},
        "official_observational_sources": {
            "provider_lat": "FERMITOOLS_OFFICIAL_FERMI_CHANNEL_PACKAGE",
            "lat_irf_name": manifest["lat_caldb"]["irf_name"],
            "irf_source_package": manifest["lat_caldb"]["package_url"],
            "irf_source_package_sha256": manifest["lat_caldb"]["package_sha256"],
            "lat_caldb_all_nine_tabular_families": manifest["lat_caldb"]["caldb_files"],
            "provider_gbm": "NASA_HEASARC_FERMI_GBM_TRIGGER_ARCHIVE",
            "gbm_all_twenty_seven_tte_cspec_rsp2_files": manifest["gbm_bursts"],
            "byte_pinning": "PER_FILE_SHA256_AND_BYTE_SIZE_FAIL_CLOSED",
            "source_binary_retention": "OFFICIAL_UPSTREAM_BINARIES_REFERENCED_BY_HASH_NOT_COMMITTED_IN_GIT",
            "preregistered_energy_bands_keV": manifest["predeclared_diagnostic_bands_keV"],
            "preregistered_prompt_windows_s": manifest["prompt_windows_s_relative_trigger"],
            "preregistered_background_window_s": manifest["background_window_relative_trigger_s"],
            "preregistered_time_bin_width_s": manifest["pulse_bin_width_s"],
            "complete_manifest": manifest},
        "actual_gbm_observations_and_response_audit": {
            "total_real_tte_events": obs["counts"]["gbm_tte_events_all_selected_detectors"],
            "real_v4_lat_roi_candidates": obs["counts"]["lat_v4_selected_roi_candidates"],
            "instrument_roles_count": obs["counts"],
            "gbm_tte_units": "SECONDS_RELATIVE_TO_EXACT_TRIGTIME",
            "gbm_rsp2_tstart_tstop_units": "FERMI_MET_SECONDS",
            "time_alignment": "TRIGTIME_FROM_TTE_CSPEC_RSP2_MATCHED_PER_DETECTOR_AND_COHORT",
            "pretrigger_background_window": "[-20,-5]_s_RELATIVE_TRIGGER_STATIONARY_RATE_PROXY_ONLY",
            "per_grb_per_detector_counts_energy_edges_rsp2_and_complete_time_series": obs["gbm_bursts"],
            "lat_irf_table_shape_and_aeff_diagnostics": obs["lat_caldb"],
            "real_lat_sample_effective_area_m2": obs["lat_real_candidate_aeff_diagnostics"],
            "specific_calibration_mismatches": [
                "130427A CSPEC_v01 channel edges versus TTE_v00 up to ~2.15%",
                "some RSP2 EBOUNDS differ from GBM TTE by detector, stored quantitatively",
                "090902B RSP2 terminal matrix extension has zero-duration interval"],
            "rsp2_matrix_spectral_forward_folding": "NOT_PERFORMED",
            "real_instrument_receipt_sha256": files["real_instrument_pulse_receipt"]["sha256"],
            "technical_qualifications": [
                "NaI/BGO channel energy overlaps are proxies, not spectral deconvolution",
                "GBM chosen detector cohort is pilot, not optimized per-burst spectral fit",
                "LAT nearest FRONT/BACK effective-area bin does not equal exposure",
                "LAT PSF/EDISP tables were parsed and numerically validated but kernels not evaluated",
                "Real LAT ROI candidates have not been assigned probabilistic GRB source memberships"]},
        "joint_identifiability_and_hypothesis_controls": {
            "structural_model": toy["source_lag_parameterization"],
            "exact_unidentifiable_reparameterization": toy["exact_reparameterization"],
            "prior_family_is_hypothetical": True,
            "source_lag_prior_sensitivity_synthetic_only": toy["source_lag_prior_sensitivity_synthetic_only"],
            "synthetic_liv_injection_unconstrained": toy["unconstrained_source_lag_true_identifiability"],
            "synthetic_source_only_mimic": toy["source_only_mimic_unconstrained"],
            "naive_zero_source_lag_false_positive": toy["synthetic_zero_source_prior_false_positive_demo"],
            "real_gbm_counts_used_to_estimate_intrinsic_GeV_prior": 0,
            "real_lat_geV_events_used_in_joint_toy_likelihood": 0,
            "population_intrinsic_lag_prior": None,
            "observed_test_statistic": None,
            "observed_p_value": None,
            "observed_liv_confidence_interval": None,
            "conditional_toy_receipt_sha256": files["source_lag_toy_receipt"]["sha256"]},
        "scientific_nonpromotion_gates": {
            **obs["scientific_gates"],
            "lorentz_invariance_violation": "NOT_DETECTED",
            "spacetime_discreteness": "NOT_DETECTED",
            "amps_quantum_firewall": "THEORETICAL_HYPOTHESIS_UNCONFIRMED_NO_OBSERVATION",
            "black_hole_boom": "NOT_DETECTED_NO_VALIDATED_FIREWALL_MATTER_TO_EXTERIOR_RADIATION_MODEL",
            "finite_energy_gamma_burst": "NOT_A_DIRECT_BLACK_HOLE_HORIZON_EXPERIMENT",
            "actual_joint_astrophysical_p_value": None,
            "physical_model_for_horizon_firewall_collision": None,
            "ordinary_calibration_and_emission_lag_explanations_excluded": False},
        "next_scientific_frontier_v6": {
            "goal": "Genuine predeclared observational joint LAT/GBM time-energy likelihood with population intrinsic lag", 
            "required_lat": ["type-resolved effective exposure including pointing and livetime", "calibrated PSF/EDISP kernel evaluation", "time-dependent diffuse and point source background", "GRB photon association and ROI selection function"],
            "required_gbm": ["time-variable background and detector optimality", "RSP2 variable-length matrix decode and time-specific spectral forward folding", "TTE-to-CSPEC/RSP2 channel-edge reconciliation under uncertainty"],
            "required_physics": ["independent GeV intrinsic emission lag constraints", "hierarchical burst spectral/timing evolution", "instrument cross-calibration systematic uncertainties"],
            "required_statistics": ["null Monte Carlo with nuisance profiling", "injection/recovery and interval coverage", "trials budget and preregistered decision threshold"],
            "must_not_claim_detection_until_complete": True},
        "verification": {
            "cosmos_repository": "Hawkar-usls/Janus-Cosmos",
            "cosmos_pull_request": args.cosmos_pr,
            "cosmos_pr_url": f"https://github.com/Hawkar-usls/Janus-Cosmos/pull/{args.cosmos_pr}",
            "cosmos_merged_main_commit": args.cosmos_merge_sha,
            "scientific_ci_workflows_passed": args.ci_passed,
            "offline_unittests_passed": args.tests_passed,
            "all_tests_failed": 0,
            "source_data_replay_status": "COMPLETED_LOCAL_WITH_SOURCE_PINNED_EXTERNAL_FITS",
            "offline_gate": "PASS_87_TESTS_AND_FROZEN_RECEIPTS",
            "source_tree_sha256": files,
            "registry_ingest_requires_integrity_and_scientific_nonpromotion": True},
        "claim_ceiling": "REAL_GBM_TTE_CSPEC_RSP2_AND_LAT_CALDB_TABLE_PROVENANCE_AND_LIGHTCURVES_ONLY_NO_OBSERVED_LIV_OR_FIREWALL_DETECTION",
        "science_gates": {
            "REAL_GBM_TTE_CSPEC_RSP2": "PASS_SOURCE_SHA256_VERIFIED",
            "REAL_LAT_P8R3_SOURCE_V3_CALDB_TABLES": "PASS_SOURCE_SHA256_GRID_SHAPE_ONLY",
            "JOINT_LAT_GBM_FORWARD_FOLDED_LRT": "BLOCKED_NO_VALIDATED_SOURCE_LAG_MODEL",
            "SOURCE_LAG_IDENTIFIABILITY": "EXACT_DEGENERACY_UNBOUNDED_BETA_PER_BURST",
            "LIV_OBSERVED": "NOT_ESTABLISHED",
            "SPACETIME_DISCRETENESS": "NOT_ESTABLISHED",
            "AMPS_FIREWALL": "OPEN_UNCONFIRMED",
            "BHB_BOOM": "NOT_ESTABLISHED"},
        "unresolved_scientific_requirements": [
            "full event-type-specific LAT PSF/EDISP and time-dependent exposure",
            "validated GBM detector response-channel edge reconciliation and spectral fit",
            "GRB population intrinsic emission lag constraints independent of propagation",
            "real LAT/GBM jointly calibrated likelihood and null-test coverage",
            "physically specified quantum firewall scattering/emission mechanism and external observable"],
        "epistemic_guard": "A calibrated instrument input is necessary but insufficient for a calibrated physics likelihood. No observational speed-of-light deviation, granular spacetime, firewall radiation, or BHB explosion is established."}
    registry["payload_sha256"] = digest(registry)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(registry, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "sha256": sha(args.output.read_bytes()),
                      "manifested_external_fits": len(manifest["lat_caldb"]["caldb_files"]) + sum(len(b["detectors"]) * 3 for b in manifest["gbm_bursts"])}))


if __name__ == "__main__":
    main()
