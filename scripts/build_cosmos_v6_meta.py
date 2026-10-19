#!/usr/bin/env python3
"""Generate long-form V6-R0 scientific meta-registry receipt from proven sources."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = {
    "gbm_rsp2_engine": "workspace/janus_cosmos_fermi_response_v6.py",
    "lat_sc_gtis_engine": "workspace/janus_cosmos_fermi_lat_diagnostics_v6.py",
    "joint_readiness_gate": "workspace/janus_cosmos_fermi_joint_v6.py",
    "official_rehydration_script": "scripts/fetch_fermi_v6_inputs.py",
    "meta_registry_builder": "scripts/build_cosmos_v6_meta.py",
    "gbm_real_forwardfold_result": "data/cosmos/COSMOS-ALI-V6-GBM-RSP2-ONOFF-SPECTRAL-2026-10-19.json",
    "lat_real_roi_exposure_result": "data/cosmos/COSMOS-ALI-V6-LAT-SOURCE-CONTROL-EXPOSURE-2026-10-19.json",
    "cross_instrument_readiness_result": "data/cosmos/COSMOS-ALI-V6-JOINT-OBSERVATION-READINESS-2026-10-19.json",
    "v6_unit_tests": "tests/test_cosmos_fermi_v6.py",
    "v6_scientific_guide": "docs/JANUS-COSMOS-LAT-GBM-JOINT-V6-R0-2026-10-19.md",
    "frozen_v5_manifest": "data/cosmos/COSMOS-ALI-V5-FERMI-LAT-GBM-CALDB-MANIFEST-2026-10-19.json",
    "frozen_v5_instrument_result": "data/cosmos/COSMOS-ALI-V5-CALDB-GBM-LAT-INSTRUMENT-PULSE-2026-10-19.json",
    "frozen_v5_source_lag_toy": "data/cosmos/COSMOS-ALI-V5-SOURCE-LAG-TOY-SENSITIVITY-2026-10-19.json",
    "frozen_v4_manifest": "data/cosmos/COSMOS-ALI-V4-FERMI-RAW-FITS-MANIFEST-2026-10-18.json",
    "frozen_v4_observation": "data/cosmos/COSMOS-ALI-V4-FERMI-RAW-EVENT-DIAGNOSTICS-2026-10-18.json",
}


def sha(x: bytes) -> str:
    return hashlib.sha256(x).hexdigest()


def canonical_digest(x: dict) -> str:
    return sha(json.dumps(x, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cosmos-merge-sha", required=True)
    ap.add_argument("--cosmos-pr", required=True, type=int)
    ap.add_argument("--ci-passed", required=True, type=int)
    ap.add_argument("--tests-passed", required=True, type=int)
    ap.add_argument("--output", required=True, type=Path)
    a = ap.parse_args()
    if len(a.cosmos_merge_sha) != 40 or a.ci_passed < 1 or a.tests_passed < 99:
        raise ValueError("V6 record requires merged Cosmos commit, passing CI and >=99 unit tests")
    files = {name: {"path": path, "sha256": sha((ROOT / path).read_bytes()),
                    "github_url": "https://github.com/Hawkar-usls/Janus-Cosmos/blob/main/" + path}
             for name, path in FILES.items()}
    gbm = json.loads((ROOT / FILES["gbm_real_forwardfold_result"]).read_text())
    lat = json.loads((ROOT / FILES["lat_real_roi_exposure_result"]).read_text())
    joint = json.loads((ROOT / FILES["cross_instrument_readiness_result"]).read_text())
    m5 = json.loads((ROOT / FILES["frozen_v5_manifest"]).read_text())
    m4 = json.loads((ROOT / FILES["frozen_v4_manifest"]).read_text())
    if m5["v4_frozen_manifest_sha256"] != files["frozen_v4_manifest"]["sha256"] or m5["v4_frozen_observed_json_sha256"] != files["frozen_v4_observation"]["sha256"]:
        raise ValueError("V4/V5 byte lineage has changed")
    for value in (gbm, lat, joint):
        copied = dict(value)
        claimed = copied.pop("result_sha256")
        if canonical_digest(copied) != claimed:
            raise ValueError("V6 result content digest invalid")
    if gbm["observational_joint_p_value"] is not None or joint["observational_liv_p_value"] is not None or lat["real_joint_observed_liv_p_value"] is not None:
        raise ValueError("unsafe observed LIV p-value promotion")
    if gbm["scientific_gates"]["OBSERVED_LAT_GBM_LIV_HYPOTHESIS_TEST"] != "BLOCKED" or joint["scientific_gates"]["OBSERVED_MULTIBURST_LIV_LRT"] != "BLOCKED_NO_FREQUENTIST_COVERAGE_OR_VALID_SOURCE_MODEL":
        raise ValueError("V6 source/instrument science gate not locked")
    if lat["counts"]["v4_lat_source_region_events_recovered_exactly"] != 175:
        raise ValueError("frozen LAT selected events changed")
    if sum(len(x["time_resolved_spectral_fit_segments"]) for x in gbm["gbm_bursts"]) != 17:
        raise ValueError("V6 time resolved spectral windows unexpectedly changed")
    six_lats = []
    for burst in m4["bursts"]:
        for role in ("event", "spacecraft"):
            meta = burst[role]
            six_lats.append({"grb": burst["grb"], "role": role, "official_url": meta["url"],
                             "sha256": meta["sha256"], "bytes": meta["bytes"]})
    receipt = {
        "artifact_uuid": "JANUS-COSMOS-ALI-FERMI-LAT-GBM-REAL-RESPONSE-FORWARDFOLD-V6R0-2026-10-19",
        "registry_schema_id": "janus.registry.physics.fermi_lat_gbm_joint_response_readiness.v6r0",
        "version": "0.1.0", "date_utc": "2026-10-19", "owner": "Hawkar-usls",
        "title": "JANUS Cosmos V6-R0 — real 17-window time-weighted GBM RSP2 forward folding, exact Poisson on/off diagnostics, LAT pointing and control annulus, source-lag nonpromotion gate",
        "experiment_id": gbm["experiment_id"],
        "classification": "OBSERVED_GBM_SPECTRAL_COUNT_PILOT_AND_LAT_IRF_EXPOSURE_PROXY_NOT_FULL_JOINT_LIV_TEST",
        "scientific_status": "REAL_RESPONSE_AND_CONTROL_DIAGNOSTICS_PASS_OBSERVED_LIV_LRT_BLOCKED",
        "research_lineage": {
            "youtube_original_context": "https://www.youtube.com/watch?v=QBpuorUxuSM",
            "video_is_not_scientific_evidence": True,
            "v5_cosmos_merged_commit": "861cf471ff1b63eee80b4bc0d86b4c5bf0251e88",
            "v5_registry_merge_commit": "9c24979b29836dea07537994ddb5b8b16972a2d0",
            "v5_cosmos_pr": "https://github.com/Hawkar-usls/Janus-Cosmos/pull/41",
            "v5_registry_pr": "https://github.com/Hawkar-usls/janus-meta-registry/pull/262",
            "historical_v4_v5_payload_sha256": {name: files[name]["sha256"] for name in ("frozen_v5_manifest", "frozen_v5_instrument_result", "frozen_v5_source_lag_toy", "frozen_v4_manifest", "frozen_v4_observation")},
            "frozen_history_unchanged": True,
            "osiris_authority": "UNMODIFIED"},
        "archive_observational_assets": {
            "provider": "NASA_FERMI_FSSC_GBM_HEASARC_AND_FERMITOOLS",
            "original_v5_manifest_containing_36_byte_pinned_irf_gbm_sources": m5,
            "six_original_v4_lat_event_sc_fits": six_lats,
            "total_rehydrated_upstream_files": 42,
            "source_sha256_integrity": "PASS_ALL_42_EXTERNALLY_REHYDRATED_USING_FROZEN_MANIFESTS",
            "v5_frozen_receipt_byte_replay_sha256": "53c4e8bacb55c64914ae263c0295aad53055d342cdac8fa95bb46e7974d38424",
            "heavy_raw_fits_not_committed_to_github": True,
            "data_kind": "REAL_BYTES_NOT_INTERPOLATED_SCIENTIFIC_SOURCE_CLAIMS"},
        "actual_gbm_time_resolved_forward_folding": {
            "onoff_count_model": "Poisson(n_on|s+b)*Poisson(m_off|tau*b), per-PHA-channel background profiled",
            "source_expected_counts": "duration*amplitude*sum_true_energy [DRM_cm2*cutoff_powerlaw_flux_integral]",
            "rsp2_response_matrix": "OGIP_PI(1/128)_PE(128)_VARIABLE_LENGTH_HEAP_DECODED",
            "ogip_F_CHAN_indexing": "ONE_BASED_RESPONSE_TO_ZERO_BASED_TTE_PHA",
            "response_time_evolution": "DURATION_OVERLAP_WEIGHTED_ACROSS_ATTITUDE_SPECIFIC_RSP2_DRMS",
            "number_of_spectral_windows": 17,
            "number_of_detectors_per_grb": 3,
            "source_spectrum": "SINGLE_CUTOFF_POWERLAW_COARSE_GRID_NOT_BAND_PLUS_EXTRA_POWERLAW",
            "source_spectral_alpha_grid": [-1.6, -1.2, -0.8, -0.4],
            "source_spectral_epeak_grid_keV": [150, 300, 600, 1200, 2400, 4800],
            "background_method": "PREBURST_-20_TO_-5_SEC_SHARED_ONOFF_NOT_VALIDATED_TIME_VARIABLE_BACKGROUND",
            "shape_residual_stress_descriptive_only": "deviance_per_fitted_channel >5 indicates likely severe single-component model mismatch",
            "significant_model_mismatch_intervals": {
                "GRB090902B": [[5, 10]], "GRB130427A": [[0, 10], [10, 20]], "GRB090510": []},
            "all_real_burst_spectral_intervals_including_detector_onoff_and_rsp2_weights": gbm["gbm_bursts"],
            "real_data_not_a_sufficient_physical_source_emission_model": True,
            "standalone_GBM_receipt_sha256": files["gbm_real_forwardfold_result"]["sha256"]},
        "actual_lat_exposure_and_control": {
            "source_roi": "5deg_SOURCE_event_class128_E>1GeV",
            "control_region": "5to10deg_ANNULUS",
            "full_v4_lat_candidate_replay": "EXACT_175_REAL_EVENTS",
            "prompt_candidates_per_burst": {"090510": 19, "090902B": 42, "130427A": 47},
            "prompt_control_annulus_count_per_burst": {"090510": 0, "090902B": 0, "130427A": 1},
            "exposure_proxy": "REAL_GTIS_SC_LIVETIME_NEAREST_FB_AEFF_AT_1_10GeV",
            "exposure_missing": ["PSF", "EDISP", "azimuth", "full livetime efficiency response", "diffuse background", "source probability"],
            "full_three_grb_control_windows_and_lat_exposure_proxies": lat["GRB_sources"],
            "standalone_LAT_receipt_sha256": files["lat_real_roi_exposure_result"]["sha256"]},
        "actual_cross_instrument_status": {
            "paired_GBM_LAT_three_bursts": joint["real_multiwavelength_diagnostic_pairs"],
            "synthetic_only_intrinsic_mimic_negative_control": joint["conditional_source_mimic_negative_control_SYNTHETIC_ONLY"],
            "GRB_GeV_emission_time_assumption": "NOT_ESTABLISHED_FROM_KEV_PEAK",
            "exact_intrinsic_source_propagation_degeneracy": "eta_plus_beta_j_is_only_identified_by_toy_if_intrinsic_beta_j_unconstrained",
            "observed_joint_liv_estimate": None,
            "observed_joint_liv_test_statistic": None,
            "observed_joint_liv_p_value": None,
            "observed_lattice_spacing_upper_bound": None,
            "standalone_cross_instrument_receipt_sha256": files["cross_instrument_readiness_result"]["sha256"]},
        "science_gates": {**joint["scientific_gates"],
            "LAT_CALDB_VERSION": "P8R3_SOURCE_V3_SHA256_VERIFIED",
            "REAL_GBM_FORWARD_FOLD": "PILOT_NOT_CALIBRATED_ASTROPHYSICAL_LRT",
            "SPECTRAL_MODEL_COMPLETE": "NOT_COMPLETE_HIGH_RESIDUALS_090902B_130427A",
            "INTRINSIC_EMISSION_LAG_PHYSICALLY_CONSTRAINED": "NO",
            "LIV_OBSERVED": "NOT_ESTABLISHED",
            "SPACETIME_DISCRETENESS": "NOT_ESTABLISHED",
            "AMPS_FIREWALL": "OPEN_UNCONFIRMED",
            "BHB_BOOM": "NOT_ESTABLISHED"},
        "next_gate_v6_r1": {
            "goal": "First validated full physical joint instrument+source likelihood, still blinded to any discovery claim",
            "required_gbm": ["Band plus hard extra-power-law component", "time-dependent background model and offsource holdout", "GBM channel-edge mismatch calibration nuisance and cross-detector normalization"],
            "required_lat": ["calibrated Pass8 SOURCE PSF/EDISP/event-type effective exposure", "true diffuse/point-source likelihood model", "photon association and source selection function"],
            "required_physics": ["independent source emission lag prior from temporal and spectral GRB evolution", "hierarchical population delay model with variation"],
            "required_validation": ["independent official-tool comparison", "null injections and profile coverage", "calibrated frequentist model checks and preregistered thresholds"],
            "allowed_observed_p_value_before_validation": None},
        "claim_ceiling": "REAL_GBM_OGIP_RESPONSE_COUNT_FOLD_AND_LAT_CONTROL_EXPOSURE_PROXY_ONLY_NO_FULL_LAT_GBM_LIV_STATISTICS_OR_BLACK_HOLE_HORIZON_DETECTION",
        "validation": {
            "cosmos_repository": "Hawkar-usls/Janus-Cosmos",
            "cosmos_merged_main_commit": a.cosmos_merge_sha,
            "cosmos_pull_request": a.cosmos_pr,
            "cosmos_pr_url": f"https://github.com/Hawkar-usls/Janus-Cosmos/pull/{a.cosmos_pr}",
            "all_unittest_passed": a.tests_passed,
            "scientific_ci_workflows_passed": a.ci_passed,
            "fail_closed_negative_tests": True,
            "source_tree_files_sha256": files,
            "observational_liv_p_value": None},
        "epistemic_guard": "RSP2 physical response folding and LAT control exposure diagnostics are real advances, but neither physical source lag nor complete LAT/GBM response has been validated. A true observational LIV p-value would be scientifically unjustified."}
    receipt["payload_sha256"] = canonical_digest(receipt)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(receipt, ensure_ascii=False, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"output": str(a.output), "sha256": sha(a.output.read_bytes()),
                      "science_status": receipt["scientific_status"]}))


if __name__ == "__main__":
    main()
