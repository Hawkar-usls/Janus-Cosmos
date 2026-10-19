#!/usr/bin/env python3
"""Generate provenance-rich V6-R1 *partial gate* meta-registry JSON.

Requires the real JANUS Cosmos merge SHA and green CI before emission.
The builder never upgrades a descriptive instrument-model improvement to
observational LIV significance. Raw Fermi binaries stay outside git.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_FILES = {
    "gbm_spectral_engine": "workspace/janus_cosmos_fermi_spectral_v6r1.py",
    "lat_irf_audit_engine": "workspace/janus_cosmos_fermi_lat_irf_v6r1.py",
    "intrinsic_lag_negative_control_engine": "workspace/janus_cosmos_fermi_source_lag_v6r1.py",
    "gbm_real_data_receipt": "data/cosmos/COSMOS-ALI-V6R1-GBM-BAND-POWERLAW-DRIFT-2026-10-19.json",
    "lat_real_irf_receipt": "data/cosmos/COSMOS-ALI-V6R1-LAT-FULL-IRF-READINESS-2026-10-19.json",
    "intrinsic_lag_synthetic_control_receipt": "data/cosmos/COSMOS-ALI-V6R1-SOURCE-LAG-TOY-INJECTION-COVERAGE-2026-10-19.json",
    "r1_offline_tests": "tests/test_cosmos_fermi_v6r1.py",
    "r1_scientific_protocol": "docs/JANUS-COSMOS-V6-R1-BAND-IRF-SOURCE-LAG-GATE-2026-10-19.md",
    "r1_meta_builder": "scripts/build_cosmos_v6r1_meta.py",
    "v5_frozen_manifest": "data/cosmos/COSMOS-ALI-V5-FERMI-LAT-GBM-CALDB-MANIFEST-2026-10-19.json",
    "v4_frozen_manifest": "data/cosmos/COSMOS-ALI-V4-FERMI-RAW-FITS-MANIFEST-2026-10-18.json",
    "v6r0_joint": "data/cosmos/COSMOS-ALI-V6-JOINT-OBSERVATION-READINESS-2026-10-19.json",
}


def sha_bytes(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def canonical_digest(value: dict) -> str:
    return sha_bytes(json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode())


def read_and_validate(relpath: str) -> dict:
    data = json.loads((ROOT / relpath).read_text())
    without = dict(data)
    expected = without.pop("result_sha256")
    if canonical_digest(without) != expected:
        raise ValueError(f"frozen result digest changed: {relpath}")
    return data


def build(cosmos_merge_sha: str, cosmos_pr: int, tests_passed: int, ci_passed: int) -> dict:
    if len(cosmos_merge_sha) != 40 or any(c not in "0123456789abcdef" for c in cosmos_merge_sha):
        raise ValueError("must provide actual merged JANUS Cosmos commit")
    if tests_passed < 112 or ci_passed < 1 or cosmos_pr <= 42:
        raise ValueError("V6-R1 requires CI+112 tests and successor Cosmos PR")
    files = {name: {"path": path, "sha256": sha_bytes((ROOT / path).read_bytes()),
                    "github_url": "https://github.com/Hawkar-usls/Janus-Cosmos/blob/main/" + path}
             for name, path in SOURCE_FILES.items()}
    gbm = read_and_validate(SOURCE_FILES["gbm_real_data_receipt"])
    lat = read_and_validate(SOURCE_FILES["lat_real_irf_receipt"])
    toy = read_and_validate(SOURCE_FILES["intrinsic_lag_synthetic_control_receipt"])
    v6 = read_and_validate(SOURCE_FILES["v6r0_joint"])
    m5 = json.loads((ROOT / SOURCE_FILES["v5_frozen_manifest"]).read_text())
    m4 = json.loads((ROOT / SOURCE_FILES["v4_frozen_manifest"]).read_text())
    if m5["v4_frozen_manifest_sha256"] != files["v4_frozen_manifest"]["sha256"]:
        raise ValueError("frozen V4/V5 lineage drifted")
    if not gbm["observed_liv_p_value"] is lat["observed_liv_p_value"] is toy["observed_multiburst_LIV_p_value"] is v6["observational_liv_p_value"] is None:
        raise ValueError("observed p-values must remain null in incomplete scientific gate")
    if sum(len(b["spectral_segments"]) for b in gbm["gbm_bursts"]) != 17 or len(lat["response_families"]) != 9:
        raise ValueError("GBM/IRF fixed cohort mismatch")
    if toy["scenarios"][2]["naive_95pct_interval_empirical_coverage"] > .1:
        raise ValueError("synthetic intrinsic-lag false-positive control failed")
    def segment(grb, lo, hi):
        matches = [s for b in gbm["gbm_bursts"] if b["grb"] == grb
                   for s in b["spectral_segments"] if s["relative_trigger_interval_s"] == [lo, hi]]
        if len(matches) != 1:
            raise ValueError("predeclared GRB spectral segment not found")
        return matches[0]
    s090902 = segment("090902B", 20., 30.)
    s130427a = segment("130427A", 20., 40.)
    s130427a_early = segment("130427A", 0., 10.)
    r1 = {
        "artifact_uuid": "JANUS-COSMOS-ALI-FERMI-LAT-GBM-BAND-IRF-LAG-COVERAGE-V6R1-2026-10-19",
        "registry_schema_id": "janus.registry.physics.fermi_lat_gbm_source_lag_v6r1_partial_gate",
        "version": "0.1.0",
        "date_utc": "2026-10-19",
        "owner": "Hawkar-usls",
        "title": "JANUS Cosmos V6-R1 real GBM Band+PL response folding and independent BGO holdout, full LAT PSF/EDISP parameter audit, source-lag identifiability and conditional synthetic coverage",
        "scientific_status": "PARTIAL_V6_R1_GBM_MULTI_COMPONENT_AND_DRIFT_PILOT_LAT_IRF_PARAMETER_AUDIT_FULL_JOINT_LRT_BLOCKED",
        "claim_ceiling": "NO_CALIBRATED_LAT_SOURCE_PROBABILITY_NO_VALIDATED_GBM_BACKGROUND_SYSTEMATICS_NO_SOURCE_LAG_PRIOR_NO_OBSERVED_LIV_OR_FIREWALL_OR_BHB_DETECTION",
        "frozen_science_lineage": {
            "v6r0_cosmos_merge_sha": "d41537f8c74c5e8b281da72059753bf105edacc9",
            "v6r0_cosmos_pr": "https://github.com/Hawkar-usls/Janus-Cosmos/pull/42",
            "v6r0_registry_merge_sha": "e41d5f634a1b1f729fc6345e9caa1b21c00385b0",
            "v6r0_registry_pr": "https://github.com/Hawkar-usls/janus-meta-registry/pull/267",
            "v6r0_joint_result_sha256": files["v6r0_joint"]["sha256"],
            "current_v6r1_cosmos_merge_sha": cosmos_merge_sha,
            "current_v6r1_cosmos_pr": "https://github.com/Hawkar-usls/Janus-Cosmos/pull/" + str(cosmos_pr),
            "V4_V5_frozen_manifest_byte_hashes_unchanged": True,
            "external_video_heuristic_not_authoritative_observational_source": "https://www.youtube.com/watch?v=QBpuorUxuSM",
        },
        "source_and_method_reference_urls": [
            "https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm/gbm_data_tools/gdt-docs/notebooks/SpectralFitting.html",
            "https://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_LAT_IRFs/IRF_PSF.html",
            "https://fermi.gsfc.nasa.gov/ssc/data/analysis/lat_irfs/irf_energy_dispersion.html",
            "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/BackgroundModels.html",
            "https://arxiv.org/abs/0909.2470",
            "https://doi.org/10.1103/PhysRevD.87.122001",
        ],
        "external_real_data_provenance": {
            "provider": "NASA_FERMI_FSSC_GBM_HEASARC_AND_FERMITOOLS_P8R3_SOURCE_V3",
            "v5_source_manifest": m5,
            "v4_lat_event_sc_fits": [
                {"grb": b["grb"], "role": key, "filename": b[key]["filename"],
                 "source_url": b[key]["url"], "source_sha256": b[key]["sha256"]}
                for b in m4["bursts"] for key in ("event", "spacecraft")],
            "number_verified_external_raw_files": 42,
            "raw_archive_bytes_committed_to_repo": False,
            "v5_manifest_sha256": files["v5_frozen_manifest"]["sha256"],
        },
        "GBM_real_count_domain_multi_component": {
            "data_kind": "REAL_OBSERVED_COUNTS_AND_REAL_OGIP_RSP2_RESPONSE",
            "spectrum": "Band plus q times power law normalized to Band photon flux at 1 MeV",
            "source_model_shape_parameters": {"alpha": [-1.4, -1., -.6], "beta": [-3., -2.4],
                                             "epeak_keV": [150., 400., 1000., 2500.],
                                             "extra_powerlaw_index": -1.8,
                                             "extra_q_at_1MeV": [0., .5, 2., 8.]},
            "preregistered_template_grid_nodes": 96,
            "response_integrator": "SIMPSON4_SUBINTERVAL_TRUE_ENERGY_X_EXACT_OGIP_RSP2_ATTITUDE_WEIGHTED",
            "source_normalization": "ONE_PER_SEGMENT_SHARED_TWO_NAI_ONE_BGO",
            "background_train": "TWO_DISJOINT_PRETRIGGER_WINDOWS_PER_DETECTOR_FIT_ONE_COMMON_LOG_LINEAR_SLOPE",
            "background_holdout": "THIRD_DISJOINT_PRETRIGGER_WINDOW_NOT_USED_TO_FIT_SLOPE",
            "detector_holdout": "TWO_NAI_SPECTRAL_TRAIN_BGO_SEPARATE_ON_SOURCE_POISSON_DEVIANCE",
            "real_data_17_segments_and_all_diagnostics": gbm["gbm_bursts"],
            "likelihood_limitations": gbm["background_and_fit_important_limitations"],
            "three_detector_band_vs_hard_component_significance": None,
            "real_burst_090902B_20to30s_descriptive_delta_deviance": s090902["three_detector_fit"]["descriptive_deviance_improvement_not_lrt"],
            "independent_090902B_BGO_deviance_band_vs_combo": [s090902["nai_train_bgo_holdout"][k] for k in ("bgo_holdout_band_deviance", "bgo_holdout_band_plus_powerlaw_deviance")],
            "independent_130427A_20to40s_BGO_deviance_band_vs_combo": [s130427a["nai_train_bgo_holdout"][k] for k in ("bgo_holdout_band_deviance", "bgo_holdout_band_plus_powerlaw_deviance")],
            "independent_130427A_0to10s_BGO_HUGE_MODEL_MISMATCH": s130427a_early["nai_train_bgo_holdout"]["bgo_holdout_band_deviance"],
            "observed_model_selection_p_value": None,
        },
        "LAT_actual_CALDB_response_audit": {
            "instrument_response_release": "P8R3_SOURCE_V3",
            "all_nine_irf_families_checked": True,
            "individual_event_type_subtables_checked": 30,
            "full_real_response_family_coefficients_and_bounds": lat["response_families"],
            "native_psf_moffat_scale_extracted": True,
            "native_energy_dispersion_mixture_scale_extracted": True,
            "azimuth_and_livetime_efficiency_parameters_audited_not_applied": True,
            "normalized_PSF_EDISP_event_likelihood": None,
            "galactic_isotropic_diffuse_model_likelihood": None,
            "full_LAT_likelihood_blockers": lat["cannot_promote_to_calibrated_lat_forward_model_without"],
        },
        "physical_intrinsic_lag_assumption_gate": {
            "exact_nonidentifiability_theorem": toy["theorem"],
            "toy_weight_units": "DIMENSIONLESS_NOT_PHYSICAL_COSMOLOGY",
            "toy_random_seed": toy["seed"],
            "toy_injections_and_coverage_2000_each": toy["scenarios"],
            "unconstrained_lag_profile_infinite_eta_family": toy["unconstrained_source_delay_exact_profile"],
            "naive_interval_coverage_with_correct_beta_null": toy["scenarios"][0]["naive_95pct_interval_empirical_coverage"],
            "naive_interval_coverage_when_beta_mimics_eta": toy["scenarios"][2]["naive_95pct_interval_empirical_coverage"],
            "externally_justified_source_GeV_emission_lag_prior": None,
            "conditional_toy_is_observed_LAT_GBM_result": False,
        },
        "validation": {
            "offline_total_unittest_passed": tests_passed,
            "cosmos_GitHub_Action_workflows_green": ci_passed,
            "source_result_sha256_verified": True,
            "source_sha256_map_and_links": files,
            "observed_joint_LIV_p_value": None,
            "observed_LIV_eta_estimate": None,
            "observed_lattice_spacing_limit": None,
            "black_hole_firewall_or_horizon_detection": False,
        },
        "science_gates": {
            "REAL_GBM_RSP2_BAND_PLUS_PL_17_SEGMENTS": "PASS_EXPERIMENTAL_MODEL_DIAGNOSTIC",
            "GBM_TIME_DEPENDENT_OFFSOURCE_MODEL": "PILOT_LOG_LINEAR_DRIFT_UNCERTAINTY_UNCALIBRATED",
            "INDEPENDENT_PRETRIGGER_BACKGROUND_HOLDOUT": "PASS_REAL_COUNT_DIAGNOSTIC",
            "INDEPENDENT_BGO_CROSS_DETECTOR_CHECK": "PASS_REAL_COUNT_DIAGNOSTIC_MIXED",
            "GBM_CROSS_DETECTOR_CALIBRATION": "BLOCKED_CHANNEL_EDGES_AND_NORMALIZATION",
            "ACTUAL_LAT_P8R3_SOURCE_V3_PSF_EDISP_TABLES": "PASS_REAL_PARAMETER_AUDIT",
            "LAT_FULL_PSF_EDISP_EFFICIENCY_DIFFUSE_SOURCE_LIKELIHOOD": "BLOCKED_NOT_EVALUATED",
            "PHYSICALLY_JUSTIFIED_SOURCE_INTRINSIC_LAG_PRIOR": "BLOCKED",
            "SYNTHETIC_LAG_INJECTION_AND_GAUSSIAN_COVERAGE": "PASS_ONLY_TOY",
            "END_TO_END_REAL_FERMI_POISSON_INJECTION_RECOVERY_COVERAGE": "NOT_DONE",
            "OBSERVED_MULTIBURST_LIV_P_VALUE": "BLOCKED_NULL",
            "LORENTZ_VIOLATION": "NOT_ESTABLISHED",
            "DISCRETE_SPACETIME": "NOT_ESTABLISHED",
            "AMPS_FIREWALL": "OPEN_UNCONFIRMED",
            "BHB_BOOM": "NOT_ESTABLISHED",
        },
        "next_gate_v6_r2_or_remaining_v6_r1": {
            "GBM": ["background extrapolation validation", "joint shared-off nuisance likelihood", "RSP2/TTE channel edge crosscalibration", "detector-specific normalization systematics", "source-model posterior predictive stress test"],
            "LAT": ["official normalized PSF and energy dispersion convolution", "phi/efficiency/livetime dependence", "diffuse models", "source catalog ROI", "Fermitools benchmark"],
            "SOURCE": ["independent conservative intrinsic lag constraints", "hierarchical burst spectral evolution", "provenance and sensitivity to source model"],
            "STATISTICS": ["end-to-end null+signal instrument injections", "coverage and false-positive calibration", "precommitted acceptance thresholds", "only afterward observed LRT and p-value"],
            "allowed_observed_p_value_now": None,
        },
    }
    r1["payload_sha256"] = canonical_digest(r1)
    return r1


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cosmos-merge-sha", required=True)
    p.add_argument("--cosmos-pr", type=int, required=True)
    p.add_argument("--tests-passed", type=int, required=True)
    p.add_argument("--ci-passed", type=int, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    result = build(a.cosmos_merge_sha, a.cosmos_pr, a.tests_passed, a.ci_passed)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n")
    print(json.dumps({"output": str(a.output), "scientific_status": result["scientific_status"],
                      "file_count": len(SOURCE_FILES), "observed_p_value": result["validation"]["observed_joint_LIV_p_value"]}))


if __name__ == "__main__":
    main()
