#!/usr/bin/env python3
"""V6 cross-instrument observed diagnostics and explicit joint-inference lock.

The only products combined are *real* per-source GBM source pulses/spectral
diagnostics and LAT ROI/spacecraft summaries; no unvalidated observed LIV
test statistic is computed. Intrinsic GeV release times remain unidentifiable.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from janus_cosmos_fermi_joint_v4 import fit_grid, toy_bursts
from janus_cosmos_simulation_audit import AuditError, ROOT, digest

SCHEMA = "janus.cosmos.fermi-lat-gbm-cross-instrument-readiness.v6"
GBM = ROOT / "data/cosmos/COSMOS-ALI-V6-GBM-RSP2-ONOFF-SPECTRAL-2026-10-19.json"
LAT = ROOT / "data/cosmos/COSMOS-ALI-V6-LAT-SOURCE-CONTROL-EXPOSURE-2026-10-19.json"
V5 = ROOT / "data/cosmos/COSMOS-ALI-V5-CALDB-GBM-LAT-INSTRUMENT-PULSE-2026-10-19.json"
V4 = ROOT / "data/cosmos/COSMOS-ALI-V4-FERMI-RAW-EVENT-DIAGNOSTICS-2026-10-18.json"
OUT = ROOT / "data/cosmos/COSMOS-ALI-V6-JOINT-OBSERVATION-READINESS-2026-10-19.json"


def run(gbm: dict, lat: dict, prior: dict, raw: dict) -> dict:
    if not gbm.get("schema", "").endswith("spectral-pilot.v6") or not lat.get("schema", "").endswith("controls.v6"):
        raise AuditError("V6 observed readiness requires both validated pilots")
    if len(gbm["gbm_bursts"]) != 3 or len(lat["GRB_sources"]) != 3:
        raise AuditError("wrong V6 burst count")
    if lat["counts"]["v4_lat_source_region_events_recovered_exactly"] != 175:
        raise AuditError("V4 source region had unexpected photon count")
    if gbm["scientific_gates"]["OBSERVED_LAT_GBM_LIV_HYPOTHESIS_TEST"] != "BLOCKED":
        raise AuditError("spectral diagnostic was improperly promoted")
    if prior["scientific_gates"]["joint_observational_liv_hypothesis_test"] != "BLOCKED_NO_VALIDATED_LIKELIHOOD":
        raise AuditError("frozen V5 observed science gate weakened")
    if raw["counts"]["diagnostic_source_region_candidates"] != 175:
        raise AuditError("frozen V4 source set modified")
    pairs = []
    for b, l, prev, r in zip(gbm["gbm_bursts"], lat["GRB_sources"], prior["gbm_bursts"], raw["sources"]):
        if len({b["grb"], l["grb"], prev["grb"], r["grb"]}) != 1:
            raise AuditError("GRB pairing mismatch")
        prompt = next(w for w in l["three_time_windows"] if w["time_window"] == "prompt")
        soft = next(d["pulse_bands"]["soft"] for d in prev["detectors"] if "soft" in d["pulse_bands"])
        off = {"GBM_NaI_soft_band_peak_relative_trigger_s": soft["peak_time_relative_trigger_s"],
               "LAT_prompt_ROI_5deg_candidates": prompt["source_circle_5deg_count"],
               "LAT_prompt_ring_5to10deg_controls": prompt["control_annulus_5to10deg_count"],
               "LAT_prompt_exposure_proxy_at_1GeV_cm2_s": prompt["nearest_FB_aeff_times_livetime_cm2_s_by_energy_MeV"]["1000"]}
        mismatches = [s for s in b["time_resolved_spectral_fit_segments"] if s["onoff_spectral_grid_fit"]["spectral_mismatch_stress_gate"] == "HIGH_RESIDUAL_MODEL_REVISE"]
        pairs.append({"grb": b["grb"], "redshift": l["redshift"], "prompt_window_s": prompt["relative_to_exact_gbm_trigger_s"],
                      **off,
                      "GBM_RSP2_time_resolved_spectral_segments": len(b["time_resolved_spectral_fit_segments"]),
                      "GBM_CPL_high_residual_segments": len(mismatches),
                      "GBM_CPL_high_residual_window_bounds": [[s["start_relative_trigger_s"], s["stop_relative_trigger_s"]] for s in mismatches],
                      "GBM_peak_not_GeV_emission_start": True,
                      "LAT_ROI_candidates_not_event_source_posteriors": True,
                      "lat_geV_minus_gbm_keV_intrinsic_emission_offset_constraint": None})
    toy = toy_bursts(seed=20261019, eta_true=0., source_beta=0.8)
    mimicked = fit_grid(toy, "known_zero")
    unbounded = fit_grid(toy, "unbounded")
    if unbounded["eta_mle"] is not None or unbounded["twice_delta_log_likelihood"] != 0.:
        raise AuditError("source-lag degeneracy must remain exact")
    out = {"schema": SCHEMA, "experiment_id": gbm["experiment_id"],
           "real_data_kind": "JOINT_READY_OBSERVATIONAL_INSTRUMENT_DIAGNOSTICS_NOT_JOINT_LIV_LIKELIHOOD",
           "real_GBM_forward_folding": "OGIP_RSP2_WITH_POISSON_ONOFF_CPL_GRID_PILOT",
           "real_LAT_forward_folding": "SC_LIVETIME_NEAREST_AEFF_EXPOSURE_PROXY_NO_PSF_OR_EDISP",
           "source_physics": "GEV_EMISSION_TIME_UNRECONSTRUCTED_GBM_SOFT_PEAK_NOT_PROPAGATION_BASELINE",
           "real_multiwavelength_diagnostic_pairs": pairs,
           "conditional_source_mimic_negative_control_SYNTHETIC_ONLY": {"injected_liv_eta": 0.,
                                                    "intrinsic_beta": 0.8,
                                                    "naive_eta_fit": mimicked,
                                                    "unconstrained_eta": unbounded},
           "unreconciled_response_systematics": ["LAT PSF and energy dispersion not included in likelihood",
               "GBM response-channel edge mismatches not propagated as nuisance", "GBM time varying background and additional source spectral components missing"],
           "conditional_hypothetical_GeV_lag_prior": None,
           "observational_liv_test_statistic": None,
           "observational_liv_p_value": None,
           "observational_lattice_spacing_confidence_interval": None,
           "firewall_emission_observation": None,
           "scientific_gates": {"ACTUAL_LAT_GBM_COHORT_AND_INSTRUMENT_PROVENANCE": "PASS",
                                "REAL_GBM_RSP2_FORWARD_FOLD": "PILOT_CPL_ONOFF_NOT_A_VALIDATED_GLOBAL_SOURCE_MODEL",
                                "REAL_LAT_FULL_IRF_LIKELIHOOD": "BLOCKED_PSF_EDISP_EXPOSURE_BACKGROUND",
                                "PHYSICAL_GE_V_SOURCE_EMISSION_PRIOR": "BLOCKED_NO_INDEPENDENT_INFORMATION",
                                "OBSERVED_MULTIBURST_LIV_LRT": "BLOCKED_NO_FREQUENTIST_COVERAGE_OR_VALID_SOURCE_MODEL",
                                "LORENTZ_VIOLATION": "NOT_ESTABLISHED",
                                "SPACETIME_DISCRETENESS": "NOT_ESTABLISHED",
                                "AMPS_FIREWALL": "OPEN_UNCONFIRMED",
                                "BHB_BOOM": "NOT_ESTABLISHED"},
           "next_gate": "V6_R1_BAND_PLUS_PL_BG_POLYNOMIAL_LAT_FULL_PSF_EDISP_FORWARD_MODEL_INJECTION_COVERAGE",
           "claim_ceiling": "NO_OBSERVED_LIV_P_VALUE_AND_NO_BHB_FIREWALL_INFERENCE"}
    out["result_sha256"] = digest(out)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gbm-result", type=Path, default=GBM)
    ap.add_argument("--lat-result", type=Path, default=LAT)
    ap.add_argument("--output", type=Path, default=OUT)
    a = ap.parse_args()
    result = run(json.loads(a.gbm_result.read_text()), json.loads(a.lat_result.read_text()),
                 json.loads(V5.read_text()), json.loads(V4.read_text()))
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"output": str(a.output), "pair_count": len(result["real_multiwavelength_diagnostic_pairs"]),
                      "observational_p_value": result["observational_liv_p_value"]}))


if __name__ == "__main__":
    main()
