#!/usr/bin/env python3
"""V5 source-lag sensitivity experiment: purely synthetic, observationally gated.

Real LAT/GBM counts are NOT supplied to the Gaussian toy likelihood. The V5
instrument audit constrains which nuisance components remain missing; it is
not legitimate to turn toy significance into an observed likelihood ratio.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from janus_cosmos_fermi_joint_v4 import fit_grid, toy_bursts
from janus_cosmos_fermi_calibration_v5 import SCHEMA as INSTRUMENT_SCHEMA
from janus_cosmos_simulation_audit import AuditError, digest

SCHEMA = "janus.cosmos.fermi-multiburst-source-lag-sensitivity.v5"


def run(instrument: dict, seed: int = 20261019) -> dict:
    if instrument.get("schema") != INSTRUMENT_SCHEMA or instrument.get("counts", {}).get("gbm_tte_files") != 9:
        raise AuditError("V5 pinned real instrument audit is required")
    if instrument["scientific_gates"]["joint_observational_liv_hypothesis_test"] != "BLOCKED_NO_VALIDATED_LIKELIHOOD":
        raise AuditError("real-data physics gate must stay closed")
    if [b["grb"] for b in instrument["gbm_bursts"]] != ["090510", "090902B", "130427A"]:
        raise AuditError("wrong V5 GRB cohort")
    # Synthetic trials are deliberately independent of the *measured* GBM
    # light curves, which do not constrain GeV intrinsic emission lags.
    signal = toy_bursts(seed, eta_true=0.8)
    mimic = toy_bursts(seed, eta_true=0.0, source_beta=0.8)
    null = toy_bursts(seed, eta_true=0.0, source_beta=0.0)
    families = []
    for sigma in (0.1, 0.5, 2.0):
        families.append({"hypothetical_intrinsic_beta_prior_sigma": sigma,
                         "synthetic_liv_injection": fit_grid(signal, "gaussian_prior", sigma),
                         "source_only_mimic": fit_grid(mimic, "gaussian_prior", sigma),
                         "zero_effect_null": fit_grid(null, "gaussian_prior", sigma),
                         "prior_is_not_derived_from_gbm_tte": True})
    free = fit_grid(signal, "unbounded")
    free_mimic = fit_grid(mimic, "unbounded")
    if free["eta_mle"] is not None or free_mimic["twice_delta_log_likelihood"] != 0:
        raise AuditError("unbounded per-burst source lag must preserve exact degeneracy")
    out = {"schema": SCHEMA, "experiment_id": instrument["experiment_id"],
           "data_kind": "THREE_BURST_SYNTHETIC_TIMING_PRIOR_SENSITIVITY_CONDITIONED_ON_REAL_INSTRUMENT_READINESS_ONLY",
           "seed": seed,
           "actual_gbm_tte_count_used_for_source_lag_prior_calibration": 0,
           "actual_lat_event_count_used_in_likelihood": 0,
           "real_instrument_data_are_in_separate_v5_audit": True,
           "source_lag_parameterization": "t_ji = mu_j + (eta + beta_j) K_j(E,z) + normal_noise",
           "exact_reparameterization": "eta_prime=eta+delta; beta_j_prime=beta_j-delta; predicted_times_invariant",
           "source_lag_prior_sensitivity_synthetic_only": families,
           "unconstrained_source_lag_true_identifiability": free,
           "source_only_mimic_unconstrained": free_mimic,
           "synthetic_zero_source_prior_false_positive_demo": fit_grid(mimic, "known_zero"),
           "independent_high_energy_intrinsic_emission_prior": None,
           "observational_liv_test_statistic": None,
           "observational_liv_p_value": None,
           "observational_liv_limit": None,
           "gates": {"SYNTHETIC_SOURCE_LAG_PRIOR_SENSITIVITY": "PASS_PRIOR_DEPENDENT_ONLY",
                     "SYNTHETIC_EXACT_NONIDENTIFIABILITY": "PASS",
                     "OBSERVED_LIV_TEST": "BLOCKED_NO_VALIDATED_RESPONSE_AND_SOURCE_PRIOR",
                     "AMPS_FIREWALL": "OPEN_UNCONFIRMED", "BHB_BOOM": "NOT_ESTABLISHED"},
           "claim_ceiling": "source lag priors change synthetic fits; real GBM pulse timings do not identify GeV intrinsic emission delays"}
    out["result_sha256"] = digest(out)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--instrument-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20261019)
    a = parser.parse_args()
    audit = json.loads(a.instrument_audit.read_text(encoding="utf-8"))
    result = run(audit, a.seed)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, sort_keys=True, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(a.output), "free_source_lag": result["unconstrained_source_lag_true_identifiability"]}))


if __name__ == "__main__":
    main()
