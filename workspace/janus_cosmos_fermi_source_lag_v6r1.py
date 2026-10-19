#!/usr/bin/env python3
"""V6-R1 *synthetic-only* source-delay identifiability/injection controls.

No synthetic lag data or Gaussian coverage check is a Fermi-LAT observation.
An unconstrained burst-specific intrinsic delay beta_j makes the universal
propagation parameter eta non-identifiable: y_j=K_j*eta+beta_j+noise.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path

from janus_cosmos_simulation_audit import AuditError, ROOT, digest

SCHEMA = "janus.cosmos.fermi-synthetic-source-lag-degeneracy-coverage-controls.v6r1"
OUT = ROOT / "data/cosmos/COSMOS-ALI-V6R1-SOURCE-LAG-TOY-INJECTION-COVERAGE-2026-10-19.json"
SEED = 20261019
N = 2000


def naive_eta_ci(y: tuple[float, ...], k: tuple[float, ...], sigma: tuple[float, ...]) -> tuple[float, float, float]:
    """Correct ONLY conditional on externally established beta_j == 0."""
    if not (len(y) == len(k) == len(sigma) >= 2 and all(x > 0 for x in sigma)):
        raise AuditError("invalid synthetic delay array")
    precision = sum(q * q / (s * s) for q, s in zip(k, sigma))
    if precision <= 0:
        raise AuditError("no synthetic propagation leverage")
    eta_hat = sum(q * z / (s * s) for z, q, s in zip(y, k, sigma)) / precision
    error = 1. / math.sqrt(precision)
    return eta_hat, eta_hat - 1.95996398454 * error, eta_hat + 1.95996398454 * error


def exact_lag_profile_degeneracy(k: tuple[float, ...], y: tuple[float, ...], sigma: tuple[float, ...], eta: float) -> dict:
    """For every eta, beta_j=y_j-K_j*eta gives identically zero residual."""
    if len(k) != len(y) or len(k) != len(sigma):
        raise AuditError("wrong lag dimensions")
    beta = tuple(z - q * eta for z, q in zip(y, k))
    chisq = sum(((z - q * eta - b) / s) ** 2 for z, q, b, s in zip(y, k, beta, sigma))
    return {"eta_candidate": eta, "nuisance_intrinsic_beta_by_burst": list(beta),
            "profile_chisquare": round(chisq, 12)}


def simulate(seed: int = SEED, n: int = N) -> dict:
    if not (1 <= n <= 100000):
        raise AuditError("invalid synthetic sample count")
    rng = random.Random(seed)
    burst_names = ("090510", "090902B", "130427A")
    # Dimensionless *illustrative* propagation lever arms, not a cosmology
    # integral or observed burst delay. Their scales define an arbitrary toy.
    k = (.8, 1.5, .3)
    sigma = (.2, .2, .3)
    scenarios = (("null_no_intrinsic_delay", 0., (0., 0., 0.)),
                 ("injected_eta_without_intrinsic_delay", .5, (0., 0., 0.)),
                 ("adversarial_intrinsic_delay_mimics_propagation", 0., (.64, 1.2, .24)))
    reports = []
    for name, true_eta, beta in scenarios:
        coverage = 0
        false_rejection_null = 0
        bias_sum = 0.
        for _ in range(n):
            y = tuple(true_eta * a + b + rng.gauss(0., s) for a, b, s in zip(k, beta, sigma))
            estimate, low, high = naive_eta_ci(y, k, sigma)
            coverage += low <= true_eta <= high
            false_rejection_null += not (low <= 0. <= high)
            bias_sum += estimate - true_eta
        reports.append({"scenario": name, "synthetic_true_eta": true_eta,
                        "synthetic_intrinsic_beta": list(beta), "pseudo_experiments": n,
                        "naive_95pct_interval_empirical_coverage": round(coverage / n, 6),
                        "naive_null_eta0_excluded_fraction": round(false_rejection_null / n, 6),
                        "mean_naive_estimator_bias": round(bias_sum / n, 6),
                        "naive_model_asserts_beta_zero": True,
                        "not_an_observed_fermi_significance": True})
    witness = (1.2, -.3, .8)
    profiles = [exact_lag_profile_degeneracy(k, witness, sigma, eta) for eta in (-1., 0., .5, 5.)]
    if not all(p["profile_chisquare"] == 0 for p in profiles):
        raise AuditError("synthetic intrinsic delay degeneracy was not exact")
    result = {"schema": SCHEMA, "seed": seed, "pseudo_experiments_each_scenario": n,
              "synthetic_dimensionless_weights_NOT_cosmology": dict(zip(burst_names, k)),
              "synthetic_observation_errors": dict(zip(burst_names, sigma)),
              "scenarios": reports,
              "unconstrained_source_delay_exact_profile": profiles,
              "theorem": "For independent beta_j per burst, replace beta_j by beta_j+(eta_old-eta_new)*K_j; all predicted delays are unchanged, for any eta_new.",
              "intrinsic_GeV_lag_prior_physically_justified": False,
              "valid_real_LAT_GBM_coverage_calibration": False,
              "injection_recovery_status": "SYNTHETIC_GAUSSIAN_TOY_ONLY_NOT_INSTRUMENT_PLUS_SOURCE_COVERAGE",
              "observed_multiburst_LIV_p_value": None,
              "observed_eta_confidence_interval": None,
              "scientific_gates": {"TOY_BETA_ETA_EXACT_DEGENERACY": "PASS_ALGEBRAIC",
                                   "TOY_INJECTION_COVERAGE": "PASS_SEE_CONDITIONAL_RESULTS",
                                   "PHYSICAL_INTRINSIC_GE_V_SOURCE_DELAY_PRIOR": "BLOCKED",
                                   "FULL_FERMI_INJECTION_RECOVERY_COVERAGE": "NOT_DONE",
                                   "OBSERVED_LIV_TEST": "BLOCKED",
                                   "BHB_FIREWALL": "OPEN_UNCONFIRMED",
                                   "BHB_BOOM": "NOT_ESTABLISHED"},
              "claim_ceiling": "SYNTHETIC_NEGATIVE_CONTROL_ONLY_NO_OBSERVATIONAL_LIV_CONSTRAINT"}
    result["result_sha256"] = digest(result)
    return result


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=OUT)
    args = p.parse_args()
    result = simulate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"output": str(args.output), "coverage": [s["naive_95pct_interval_empirical_coverage"] for s in result["scenarios"]],
                      "observed_p": result["observed_multiburst_LIV_p_value"]}))


if __name__ == "__main__":
    main()
