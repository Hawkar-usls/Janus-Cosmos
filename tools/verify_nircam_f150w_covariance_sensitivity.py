#!/usr/bin/env python3
"""Reproduce the analytic (NOT injection-calibrated) NIRCam sensitivity grid.

Only reads committed summary receipts. No original FITS, Stage 2 JWST
spectroscopy, empirical PSF recovery, cross-epoch covariance, or source
classification is computed here.
"""
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / "domains" / "space"
x = json.loads((D / "NIRCAM_F150W_COVARIANCE_SENSITIVITY_2026-10-18.json").read_text())
old = json.loads((D / "NIRCAM_F150W_ASTROMETRY_RECEIPT_2026-10-16.json").read_text())

assert x["schema"] == "janus-cosmos-nircam-f150w-analytic-covariance-sensitivity-v1"
assert x["status"] == "SCENARIO_SENSITIVITY_ONLY_NOT_CALIBRATED_INJECTION_RECOVERY"
assert len(old["epochs"]) == 4
epochs = old["epochs"]
assert len({r["product"] for r in epochs}) == 4
assert [r["date"][:10] for r in epochs] == x["source_four_dates_utc"]
ys = [r["target_flux_aperture_njy"] for r in epochs]
sigmas = [r["one_sigma_controls_njy"] for r in epochs]
assert ys == x["measured_aperture_njy"]
assert sigmas == x["empirical_one_sigma_njy"]
assert (old["source_ra"], old["source_dec"]) == (x["source_ra_deg"], x["source_dec_deg"])
assert all(math.isfinite(a) and b > 0 for a, b in zip(ys, sigmas))

def near(actual, expected, tol=1e-8):
    assert math.isclose(actual, expected, rel_tol=0, abs_tol=tol), (actual, expected)

rawweights = [1 / (s * s) for s in sigmas]
weights = [v / sum(rawweights) for v in rawweights]
mu = sum(w * y for w, y in zip(weights, ys))
v_ind = sum((w * s) ** 2 for w, s in zip(weights, sigmas))
v_max = sum(w * s for w, s in zip(weights, sigmas)) ** 2
line_sigma = (x["published_line_flux_sigma_cgs"] /
              x["published_line_flux_cgs"]) * x["published_line_f150w_equivalent_total_njy"]
near(mu, x["weighted_mean_njy"])
near(math.sqrt(v_ind), x["nominal_independent_mean_sigma_njy"])
near(math.sqrt(v_max), x["mean_sigma_if_all_pairwise_rho_one_njy"])
near(line_sigma, x["propagated_spectral_to_f150w_sigma_njy"])
for calculated, recorded in zip(weights, x["weights_nominal_independent"]):
    near(calculated, recorded)

assert len(x["scenarios"]) >= 6
ids = set()
for scene in x["scenarios"]:
    assert scene["id"] not in ids
    ids.add(scene["id"])
    f = scene["fraction"]
    rho = scene["rho"]
    systematic = scene["additive_systematic_njy"]
    assert 0 < f <= 1 and 0 <= rho <= 1 and systematic >= 0
    expected = f * x["published_line_f150w_equivalent_total_njy"]
    variance = v_ind + rho * (v_max - v_ind) + (f * line_sigma) ** 2 + systematic ** 2
    error = math.sqrt(variance)
    near(scene["predicted_aperture_njy"], expected)
    near(scene["combined_error_njy"], error)
    near(scene["conditional_residual_over_combined_error"], (expected - mu) / error)

guardrails = set(x["distinguish"])
assert {"PSF_INJECTION_NOT_RUN", "EPOCH_COVARIANCE_NOT_ESTIMATED",
        "SOURCE_CLASS_UNRESOLVED", "REDSHIFT_UNRESOLVED",
        "DIRECTIONAL_TRAIL_NOT_DETECTED"}.issubset(guardrails)
assert any("not measured" in s or "not estimated" in s for s in x["assumptions"])

print("PASS: four independent dated receipts reproduce analytic conditional grid")
print("F150W_APERTURE_MEAN_NJY={:.6f}".format(mu))
print("F150W_INDEPENDENT_SIGMA_NJY={:.6f}".format(math.sqrt(v_ind)))
print("F150W_COVARIANCE_AND_PSF_INJECTION=NOT_MEASURED")
print("INTRINSIC_VARIABILITY=UNCONFIRMED")
print("SOURCE_CLASS_AND_DISTANCE=UNRESOLVED")
