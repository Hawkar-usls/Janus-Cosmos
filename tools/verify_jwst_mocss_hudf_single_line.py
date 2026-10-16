#!/usr/bin/env python3
"""Finite, source-bound consistency checks for the JWST MOCSS COSMOS case.

This validates arithmetic and anti-inference gates, not MAST sensor data or
an astrophysical source identity. Run beside cosmos_case.json.
"""
import datetime as dt
import json
from pathlib import Path


def run():
    # In Janus-Cosmos this verifier lives in tools/ and the source-bound
    # dataset lives in domains/space/. The fallback is for an isolated replay.
    repo_root = Path(__file__).resolve().parents[1]
    p = repo_root / "domains" / "space" / "JWST_MOCSS_HUDF_SINGLE_LINE_COSMOS_2026-10-16.json"
    if not p.exists():
        p = Path(__file__).with_name("cosmos_case.json")
    case = json.loads(p.read_text(encoding="utf-8"))
    obs = case["derived_checks"]["lambda_observed_angstrom"]
    for line in case["derived_checks"]["candidate_lines"]:
        z = obs / line["rest_vacuum_angstrom"] - 1.0
        assert abs(z - line["expected_redshift_rounded"]) < 0.011, (line["id"], z)
        # Persistent-source rest-EW lower bounds are rounded values in Table 4.
        rest_ew = 7620 / (1 + z)
        assert abs(rest_ew - line["persistent_source_rest_ew_lower_angstrom"]) < 15, line["id"]
    first = dt.date.fromisoformat(case["measured_reported"]["spectroscopy"]["first_epoch"][0])
    second = dt.date.fromisoformat(case["measured_reported"]["spectroscopy"]["second_epoch"][0])
    assert (second - first).days == 88
    factor = case["measured_reported"]["line_flux_erg_s_cm2"] / case["derived_checks"]["imaging_nondetection_equivalent_line_flux_erg_s_cm2"]
    assert abs(factor - case["derived_checks"]["cross_epoch_flux_vs_imaging_limit_factor"]) < 0.01
    assert case["gates"]["PHYSICAL_SOURCE_ID"] == "UNRESOLVED"
    assert case["gates"]["RAW_JWST_REPLAY"] == "INSUFFICIENT_DATA"
    assert case["training_policy"]["gradients_or_model_weights_changed"] is False
    assert case["world_truth"] is False
    print("PASS: four redshift/EW checks, 88-day baseline, flux factor, and authority gates")


if __name__ == "__main__":
    run()
