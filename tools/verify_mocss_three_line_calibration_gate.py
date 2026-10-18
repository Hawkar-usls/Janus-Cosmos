#!/usr/bin/env python3
"""Keep the calibrated three-window line search falsifiable and non-authoritative."""
import json
import math
from pathlib import Path

D=Path(__file__).resolve().parents[1]/"domains"/"space"
p=json.loads((D/"MOCSS_THREE_LINE_CALIBRATION_GATE_2026-10-18.json").read_text())
assert p["schema"]=="janus-cosmos-mocss-calibrated-three-line-gate-v1"
assert p["science_status"]=="PRE_REGISTERED_REQUIRED_CALIBRATION_NO_NEW_LINE_FLUX_MEASURED"
assert p["known_spectroscopy"]["n_rate_input_detector_frames"]==32
assert p["known_spectroscopy"]["n_science_spectral_exposures"]==16
assert p["calibration_requirements"]["pipeline_author_method_skipped"]==["master_background","pathloss","barshadow"]
assert p["calibration_requirements"]["custom_full_row_background_for_open_msa_REQUIRED"]
assert p["calibration_requirements"]["synthetic_line_injection_recovery_per_group_and_wavelength_required"]
lines=p["three_line_hypotheses"]
assert abs(lines["oiii_4960"]["observed_lambda_A"]-15850.7932988)<0.001
assert abs(lines["pagamma"]["observed_lambda_A"]-13656.7483309)<0.001
assert abs(lines["oiii_4960"]["conditional_companion_flux_cgs"]-2.5e-18/3)<1e-28
assert abs(lines["pagamma"]["conditional_pagamma_unreddened_flux_cgs"]-2.5e-18/1.6)<1e-28
assert lines["oiii_4960"]["article_tested_incorrect_wavelength_A"]==15864
for key in ("stage2_16004_independent_photometry",
            "stage2_15850p793_companion_measurement",
            "stage2_13656p748_pagamma_measurement",
            "mocss_throughput_injection"):
    assert p["gate_states"][key]=="NOT_PERFORMED"
assert p["gate_states"]["line_ratio_redshift"]=="OPEN"
assert p["gate_states"]["directed_gas_trail"]=="NOT_DETECTED"
for v in p["known_spectroscopy"]["author_reported_per_group"].values():
    assert v["flux_1e18_cgs"]>0 and v["one_sigma_1e18_cgs"]>0
assert p["known_spectroscopy"]["author_fluxes_not_independently_remeasured"]
print("PASS: corrected 3-line positions, predicted ratios, 16-exposure and strict not-performed gates")
print("STAGE2_NIRSPEC_INDEPENDENT_SPECTRAL_FLUX=NOT_MEASURED")
