#!/usr/bin/env python3
"""Immutable diagnostic four-group x1d receipt guard. Never promotes sky flux."""
from pathlib import Path
import json,re
p=Path(__file__).resolve().parents[1]/"domains"/"space"
x=json.loads((p/"MOCSS_FOUR_CORRECT_ROWS_DIAGNOSTIC_X1D_2026-10-19.json").read_text())
g=json.loads((p/"MOCSS_FOUR_GROUP_CORRECT_ROW_AXIS_GEOMETRY_2026-10-19.json").read_text())
w=json.loads((p/"MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json").read_text())
assert x["state"]=="FOUR_REAL_ARCHIVAL_DIAGNOSTIC_X1D_CORRECT_THREE_ROW_MSA_GEOMETRY_NOT_SCIENCE_GRADE"
assert [q["group"] for q in x["groups"]]==list("ABCD")
assert x["provenance"]["run_id"]==39223778736
for i,entry in enumerate(x["groups"]):
    corresponding=g["groups"][i]
    idx=entry["exp"]
    item=w["exposures"][idx-1]
    assert entry["group"]==corresponding["group"]==item["group"]
    assert entry["shutter_rows"]==corresponding["rows"]
    assert entry["column"]==corresponding["column"]
    assert item["detectors"][item["expected_science_detector_from_WCS"]]["selected_virtual_shutter"]["row"]==entry["shutter_rows"][1]
    assert re.fullmatch("[a-f0-9]{64}",entry["fits_sha256"])
    assert len(entry["diagnostic_formal_ratio"])==3
    assert len(entry["diagnostic_chi2_per_dof"])==3
assert x["groups"][2]["diagnostic_chi2_per_dof"][0]>1000
assert x["all_formal_ratios_diagnostic_not_detection_significance"]
for k in ("full_16_exposure_4_dither_stacks_reproduced","author_row_sky_subtraction_reproduced",
          "three_shutter_real_optical_throughput_calibrated","physical_psf_detector_injection_recovery_run"):
    assert x[k] is False
for k in ("spectral_source_fluxes_or_upper_limits","companion_OIII_detected","companion_Pagamma_detected",
          "physical_source_type","redshift","directional_gas_trail"):
    assert x[k] is None
print("PASS: four SHA-pinned real diagnostic source-aware x1d with 3 dispersion rows; science flux OPEN")
