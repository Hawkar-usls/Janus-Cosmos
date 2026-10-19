#!/usr/bin/env python3
"""Check original archive provenance of 16 RATE/2 MSA/4 sky-proxy WCS stacks.

Strictly forbid interpreting median-row estimator as illuminated MSA sky model,
nor 2D WCS-aligned patch injection as original native detector injection.
"""
import json,re
from pathlib import Path
d=Path(__file__).resolve().parents[1]/"domains"/"space"
r=json.loads((d/"MOCSS_16_REAL_RATE_FOUR_STACKS_ROW_SKY_PROXY_2026-10-19.json").read_text())
o=json.loads((d/"MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json").read_text())
assert r["status"].startswith("16_ORIGINAL_RATE_2_ORIGINAL_MSA_SHA_VERIFIED")
assert r["provenance"]["number_original_RATE"]==len(o["exposures"])==16
assert r["provenance"]["number_MSA_reference"]==2
assert r["provenance"]["original_SHA_verified"] is True
assert re.fullmatch(r"[0-9a-f]{64}",r["provenance"]["original_json_sha256"])
assert re.fullmatch(r"[0-9a-f]{64}",r["provenance"]["wcs_stack_npz_sha256"])
for group,indices,det in [('A',[1,2,3,4],'NRS2'),('B',[5,6,7,8],'NRS2'),
                          ('C',[9,10,11,12],'NRS1'),('D',[13,14,15,16],'NRS1')]:
    line=r["groups"][group]
    assert line["indices"]==indices and line["detector"]==det
    assert line["all_three_center_coverage"]==4
    assert set(line["line_center_relative_robust_units_NOT_SNR"])=={'primary','OIII_4960','Pa_gamma'}
    for i in indices:
        assert o["exposures"][i-1]["group"]==group
g=r["experiment_limits"]
assert g["valid_original_MSA_open_closed_shutter_config_counted"] is True
for key in ("physical_open_closed_row_to_detector_sky_calibration_completed",
            "authors_full_row_illumination_corrected_sky_reproduced",
            "full_calibrated_2D_spectrum_from_stacked_images",
            "original_RATE_native_2D_detector_injection_in_this_artifact",
            "independent_detector_covariance_measured","stpsf_dispersed_detector_PSF_used"):
    assert g[key] is False
for key in ("line_flux_absolute_cgs","line_identification","redshift","gas_trail"):
    assert g[key] is None
print("PASS: archival 16 RATE 2 MSA SHA and 12 source windows; robust sky estimator ONLY")
print("ORIGINAL_MSA_ROW_TO_DETECTOR_SKY=NOT_CALIBRATED")
