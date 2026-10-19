#!/usr/bin/env python3
"""Strict audit of actual 16 archived RAW RATE source-window stacks and 24
parametric source injections into native 2D pixels, with sky recomputation.

Green means a reproducible engineering test, NEVER physical STPSF flux.
"""
import json, math, re
from pathlib import Path
d=Path(__file__).resolve().parents[1]/"domains"/"space"
x=json.loads((d/"MOCSS_FULL16_NATIVE_DETECTOR_2D_GAUSSIAN_INJECTION_2026-10-19.json").read_text())
ref=json.loads((d/"MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json").read_text())
prev=json.loads((d/"MOCSS_16_REAL_RATE_FOUR_STACKS_ROW_SKY_PROXY_2026-10-19.json").read_text())
assert x["status"].startswith("ALL_16_ORIGINAL_RATE_2_MSA_SHA_PASS")
assert x["inputs"]["original_rate_count"]==len(ref["exposures"])==16
assert x["inputs"]["original_MSA_count"]==2
assert x["inputs"]["hashes_verified"] is True
assert x["provenance"]["run_id"]==39231058856
assert re.fullmatch("[0-9a-f]{64}",x["provenance"]["full_json_sha256"])
assert x["provenance"]["stack_arrays_npz_sha256"]==prev["provenance"]["wcs_stack_npz_sha256"]
assert x["native_detector_trial"]["injection_done_into_original_2D_RATE_before_full_row_local_sky_recompute"]
assert not x["native_detector_trial"]["source_model_physical_STPSF"]
assert x["native_detector_trial"]["sky_row_reestimation_absolute_change_DN_per_s"]==0
count=0
for group in "ABCD":
    vals=x["per_group_response_ratio_by_wavelength"][group]
    assert len(x["inputs"]["group_four_dithers"][group])==4
    assert set(vals)=={"primary","OIII_4960","Pa_gamma"}
    for ratio in vals.values():
        assert len(ratio)==2
        assert .9<ratio[0]<1.1
        assert 1.1<ratio[1]<1.3
        count+=len(ratio)
assert count==24
limits=x["science_ceiling"]
assert limits["median_per_full_detector_row_4mad_with_mask_is_proxy_not_shutter_illumination"] is True
for field in ("full_open_closed_MSA_to_detector_pixel_mapping_done",
              "source_psf_model_matches_calibrated_observation",
              "stage2_calibrated_four_dither_physical_flux",
              "independent_Gaussian_pvalue_from_native_trials"):
    assert limits[field] is False
for field in ("absolute_flux_cgs","source_class","redshift","directed_gas_trail"):
    assert limits[field] is None
print("PASS: native 2D original rate pixel injection, 24 source-window recovery controls, no physical calibration")
