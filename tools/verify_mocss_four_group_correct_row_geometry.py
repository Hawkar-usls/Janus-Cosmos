#!/usr/bin/env python3
"""Reproducible full four-anchor MSA physical-axis and SHA checkpoint gate."""
import json,math,re
from pathlib import Path
base=Path(__file__).resolve().parents[1]/"domains"/"space"
x=json.loads((base/"MOCSS_FOUR_GROUP_CORRECT_ROW_AXIS_GEOMETRY_2026-10-19.json").read_text())
w=json.loads((base/"MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json").read_text())
assert x["status"].startswith("FOUR_ORIGINAL_ARCHIVAL_RATE_MSA_INPUTS_SHA_PASS")
assert x["rows_dispersion_direction"] and x["columns_crossdispersion"]
assert x["workflow_run_id"]==39223778736
assert not x["source_in_original_MPT"] and not x["original_archive_modified"]
assert len(x["groups"])==4
for group in x["groups"]:
    idx=group["exposure_index"]
    src=w["exposures"][idx-1]
    det=src["expected_science_detector_from_WCS"]
    p=src["detectors"][det]
    assert src["group"]==group["group"] and group["detector"]==det
    assert p["sha256"]==group["original_rate_sha256"]
    assert p["selected_virtual_shutter"]["row"]==group["source_row"]
    assert p["selected_virtual_shutter"]["column"]==group["column"]
    assert p["selected_virtual_shutter"]["quadrant"]==group["quadrant"]
    assert group["rows"]==[group["source_row"]-1,group["source_row"],group["source_row"]+1]
    for key in ("original_msa_sha256","test_msa_sha256"):
        assert re.fullmatch("[0-9a-f]{64}",group[key])
    assert group["original_msa_sha256"]!=group["test_msa_sha256"]
for h in x["source_json_sha256"].values():assert re.fullmatch("[0-9a-f]{64}",h)
for key in ("source_shutter_PSF_throughput_verified","full_16_exposure_four_stacks_done",
            "custom_full_row_background_science_valid","physical_injection_recovery_done"):
    assert x[key] is False
for key in ("absolute_flux_16004_cgs","companion_OIII_15850_cgs","companion_Pagamma_13656_cgs","source_class","redshift","directed_trail"):
    assert x[key] is None
assert x["past_three_column_x1d_flux_assumptions"]=="SUPERSEDED_AND_INVALID_FOR_DISPERSION_AXIS"
print("PASS: 4 original rate and MSA SHA, 3 consecutive dispersive rows, no invalid physical source claims")
