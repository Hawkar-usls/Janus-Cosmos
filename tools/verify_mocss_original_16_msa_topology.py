#!/usr/bin/env python3
"""Guard actual 16 archival MSA quadrant rasters against overstated science claims."""
import json,hashlib,re
from pathlib import Path
d=Path(__file__).resolve().parents[1]/"domains"/"space"
x=json.loads((d/"MOCSS_ORIGINAL_MSA_16_EXPOSURE_TOPOLOGY_2026-10-19.json").read_text())
p=json.loads((d/"MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json").read_text())
assert x["science_status"]=="ORIGINAL_FITS_MSA_APERTURE_TOPOLOGY_VERIFIED_ALL16_DETECTOR_ILLUMINATION_NOT_YET_PROJECTED"
assert x["upstream"]["workflow_id"]==39237389277
assert len(x["per_exposure"])==len(p["exposures"])==16
assert x["original_shutter_mask"]["total_quadrant_masks"]==64
for field in ("raw_json_sha256","raw_rasters_npz_sha256"):
    assert re.fullmatch("[0-9a-f]{64}",x["upstream"][field])
assert x["original_shutter_mask"]["near_source_center_open_all16"]
assert x["original_shutter_mask"]["near_source_dispersion_neighbor_rows_all48"]
assert .90<x["original_shutter_mask"]["global_open_fraction_min"]<x["original_shutter_mask"]["global_open_fraction_max"]<1.
assert 0.9<x["original_shutter_mask"]["shutter_row_near_source_open_fraction_min"]<=x["original_shutter_mask"]["shutter_row_near_source_open_fraction_max"]<1
for rec,reference in zip(x["per_exposure"],p["exposures"]):
    slit=reference["detectors"][reference["expected_science_detector_from_WCS"]]["selected_virtual_shutter"]
    assert rec["index"]==reference["exposure_index"]
    assert rec["group"]==reference["group"]
    assert rec["original_MSAMETID_header"]==reference["source_msametid_unwrapped"]
    assert rec["PATT_NUM_header"]==1
    assert rec["nearest_MSA_quadrant"]==slit["quadrant"]
    assert rec["nearest_MSA_row"]==slit["row"]
    assert rec["nearest_MSA_column"]==slit["column"]
    assert rec["nearest_center_open"] and rec["three_consecutive_MSA_rows_at_fixed_column_open"]
    assert 0<rec["nearest_row_open_count_of_171"]<=171
    assert rec["original_msa_entries_after_both_filters"]>200000
g=x["calibration_gates"]
assert g["physical_MSA_plane_binary_open_closed_mask"].startswith("PASS")
assert g["correct_original_MSA_metadata_plus_dither_headers"]=="PASS_16"
for key in ("wavelength_dependent_full_detector_pixel_illumination",
            "closed_shutter_leakage_model",
            "full_2D_dispersed_four_stack_mocss_calibration",
            "physically_dispersed_STPSF_source_injection"):
    assert g[key]=="NOT_PERFORMED"
assert g["independent_absolute_line_flux_A_to_D"]=="NOT_MEASURED"
assert g["source_identity"]==g["redshift"]=="OPEN"
print("PASS: original MSA 16 FITS headers, 64 quadrant masks and source-adjacent 48/48 OPEN; detector illumination remains UNSOLVED")
