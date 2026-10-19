#!/usr/bin/env python3
"""Independent frozen original-open-shutter CRDS 4-anchor WCS pass only.

This gate validates archived numeric detector positions against the original
MSA topology and forbids interpreting finite slit WCS as photon flux.
"""
import json,re,math
from pathlib import Path
root=Path(__file__).resolve().parents[1]/"domains"/"space"
x=json.loads((root/"MOCSS_CRDS_GWCS_FOUR_ANCHOR_OPEN_SHUTTER_FOOTPRINT_2026-10-19.json").read_text())
e=json.loads((root/"MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json").read_text())
g=json.loads((root/"MOCSS_CRDS_WCS_FORWARD_ILLUMINATION_GATE_2026-10-19.json").read_text())
assert x["science_status"]=="12_ORIGINAL_OPEN_SHUTTER_CRDS_GWCS_36_FINITE_ON_DETECTOR_SAMPLE_GEOMETRIES_ONLY"
p=x["provenance"]
assert p["num_original_exposures"]==4 and p["num_actual_MSA_open_shutters"]==12
assert p["num_three_wavelength_center_tests"]==36 and p["action_run_id"]==39243207739
assert re.fullmatch("[0-9a-f]{64}",p["artifact_json_sha256"])
assert x["results"]["official_assign_wcs"]=="COMPLETE_ALL4"
assert x["results"]["finite_three_lambda_center_projection"]=="PASS_36"
assert x["results"]["all_center_coordinates_on_detector"]=="PASS_36"
for group,idx in zip("ABCD",[1,8,10,13]):
    m=x["sampled_source_nearest_shutter_detector_xy_pix"][group]
    a=e["exposures"][idx-1]
    d=a["expected_science_detector_from_WCS"]
    v=a["detectors"][d]["selected_virtual_shutter"]
    assert m["exposure"]==idx and m["detector"]==d.upper()
    assert (m["quadrant"],m["msa_row"],m["msa_col"])==(v["quadrant"],v["row"],v["column"])
    for key in x["wavelength_A"]:
        u,vv=m[key]
        assert all(math.isfinite(z) and 0<=z<2048 for z in (u,vv))
assert x["method"]["source_original_MPT_identification"] is False
assert g["real_CRDS_GWCS_pilot"]["official_WCS_detector_footprints_measured"].startswith("PASS_12")
assert g["real_CRDS_GWCS_pilot"]["numeric_slit_id_error_corrected"] is True
c=x["physical_calibration_gates"]
for key in ("spectral_pixel_illumination_weighting_for_all_original_open_shutters",
            "real_closed_shutter_leakage_scattered_sky","all_four_full_stage2_science_coadds",
            "dispersed_STPSF_in_original_detector_and_physical_flux_injection"):
    assert c[key] is False
for key in ("line_16004_absolute_cgs","line_15850_OIII_confirmed",
            "line_13656_Pagamma_confirmed","source_identity","redshift","directed_gas_trail"):
    assert c[key] is None
assert g["scientific_status"]["full_physical_calibration"]=="OPEN"
print("PASS: 36 actual original MSA-neighbor JWST CRDS GWCS detector positions; physical throughput still OPEN")
