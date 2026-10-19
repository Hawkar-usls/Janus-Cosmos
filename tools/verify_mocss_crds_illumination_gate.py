#!/usr/bin/env python3
"""JANUS COSMOS MOCSS: original MSA aperture is not detector illumination.

Evidence and missing physics checks. Never promote a WCS envelope or STPSF
monochromatic PSF ratio to a calibrated cgs source flux or detection.
"""
import json
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"domains"/"space"
g=json.loads((p/"MOCSS_CRDS_WCS_FORWARD_ILLUMINATION_GATE_2026-10-19.json").read_text())
e=json.loads((p/"MOCSS_ORIGINAL_MSA_16_EXPOSURE_TOPOLOGY_2026-10-19.json").read_text())
assert g["schema"]=="janus-mocss-physical-detector-illumination-forward-operator-gate-v1"
assert g["real_CRDS_GWCS_pilot"]["input_original_exposure_indices"]==[1,8,10,13]
assert g["real_CRDS_GWCS_pilot"]["official_WCS_detector_footprints_measured"].startswith("PASS_12_ORIGINAL_OPEN_SHUTTERS_36_FINITE")
assert g["real_CRDS_GWCS_pilot"]["full_WCS_all_shutters_all_exposures"] is False
assert len(g["physical_forward_model_requirements"])==8
assert g["physical_forward_model_requirements"][0]["status"]=="ORIGINAL_COMMAND_APERTURE_MASK_VERIFIED"
for s in g["physical_forward_model_requirements"][2:]:
    assert s["status"] in ("NOT_CONSTRUCTED","NOT_PERFORMED","ONLY_GAUSSIAN_RATE_WINDOW_DIAGNOSTIC_EXISTS")
assert g["invariants"]["source_is_original_MPT_catalog_source"] is False
assert all(v is True for k,v in g["invariants"].items() if k != "source_is_original_MPT_catalog_source")
assert e["calibration_gates"]["physical_MSA_plane_binary_open_closed_mask"].startswith("PASS")
assert e["calibration_gates"]["wavelength_dependent_full_detector_pixel_illumination"]=="NOT_PERFORMED"
assert g["scientific_status"]["full_physical_calibration"]=="OPEN"
for key in ("full_source_flux_at_16004","OIII_15850_detected","Pa_gamma_13656_detected",
            "source_redshift","source_class","gas_trail"):
    assert g["scientific_status"][key] is None
print("PASS: original shutter plane and 4 anchor 36 CRDS geometries, physical forward illumination/line flux OPEN")
