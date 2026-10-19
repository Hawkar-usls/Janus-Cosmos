#!/usr/bin/env python3
"""Full16 geometric CRDS/3 original open slit x 3 wavelength receipt validation.

Uses independent original WCS receipt as source-index authority. Exact
input evidence is retained in four SHA-256-identified GitHub artifacts.
No inferred photon throughput or line flux is authorized.
"""
import json,re,math
from pathlib import Path
base=Path(__file__).resolve().parents[1]/"domains"/"space"
x=json.loads((base/"MOCSS_CRDS_GWCS_ALL16_OPEN_NEIGHBOR_FOOTPRINT_2026-10-19.json").read_text())
w=json.loads((base/"MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json").read_text())
gate=json.loads((base/"MOCSS_CRDS_WCS_FORWARD_ILLUMINATION_GATE_2026-10-19.json").read_text())
assert x["status"]=="ALL16_ORIGINAL_SHA_VERIFIED_48_OPEN_NEIGHBOR_SLITS_144_CRDS_WCS_CENTERS_ON_DETECTOR"
assert x["n_rate"]==len(w["exposures"])==16
assert x["n_original_MSA"]==2 and x["n_open_shutter_slits"]==48
assert x["n_three_line_center_samples"]==144
assert x["upstream"]["run_id"]==39243558142
assert set(x["wavelength_angstrom"])=={"primary","OIII_4960","Pa_gamma"}
counter=[]
for group in "ABCD":
    v=x["groups"][group]
    assert v["detector"]==("NRS2" if group in "AB" else "NRS1")
    assert len(v["exps"])==len(v["xy"])==4
    assert re.fullmatch("[a-f0-9]{64}",v["sha"])
    for i,(index,xy) in enumerate(zip(v["exps"],v["xy"])):
        rec=w["exposures"][index-1]
        assert rec["group"]==group
        assert rec["expected_science_detector_from_WCS"].upper()==v["detector"]
        assert len(xy)==2 and all(math.isfinite(z) and 0<=z<2048 for z in xy)
        counter.append(index)
assert sorted(counter)==list(range(1,17))
g=x["science_gates"]
assert g["CRDS_GWCS_nearest_open_shutter_geometry"]=="PASS_144"
for key in ("all_original_open_shutter_detector_illumination_map",
            "dispersed_STPSF_PSF_throughput","custom_projected_sky_and_leakage"):
    assert g[key]=="NOT_CONSTRUCTED"
assert g["calibrated_four_stack_stage2"]=="NOT_EXECUTED"
assert g["physical_flux_injected_detector_source"]=="NOT_EXECUTED"
for key in ("source_flux_cgs","oiii_detected","pagamma_detected","source_class","redshift","gas_trail"):
    assert g[key] is None
assert gate["real_CRDS_GWCS_pilot"]["all16_near_source_3_original_open_shutter_3_lambda"]=="PASS_144_CENTERS_OFFICIAL_CRDS_GWCS"
assert gate["scientific_status"]["full_physical_calibration"]=="OPEN"
print("PASS: 144 independently archived official JWST CRDS WCS slit center geometries, physical flux remains OPEN")
