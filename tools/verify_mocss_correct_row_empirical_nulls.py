#!/usr/bin/env python3
"""Strict real-corrected-dispersion-axis MOCSS null gate; no source promotion."""
import json,math,re
from pathlib import Path
d=Path(__file__).resolve().parents[1]/"domains"/"space"
x=json.loads((d/"MOCSS_CORRECT_ROW_FOUR_X1D_EMPIRICAL_NULLS_2026-10-19.json").read_text())
y=json.loads((d/"MOCSS_FOUR_CORRECT_ROWS_DIAGNOSTIC_X1D_2026-10-19.json").read_text())
assert x["schema"]=="janus-mocss-correct-three-row-empirical-offline-null-v1"
assert x["run"]["source_four_corrected_x1d_run"]==y["provenance"]["run_id"]==39223778736
assert re.fullmatch("[a-f0-9]{64}",x["run"]["raw_artifact_json_sha256"])
assert len(x["groups"])==4 and [z["group"] for z in x["groups"]]==list("ABCD")
assert x["groups"][2]["primary_reduced_chi2"]>1000
assert x["groups"][2]["primary_fraction_formal_abs3"]>.75
for g in x["groups"]:
    assert g["primary_offline_controls"]>=15
    assert 0<=g["primary_fraction_formal_abs3"]<=1
    assert g["primary_reduced_chi2"]>0
assert x["nonindependent_null_centers"] is True
assert x["source_info_artificial"] is True
for k in ["PSF_bar_loss_calibrated","full_row_sky_background_validated","full_four_dither_stacks","two_dimensional_optical_line_injection"]:
    assert x[k] is False
for k in ["calibrated_flux_or_upper_limit","detected_companion_OIII","detected_companion_Pagamma","source_class","redshift","directed_gas_trail"]:
    assert x[k] is None
print("PASS: four real corrected-axis x1d formal-null spikes and spectroscopic scientific ceiling")
