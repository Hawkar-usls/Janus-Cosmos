#!/usr/bin/env python3
"""Proof-of-observation gate for 4 real diagnostic JWST x1d FITS.

This verifies numerical summaries from an Actions archival replay; it does not
turn wrong-axis synthetic MSA diagnostics into valid MOCSS spectroscopy.
"""
import json
import re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/"domains"/"space"
d=json.loads((ROOT/"MOCSS_FOUR_X1D_OFFLINE_NULL_STRESS_2026-10-19.json").read_text())
original=json.loads((ROOT/"MOCSS_DISPERSION_THREE_ROW_CORRECTION_2026-10-19.json").read_text())
assert d["schema"]=="janus-mocss-four-group-real-x1d-controls-and-axis-correction-v1"
assert len(d["data"])==4
assert [x["group"] for x in d["data"]]==list("ABCD")
assert re.fullmatch("[a-f0-9]{64}",d["evidence"]["control_full_json_sha256"])
assert d["evidence"]["data_reprocessed_SCIENCE_VALID"] is False
assert d["evidence"]["control_run_id"]==39223778656
for item in d["data"]:
    assert re.fullmatch("[a-f0-9]{64}",item["x1d_sha256"])
    assert set(item["windows"])=={"primary_16004","OIII_15850","Pa_gamma_13656"}
    for sub in item["windows"].values():
        assert sub["control_n"]>=15
        assert 0 < sub["control_fraction_abs3"] <= 1
        assert sub["max_abs_formal_control"]>3
        assert sub["reduced_chi2"]>0
assert d["data"][2]["windows"]["primary_16004"]["reduced_chi2"]>1000
assert d["data"][3]["windows"]["OIII_15850"]["formal_ratio"]<-10
assert original["corrected_geometry"]["physical_dispersion_axis"]=="MSA_SHUTTER_ROW"
assert original["corrected_geometry"]["cross_dispersion_axis"]=="MSA_SHUTTER_COLUMN"
for item in original["corrected_geometry"]["per_group_sample_geometry"]:
    assert item["rows"][1]-item["rows"][0]==1
    assert item["rows"][2]-item["rows"][1]==1
assert not original["previous_diagnostic_geometries"]["physical_throughput_reproduction_valid"]
assert d["one_dimensional_mismatch_injections"]["response_fraction_is_NOT_optical_throughput"]
assert d["gate_states"]["old_three_column_throughput_geometry"]=="INVALID_FOR_MOCSS_DISPERSION_LOSS"
assert d["gate_states"]["science_grade_three_line_flux"]=="NOT_MEASURED"
assert d["gate_states"]["redshift"]=="OPEN"
print("PASS: 4 real x1d SHA-pinned diagnostic controls and physically correct MSA row axis")
print("PRIOR_THREE_COLUMN_MODEL=NOT_SCIENCE_VALID; THREE_LINE_PHYSICAL_FLUX=NOT_MEASURED")
