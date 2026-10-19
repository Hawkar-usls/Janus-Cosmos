#!/usr/bin/env python3
"""Verify both Stage2 diagnostic x1d archival receipts without physical promotion."""
import json,re
from pathlib import Path
d=Path(__file__).resolve().parents[1]/"domains"/"space"
x=json.loads((d/"MOCSS_BD_DIAGNOSTIC_X1D_SUCCESS_2026-10-19.json").read_text())
w=json.loads((d/"MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json").read_text())
assert x["status"]=="DIAGNOSTIC_X1D_CREATED_FOR_TWO_ARCHIVAL_ANGLES_NOT_CALIBRATED_ASTROPHYSICAL_FLUX"
assert x["synthetic_source_info_not_original_mpt"] and x["source_exposures"]==[8,13]
assert len(x["results"])==2
for run in x["results"]:
    assert re.fullmatch("[a-f0-9]{64}",run["artifact_json_sha256"])
    assert set(run["groups"])=={"B","D"}
    for k,v in run["groups"].items():
        assert v["x1d_created"] is True
        assert v["exposure_index"] in [8,13]
        s=w["exposures"][v["exposure_index"]-1]
        assert s["group"]==k
        assert s["expected_science_detector_from_WCS"].upper()==v["detector"]
        assert re.fullmatch("[a-f0-9]{64}",v["x1d_fits_sha256"])
        assert v["x1d_fits_bytes"]>200000
        if "source_fraction_xy" in v:
            assert all(0<a<1 for a in v["source_fraction_xy"])
assert x["manual_extract_1d_number_per_group"]==1
assert {"assign_wcs","extract_2d","wavecorr","flat_field","photom"}.issubset(x["reported_completed_steps"])
for k in ["A_B_C_D_valid_MOCSS_science_spectra","source_identity","redshift"]:
    assert x["science_gates"][k]=="OPEN"
for k in ["line_16004_absolute_flux","line_15850p793_OIII","line_13656p748_Pagamma"]:
    assert x["science_gates"][k]=="NOT_MEASURED"
assert x["science_gates"]["PSF_injection_physical_flux_recovery"]=="NOT_DONE"
print("PASS: B/D diagnostic x1d real archive receipts, no calibrated astrophysical flux")
