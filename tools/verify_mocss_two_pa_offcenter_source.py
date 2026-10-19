#!/usr/bin/env python3
"""Verify two-orientation source offsets from original independent JWST sky WCS.

The source_info record in temporary MSA is a reconstruction, NOT proof that
the MPT planned a real target, and cannot be used as science flux.
"""
import json, math, re
from pathlib import Path
base=Path(__file__).resolve().parents[1]/"domains"/"space"
x=json.loads((base/"MOCSS_TWO_PA_OFFCENTER_SOURCE_INFO_2026-10-18.json").read_text())
original=json.loads((base/"MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json").read_text())
assert x["status"].startswith("ORIGINAL_RATE_AND_MSA_SHA_VERIFIED")
assert x["provenance"]["crds"]=="jwst_1413.pmap"
assert re.fullmatch(r"[a-f0-9]{64}",x["provenance"]["raw_json_sha256"])
assert len(x["frames"])==2
for g in x["frames"]:
    r=original["exposures"][g["index"]-1]
    det=r["expected_science_detector_from_WCS"]
    assert det==g["detector"] and r["group"]==g["group"]
    item=r["detectors"][det]
    assert item["sha256"]==g["original_rate_sha256"]
    assert g["msa_columns"]==[g["primary_column"]-1,g["primary_column"],g["primary_column"]+1]
    assert item["selected_virtual_shutter"]["row"]==g["msa_row"]
    assert item["selected_virtual_shutter"]["quadrant"]==g["msa_quadrant"]
    ref=item["sky_registered_trace"]["16004A"]
    for index,coord in enumerate(("slit_x","slit_y")):
        assert math.isclose(g["wcs_normalized_source_slit_xy"][index],ref[coord],abs_tol=1e-10)
        assert math.isclose(g["synthetic_source_info_fraction_xy"][index],.5+ref[coord],abs_tol=1e-10)
    assert all(0 < v < 1 for v in g["synthetic_source_info_fraction_xy"])
    assert g["test_msa_sha256"]!=g["original_msa_sha256"]
assert x["source_marked_from_original_MPT"] is False
assert x["scientific_calibrated_source_x1d"] is False
assert x["shutter_transmission_calibrated"] is False
assert x["source_redshift"] is None and x["source_class"] is None and x["gas_trail"] is None
print("PASS: two original archived NIRSpec source-offset geometries; source class, flux and redshift remain OPEN")
