#!/usr/bin/env python3
"""Offline gate: source-aware TEST geometry grounded in independent archived WCS.

Never promotes synthetic SOURCE_INFO into an original MPT source or flux.
"""
import json,math,re
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"domains"/"space"
x=json.loads((p/"MOCSS_THREE_SHUTTER_SOURCE_AWARE_GEOMETRY_2026-10-18.json").read_text())
w=json.loads((p/"MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json").read_text())
assert x["schema"]=="janus-mocss-three-shutter-source-aware-geometry-receipt-v1"
assert x["status"]=="REAL_RATE_AND_MSA_SHA_VERIFIED_SYNTHETIC_SOURCE_INFO_GEOMETRY_PASS_NOT_SCIENCE"
assert len(x["frames"])==2 and [r["group"] for r in x["frames"]]==["B","D"]
assert re.fullmatch("[a-f0-9]{64}",x["provenance"]["full_artifact_raw_json_sha256"])
assert x["source_info_change"]["orig_mpt_catalog_source"] is False
assert x["source_info_change"]["original_rate_and_msa_files_mutated"] is False
assert x["source_info_change"]["source_position_fractional_inside_shutter"]==[0.5,0.5]
assert x["source_info_change"]["one_marked_primary_per_slitlet"] is True
for rec in x["frames"]:
    source=w["exposures"][rec["exposure_index"]-1]
    det=source["expected_science_detector_from_WCS"]
    assert source["group"]==rec["group"] and det==rec["detector"]
    v=source["detectors"][det]
    assert v["sha256"]==rec["original_rate_sha256"]
    assert v["selected_virtual_shutter"]["row"]==rec["msa_row"]
    assert v["selected_virtual_shutter"]["column"]==rec["primary_column"]
    assert v["selected_virtual_shutter"]["quadrant"]==rec["msa_quadrant"]
    assert rec["columns"]==[rec["primary_column"]-1,rec["primary_column"],rec["primary_column"]+1]
    assert len(rec["original_msa_sha256"])==64 and len(rec["test_msa_sha256"])==64
    assert rec["original_msa_sha256"]!=rec["test_msa_sha256"]
    ref=v["sky_registered_trace"]["16004A"]
    for a,b in zip(rec["normalized_original_slit_xy"],[ref["slit_x"],ref["slit_y"]]):
        assert math.isclose(a,b,abs_tol=1e-9)
for k,v in x["gate_states"].items():
    if k!="three_shutter_geometry_source_info":
        assert v in {"NOT_PERFORMED","NOT_CLAIMED","OPEN","NOT_DETECTED"}
assert x["gate_states"]["three_shutter_geometry_source_info"]=="PASS_BOTH_ORIENTATIONS"
print("PASS: 2 real archival geometries, source MSA only synthetic, science x1d unresolved")
