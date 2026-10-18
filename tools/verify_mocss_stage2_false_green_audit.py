#!/usr/bin/env python3
"""Frozen false-green MOCSS Stage2 receipt auditor.

The previous Action returned success despite in-receipt AttributeError. This
verifier checks the failure provenance to ensure later training cannot count
it as a calibrated NIRSpec x1d nor a negative line detection.
"""
import json
import re
from pathlib import Path

base=Path(__file__).resolve().parents[1] / "domains" / "space"
x=json.loads((base/"MOCSS_STAGE2_LIST_RETURN_FAILURE_AUDIT_2026-10-18.json").read_text())
w=json.loads((base/"MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json").read_text())
assert x["schema"]=="janus-mocss-spec2-feasibility-failure-audit-v1"
assert x["science_status"].startswith("PIPELINE_RETURN_TYPE_DIAGNOSTIC_ERROR")
old=x["previous_run"]
assert old["run_id"]==39143330154
assert old["workflow_conclusion"]=="success" and old["workflow_success_misleading"]
assert old["scientific_negative_result_established"] is False
assert old["error_occurred_after_pipeline_return"] is True
assert re.fullmatch(r"[0-9a-f]{64}",old["artifact_json_sha256"])
for rec in old["archive_exposures"]:
    exp=w["exposures"][rec["exposure_index"]-1]
    assert exp["group"]==rec["group"]
    det=exp["expected_science_detector_from_WCS"]
    assert rec["rate_sha256"]==exp["detectors"][det]["sha256"]
    assert rec["returned_model_type"]=="list"
    assert rec["output_calibrated_x1d"] is False
    assert "AttributeError" in rec["observed_error"]
for value in x["remaining_calibration"].values():
    assert isinstance(value,bool)
assert x["remaining_calibration"]["quantitative_real_line_flux_by_group"] is False
assert x["remaining_calibration"]["oiii_companion_15850p793"] is False
assert x["remaining_calibration"]["pagamma_companion_13656p748"] is False
assert x["no_claims"]["source_identity"] is None and x["no_claims"]["redshift"] is None
print("PASS: false-green Stage2 status not a calibrated source x1d or negative spectral measurement")
