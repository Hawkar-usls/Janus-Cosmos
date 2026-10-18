#!/usr/bin/env python3
"""Cross-check four original archive RATE fingerprints and procedural injections.

The fitted parametric Gaussian is exactly the injected Gaussian. This proves
only that the local estimator responds linearly; never interpret it as
astrophysical PSF completeness, calibrated significance or real line flux.
"""
import json
import math
import re
from pathlib import Path

root=Path(__file__).resolve().parents[1]/"domains"/"space"
rec=json.loads((root/"MOCSS_THREE_LINE_REAL_RATE_INJECTION_2026-10-18.json").read_text())
wcs=json.loads((root/"MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json").read_text())
assert rec["schema"]=="janus-mocss-real-rate-injection-summary-v1"
assert rec["status"]=="FOUR_ORIGINAL_RATE_EXPOSURES_GAUSSIAN_DN_S_INJECTION_PASS_NO_FLUX_AUTHORITY"
assert rec["provenance"]["workflow_conclusion"]=="success"
assert re.fullmatch("[0-9a-f]{64}",rec["provenance"]["artifact_original_json_sha256"])
assert list(rec["input_exposure_index_by_group"].values())==[1,8,10,13]
assert len(wcs["exposures"])==16
for g,idx in rec["input_exposure_index_by_group"].items():
    exp=wcs["exposures"][idx-1]
    assert exp["group"]==g
    det=exp["expected_science_detector_from_WCS"]
    assert exp["detectors"][det]["sha256"]==rec["provenance"]["original_file_sha256"][g]
    assert len(rec["raw_template_peak_scores_dn_per_s"][g])==3
    assert len(rec["local_offline_template_mad_dn_per_s"][g])==3
    assert all(math.isfinite(x) for x in rec["raw_template_peak_scores_dn_per_s"][g].values())
    assert all(x>0 for x in rec["local_offline_template_mad_dn_per_s"][g].values())
assert rec["injected_peak_dn_per_s"]==[0.01,0.04,0.12]
assert math.isclose(rec["independently_fitted_injection_change_ratio_to_known_input_approx"],1.0)
assert set(rec["observed_wavelength_A"])=={"primary","oiii_companion","pagamma"}
assert math.isclose(rec["observed_wavelength_A"]["oiii_companion"],16004*4960.295/5008.239,abs_tol=0.005)
assert math.isclose(rec["observed_wavelength_A"]["pagamma"],16004*10941.1/12821.6,abs_tol=0.005)
assert rec["gate"]["synthetic_gaussian_raw_rate_injection"]=="PASS_PROCEDURAL"
assert rec["gate"]["three_line_calibrated_fluxes"]=="NOT_PERFORMED"
assert rec["gate"]["line_identity"]=="OPEN"
assert rec["gate"]["physical_source"]=="OPEN"
assert rec["gate"]["redshift"]=="OPEN"
assert "NOT a completeness test" in rec["interpretation"]
print("PASS: real four FITS fingerprints and 12 three-line geometric Gaussian-injection controls")
print("RATE_INJECTION_ALGORITHM_LINEARITY=PASS_PROCEDURAL")
print("NIRSPEC_CALIBRATED_LINE_FLUX=NOT_ESTABLISHED")
