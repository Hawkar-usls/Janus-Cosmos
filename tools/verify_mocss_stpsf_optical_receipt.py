#!/usr/bin/env python3
"""Audit archive-grounded STPSF finite-FOV ratios, NOT science MOCSS throughput."""
import json,re,math
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"domains"/"space"/"MOCSS_STPSF_OPTICAL_THREE_WAVELENGTHS_2026-10-19.json"
r=json.loads(p.read_text())
assert r["status"]=="OFFICIAL_STPSF_MONOCHROMATIC_OPTICAL_MASK_FOV_RATIOS_MEASURED_NOT_MOCSS_SCIENCE_THROUGHPUT"
assert r["provenance"]["stpsf_version"]=="2.1.0"
assert r["provenance"]["run_id"]==39229599393
assert re.fullmatch("[0-9a-f]{64}",r["provenance"]["artifact_json_sha256"])
assert set(r["results"])=={"primary","OIII_4960","Pa_gamma"}
for l in r["results"].values():
    for prop in ("center","offcenter"):
        assert 0<l[prop+"_fraction"]<1.0
        assert math.isclose(l[prop+"_fraction"]+l[prop+"_loss"],1.,abs_tol=1e-12)
    assert 0.01<l["center_loss"]<.20
assert math.isclose(r["author_reported_bar_only_loss_at_1p6_um_about"],.033)
assert r["source_class"] is None and r["redshift"] is None
assert "NO_DIRECT_CONFLICT_TEST" in r["model_vs_author_mismatch"]
print("PASS: official STPSF 3 lambda two-position optical PSF sums; science-grade MOCSS throughput UNCALIBRATED")
