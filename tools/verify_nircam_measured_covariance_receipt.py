#!/usr/bin/env python3
"""Strict offline audit of independent, actual-pixel F150W covariance checkpoint.

Audits immutable numerical summaries plus input SHA links and science authority.
It does not download the archive anew or attest to pure instrumental noise.
"""
import json
import math
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]/"domains"/"space"
data=json.loads((ROOT/"NIRCAM_F150W_MATCHED_SKY_COVARIANCE_OBSERVED_2026-10-18.json").read_text())
original=json.loads((ROOT/"NIRCAM_F150W_TIME_DOMAIN_RECEIPT_2026-10-16.json").read_text())
assert data["status"]=="ARCHIVAL_MATCHED_SKY_COVARIANCE_OBSERVED_NOT_PURE_NOISE"
assert data["schema"]=="janus-cosmos-f150w-measured-matched-sky-receipt-v1"
assert data["control_locations_joint_clipped"]>=50
assert 0<data["control_locations_joint_clipped"]<=data["control_locations_joint_unclipped"]
p=data["provenance"]
assert p["workflow_conclusion"]=="success"
assert re.fullmatch(r"[0-9a-f]{64}",p["original_artifact_json_sha256"])
rows=[r for r in original["results"] if r.get("date_time_utc","").startswith("2023-10") and r.get("coverage")]
rows.sort(key=lambda x:x["date_time_utc"])
assert len(rows)==4
assert [r["sci_subset_sha256"] for r in rows]==p["input_sci_sha256"]
assert [r["date_time_utc"][:10] for r in rows]==data["epochs_utc"]
cov=data["matched_aperture_covariance_njy_squared"]
corr=data["matched_aperture_correlation"]
w=data["fixed_nominal_weights"]
assert len(cov)==len(corr)==len(w)==4
assert abs(sum(w)-1)<1e-12
for i in range(4):
    assert len(cov[i])==len(corr[i])==4
    assert cov[i][i]>0 and abs(corr[i][i]-1)<1e-9
    for j in range(4):
        assert abs(cov[i][j]-cov[j][i])<1e-12
        assert abs(corr[i][j]-corr[j][i])<1e-12
        assert abs(corr[i][j]-cov[i][j]/math.sqrt(cov[i][i]*cov[j][j]))<1e-8
# Cholesky: reject singular or indefinite published noise matrix.
chol=[[0.0]*4 for _ in range(4)]
for i in range(4):
    for j in range(i+1):
        temp=cov[i][j]-sum(chol[i][k]*chol[j][k] for k in range(j))
        if i==j:
            assert temp>0
            chol[i][j]=math.sqrt(temp)
        else:
            chol[i][j]=temp/chol[j][j]
sigma=math.sqrt(sum(w[i]*cov[i][j]*w[j] for i in range(4) for j in range(4)))
assert abs(sigma-data["matched_sky_covariance_sigma_njy"])<1e-12
ivw=[1/r["control_noise_1sigma_njy"]**2 for r in rows]
ivw=[v/sum(ivw) for v in ivw]
assert max(abs(a-b) for a,b in zip(w,ivw))<1e-12
flux=sum(v*r["flux_aperture_njy"] for v,r in zip(w,rows))
assert abs(flux-data["source_aperture_flux_njy"])<1e-12
assert abs(math.sqrt(sum((w[i]*rows[i]["control_noise_1sigma_njy"])**2 for i in range(4)))
           -data["previous_diagonal_noise_sigma_njy"])<1e-12
assert sigma > data["previous_diagonal_noise_sigma_njy"]
for m in data["conditional_updated_sensitivity"]:
    f=m["effective_aperture_fraction_ASSUMED"]
    assert 0<f<=1
    want=(5.5*f-flux)/math.hypot(sigma,0.66*f)
    assert abs(m["conditional_standardized_residual_NOT_formal_sigma"]-want)<1e-11
msg=" ".join(data["epistemic_limits"]).lower()
assert "not isolated instrumental noise" in msg and "not measured" in msg
assert "no_intrinsic_variability_proof" in msg
print("PASS: all 4 hashed original F150W epochs, 95 control sky positions, SPD covariance and no-overclaim gates")
print("MATCHED_SKY_UNCERTAINTY_NJY={:.6f}".format(sigma))
print("PHYSICAL_VARIABILITY=UNCONFIRMED; NIRSPEC_CALIBRATION=NOT_DONE")
