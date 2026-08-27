#!/usr/bin/env python3
"""JANUS COSMOS third-generation catalog-epoch identity kernel.

Tests whether independent catalog positions are statistically compatible with one
Gaia astrometric worldline at their recorded/reference epochs. Shared Gaia
uncertainty is carried across epochs in one joint block covariance.

Scientific boundaries: release years are never observation epochs; missing
covariance is never silently zeroed; Mahalanobis compatibility is not identity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass, asdict
from typing import Any, Dict, Mapping, Sequence

import numpy as np

from janus_astrometric_worldlines import AstrometricState, covariance_matrix

SCHEMA = "janus.cosmos.catalog-epoch-identity.v1"
GAIA_COMMON_EPOCH_JYEAR = 2016.0
TAIL_ALPHA = 0.0027


def _sha(obj: object) -> str:
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    return hashlib.sha256(raw).hexdigest()


def mjd_to_jyear(mjd: float) -> float:
    return 2000.0 + (float(mjd) - 51544.5) / 365.25


def jd_to_mjd(jd: float) -> float:
    x = float(jd)
    if x > 2_000_000.0:
        return x - 2_400_000.5
    if 40_000.0 < x < 100_000.0:
        return x
    raise ValueError(f"UNRECOGNIZED_JD_SCALE:{x}")


def tangent_offset_mas(ref_ra_deg: float, ref_dec_deg: float,
                       ra_deg: float, dec_deg: float) -> np.ndarray:
    dra = (float(ra_deg) - float(ref_ra_deg) + 180.0) % 360.0 - 180.0
    return np.array([
        dra * math.cos(math.radians(ref_dec_deg)) * 3_600_000.0,
        (float(dec_deg) - float(ref_dec_deg)) * 3_600_000.0,
    ], dtype=float)


def ellipse_covariance_mas2(major_arcsec: float, minor_arcsec: float,
                            pa_deg_east_of_north: float) -> np.ndarray:
    a = float(major_arcsec) * 1000.0
    b = float(minor_arcsec) * 1000.0
    t = math.radians(float(pa_deg_east_of_north))
    u = np.array([math.sin(t), math.cos(t)])
    v = np.array([math.cos(t), -math.sin(t)])
    out = a*a*np.outer(u, u) + b*b*np.outer(v, v)
    _assert_psd(out, "ERROR_ELLIPSE")
    return out


def covariance_from_sigmas_corr_mas2(sx_arcsec: float, sy_arcsec: float,
                                     rho: float) -> np.ndarray:
    sx, sy, r = float(sx_arcsec)*1000.0, float(sy_arcsec)*1000.0, float(rho)
    if not -1.0 <= r <= 1.0:
        raise ValueError("CORRELATION_RANGE")
    out = np.array([[sx*sx, r*sx*sy], [r*sx*sy, sy*sy]], dtype=float)
    _assert_psd(out, "SIGMA_CORR")
    return out


def covariance_from_cosigma_mas2(sx_arcsec: float, sy_arcsec: float,
                                 cosigma_arcsec: float) -> np.ndarray:
    sx, sy, co = float(sx_arcsec)*1000.0, float(sy_arcsec)*1000.0, float(cosigma_arcsec)*1000.0
    out = np.array([[sx*sx, co*abs(co)], [co*abs(co), sy*sy]], dtype=float)
    _assert_psd(out, "COSIGMA")
    return out


def conservative_unknown_corr_bound_mas2(sx_arcsec: float, sy_arcsec: float,
                                         systematic_floor_mas: float = 0.0) -> np.ndarray:
    sx2 = (float(sx_arcsec)*1000.0)**2 + float(systematic_floor_mas)**2
    sy2 = (float(sy_arcsec)*1000.0)**2 + float(systematic_floor_mas)**2
    # Any 2x2 covariance with these diagonal variances has lambda_max <= trace.
    # trace*I is therefore a conservative Loewner upper bound without inventing rho.
    return np.eye(2) * (sx2 + sy2)


def _assert_psd(cov: np.ndarray, label: str) -> None:
    c = np.asarray(cov, dtype=float)
    if c.ndim != 2 or c.shape[0] != c.shape[1] or not np.all(np.isfinite(c)):
        raise ValueError(f"{label}:BAD_COVARIANCE")
    eig = np.linalg.eigvalsh((c+c.T)/2.0)
    if float(eig.min()) < -1e-7*max(1.0, float(np.abs(eig).max())):
        raise ValueError(f"{label}:COVARIANCE_NOT_PSD")


@dataclass(frozen=True)
class CatalogMeasurement:
    catalog: str
    source_id: str
    ra_deg: float
    dec_deg: float
    epoch_mjd: float
    covariance_mas2: Sequence[Sequence[float]]
    covariance_status: str
    epoch_status: str
    position_semantics: str
    provenance: Mapping[str, Any]

    @property
    def epoch_jyear(self) -> float:
        return mjd_to_jyear(self.epoch_mjd)

    def validate(self) -> None:
        for name in ("ra_deg", "dec_deg", "epoch_mjd"):
            if not math.isfinite(float(getattr(self, name))):
                raise ValueError(f"NONFINITE:{name}")
        if not -90.0 <= self.dec_deg <= 90.0:
            raise ValueError("DEC_RANGE")
        if "RELEASE_YEAR" in self.epoch_status:
            raise ValueError("RELEASE_YEAR_FORBIDDEN")
        _assert_psd(np.asarray(self.covariance_mas2, dtype=float), self.catalog)


def gaia_cov_native(state: AstrometricState) -> np.ndarray:
    cov, status = covariance_matrix(state)
    if status != "FULL_CATALOG_COVARIANCE":
        raise ValueError("FULL_GAIA_COVARIANCE_REQUIRED")
    out = np.asarray(cov, dtype=float)
    _assert_psd(out, "GAIA_5X5")
    return out


def design_matrix(dt_years: float) -> np.ndarray:
    # Gaia native order: [ra*, dec, parallax, pmra*, pmdec]. Catalog reference
    # positions are compared in their ICRS/J2000 frame; no apparent parallax shift
    # is invented without observer/time parallax factors.
    dt = float(dt_years)
    return np.array([[1,0,0,dt,0], [0,1,0,0,dt]], dtype=float)


def predicted_gaia_offset_mas(state: AstrometricState, epoch_jyear: float) -> np.ndarray:
    dt = float(epoch_jyear) - float(state.ref_epoch_jyear)
    return np.array([state.pmra_masyr*dt, state.pmdec_masyr*dt], dtype=float)


def residual_at_catalog_epoch(state: AstrometricState, m: CatalogMeasurement) -> np.ndarray:
    observed = tangent_offset_mas(state.ra_deg, state.dec_deg, m.ra_deg, m.dec_deg)
    return observed - predicted_gaia_offset_mas(state, m.epoch_jyear)


def transport_catalog_to_common_epoch(state: AstrometricState, m: CatalogMeasurement,
                                      common_epoch_jyear: float = GAIA_COMMON_EPOCH_JYEAR) -> Dict[str, float | bool]:
    # This is conditional transport under the identity hypothesis, not an
    # independently measured proper motion of the external catalog source.
    x = tangent_offset_mas(state.ra_deg, state.dec_deg, m.ra_deg, m.dec_deg)
    dt = float(common_epoch_jyear) - m.epoch_jyear
    moved = x + np.array([state.pmra_masyr, state.pmdec_masyr])*dt
    return {
        "x_east_mas_relative_to_gaia_ref_position": float(moved[0]),
        "y_north_mas_relative_to_gaia_ref_position": float(moved[1]),
        "common_epoch_jyear": float(common_epoch_jyear),
        "transport_is_conditional_on_gaia_identity_hypothesis": True,
    }


def chi2_survival_even_dof(chi2: float, dof: int) -> float:
    if dof <= 0 or dof % 2:
        raise ValueError("ONLY_POSITIVE_EVEN_DOF_SUPPORTED")
    x = max(0.0, float(chi2))/2.0
    n = dof//2
    return math.exp(-x)*sum((x**k)/math.factorial(k) for k in range(n))


def mahalanobis(residual: np.ndarray, covariance: np.ndarray) -> float:
    r, c = np.asarray(residual, dtype=float), np.asarray(covariance, dtype=float)
    _assert_psd(c, "MAHALANOBIS")
    inv = np.linalg.pinv(c, hermitian=True, rcond=1e-12)
    return float(r.T @ inv @ r)


def evaluate_identity(state: AstrometricState, measurements: Sequence[CatalogMeasurement],
                      common_epoch_jyear: float = GAIA_COMMON_EPOCH_JYEAR,
                      tail_alpha: float = TAIL_ALPHA) -> Dict[str, Any]:
    state.validate()
    if not measurements:
        return {"schema": SCHEMA, "status": "I_DO_NOT_KNOW", "reason": "NO_INDEPENDENT_CATALOG_MEASUREMENTS"}
    for m in measurements:
        m.validate()

    g5 = gaia_cov_native(state)
    pairwise: Dict[str, Any] = {}
    residuals, hs, cats, order = [], [], [], []
    for m in measurements:
        h = design_matrix(m.epoch_jyear-state.ref_epoch_jyear)
        cg = h @ g5 @ h.T
        cc = np.asarray(m.covariance_mas2, dtype=float)
        ctot = cg + cc
        r = residual_at_catalog_epoch(state, m)
        x2 = mahalanobis(r, ctot)
        p = chi2_survival_even_dof(x2, 2)
        pairwise[m.catalog] = {
            "source_id": m.source_id,
            "catalog_epoch_mjd": m.epoch_mjd,
            "catalog_epoch_jyear": m.epoch_jyear,
            "epoch_status": m.epoch_status,
            "position_semantics": m.position_semantics,
            "residual_east_north_mas": [float(r[0]), float(r[1])],
            "residual_norm_mas": float(np.linalg.norm(r)),
            "gaia_prediction_covariance_mas2": cg.tolist(),
            "catalog_position_covariance_mas2": cc.tolist(),
            "catalog_covariance_status": m.covariance_status,
            "total_covariance_mas2": ctot.tolist(),
            "mahalanobis_chi2": x2,
            "dof": 2,
            "tail_probability": p,
            "tail_alpha_preregistered": tail_alpha,
            "compatibility_status": "NOT_REJECTED_AT_PREREGISTERED_TAIL" if p >= tail_alpha else "REJECTED_AT_PREREGISTERED_TAIL",
            "transported_to_common_epoch": transport_catalog_to_common_epoch(state, m, common_epoch_jyear),
            "provenance": dict(m.provenance),
        }
        residuals.append(r); hs.append(h); cats.append(cc); order.append(m.catalog)

    # Joint GLS: off-diagonal blocks carry the shared Gaia uncertainty between
    # catalog epochs; only each catalog's own measurement covariance is diagonal.
    n = len(measurements)
    rjoint = np.concatenate(residuals)
    cjoint = np.zeros((2*n, 2*n), dtype=float)
    for i in range(n):
        for j in range(n):
            block = hs[i] @ g5 @ hs[j].T
            if i == j:
                block = block + cats[i]
            cjoint[2*i:2*i+2, 2*j:2*j+2] = block
    _assert_psd(cjoint, "JOINT_GLS")
    x2_joint = mahalanobis(rjoint, cjoint)
    dof = 2*n
    p_joint = chi2_survival_even_dof(x2_joint, dof)
    conservative = [m.catalog for m in measurements if "CONSERVATIVE" in m.covariance_status or "BOUND" in m.covariance_status]

    out: Dict[str, Any] = {
        "schema": SCHEMA,
        "formula": "RESPICIENS_ET_PROSPICIENS_GEN3",
        "common_epoch_jyear": common_epoch_jyear,
        "gaia_source": {
            "source_id": state.source_id, "catalog": state.catalog,
            "ra_deg": state.ra_deg, "dec_deg": state.dec_deg,
            "reference_epoch_jyear": state.ref_epoch_jyear,
            "pmra_masyr": state.pmra_masyr, "pmdec_masyr": state.pmdec_masyr,
            "parallax_mas": state.parallax_mas,
        },
        "catalog_order": order,
        "pairwise": pairwise,
        "joint_gls": {
            "residual_vector_east_north_mas": [float(x) for x in rjoint],
            "block_covariance_mas2": cjoint.tolist(),
            "shared_gaia_cross_epoch_covariance_included": True,
            "mahalanobis_chi2": x2_joint,
            "dof": dof,
            "tail_probability": p_joint,
            "tail_alpha_preregistered": tail_alpha,
            "compatibility_status": "NOT_REJECTED_AT_PREREGISTERED_TAIL" if p_joint >= tail_alpha else "REJECTED_AT_PREREGISTERED_TAIL",
            "catalogs_using_conservative_covariance_upper_bound": conservative,
        },
        "epistemic_firewall": {
            "compatibility_is_identity_proof": False,
            "cross_catalog_agreement_is_love_edem_identity": False,
            "release_year_is_observation_epoch": False,
            "missing_covariance_replaced_with_zero": False,
            "shared_gaia_uncertainty_double_counted": False,
            "common_epoch_transport_is_independent_motion_measurement": False,
            "negative_result_is_valid": True,
            "i_do_not_know_is_valid": True,
        },
        "claim_ceiling": "INDEPENDENT_CATALOG_COMPATIBILITY_WITH_A_GAIA_STELLAR_WORLDLINE_ONLY__NOT_LOVE_EDEM_IDENTITY_OR_ANOMALY",
    }
    out["input_sha256"] = _sha({
        "gaia": asdict(state), "measurements": [asdict(m) for m in measurements],
        "common_epoch_jyear": common_epoch_jyear, "tail_alpha": tail_alpha,
    })
    out["result_sha256"] = _sha(out)
    return out


def self_test() -> None:
    keys = ["ra_dec","ra_parallax","ra_pmra","ra_pmdec","dec_parallax","dec_pmra","dec_pmdec","parallax_pmra","parallax_pmdec","pmra_pmdec"]
    s = AstrometricState(
        source_id="synthetic", catalog="GAIA_DR3", ra_deg=10.0, dec_deg=20.0,
        ref_epoch_jyear=2016.0, parallax_mas=1.0, pmra_masyr=10.0, pmdec_masyr=-5.0,
        ra_error_mas=.1, dec_error_mas=.1, parallax_error_mas=.1,
        pmra_error_masyr=.05, pmdec_error_masyr=.05,
        correlations={k:0.0 for k in keys},
    )
    ep = 2010.0; dt = ep-2016.0
    x, y = s.pmra_masyr*dt, s.pmdec_masyr*dt
    ra = s.ra_deg + x/(3_600_000.0*math.cos(math.radians(s.dec_deg)))
    dec = s.dec_deg + y/3_600_000.0
    m = CatalogMeasurement(
        catalog="SYNTH", source_id="x", ra_deg=ra, dec_deg=dec,
        epoch_mjd=51544.5+(ep-2000.0)*365.25,
        covariance_mas2=[[100.0,0.0],[0.0,100.0]],
        covariance_status="MEASURED_FULL_2D_COVARIANCE",
        epoch_status="EXACT_REFERENCE_EPOCH", position_semantics="SYNTHETIC",
        provenance={"test":True},
    )
    a=evaluate_identity(s,[m]); b=evaluate_identity(s,[m])
    assert a["result_sha256"]==b["result_sha256"]
    assert a["pairwise"]["SYNTH"]["mahalanobis_chi2"] < 1e-8
    assert a["joint_gls"]["shared_gaia_cross_epoch_covariance_included"] is True
    assert a["epistemic_firewall"]["compatibility_is_identity_proof"] is False
    print("JANUS_CATALOG_EPOCH_IDENTITY_SELF_TEST=PASS")


def main() -> int:
    p=argparse.ArgumentParser(); p.add_argument("--self-test",action="store_true"); a=p.parse_args()
    if a.self_test: self_test(); return 0
    p.error("use this kernel through a target-specific runner"); return 2


if __name__=="__main__":
    raise SystemExit(main())
