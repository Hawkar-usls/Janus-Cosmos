#!/usr/bin/env python3
"""
JANUS COSMOS epoch-aware astrometric worldlines.

Scientific scope:
- propagate catalog astrometry in a local tangent plane across a bounded epoch window;
- sample measured astrometric uncertainty, including a catalog covariance matrix when available;
- use parallax as a conservative annual-displacement envelope unless observer-specific
  parallax factors are explicitly supplied;
- never treat propagation as observation or object identity.

This module is for astronomical catalog cross-checking only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import statistics
from dataclasses import dataclass, asdict
from typing import Dict, List, Mapping, Sequence

SCHEMA = "janus.cosmos.astrometric-worldlines.v1"
GAIA_DR3_REF_EPOCH_JYEAR = 2016.0
DEFAULT_EPOCH_START = 1900.0
DEFAULT_EPOCH_END = 2100.0
DEFAULT_GATES_ARCSEC = (1.5, 3.0)
PARAM_ORDER = ("ra", "dec", "parallax", "pmra", "pmdec")


def _sha(obj: object) -> str:
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    return hashlib.sha256(raw).hexdigest()


def _quantile(values: Sequence[float], q: float) -> float:
    xs = sorted(float(x) for x in values)
    if not xs:
        raise ValueError("empty quantile")
    if len(xs) == 1:
        return xs[0]
    p = (len(xs) - 1) * q
    lo, hi = int(math.floor(p)), int(math.ceil(p))
    if lo == hi:
        return xs[lo]
    f = p - lo
    return xs[lo] * (1 - f) + xs[hi] * f


def _stats(xs: Sequence[float]) -> Dict[str, float]:
    return {
        "mean": statistics.fmean(xs),
        "std": statistics.pstdev(xs),
        "q05": _quantile(xs, .05),
        "q50": _quantile(xs, .50),
        "q95": _quantile(xs, .95),
        "min": min(xs),
        "max": max(xs),
    }


@dataclass(frozen=True)
class AstrometricState:
    source_id: str
    catalog: str
    ra_deg: float
    dec_deg: float
    ref_epoch_jyear: float
    parallax_mas: float
    pmra_masyr: float
    pmdec_masyr: float
    ra_error_mas: float
    dec_error_mas: float
    parallax_error_mas: float
    pmra_error_masyr: float
    pmdec_error_masyr: float
    correlations: Mapping[str, float] | None = None

    def validate(self) -> None:
        for name, value in asdict(self).items():
            if name in {"source_id", "catalog", "correlations"}:
                continue
            if not math.isfinite(float(value)):
                raise ValueError(f"NONFINITE:{name}")
        for name in ("ra_error_mas", "dec_error_mas", "parallax_error_mas",
                     "pmra_error_masyr", "pmdec_error_masyr"):
            if getattr(self, name) < 0:
                raise ValueError(f"NEGATIVE_ERROR:{name}")
        if not (-90.0 <= self.dec_deg <= 90.0):
            raise ValueError("DEC_RANGE")


CORR_KEYS = {
    ("ra", "dec"): "ra_dec",
    ("ra", "parallax"): "ra_parallax",
    ("ra", "pmra"): "ra_pmra",
    ("ra", "pmdec"): "ra_pmdec",
    ("dec", "parallax"): "dec_parallax",
    ("dec", "pmra"): "dec_pmra",
    ("dec", "pmdec"): "dec_pmdec",
    ("parallax", "pmra"): "parallax_pmra",
    ("parallax", "pmdec"): "parallax_pmdec",
    ("pmra", "pmdec"): "pmra_pmdec",
}


def covariance_matrix(state: AstrometricState) -> tuple[List[List[float]], str]:
    state.validate()
    sigmas = [
        state.ra_error_mas,
        state.dec_error_mas,
        state.parallax_error_mas,
        state.pmra_error_masyr,
        state.pmdec_error_masyr,
    ]
    corr = dict(state.correlations or {})
    full = all(k in corr for k in CORR_KEYS.values())
    out = [[0.0] * 5 for _ in range(5)]
    for i, s in enumerate(sigmas):
        out[i][i] = s * s
    for i in range(5):
        for j in range(i + 1, 5):
            key = CORR_KEYS[(PARAM_ORDER[i], PARAM_ORDER[j])]
            rho = float(corr.get(key, 0.0))
            if not -1.0 <= rho <= 1.0:
                raise ValueError(f"CORRELATION_RANGE:{key}")
            out[i][j] = out[j][i] = rho * sigmas[i] * sigmas[j]
    return out, ("FULL_CATALOG_COVARIANCE" if full else "DIAGONAL_OR_PARTIAL_COVARIANCE")


def _cholesky_psd(a: Sequence[Sequence[float]]) -> List[List[float]]:
    n = len(a)
    l = [[0.0] * n for _ in range(n)]
    scale = max(1.0, max(abs(float(a[i][i])) for i in range(n)))
    tol = 1e-10 * scale
    for i in range(n):
        for j in range(i + 1):
            s = float(a[i][j]) - sum(l[i][k] * l[j][k] for k in range(j))
            if i == j:
                if s < -tol:
                    raise ValueError("COVARIANCE_NOT_PSD")
                l[i][j] = math.sqrt(max(0.0, s))
            else:
                l[i][j] = 0.0 if l[j][j] == 0.0 else s / l[j][j]
    return l


def _sample_state(state: AstrometricState, rng: random.Random) -> tuple[float, float, float, float, float]:
    cov, _ = covariance_matrix(state)
    l = _cholesky_psd(cov)
    z = [rng.gauss(0.0, 1.0) for _ in range(5)]
    d = [sum(l[i][j] * z[j] for j in range(i + 1)) for i in range(5)]
    return (
        d[0],
        d[1],
        state.parallax_mas + d[2],
        state.pmra_masyr + d[3],
        state.pmdec_masyr + d[4],
    )


def tangent_offset_mas(target_ra_deg: float, target_dec_deg: float,
                       source_ra_deg: float, source_dec_deg: float) -> tuple[float, float]:
    dra = (source_ra_deg - target_ra_deg + 180.0) % 360.0 - 180.0
    x = dra * math.cos(math.radians(target_dec_deg)) * 3_600_000.0
    y = (source_dec_deg - target_dec_deg) * 3_600_000.0
    return x, y


def min_linear_distance_in_window(
    x0_mas: float, y0_mas: float, vx_masyr: float, vy_masyr: float,
    ref_epoch: float, epoch_start: float, epoch_end: float,
) -> tuple[float, float]:
    if epoch_end < epoch_start:
        raise ValueError("EPOCH_WINDOW_REVERSED")
    vv = vx_masyr * vx_masyr + vy_masyr * vy_masyr
    if vv == 0.0:
        epoch = min(max(ref_epoch, epoch_start), epoch_end)
    else:
        dt_star = -(x0_mas * vx_masyr + y0_mas * vy_masyr) / vv
        epoch = min(max(ref_epoch + dt_star, epoch_start), epoch_end)
    dt = epoch - ref_epoch
    x = x0_mas + vx_masyr * dt
    y = y0_mas + vy_masyr * dt
    return math.hypot(x, y), epoch


def unconstrained_linear_closest_approach(
    x0_mas: float, y0_mas: float, vx_masyr: float, vy_masyr: float,
    ref_epoch: float,
) -> tuple[float, float]:
    vv = vx_masyr * vx_masyr + vy_masyr * vy_masyr
    if vv == 0:
        return math.hypot(x0_mas, y0_mas), ref_epoch
    dt = -(x0_mas * vx_masyr + y0_mas * vy_masyr) / vv
    return math.hypot(x0_mas + vx_masyr * dt, y0_mas + vy_masyr * dt), ref_epoch + dt


def evaluate_worldline(
    state: AstrometricState,
    target_ra_deg: float,
    target_dec_deg: float,
    samples: int = 8192,
    seed: int = 16016,
    epoch_start: float = DEFAULT_EPOCH_START,
    epoch_end: float = DEFAULT_EPOCH_END,
    gates_arcsec: Sequence[float] = DEFAULT_GATES_ARCSEC,
) -> Dict[str, object]:
    state.validate()
    if samples < 2:
        raise ValueError("SAMPLES_LT_2")
    x_nom, y_nom = tangent_offset_mas(target_ra_deg, target_dec_deg, state.ra_deg, state.dec_deg)
    nom_min_mas, nom_epoch = min_linear_distance_in_window(
        x_nom, y_nom, state.pmra_masyr, state.pmdec_masyr,
        state.ref_epoch_jyear, epoch_start, epoch_end,
    )
    uncon_mas, uncon_epoch = unconstrained_linear_closest_approach(
        x_nom, y_nom, state.pmra_masyr, state.pmdec_masyr, state.ref_epoch_jyear
    )
    cov, cov_status = covariance_matrix(state)
    rng = random.Random(seed)
    min_seps: List[float] = []
    parallax_lower_bounds: List[float] = []
    closest_epochs: List[float] = []
    for _ in range(samples):
        dx, dy, plx, pmra, pmdec = _sample_state(state, rng)
        min_mas, ep = min_linear_distance_in_window(
            x_nom + dx, y_nom + dy, pmra, pmdec,
            state.ref_epoch_jyear, epoch_start, epoch_end,
        )
        sep_as = min_mas / 1000.0
        lower_as = max(0.0, sep_as - abs(plx) / 1000.0)
        min_seps.append(sep_as)
        parallax_lower_bounds.append(lower_as)
        closest_epochs.append(ep)

    gate_counts = {
        str(float(g)): sum(1 for x in parallax_lower_bounds if x <= float(g))
        for g in gates_arcsec
    }
    result = {
        "schema": SCHEMA,
        "formula": "RESPICIENS_ET_PROSPICIENS",
        "source": {
            "source_id": state.source_id,
            "catalog": state.catalog,
            "ra_deg": state.ra_deg,
            "dec_deg": state.dec_deg,
            "ref_epoch_jyear": state.ref_epoch_jyear,
            "parallax_mas": state.parallax_mas,
            "pmra_masyr": state.pmra_masyr,
            "pmdec_masyr": state.pmdec_masyr,
        },
        "target": {"ra_deg": target_ra_deg, "dec_deg": target_dec_deg, "frame": "ICRS"},
        "epoch_window_jyear": [epoch_start, epoch_end],
        "samples": samples,
        "seed": seed,
        "covariance_status": cov_status,
        "covariance_parameter_order": [
            "ra_offset_mas", "dec_offset_mas", "parallax_mas", "pmra_masyr", "pmdec_masyr"
        ],
        "covariance_parameter_units": ["mas", "mas", "mas", "mas/yr", "mas/yr"],
        "measurement_uncertainty": {
            "ra_error_mas": state.ra_error_mas,
            "dec_error_mas": state.dec_error_mas,
            "parallax_error_mas": state.parallax_error_mas,
            "pmra_error_masyr": state.pmra_error_masyr,
            "pmdec_error_masyr": state.pmdec_error_masyr,
            "correlations": dict(state.correlations or {}),
        },
        "covariance_matrix_native_units": cov,
        "nominal": {
            "separation_at_ref_epoch_arcsec": math.hypot(x_nom, y_nom) / 1000.0,
            "minimum_linear_separation_in_window_arcsec": nom_min_mas / 1000.0,
            "closest_epoch_in_window_jyear": nom_epoch,
            "unconstrained_linear_closest_approach_arcsec": uncon_mas / 1000.0,
            "unconstrained_linear_closest_epoch_jyear": uncon_epoch,
        },
        "sampled_minimum_linear_separation_arcsec": _stats(min_seps),
        "sampled_parallax_envelope_lower_bound_arcsec": _stats(parallax_lower_bounds),
        "sampled_closest_epoch_jyear": _stats(closest_epochs),
        "gate_counts_after_conservative_parallax_envelope": gate_counts,
        "gate_fractions_after_conservative_parallax_envelope": {
            k: v / samples for k, v in gate_counts.items()
        },
        "parallax_policy": {
            "mode": "CONSERVATIVE_AMPLITUDE_ENVELOPE",
            "observer_specific_parallax_factors_used": False,
            "meaning": "A zero hit after subtracting the full sampled parallax amplitude excludes rescue by annual parallax within this linearized model.",
        },
        "epistemic_firewall": {
            "propagation_is_observation": False,
            "worldline_is_identity_proof": False,
            "catalog_covariance_must_not_be_invented": True,
            "missing_catalog_epoch_blocks_cross_catalog_epoch_match": True,
            "negative_result_is_valid": True,
        },
        "claim_ceiling": "BOUNDED_EPOCH_LINEAR_ASTROMETRIC_ROBUSTNESS_ONLY__NOT_OBJECT_IDENTITY",
    }
    result["input_sha256"] = _sha({
        "state": asdict(state), "target": result["target"], "epoch_window": result["epoch_window_jyear"],
        "samples": samples, "seed": seed, "gates": list(gates_arcsec),
    })
    result["result_sha256"] = _sha(result)
    return result


def self_test() -> None:
    corr = {k: 0.0 for k in CORR_KEYS.values()}
    s = AstrometricState(
        source_id="synthetic", catalog="GAIA_DR3",
        ra_deg=10.01, dec_deg=20.0, ref_epoch_jyear=2016.0,
        parallax_mas=1.0, pmra_masyr=-10.0, pmdec_masyr=0.0,
        ra_error_mas=.1, dec_error_mas=.1, parallax_error_mas=.1,
        pmra_error_masyr=.05, pmdec_error_masyr=.05, correlations=corr,
    )
    a = evaluate_worldline(s, 10.0, 20.0, samples=128, seed=42, epoch_start=2000, epoch_end=2030)
    b = evaluate_worldline(s, 10.0, 20.0, samples=128, seed=42, epoch_start=2000, epoch_end=2030)
    assert a["result_sha256"] == b["result_sha256"]
    assert a["covariance_status"] == "FULL_CATALOG_COVARIANCE"
    assert a["covariance_parameter_units"] == ["mas", "mas", "mas", "mas/yr", "mas/yr"]
    assert len(a["measurement_uncertainty"]["correlations"]) == 10
    assert a["epistemic_firewall"]["propagation_is_observation"] is False
    assert a["sampled_parallax_envelope_lower_bound_arcsec"]["min"] >= 0.0
    print("JANUS_ASTROMETRIC_WORLDLINES_SELF_TEST=PASS")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        self_test()
        return 0
    p.error("use this module through a target-specific runner")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
