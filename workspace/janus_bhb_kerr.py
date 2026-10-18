#!/usr/bin/env python3
"""BHB-R2: equatorial Kerr test-particle collisions and ZAMO photon escape.

Geometrized G=M=c=1 units; r is in gravitational radii GM/c².
Each photon is an independently specified *hypothetical* ZAMO-frame emission,
NOT a collision product. All calculations are exterior r>r_+ and equatorial.
No backreaction, stochastic emission, flux, detector, or astrophysical event.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any

from janus_bhb import BHBError, number, show, DEFAULT_MASS


SCHEMA = "janus.cosmos.bhb.kerr-equatorial.v1"
PRECISION = 90
ROOT = Path(__file__).resolve().parents[1]
EPS = Decimal("1e-58")


def minimum_massive_q(r: Decimal, a: Decimal, e: Decimal, ell: Decimal) -> Decimal:
    """Exact minimum of Q=R_timelike/r for all radii x>=r, e>=1.

    Q=(e²-1)x³+2x²+[a²(e²-1)-ell²]x+2(ell-ae)².
    Evaluating the positive derivative root certifies no radial forbidden gap.
    """
    c3 = e * e - 1
    c1 = a * a * c3 - ell * ell

    def q(x: Decimal) -> Decimal:
        return ((c3 * x + 2) * x + c1) * x + 2 * (ell - a * e) ** 2

    candidates = [q(r)]
    if c3 == 0:
        vertex = ell * ell / 4
        if vertex > r:
            candidates.append(q(vertex))
    else:
        disc = 16 - 12 * c3 * c1
        if disc > 0:
            vertex = (-4 + disc.sqrt()) / (6 * c3)
            if vertex > r:
                candidates.append(q(vertex))
    return min(candidates)


def minimum_photon_n(r: Decimal, a: Decimal, b: Decimal) -> Decimal:
    """Exact minimum of N=R_null/r for all x>=r, E_inf normalized to 1.

    N=x³+(a²-b²)x+2(b-a)². Positive minimum => outward ray reaches infinity.
    """
    def n(x: Decimal) -> Decimal:
        return x**3 + (a * a - b * b) * x + 2 * (b - a) ** 2

    candidates = [n(r)]
    if b * b > a * a:
        vertex = ((b * b - a * a) / 3).sqrt()
        if vertex > r:
            candidates.append(n(vertex))
    return min(candidates)


def geometry(r: Decimal, a: Decimal) -> dict[str, Decimal]:
    if abs(a) > 1:
        raise BHBError("spin a/M must satisfy -1 <= a/M <= 1")
    rp = 1 + (1 - a * a).sqrt()
    if r <= rp:
        raise BHBError("r/M must be strictly outside Kerr horizon r_+/M")
    delta = r * r - 2 * r + a * a
    area = (r * r + a * a) ** 2 - a * a * delta
    gphi = area / (r * r)
    grr = r * r / delta
    lapse = (delta * r * r / area).sqrt()
    omega = 2 * a * r / area
    return {"r_plus": rp, "delta": delta, "area": area, "gphi": gphi,
            "grr": grr, "lapse": lapse, "omega": omega,
            "gtt": -(1 - 2 / r), "gtphi": -2 * a / r}


def timelike(raw: Any, g: dict[str, Decimal], r: Decimal, a: Decimal, idx: int) -> dict[str, Any]:
    if not isinstance(raw, dict):
        raise BHBError(f"particle_{idx} must be object")
    mass = number(raw.get("rest_energy_GeV", DEFAULT_MASS), "rest_energy_GeV")
    e = number(raw.get("specific_energy", "1"), "specific_energy")
    ell = number(raw.get("specific_Lz_over_M", "0"), "specific_Lz_over_M")
    direction = raw.get("direction", "in")
    origin = raw.get("origin", "infinity")
    if mass <= 0 or e <= 0 or direction not in ("in", "out"):
        raise BHBError("positive mass/energy and in/out direction required")
    if origin not in ("infinity", "local_production", "unverified", "white_hole"):
        raise BHBError("invalid particle origin")
    p = e * (r * r + a * a) - a * ell
    big_r = p * p - g["delta"] * (r * r + (ell - a * e) ** 2)
    x = e - g["omega"] * ell
    local = big_r >= 0 and x > 0
    # Future-directed from infinity is a sufficient model-specific provenance
    # witness only for inward E>=1 geodesics with no forbidden radial gap.
    path_min = minimum_massive_q(r, a, e, ell) if e >= 1 else None
    if not local:
        status = "REJECTED_AT_COLLISION_EVENT"
    elif origin == "infinity" and direction == "in" and e >= 1 and path_min > 0:
        status = "INFINITY_INFALL_GEODESIC_ADMISSIBLE"
    elif origin == "infinity" and direction == "in" and e >= 1 and path_min <= 0:
        status = "REJECTED_FORBIDDEN_RADIAL_BARRIER_FROM_INFINITY"
    else:
        status = "PROVENANCE_NOT_ESTABLISHED"
    result = {"rest_energy_GeV": show(mass), "specific_energy": show(e),
              "specific_Lz_over_M": show(ell), "direction": direction,
              "origin": origin, "radial_R_over_M4": show(big_r),
              "radial_path_min_Q": show(path_min) if path_min is not None else None,
              "trajectory_gate": status}
    if local:
        ur = (1 if direction == "out" else -1) * big_r.sqrt() / (r * r)
        ut = (a * (ell - a * e) + (r * r + a * a) * p / g["delta"]) / (r * r)
        uphi = (ell - a * e + a * p / g["delta"]) / (r * r)
        gamma = x / g["lapse"]
        vr = g["grr"].sqrt() * ur / gamma
        vp = ell / (g["gphi"].sqrt() * gamma)
        norm = gamma * gamma * (1 - vr * vr - vp * vp)
        if abs(norm - 1) > EPS:
            raise BHBError("ZAMO timelike norm failed")
        result.update({"zamo_gamma": show(gamma), "zamo_vr_over_c": show(vr),
                       "zamo_vphi_over_c": show(vp)})
        result["_v"] = (ut, ur, uphi, gamma, vr, vp)
        result["_m"] = mass
    return result


def collision(g: dict[str, Decimal], p1: dict, p2: dict) -> dict[str, Any]:
    if "_v" not in p1 or "_v" not in p2:
        return {"local_collision_gate": "REJECTED_NONREAL_OR_PAST_DIRECTED_WORLDLINE",
                "E_cm_GeV": None, "tetrad_cross_check": "NOT_APPLICABLE"}
    t1, r1, f1, gam1, vr1, vp1 = p1["_v"]
    t2, r2, f2, gam2, vr2, vp2 = p2["_v"]
    rel = gam1 * gam2 * (1 - vr1 * vr2 - vp1 * vp2)
    coord_dot = (g["gtt"] * t1 * t2 + g["gtphi"] * (t1 * f2 + f1 * t2)
                 + g["gphi"] * f1 * f2 + g["grr"] * r1 * r2)
    residual = abs(rel + coord_dot)
    if residual > EPS * max(Decimal(1), rel):
        raise BHBError("coordinate versus ZAMO invariant failed")
    if rel < 1 - EPS:
        raise BHBError("relative gamma less than 1")
    if rel < 1:
        rel = Decimal(1)
    m1, m2 = p1["_m"], p2["_m"]
    ecm = (m1 * m1 + m2 * m2 + 2 * m1 * m2 * rel).sqrt()
    return {"local_collision_gate": "PASS_LOCAL_KERR_GEODESIC_INVARIANT",
            "relative_gamma": show(rel), "E_cm_GeV": show(ecm),
            "tetrad_cross_check": "PASS", "coordinate_tetrad_residual": show(residual, 8)}


def photon(g: dict[str, Decimal], r: Decimal, a: Decimal, raw: Any) -> dict[str, Any]:
    if not isinstance(raw, dict):
        raise BHBError("photon must specify ZAMO local emission direction")
    if raw.get("source") != "hypothetical_zamo_emission":
        raise BHBError("photon source must be hypothetical_zamo_emission")
    nphi = number(raw.get("azimuth_sine"), "photon.azimuth_sine")
    if abs(nphi) >= 1:
        raise BHBError("outward radial cosine must be >0: |azimuth_sine|<1")
    nr = (1 - nphi * nphi).sqrt()
    sqrt_gphi = g["gphi"].sqrt()
    shift = g["lapse"] + g["omega"] * sqrt_gphi * nphi
    if shift <= 0:
        return {"escape_gate": "NO_POSITIVE_INFINITY_ENERGY", "E_inf_over_E_zamo": show(shift)}
    b = sqrt_gphi * nphi / shift
    min_n = minimum_photon_n(r, a, b)
    rnorm = (r * r + a * a - a * b) ** 2 - g["delta"] * (b - a) ** 2
    # At the emission point the null radial potential must agree with the
    # observer's orthonormal radial momentum, independent of polynomial test.
    null_residual = abs(rnorm * shift * shift - g["delta"] * r * r * nr * nr)
    if null_residual > EPS * max(Decimal(1), rnorm * shift * shift):
        raise BHBError("null radial ZAMO frame cross-check failed")
    return {"source": "HYPOTHETICAL_NOT_DERIVED_FROM_COLLISION",
            "azimuth_sine": show(nphi), "outward_radial_cosine": show(nr),
            "E_inf_over_E_zamo": show(shift), "Lz_over_M_Einf": show(b),
            "min_null_radial_N_to_infinity": show(min_n),
            "null_frame_cross_check": "PASS",
            "escape_gate": ("EQUATORIAL_NULL_GEODESIC_REACHES_INFINITY"
                            if min_n > 0 else "BLOCKED_OR_CRITICAL_RADIAL_POTENTIAL"),
            "collision_product_gate": "NOT_MODELED_NO_EMISSION_SPECTRUM_OR_YIELD"}


def evaluate_case(raw: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(raw, dict) or raw.get("spacetime", "kerr") != "kerr":
        raise BHBError("R2 requires Kerr spacetime object")
    with localcontext() as ctx:
        ctx.prec = PRECISION
        a = number(raw.get("spin_a_over_M"), "spin_a_over_M")
        r = number(raw.get("r_over_M"), "r_over_M")
        g = geometry(r, a)
        p1 = timelike(raw.get("particle_1"), g, r, a, 1)
        p2 = timelike(raw.get("particle_2"), g, r, a, 2)
        coll = collision(g, p1, p2)
        ray = photon(g, r, a, raw.get("photon"))
        for p in (p1, p2):
            p.pop("_v", None)
            p.pop("_m", None)
        source = "PASS" if all(p["trajectory_gate"] == "INFINITY_INFALL_GEODESIC_ADMISSIBLE" for p in (p1,p2)) else "NOT_ESTABLISHED_OR_REJECTED"
        row = {"schema": SCHEMA, "scenario_id": raw.get("scenario_id", "BHB-R2-CASE"),
               "model_class": "EQUATORIAL_STATIONARY_KERR_TEST_GEODESICS_ONLY",
               "spin_a_over_M": show(a), "r_over_M": show(r), "outer_horizon_r_plus_over_M": show(g["r_plus"]),
               "zamo_lapse": show(g["lapse"]), "frame_drag_omega_M": show(g["omega"]),
               "particle_1": p1, "particle_2": p2,
               "trajectory_pair_gate": source, "collision": coll, "test_photon": ray,
               "escaped_collision_signal_gate": "OPEN_NO_COLLISION_PRODUCT_FOUR_MOMENTUM",
               "real_bhb_detection_gate": "NOT_ESTABLISHED_NO_OBSERVATIONS",
               "claim_ceiling": "LOCAL_KINEMATICS_AND_HYPOTHETICAL_RAY_ESCAPE_ONLY"}
        row["result_sha256"] = hashlib.sha256(json.dumps(row,sort_keys=True,separators=(",", ":")).encode()).hexdigest()
        return row


def evaluate_batch(raw: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(raw, dict) or not isinstance(raw.get("cases"),list) or not 1 <= len(raw["cases"]) <= 1000:
        raise BHBError("batch requires 1..1000 cases")
    rows = [evaluate_case(c) for c in raw["cases"]]
    out = {"schema": "janus.cosmos.bhb.kerr-batch.v1", "experiment_id": raw.get("experiment_id","BHB-COSMOS-R2"),
           "status": "REPRODUCIBLE_MODEL_ONLY", "precision_digits": PRECISION,
           "cases": rows,
           "gates": {"KERR_EQUATORIAL_INVARIANT": "PASS" if all(r["collision"]["tetrad_cross_check"] in ("PASS","NOT_APPLICABLE") for r in rows) else "FAIL",
                     "RADIAL_TRAJECTORY_ACCESS": "MIXED_CHECK_PER_CASE",
                     "HYPOTHETICAL_EQUATORIAL_PHOTON_TRANSPORT": "CHECKED_PER_CASE",
                     "ESCAPED_COLLISION_SIGNAL": "OPEN_NO_PRODUCT_MODEL",
                     "REAL_BHB_DETECTION": "NOT_ESTABLISHED"},
           "sources": ["arXiv:0909.0169", "arXiv:0911.3363", "arXiv:2104.07889"]}
    out["result_sha256"] = hashlib.sha256(json.dumps(out,sort_keys=True,separators=(",", ":")).encode()).hexdigest()
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--batch", type=Path, required=True)
    ap.add_argument("--output", type=Path, help="write reproducible JSON receipt")
    args = ap.parse_args(argv)
    try:
        doc = evaluate_batch(json.loads(args.batch.read_text(encoding="utf-8")))
        serialized = json.dumps(doc, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
        if args.output:
            args.output.write_text(serialized, encoding="utf-8")
        else:
            sys.stdout.write(serialized)
        return 0
    except (BHBError, ValueError, OSError, KeyError) as exc:
        print(f"BHB-R2-ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
