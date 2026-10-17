#!/usr/bin/env python3
"""JANUS Cosmos BHB: provenance-gated Schwarzschild collision kinematics.

Computes local invariant centre-of-mass energy for two *radial* massive
test-particle geodesics at the same event, outside a nonrotating horizon.
All dimensionless radii are r/r_s, conserved specific energies are E/m,
and particle rest energies are in GeV. No assumption about radiation
production, escaping luminosity, merger rates, or Kerr geometry is made.

The radial geodesic equations and orthonormal static tetrad give
 gamma_i = e_i/sqrt(f),  beta_i = sign_i sqrt(1-f/e_i**2),
 f=1-r_s/r, and E_cm**2=m1**2+m2**2+2*m1*m2*gamma_relative.
An algebraically rationalized expression prevents catastrophic
cancellation for particles travelling in the same radial direction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from decimal import Decimal, InvalidOperation, localcontext
from pathlib import Path
from typing import Any


SCHEMA = "janus.cosmos.bhb.schwarzschild-radial.v1"
ROOT = Path(__file__).resolve().parents[1]
FROZEN_UPSTREAM = ROOT / "data/bhb/BHB-R0-UPSTREAM-2026-10-17.json"
FROZEN_SHA256 = "3fa8b5027edd79a9a9e5d1d7561d7ca029debc7f6656142a34167842ef76f067"
UPSTREAM_GIT_COMMIT = "5b41a09c0e7885a649610a9e3b3ffe2a1c4778f2"
UPSTREAM_BLOB_SHA = "10ca85cab017a936ef5d3511d932dcc85cbb7696"
DEFAULT_MASS = "0.938272"
DECIMAL_PRECISION = 90
VALID_DIRECTIONS = {"in": -1, "out": 1}
VALID_ORIGINS = {"infinity", "white_hole", "local_production", "unverified"}
VALID_CONTEXTS = {"collapsing_black_hole", "eternal_schwarzschild"}


class BHBError(ValueError):
    """Invalid input, malformed provenance, or failed evidence replay."""


def number(value: Any, label: str) -> Decimal:
    if isinstance(value, bool) or value is None:
        raise BHBError(f"{label} must be a finite decimal number")
    try:
        n = Decimal(str(value))
    except (InvalidOperation, ValueError) as exc:
        raise BHBError(f"{label} must be a finite decimal number") from exc
    if not n.is_finite():
        raise BHBError(f"{label} must be finite")
    return n


def show(value: Decimal, sig: int = 24) -> str:
    return format(value, f".{sig}g")


def particle(raw: Any, f: Decimal, index: int) -> dict[str, Any]:
    if not isinstance(raw, dict):
        raise BHBError(f"particle_{index} must be an object")
    mass = number(raw.get("rest_energy_GeV", DEFAULT_MASS), f"particle_{index}.rest_energy_GeV")
    e = number(raw.get("specific_energy", "1"), f"particle_{index}.specific_energy")
    direction = raw.get("direction", "in")
    origin = raw.get("origin", "unverified")
    if mass <= 0 or e <= 0 or e * e < f:
        raise BHBError(f"particle_{index}: positive mass, e>0, e^2>=f required")
    if direction not in VALID_DIRECTIONS or origin not in VALID_ORIGINS:
        raise BHBError(f"particle_{index}: invalid direction or origin")
    beta_abs = (1 - f / (e * e)).sqrt()
    beta = beta_abs * VALID_DIRECTIONS[direction]
    return {
        "m": mass, "e": e, "beta": beta, "beta_abs": beta_abs,
        "direction": direction, "origin": origin,
        "root": (e * e - f).sqrt(),
    }


def source_gate(context: str, p1: dict, p2: dict) -> tuple[str, list[str]]:
    reasons: list[str] = []
    parts = (p1, p2)
    if any(p["origin"] == "infinity" and p["e"] < 1 for p in parts):
        reasons.append("INFINITY_ORIGIN_REQUIRES_SPECIFIC_ENERGY_AT_LEAST_ONE")
    if any(p["direction"] == "out" and p["origin"] == "infinity" for p in parts):
        reasons.append("RADIAL_OUTGOING_GEODESIC_CANNOT_ARRIVE_FROM_INFINITY_NEAR_HORIZON")
    if any(p["direction"] == "out" and p["origin"] == "unverified" for p in parts):
        reasons.append("OUTGOING_SOURCE_UNVERIFIED")
    if any(p["direction"] == "out" and p["origin"] == "white_hole" for p in parts):
        if context == "collapsing_black_hole":
            reasons.append("WHITE_HOLE_WORLDLINE_NOT_IN_COLLAPSING_BH_SPACETIME")
    if any(p["origin"] == "local_production" for p in parts):
        reasons.append("LOCAL_PRODUCTION_REQUIRES_INDEPENDENT_ENERGY_BUDGET_AND_CAUSAL_SOURCE")
    if any(p["origin"] == "unverified" for p in parts):
        reasons.append("PARTICLE_ORIGIN_UNVERIFIED")
    if reasons:
        return "NOT_ESTABLISHED_OR_INCOMPATIBLE", sorted(set(reasons))
    if any(p["direction"] == "out" for p in parts):
        return "CONDITIONAL_ETERNAL_WHITE_HOLE_ONLY", ["NOT_AN_ASTROPHYSICAL_COLLAPSING_BLACK_HOLE_PREDICTION"]
    if any(p["origin"] == "white_hole" for p in parts):
        return "NOT_ESTABLISHED_OR_INCOMPATIBLE", ["WHITE_HOLE_ORIGIN_FOR_INGOING_STREAM_NOT_VALIDATED"]
    return "IDEALIZED_GEODESIC_KINEMATICS_ADMISSIBLE", ["JOINT_SPACETIME_EVENT_AND_COLLISION_RATE_NOT_MODELED"]


def evaluate_case(raw: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(raw, dict):
        raise BHBError("case must be an object")
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PRECISION
        geometry = raw.get("spacetime", "schwarzschild")
        context = raw.get("physical_context", "collapsing_black_hole")
        if geometry != "schwarzschild" or context not in VALID_CONTEXTS:
            raise BHBError("only Schwarzschild and known physical_context are supported")
        x = number(raw.get("r_over_rs"), "r_over_rs")
        if x <= 1:
            raise BHBError("r/r_s must be > 1; static tetrad undefined at horizon")
        f = (x - 1) / x
        p1 = particle(raw.get("particle_1"), f, 1)
        p2 = particle(raw.get("particle_2"), f, 2)
        prod = p1["e"] * p2["e"]
        roots = p1["root"] * p2["root"]
        same = p1["direction"] == p2["direction"]
        # Rationalized same-direction relative Lorentz factor; no 0/0 at f>0.
        relative_gamma = ((p1["e"]**2 + p2["e"]**2 - f) / (prod + roots)
                          if same else (prod + roots) / f)
        if relative_gamma < 1:
            if 1 - relative_gamma < Decimal("1e-75"):
                relative_gamma = Decimal(1)
            else:
                raise BHBError("relative gamma below 1: invalid physical configuration")
        m1, m2 = p1["m"], p2["m"]
        ecm2 = m1**2 + m2**2 + 2*m1*m2*relative_gamma
        ecm = ecm2.sqrt()
        # Independent 1+1 Minkowski invariant in local static tetrad.
        gamma1 = p1["e"] / f.sqrt()
        gamma2 = p2["e"] / f.sqrt()
        energy = m1*gamma1 + m2*gamma2
        momentum = m1*gamma1*p1["beta"] + m2*gamma2*p2["beta"]
        tetrad_ecm2 = energy*energy - momentum*momentum
        err = abs(tetrad_ecm2 - ecm2)
        if err > Decimal("1e-54") * max(Decimal(1), ecm2):
            raise BHBError("local tetrad invariant cross-check failed")
        gate, reasons = source_gate(context, p1, p2)
        row = {
            "scenario_id": raw.get("scenario_id", "BHB-CASE"),
            "schema": SCHEMA,
            "model_class": "MODEL_ONLY_SCHWARZSCHILD_TEST_PARTICLE_GEODESICS",
            "physical_context": context,
            "r_over_rs": str(x),
            "lapse_squared": show(f),
            "particle_1": {"direction": p1["direction"], "origin": p1["origin"], "rest_energy_GeV": show(m1), "specific_energy": show(p1["e"]), "beta_over_c_signed": show(p1["beta"]), "one_minus_abs_beta": show(1-p1["beta_abs"])},
            "particle_2": {"direction": p2["direction"], "origin": p2["origin"], "rest_energy_GeV": show(m2), "specific_energy": show(p2["e"]), "beta_over_c_signed": show(p2["beta"]), "one_minus_abs_beta": show(1-p2["beta_abs"])},
            "relative_gamma": show(relative_gamma),
            "E_cm_GeV": show(ecm),
            "E_cm_squared_GeV2": show(ecm2),
            "kinematics_cross_check": "PASS",
            "tetrad_invariant_absolute_residual_GeV2": show(err, 8),
            "trajectory_source_gate": gate,
            "trajectory_source_reasons": reasons,
            "escaped_energy_gate": "OPEN_NO_RAY_TRACING_OR_EMISSION_MODEL",
            "observation_gate": "NO_OBSERVATION_INPUT",
            "claim_ceiling": "LOCAL_KINEMATIC_INVARIANT_ONLY__NO_EXPLOSION_OR_OBSERVABLE_ENERGY_CLAIM",
        }
        row["result_sha256"] = hashlib.sha256(json.dumps(row, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()
        return row


def evaluate_batch(document: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(document, dict) or not isinstance(document.get("cases"), list):
        raise BHBError("batch must be an object with cases list")
    if not document["cases"] or len(document["cases"]) > 1000:
        raise BHBError("cases must contain 1..1000 entries")
    rows = [evaluate_case(case) for case in document["cases"]]
    result = {
        "schema": "janus.cosmos.bhb.batch.v1",
        "experiment_id": document.get("experiment_id", "BHB-COSMOS-R1"),
        "status": "REPRODUCIBLE_MODEL_ONLY",
        "model_precision_digits": DECIMAL_PRECISION,
        "cases": rows,
        "gates": {
            "SCHWARZSCHILD_LOCAL_INVARIANT": "PASS" if all(r["kinematics_cross_check"] == "PASS" for r in rows) else "FAIL",
            "ALL_TRAJECTORIES_ADMISSIBLE": "PASS_CONDITIONAL" if all(r["trajectory_source_gate"] == "IDEALIZED_GEODESIC_KINEMATICS_ADMISSIBLE" for r in rows) else "NOT_ESTABLISHED",
            "ESCAPING_SIGNAL": "OPEN",
            "KERR_ROTATION": "NOT_MODELED",
            "REAL_BHB_DETECTION": "NOT_ESTABLISHED",
        },
        "upstream": {"repo": "Hawkar-usls/janus-meta-registry", "commit": UPSTREAM_GIT_COMMIT, "blob_sha": UPSTREAM_BLOB_SHA},
    }
    result["result_sha256"] = hashlib.sha256(json.dumps(result, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()
    return result


def replay_upstream(source_path: Path = FROZEN_UPSTREAM) -> dict[str, Any]:
    raw = source_path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != FROZEN_SHA256:
        raise BHBError("UPSTREAM_FIXTURE_SHA256_MISMATCH")
    upstream = json.loads(raw)
    if upstream.get("artifact_uuid") != "BHB-BLACK-HOLE-BOOM-R0-2026-10-17-v1.0":
        raise BHBError("UPSTREAM_ARTIFACT_ID_MISMATCH")
    if upstream.get("provenance", {}).get("source_code_sha256") != "cfdcf3a1ed88a7f0a6a59bd00e818fb1a18fe9a447a78ee64f5988aa59158d84":
        raise BHBError("UPSTREAM_CODE_ID_MISMATCH")
    checked = []
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PRECISION
        for source in upstream.get("rows", []):
            x = source["r_over_rs"]
            p_in = {"direction": "in", "origin": "infinity"}
            p_out = {"direction": "out", "origin": "white_hole"}
            same = evaluate_case({"scenario_id": "R0-CONTROL", "r_over_rs": x, "particle_1": p_in, "particle_2": p_in})
            opposite = evaluate_case({"scenario_id": "R0-COUNTERFACTUAL", "r_over_rs": x, "physical_context": "eternal_schwarzschild", "particle_1": p_in, "particle_2": p_out})
            for key, value in [("infall_infall_Ecm_GeV", same["E_cm_GeV"]), ("infall_outgoing_Ecm_GeV", opposite["E_cm_GeV"])]:
                expected = number(source[key], key)
                observed = number(value, key)
                if abs(observed-expected) > Decimal("5e-11") * max(expected, Decimal(1)):
                    raise BHBError(f"UPSTREAM_REPLAY_NUMERIC_MISMATCH_{x}_{key}")
            if same["trajectory_source_gate"] != "IDEALIZED_GEODESIC_KINEMATICS_ADMISSIBLE":
                raise BHBError("R0_CONTROL_PROVENANCE_CHANGED")
            if opposite["trajectory_source_gate"] != "CONDITIONAL_ETERNAL_WHITE_HOLE_ONLY":
                raise BHBError("R0_WHITE_HOLE_FIREWALL_CHANGED")
            checked.append(x)
    if len(checked) != 7:
        raise BHBError("UPSTREAM_EXPECTED_SEVEN_ROWS")
    return {
        "schema": "janus.cosmos.bhb.replay-receipt.v1",
        "status": "PASS_FROZEN_R0_INDEPENDENT_REPLAY",
        "source_repo": "Hawkar-usls/janus-meta-registry",
        "source_commit": UPSTREAM_GIT_COMMIT,
        "source_git_blob_sha": UPSTREAM_BLOB_SHA,
        "source_sha256": digest,
        "checked_rows": checked,
        "checked_scenarios": 2 * len(checked),
        "precision_digits": DECIMAL_PRECISION,
        "claim_ceiling": "REPLAY_OF_ANALYTIC_MODEL_NOT_PHYSICAL_DETECTION",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--batch", type=Path, help="JSON file containing {'cases': [...]} to evaluate")
    group.add_argument("--replay-upstream", action="store_true", help="replay the pinned Meta Registry BHB R0 fixture")
    group.add_argument("--self-test", action="store_true", help="run built-in frozen replay")
    args = parser.parse_args(argv)
    try:
        if args.batch:
            payload = json.loads(args.batch.read_text(encoding="utf-8"))
            result = evaluate_batch(payload)
        else:
            result = replay_upstream()
            if args.self_test:
                result["self_test"] = "PASS"
        print(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True))
        return 0
    except (BHBError, OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
        print(f"BHB_FAIL_CLOSED: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
