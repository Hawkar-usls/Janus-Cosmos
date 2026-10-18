#!/usr/bin/env python3
"""BHB-FW1: AMPS quantum-horizon hypothesis overlay (not a firewall simulator).

This fail-closed epistemic layer binds an unchanged BHB-R2 Kerr receipt to a
machine-readable account of the conditional AMPS information paradox.
It produces no quantum stress-energy, scattering probability or emissions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from decimal import localcontext
from pathlib import Path
from typing import Any

from janus_bhb import BHBError, number, show
from janus_bhb_kerr import PRECISION, ROOT, geometry


FW_SCHEMA = "janus.cosmos.bhb.quantum-firewall-hypotheses.v1"
RESULT_SCHEMA = "janus.cosmos.bhb.firewall-overlay.v1"
DEFAULT_HYPOTHESES = ROOT / "data/bhb/BHB-FIREWALL-FW1-HYPOTHESES-2026-10-18.json"
DEFAULT_R2 = ROOT / "data/bhb/BHB-COSMOS-R2-RESULT-2026-10-18.json"
PREMISES = (
    "UNITARY_PURE_HAWKING_RADIATION",
    "SEMICLASSICAL_LOW_ENERGY_QFT_OUTSIDE_STRETCHED_HORIZON",
    "SMOOTH_HORIZON_NO_DRAMA_FOR_INFALLER",
)


def digest(document: Any) -> str:
    return hashlib.sha256(json.dumps(document, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()).hexdigest()


def evaluate_overlay(spec: dict[str, Any], r2: dict[str, Any], r2_bytes: bytes) -> dict[str, Any]:
    if not isinstance(spec, dict) or not isinstance(r2, dict):
        raise BHBError("spec and R2 must be JSON objects")
    if spec.get("schema") != FW_SCHEMA or spec.get("scientific_status") != "OPEN_UNCONFIRMED":
        raise BHBError("AMPS hypothesis must remain OPEN_UNCONFIRMED")
    if spec.get("physical_firewall_presence") is not None:
        raise BHBError("firewall presence must be null: no physical evidence gate")
    if spec.get("premises") != list(PREMISES):
        raise BHBError("AMPS three-premise contract was changed")
    if spec.get("argued_joint_compatibility") != "CONTRADICTORY_UNDER_AMPS_ARGUMENT":
        raise BHBError("invalid AMPS conditional logic")
    if spec.get("applicability_conditions") != ["EVAPORATING_BLACK_HOLE", "POST_PAGE_TIME_STATE", "ENTANGLEMENT_OF_EARLY_AND_LATE_RADIATION"]:
        raise BHBError("AMPS applicability prerequisites are incomplete")
    for lane in ("CLASSICAL_SMOOTH_HORIZON", "AMPS_FIREWALL_CANDIDATE", "NONLOCAL_OR_MICROSTATE_RESOLUTION"):
        if lane not in spec.get("competing_hypotheses", {}):
            raise BHBError(f"missing competing hypothesis: {lane}")
    if spec["competing_hypotheses"].get("AMPS_FIREWALL_CANDIDATE") != "UNCONFIRMED":
        raise BHBError("physical AMPS firewall cannot be promoted by configuration")
    if hashlib.sha256(r2_bytes).hexdigest() != spec.get("r2_receipt_sha256"):
        raise BHBError("R2 frozen receipt file SHA256 mismatch")
    if r2.get("schema") != "janus.cosmos.bhb.kerr-batch.v1" or r2.get("result_sha256") != spec.get("r2_internal_sha256"):
        raise BHBError("R2 frozen internal receipt or schema mismatch")
    if digest({k: v for k,v in r2.items() if k != "result_sha256"}) != r2["result_sha256"]:
        raise BHBError("R2 contents altered without recalculated internal hash")
    if r2.get("gates", {}).get("REAL_BHB_DETECTION") != "NOT_ESTABLISHED":
        raise BHBError("R2 observation claim cannot be promoted")
    if len(r2.get("cases", [])) != 7:
        raise BHBError("expected seven frozen R2 scenarios")

    rows = []
    with localcontext() as ctx:
        ctx.prec = PRECISION
        for case in r2["cases"]:
            if case.get("real_bhb_detection_gate") != "NOT_ESTABLISHED_NO_OBSERVATIONS":
                raise BHBError("R2 case observation gate mismatch")
            # This is only coordinate separation from the horizon, NOT
            # Planck distance, firewall thickness, stress tensor or cut-off.
            a = number(case["spin_a_over_M"], "spin")
            r = number(case["r_over_M"], "radius")
            rp = geometry(r, a)["r_plus"]
            row = {
                "scenario_id": case["scenario_id"],
                "r_over_M": show(r),
                "outer_horizon_r_plus_over_M": show(rp),
                "coordinate_gap_r_minus_rplus_over_M": show(r - rp),
                "coordinate_gap_is_not_proper_distance_or_firewall_thickness": True,
                "classical_R2_local_E_cm_GeV": case["collision"]["E_cm_GeV"],
                "classical_R2_trajectory_gate": case["trajectory_pair_gate"],
                "classical_R2_independent_test_photon_gate": case["test_photon"]["escape_gate"],
                "classical_R2_applicability": "CONDITIONAL_SMOOTH_STATIONARY_KERR_BACKGROUND",
                "black_hole_evaporation_state": "UNKNOWN_NO_AGE_INPUT",
                "page_time_relation": "UNKNOWN_NO_ENTANGLEMENT_DATA",
                "firewall_presence": "UNDETERMINED_NOT_MEASURED",
                "quantum_firewall_collision_E_cm_GeV": None,
                "quantum_firewall_energy_deposition_GeV": None,
                "quantum_firewall_photon_spectrum": None,
                "quantum_firewall_escaping_luminosity": None,
                "amplified_collider_boom_from_firewall": "NOT_INFERABLE",
                "observational_gate": "NO_FIREWALL_DETECTION_INPUT",
            }
            rows.append(row)

    out = {
        "schema": RESULT_SCHEMA,
        "experiment_id": spec["experiment_id"],
        "source_r2_receipt_sha256": spec["r2_receipt_sha256"],
        "source_r2_internal_sha256": spec["r2_internal_sha256"],
        "quantum_hypothesis": "AMPS_BLACK_HOLE_FIREWALL",
        "scientific_status": "UNRESOLVED_CONDITIONAL_INFORMATION_PARADOX",
        "premises": list(PREMISES),
        "hypothesis_comparison": spec["competing_hypotheses"],
        "cases": rows,
        "gates": {
            "CLASSICAL_R2_RECEIPT_INTEGRITY": "PASS",
            "AMPS_LOGICAL_PREMISES": "REGISTERED_CONDITIONAL_CONFLICT",
            "FIREWALL_PHYSICAL_EXISTENCE": "OPEN_UNCONFIRMED",
            "BLACK_HOLE_AGE_POST_PAGE_TIME": "UNKNOWN_NO_QUANTUM_STATE_DATA",
            "QUANTUM_FIREWALL_DYNAMICS": "NOT_MODELED",
            "ESCAPED_COLLISION_SIGNAL": "OPEN_NO_PRODUCT_MODEL",
            "REAL_BHB_OR_FIREWALL_DETECTION": "NOT_ESTABLISHED",
        },
        "claim_ceiling": "HYPOTHESIS_RECORD_ONLY_NO_FIREWALL_PREDICTION_OR_DETECTED_EXPLOSION",
    }
    out["result_sha256"] = digest(out)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hypotheses", type=Path, default=DEFAULT_HYPOTHESES)
    ap.add_argument("--r2", type=Path, default=DEFAULT_R2)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args(argv)
    try:
        raw = args.r2.read_bytes()
        spec = json.loads(args.hypotheses.read_text(encoding="utf-8"))
        result = evaluate_overlay(spec, json.loads(raw), raw)
        encoded = json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
        if args.output:
            args.output.write_text(encoded, encoding="utf-8")
        else:
            sys.stdout.write(encoded)
        return 0
    except (BHBError, ValueError, OSError, KeyError) as exc:
        print(f"BHB-FW1-ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
