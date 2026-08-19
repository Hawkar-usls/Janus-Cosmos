#!/usr/bin/env python3
"""OSIRIS v3 spiral compute runtime.

Modern software mechanic inspired by the frozen OSIRIS semantic specification:

    ORIGIN_n -> COMPUTE/EXPERIENCE -> VERIFIED_RETURN -> ORIGIN_PRIME_(n+1)

Unlike the earlier passive ribbon-state receipt, ORIGIN_PRIME experience may now
participate in the next computation. Reuse is deliberately narrow and
proof-carrying:

* a cached SAT assignment is useful only after it satisfies the exact current
  canonical CNF again;
* a cached separator route is useful only after the separator and component
  partition are revalidated on the exact current residual CNF;
* cached UNSAT is never a verdict shortcut. Current boundary/component closure
  still executes.

Memory can alter search/repeated-discovery cost. It cannot manufacture
SAT/UNSAT authority. P_VS_NP remains OPEN.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterable, Mapping

ROOT = Path(__file__).resolve().parents[1]
DIRECT = ROOT / "experiments" / "direct"
if str(DIRECT) not in sys.path:
    sys.path.insert(0, str(DIRECT))

import s_phallus_h_gate_2_bounded_k_scaling_holdout_budget_guard as guarded  # noqa: E402
from janus_c025_core import canonical_cnf, cnf_hash, satisfies  # noqa: E402

guarded.install_budget_guard()
gate = guarded.gate

CONTRACT_PATH = Path(__file__).with_name("OSIRIS_V3_SPIRAL_COMPUTE_FROZEN_CONTRACT.json")
SEMANTIC_PATH = Path(__file__).with_name("OSIRIS_ORIGIN_PRIME_SPIRAL_SEMANTIC_SPEC_v1.0.json")
EXPECTED_SEMANTIC_BLOB = "b3163013ff7b0c22a86dff732092374a191e5a2a"
COSMOS_PARENT_SHA = "c77f920d764229efb6932bc4ea522a4ec0342c64"
FUNDAMENTUM_PIN = "00b09d778a6d57fc7f905df9b6235fb30e29c5a3"
STORE_SCHEMA = "janus.cosmos.osiris_spiral_state_store.v1"
STATE_SCHEMA = "janus.cosmos.osiris_origin_prime_state.v1"
EXPERIENCE_SCHEMA = "janus.cosmos.osiris_formula_experience.v1"
RUNTIME_ID = "OSIRIS-V3-ORIGIN-PRIME-SPIRAL-COMPUTE-2026-08-19-v1"
SPIRAL_SAT_LANE = "OSIRIS_V3_SPIRAL_REVERIFIED_SAT_WITNESS"
SPIRAL_ROUTE_LANE = "OSIRIS_V3_SPIRAL_REVALIDATED_SEPARATOR_ROUTE"


def stable_bytes(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def digest(value: Any) -> str:
    return hashlib.sha256(stable_bytes(value)).hexdigest()


def git_blob_sha(path: Path) -> str:
    raw = path.read_bytes()
    header = f"blob {len(raw)}\0".encode("ascii")
    return hashlib.sha1(header + raw).hexdigest()


def mutated_hex(value: str) -> str:
    if not value:
        return value
    return ("0" if value[0] != "0" else "1") + value[1:]


def provider_identity() -> dict[str, Any]:
    return {
        "cosmos_parent_sha": COSMOS_PARENT_SHA,
        "fundamentum_pin": FUNDAMENTUM_PIN,
        "parent_gate": "S𓂸ḥ/2",
        "semantic_spec_blob": EXPECTED_SEMANTIC_BLOB,
    }


def verify_provider_checkout() -> bool:
    try:
        actual = subprocess.check_output(
            ["git", "-C", str(ROOT / "vendor" / "Janus-Fundamentum"), "rev-parse", "HEAD"],
            text=True,
        ).strip()
    except Exception:
        return False
    return actual == FUNDAMENTUM_PIN


def load_frozen_inputs() -> tuple[dict[str, Any], dict[str, Any]]:
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    semantic = json.loads(SEMANTIC_PATH.read_text(encoding="utf-8"))
    assert contract["status"] == "FROZEN_BEFORE_IMPLEMENTATION_AND_RUN"
    assert contract["parent"]["sha"] == COSMOS_PARENT_SHA
    assert contract["parent"]["fundamentum_pin"] == FUNDAMENTUM_PIN
    assert contract["semantic_spec"]["git_blob_sha"] == EXPECTED_SEMANTIC_BLOB
    assert git_blob_sha(SEMANTIC_PATH) == EXPECTED_SEMANTIC_BLOB
    assert semantic["status"] == "FROZEN_BEFORE_SPIRAL_COMPUTE_IMPLEMENTATION"
    assert semantic["canonical_path"] == ["ORIGIN", "EXPERIENCE", "RETURN", "ORIGIN_PRIME"]
    assert semantic["state_laws"]["NEXT_GENERATION_MAY_USE_VERIFIED_EXPERIENCE"] is True
    assert semantic["experience_authority_firewall"]["MEMORY_NE_TRUTH"] is True
    return contract, semantic


def committed(core: Mapping[str, Any], field: str) -> dict[str, Any]:
    out = dict(core)
    out[field] = digest(core)
    return out


def genesis_state() -> dict[str, Any]:
    core = {
        "schema": STATE_SCHEMA,
        "state_type": "ORIGIN",
        "generation": 0,
        "previous_state_commitment": None,
        "position_commitment": None,
        "experience_commitment": None,
        "path_history_digest": None,
        "return_commitment": None,
        "provider": provider_identity(),
    }
    return committed(core, "state_commitment")


def verify_state_commitment(state: Mapping[str, Any]) -> bool:
    if not isinstance(state, Mapping) or not isinstance(state.get("state_commitment"), str):
        return False
    core = dict(state)
    got = core.pop("state_commitment", None)
    return got == digest(core)


def experience_key(formula_hash: str, budget: int) -> str:
    return f"{formula_hash}:{int(budget)}:{FUNDAMENTUM_PIN}"


def empty_store() -> dict[str, Any]:
    origin = genesis_state()
    return {
        "schema": STORE_SCHEMA,
        "provider": provider_identity(),
        "origin_state": origin,
        "current_state": origin,
        "history": [],
        "experiences": {},
    }


def verify_store(store: Mapping[str, Any]) -> tuple[bool, str | None]:
    if not isinstance(store, Mapping) or store.get("schema") != STORE_SCHEMA:
        return False, "STORE_SCHEMA_MISMATCH"
    if store.get("provider") != provider_identity():
        return False, "STORE_PROVIDER_MISMATCH"
    origin = store.get("origin_state")
    current = store.get("current_state")
    history = store.get("history")
    experiences = store.get("experiences")
    if not isinstance(origin, Mapping) or not verify_state_commitment(origin):
        return False, "ORIGIN_COMMITMENT_INVALID"
    if origin.get("state_type") != "ORIGIN" or origin.get("generation") != 0:
        return False, "ORIGIN_TYPE_OR_GENERATION_INVALID"
    if not isinstance(history, list) or not isinstance(experiences, Mapping):
        return False, "STORE_COLLECTION_INVALID"
    previous = origin["state_commitment"]
    expected_generation = 1
    for state in history:
        if not verify_state_commitment(state):
            return False, "HISTORY_STATE_COMMITMENT_INVALID"
        if state.get("state_type") != "ORIGIN_PRIME":
            return False, "HISTORY_STATE_TYPE_INVALID"
        if state.get("generation") != expected_generation:
            return False, "HISTORY_GENERATION_GAP"
        if state.get("previous_state_commitment") != previous:
            return False, "HISTORY_LINEAGE_BROKEN"
        previous = state["state_commitment"]
        expected_generation += 1
    if history:
        if current != history[-1]:
            return False, "CURRENT_STATE_NOT_HISTORY_TIP"
    else:
        if current != origin:
            return False, "CURRENT_STATE_NOT_ORIGIN"
    for key, record in experiences.items():
        if not isinstance(record, Mapping) or not verify_experience_commitment(record):
            return False, f"EXPERIENCE_COMMITMENT_INVALID:{key}"
    return True, None


def load_store(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    if not path.exists():
        return empty_store(), {"memory_status": "NEW_ORIGIN", "rejected_memory_digest": None, "reason": None}
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return empty_store(), {
            "memory_status": "REJECTED_FALLBACK_FRESH_ORIGIN",
            "rejected_memory_digest": hashlib.sha256(path.read_bytes()).hexdigest(),
            "reason": "STORE_PARSE_FAILURE",
        }
    ok, reason = verify_store(raw)
    if ok:
        return dict(raw), {"memory_status": "VALID", "rejected_memory_digest": None, "reason": None}
    return empty_store(), {
        "memory_status": "REJECTED_FALLBACK_FRESH_ORIGIN",
        "rejected_memory_digest": digest(raw),
        "reason": reason,
    }


def save_store(path: Path, store: Mapping[str, Any]) -> None:
    ok, reason = verify_store(store)
    if not ok:
        raise ValueError(f"refusing to save invalid spiral store: {reason}")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(store, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def assignment_to_json(assignment: Mapping[int, bool] | None) -> dict[str, bool] | None:
    if assignment is None:
        return None
    return {str(int(k)): bool(v) for k, v in sorted(assignment.items())}


def assignment_from_json(assignment: Mapping[str, Any] | None) -> dict[int, bool] | None:
    if assignment is None:
        return None
    return {int(k): bool(v) for k, v in assignment.items()}


def make_experience(formula: Any, budget: int, solved: Mapping[str, Any], generation: int) -> dict[str, Any] | None:
    formula = canonical_cnf(formula)
    if not solved.get("authorized") or solved.get("status") not in {"SAT", "UNSAT"}:
        return None
    residual = guarded.residual_after_parent_preprocessing(formula)
    primary = solved.get("primary") or {}
    detector = primary.get("detector") or {}
    separator = detector.get("separator")
    components = detector.get("components") or []
    route_reusable = False
    if separator:
        claim = gate.verify_separator_claim(residual, separator, components)
        route_reusable = bool(claim.get("passed"))
    sat_assignment = assignment_to_json(solved.get("assignment")) if solved.get("status") == "SAT" else None
    if sat_assignment is not None:
        if not satisfies(formula, assignment_from_json(sat_assignment) or {}):
            raise AssertionError("parent SAT assignment failed independent current-formula check")
    core = {
        "schema": EXPERIENCE_SCHEMA,
        "created_generation": int(generation),
        "formula_hash": cnf_hash(formula),
        "residual_formula_hash": cnf_hash(residual),
        "budget": int(budget),
        "provider": provider_identity(),
        "prior_status": solved.get("status"),
        "prior_authorized": bool(solved.get("authorized")),
        "prior_lane": solved.get("primary_lane"),
        "found_k": primary.get("found_k"),
        "separator": list(separator) if separator else None,
        "components": [list(c) for c in components],
        "prior_separator_certificate": detector.get("certificate_sha256"),
        "prior_minimality_proof": detector.get("minimality_proof"),
        "route_reusable_after_revalidation": route_reusable,
        "sat_assignment": sat_assignment,
        "unsat_memory_is_verdict_shortcut": False,
    }
    return committed(core, "experience_commitment")


def verify_experience_commitment(record: Mapping[str, Any]) -> bool:
    if record.get("schema") != EXPERIENCE_SCHEMA or not isinstance(record.get("experience_commitment"), str):
        return False
    core = dict(record)
    got = core.pop("experience_commitment", None)
    return got == digest(core)


def validate_experience_for_formula(record: Mapping[str, Any], formula: Any, budget: int) -> tuple[bool, str | None, dict[str, Any]]:
    formula = canonical_cnf(formula)
    residual = guarded.residual_after_parent_preprocessing(formula)
    if not verify_experience_commitment(record):
        return False, "EXPERIENCE_COMMITMENT_INVALID", {}
    if record.get("provider") != provider_identity():
        return False, "EXPERIENCE_PROVIDER_MISMATCH", {}
    if record.get("formula_hash") != cnf_hash(formula) or record.get("budget") != int(budget):
        return False, "EXPERIENCE_FORMULA_OR_BUDGET_MISMATCH", {}
    if record.get("residual_formula_hash") != cnf_hash(residual):
        return False, "EXPERIENCE_RESIDUAL_HASH_MISMATCH", {}

    details: dict[str, Any] = {"sat_witness_reverified": False, "separator_revalidated": False}
    assignment = assignment_from_json(record.get("sat_assignment"))
    if assignment is not None:
        if not satisfies(formula, assignment):
            return False, "CACHED_SAT_WITNESS_REJECTED", details
        details["sat_witness_reverified"] = True

    separator = record.get("separator")
    if separator:
        claim = gate.verify_separator_claim(residual, separator, record.get("components") or [])
        details["separator_revalidation"] = claim
        details["separator_revalidated"] = bool(claim.get("passed"))
        if not claim.get("passed"):
            return False, "CACHED_SEPARATOR_REVALIDATION_FAILED", details
    return True, None, details


def reusable_detection(record: Mapping[str, Any], residual: Any) -> dict[str, Any]:
    residual = canonical_cnf(residual)
    claim = gate.verify_separator_claim(residual, record["separator"], record["components"])
    if not claim.get("passed"):
        raise ValueError("separator route cannot be reused without current-formula revalidation")
    k = int(record["found_k"])
    metrics = {
        "graph_clause_visits": int(claim.get("graph_clause_visits", 0)),
        "graph_pair_edge_attempts": int(claim.get("graph_pair_edge_attempts", 0)),
        "base_component_flood_rounds": 0,
        "candidate_counts_by_k": {},
        "connectivity_checks_by_k": {},
        "component_flood_rounds_by_k": {},
        "separator_clause_partition_checks_by_k": {str(k): int(claim.get("separator_clause_partition_checks", 0))},
        "verified_separator_counts_by_k": {str(k): 1},
        "cumulative_candidate_count": 0,
        "experience_revalidation_only": True,
    }
    core = {
        "formula_hash": cnf_hash(residual),
        "matched": True,
        "reason": None,
        "kmax": gate.K_MAX,
        "found_k": k,
        "audited_k": [],
        "separator": list(record["separator"]),
        "components": [list(c) for c in record["components"]],
        "metrics": metrics,
        "minimality_proof": {
            "fresh_minimality_recomputed": False,
            "prior_generation_minimality_recorded": bool(record.get("prior_minimality_proof")),
            "current_separator_revalidated": True,
            "all_candidates_below_and_at_found_k_exhausted": False,
        },
        "experience_source_commitment": record["experience_commitment"],
        "reuse_class": "FORMULA_BOUND_REVALIDATED_SEPARATOR_ROUTE",
    }
    return {**core, "certificate_sha256": gate.digest(core)}


@contextmanager
def install_experience_detector(record: Mapping[str, Any] | None):
    original = gate.detect_min_separator
    telemetry = {"root_route_reused": False, "reuse_calls": 0}
    target_hash = record.get("residual_formula_hash") if record else None

    def detector(formula, kmax=gate.K_MAX):
        current_hash = cnf_hash(canonical_cnf(formula))
        if record and current_hash == target_hash:
            candidate = reusable_detection(record, formula)
            telemetry["root_route_reused"] = True
            telemetry["reuse_calls"] += 1
            return candidate
        return original(formula, kmax)

    gate.detect_min_separator = detector
    try:
        yield telemetry
    finally:
        gate.detect_min_separator = original


def root_discovery_candidates(solved: Mapping[str, Any]) -> int:
    primary = solved.get("primary") or {}
    detector = primary.get("detector") or {}
    metrics = detector.get("metrics") or {}
    return int(metrics.get("cumulative_candidate_count", 0) or 0)


def technical_position(formula: Any, budget: int) -> dict[str, Any]:
    formula = canonical_cnf(formula)
    core = {
        "formula_hash": cnf_hash(formula),
        "budget": int(budget),
        "provider": provider_identity(),
    }
    return {**core, "position_commitment": digest(core)}


def append_origin_prime(
    store: dict[str, Any],
    position: Mapping[str, Any],
    experience: Mapping[str, Any] | None,
    result: Mapping[str, Any],
    reuse_mode: str,
    memory_event: Mapping[str, Any],
    semantic: Mapping[str, Any],
) -> dict[str, Any]:
    previous = store["current_state"]
    generation = int(previous["generation"]) + 1
    path_event = {
        "semantic_path": semantic["canonical_path"],
        "generation": generation,
        "formula_hash": position["formula_hash"],
        "reuse_mode": reuse_mode,
        "experience_commitment": experience.get("experience_commitment") if experience else None,
        "technical_status": result.get("status"),
        "authorized": bool(result.get("authorized")),
        "memory_status": memory_event.get("memory_status"),
    }
    path_history_digest = digest(path_event)
    return_core = {
        "formula_hash": position["formula_hash"],
        "position_commitment": position["position_commitment"],
        "technical_status": result.get("status"),
        "authorized": bool(result.get("authorized")),
        "sat_assignment_verified": bool(result.get("assignment_verified")) if result.get("status") == "SAT" else None,
        "unsat_current_closure_executed": bool(result.get("spiral", {}).get("unsat_current_closure_executed")),
        "memory_created_verdict": False,
    }
    return_commitment = digest(return_core)
    state_core = {
        "schema": STATE_SCHEMA,
        "state_type": "ORIGIN_PRIME",
        "generation": generation,
        "previous_state_commitment": previous["state_commitment"],
        "position_commitment": position["position_commitment"],
        "experience_commitment": experience.get("experience_commitment") if experience else None,
        "path_history_digest": path_history_digest,
        "return_commitment": return_commitment,
        "provider": provider_identity(),
    }
    state = committed(state_core, "state_commitment")
    store["history"].append(state)
    store["current_state"] = state
    return state


def solve_spiral(formula: Any, budget: int, state_path: Path) -> dict[str, Any]:
    contract, semantic = load_frozen_inputs()
    formula = canonical_cnf(formula)
    position = technical_position(formula, budget)
    store, memory_event = load_store(state_path)
    prior_state = copy.deepcopy(store["current_state"])
    key = experience_key(position["formula_hash"], budget)
    prior_experience = store["experiences"].get(key)
    experience_valid = False
    experience_reason = None
    validation_details: dict[str, Any] = {}
    if prior_experience is not None:
        experience_valid, experience_reason, validation_details = validate_experience_for_formula(prior_experience, formula, budget)

    # SAT memory is a proof hint only after the witness is rechecked on the exact current CNF.
    if experience_valid and prior_experience.get("prior_status") == "SAT" and validation_details.get("sat_witness_reverified"):
        assignment = assignment_from_json(prior_experience.get("sat_assignment")) or {}
        solved: dict[str, Any] = {
            "status": "SAT",
            "authorized": True,
            "assignment": assignment,
            "assignment_verified": satisfies(formula, assignment),
            "primary_lane": SPIRAL_SAT_LANE,
            "primary": {
                "lane": SPIRAL_SAT_LANE,
                "formula_hash": position["formula_hash"],
                "found_k": prior_experience.get("found_k"),
                "experience_commitment": prior_experience["experience_commitment"],
            },
            "cost_vector": {
                "residual_states": 0,
                "transition_checks": 0,
                "spiral_witness_clause_check": len(formula),
                "heterogeneous_units_not_summed_as_runtime": True,
            },
            "spiral": {
                "experience_reused": True,
                "reuse_mode": "REVERIFIED_SAT_WITNESS",
                "base_solver_invoked": False,
                "root_min_k_discovery_candidates": 0,
                "sat_witness_reverified": True,
                "unsat_current_closure_executed": False,
                "memory_created_verdict": False,
            },
        }
        current_experience = prior_experience
    else:
        route_record = None
        if experience_valid and prior_experience.get("route_reusable_after_revalidation") and prior_experience.get("separator"):
            route_record = prior_experience
        with install_experience_detector(route_record) as reuse_telemetry:
            solved = gate.solve_with_gate(formula, int(budget))
        solved = dict(solved)
        solved["spiral"] = {
            "experience_reused": bool(reuse_telemetry["root_route_reused"]),
            "reuse_mode": "REVALIDATED_SEPARATOR_ROUTE" if reuse_telemetry["root_route_reused"] else "FULL_DISCOVERY",
            "base_solver_invoked": True,
            "root_min_k_discovery_candidates": root_discovery_candidates(solved),
            "separator_revalidated": bool(validation_details.get("separator_revalidated")) if route_record else False,
            "fresh_minimality_recomputed": not bool(reuse_telemetry["root_route_reused"]),
            "unsat_current_closure_executed": bool(solved.get("status") == "UNSAT" and solved.get("authorized")),
            "memory_created_verdict": False,
        }
        new_experience = make_experience(formula, budget, solved, int(prior_state["generation"]) + 1)
        if new_experience is not None:
            store["experiences"][key] = new_experience
            current_experience = new_experience
        else:
            current_experience = prior_experience if experience_valid else None

    reuse_mode = solved["spiral"]["reuse_mode"]
    new_state = append_origin_prime(store, position, current_experience, solved, reuse_mode, memory_event, semantic)
    save_store(state_path, store)
    output = {
        "artifact_id": RUNTIME_ID,
        "status": solved.get("status"),
        "authorized": bool(solved.get("authorized")),
        "technical": solved,
        "position": position,
        "state_transition": {
            "from": prior_state,
            "to": new_state,
            "position_same_as_previous": prior_state.get("position_commitment") == new_state.get("position_commitment") if prior_state.get("position_commitment") else None,
            "state_changed": prior_state.get("state_commitment") != new_state.get("state_commitment"),
            "generation_advanced": new_state["generation"] == int(prior_state["generation"]) + 1,
        },
        "memory": {
            **memory_event,
            "experience_found": prior_experience is not None,
            "experience_valid": experience_valid,
            "experience_rejection_reason": experience_reason,
            "validation_details": validation_details,
            "active_experience_commitment": current_experience.get("experience_commitment") if current_experience else None,
        },
        "semantic_binding": {
            "semantic_spec_blob": EXPECTED_SEMANTIC_BLOB,
            "path": semantic["canonical_path"],
            "RETURN_NE_RESET": True,
            "STATE_RETURN_NE_ORIGINAL_STATE": True,
        },
        "scientific_boundary": contract["scientific_boundary"],
    }
    output["integrity_sha256"] = digest(output)
    return output


def fixture_formula(spec: Mapping[str, Any]):
    p = spec["params"]
    return gate.separator_clique_dumbbell(int(p["k"]), int(p["left_size"]), int(p["right_size"]), bool(p["conflict"]))


def _rewrite_experience(path: Path, key: str, mutate) -> None:
    store = json.loads(path.read_text(encoding="utf-8"))
    record = copy.deepcopy(store["experiences"][key])
    mutate(record)
    store["experiences"][key] = record
    path.write_text(json.dumps(store, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _recommit_experience(record: dict[str, Any]) -> None:
    record.pop("experience_commitment", None)
    record["experience_commitment"] = digest(record)


def self_test() -> dict[str, Any]:
    contract, _semantic = load_frozen_inputs()
    assert verify_provider_checkout(), "exact Fundamentum provider pin is not checked out"
    specs = {row["id"]: row for row in contract["frozen_fixtures"]}
    sat_formula = fixture_formula(specs["SPIRAL_K4_SAT"])
    unsat_formula = fixture_formula(specs["SPIRAL_K4_UNSAT"])
    mismatch_formula = fixture_formula(specs["SPIRAL_MISMATCH_K3_SAT"])
    assert cnf_hash(sat_formula) == specs["SPIRAL_K4_SAT"]["formula_sha256"]
    assert cnf_hash(unsat_formula) == specs["SPIRAL_K4_UNSAT"]["formula_sha256"]
    assert cnf_hash(mismatch_formula) == specs["SPIRAL_MISMATCH_K3_SAT"]["formula_sha256"]

    with tempfile.TemporaryDirectory(prefix="osiris-spiral-") as td:
        td = Path(td)
        sat_state = td / "sat_state.json"
        unsat_state = td / "unsat_state.json"
        mismatch_state = td / "mismatch_state.json"

        sat1 = solve_spiral(sat_formula, 50000, sat_state)
        sat2 = solve_spiral(sat_formula, 50000, sat_state)
        unsat1 = solve_spiral(unsat_formula, 50000, unsat_state)
        unsat2 = solve_spiral(unsat_formula, 50000, unsat_state)
        mismatch1 = solve_spiral(sat_formula, 50000, mismatch_state)
        mismatch2 = solve_spiral(mismatch_formula, 50000, mismatch_state)

        sat_first_cost = int(sat1["technical"]["spiral"]["root_min_k_discovery_candidates"])
        sat_second_cost = int(sat2["technical"]["spiral"]["root_min_k_discovery_candidates"])
        unsat_first_cost = int(unsat1["technical"]["spiral"]["root_min_k_discovery_candidates"])
        unsat_second_cost = int(unsat2["technical"]["spiral"]["root_min_k_discovery_candidates"])

        # Negative controls operate on committed memory, including attacks that
        # recompute the outer experience hash after semantic mutation.
        negatives: dict[str, bool] = {}

        p = td / "neg_commit.json"
        copy.copy(sat_state)
        p.write_bytes(sat_state.read_bytes())
        store = json.loads(p.read_text(encoding="utf-8"))
        key = next(iter(store["experiences"]))
        store["experiences"][key]["experience_commitment"] = mutated_hex(store["experiences"][key]["experience_commitment"])
        p.write_text(json.dumps(store, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        neg = solve_spiral(sat_formula, 50000, p)
        negatives["experience_commitment_bitflip_rejected"] = not neg["technical"]["spiral"]["experience_reused"] and neg["memory"]["memory_status"] == "REJECTED_FALLBACK_FRESH_ORIGIN"

        base_store = json.loads(sat_state.read_text(encoding="utf-8"))
        base_key = next(iter(base_store["experiences"]))
        base_exp = base_store["experiences"][base_key]

        swapped = copy.deepcopy(base_exp)
        swapped["formula_hash"] = cnf_hash(mismatch_formula)
        _recommit_experience(swapped)
        ok, reason, _ = validate_experience_for_formula(swapped, sat_formula, 50000)
        negatives["formula_hash_swap_rejected"] = (not ok and reason == "EXPERIENCE_FORMULA_OR_BUDGET_MISMATCH")

        provider_bad = copy.deepcopy(base_exp)
        provider_bad["provider"] = dict(provider_bad["provider"])
        provider_bad["provider"]["fundamentum_pin"] = mutated_hex(FUNDAMENTUM_PIN)
        _recommit_experience(provider_bad)
        ok, reason, _ = validate_experience_for_formula(provider_bad, sat_formula, 50000)
        negatives["provider_pin_mismatch_rejected"] = (not ok and reason == "EXPERIENCE_PROVIDER_MISMATCH")

        separator_bad = copy.deepcopy(base_exp)
        separator_bad["separator"] = list(separator_bad["separator"])
        separator_bad["separator"][0] = 12
        _recommit_experience(separator_bad)
        ok, reason, _ = validate_experience_for_formula(separator_bad, sat_formula, 50000)
        negatives["separator_mutation_rejected"] = (not ok and reason == "CACHED_SEPARATOR_REVALIDATION_FAILED")

        components_bad = copy.deepcopy(base_exp)
        components_bad["components"] = [sum((list(c) for c in components_bad["components"]), [])]
        _recommit_experience(components_bad)
        ok, reason, _ = validate_experience_for_formula(components_bad, sat_formula, 50000)
        negatives["component_partition_mutation_rejected"] = (not ok and reason == "CACHED_SEPARATOR_REVALIDATION_FAILED")

        assignment_bad = copy.deepcopy(base_exp)
        amap = dict(assignment_bad["sat_assignment"])
        first_var = sorted(amap, key=int)[0]
        amap[first_var] = not bool(amap[first_var])
        assignment_bad["sat_assignment"] = amap
        _recommit_experience(assignment_bad)
        ok, reason, _ = validate_experience_for_formula(assignment_bad, sat_formula, 50000)
        negatives["cached_sat_assignment_mutation_rejected"] = (not ok and reason == "CACHED_SAT_WITNESS_REJECTED")

        negatives["cached_unsat_not_used_as_verdict_shortcut"] = bool(
            unsat2["technical"]["spiral"]["base_solver_invoked"]
            and unsat2["technical"]["spiral"]["unsat_current_closure_executed"]
            and unsat2["technical"]["status"] == "UNSAT"
        )

        forced = json.loads(sat_state.read_text(encoding="utf-8"))
        forced_state = dict(forced["current_state"])
        forced_state["state_type"] = "ORIGIN"
        forced_state.pop("state_commitment", None)
        forced_state["state_commitment"] = digest(forced_state)
        forced["current_state"] = forced_state
        forced_path = td / "forced_origin.json"
        forced_path.write_text(json.dumps(forced, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        neg = solve_spiral(sat_formula, 50000, forced_path)
        negatives["forced_origin_reuse_rejected"] = neg["memory"]["memory_status"] == "REJECTED_FALLBACK_FRESH_ORIGIN"

        gates = {
            "sat_first_exact": sat1["status"] == "SAT" and sat1["authorized"],
            "sat_second_exact": sat2["status"] == "SAT" and sat2["authorized"],
            "sat_second_reuses_reverified_witness": sat2["technical"]["spiral"]["reuse_mode"] == "REVERIFIED_SAT_WITNESS" and sat2["technical"]["assignment_verified"],
            "sat_state_advances_without_position_reset": sat2["state_transition"]["generation_advanced"] and sat2["state_transition"]["state_changed"] and sat2["state_transition"]["position_same_as_previous"] is True,
            "sat_reuse_reduces_root_discovery": sat_first_cost > sat_second_cost,
            "unsat_first_exact": unsat1["status"] == "UNSAT" and unsat1["authorized"],
            "unsat_second_exact": unsat2["status"] == "UNSAT" and unsat2["authorized"],
            "unsat_second_reuses_only_revalidated_route": unsat2["technical"]["spiral"]["reuse_mode"] == "REVALIDATED_SEPARATOR_ROUTE" and unsat2["technical"]["spiral"]["separator_revalidated"],
            "unsat_current_closure_still_executes": unsat2["technical"]["spiral"]["unsat_current_closure_executed"] and unsat2["technical"]["spiral"]["base_solver_invoked"],
            "unsat_reuse_reduces_root_discovery": unsat_first_cost > unsat_second_cost,
            "formula_mismatch_does_not_reuse_foreign_experience": mismatch1["status"] == "SAT" and mismatch2["technical"]["spiral"]["experience_reused"] is False and mismatch2["memory"]["experience_found"] is False,
            "all_negative_controls_reject": all(negatives.values()),
            "P_VS_NP_OPEN": contract["scientific_boundary"]["P_VS_NP"] == "OPEN",
        }
        passed = all(gates.values())
        result = {
            "artifact_id": "OSIRIS-V3-SPIRAL-COMPUTE-SELF-TEST-2026-08-19-v1",
            "status": "PASS_KEEP_OSIRIS_V3_ORIGIN_PRIME_SPIRAL_COMPUTE" if passed else "STOP_OSIRIS_V3_SPIRAL_COMPUTE_GATE_FAILURE",
            "provider": provider_identity(),
            "semantic_spec_blob": EXPECTED_SEMANTIC_BLOB,
            "cost_observation": {
                "SAT_first_root_discovery_candidates": sat_first_cost,
                "SAT_second_root_discovery_candidates": sat_second_cost,
                "UNSAT_first_root_discovery_candidates": unsat_first_cost,
                "UNSAT_second_root_discovery_candidates": unsat_second_cost,
                "heterogeneous_units_not_summed": True,
                "claim": "repeated exact-input discovery reuse only; not a general complexity result",
            },
            "state_observation": {
                "SAT_generation_1": sat1["state_transition"]["to"]["state_commitment"],
                "SAT_generation_2": sat2["state_transition"]["to"]["state_commitment"],
                "same_position_generation_2": sat2["state_transition"]["position_same_as_previous"],
                "state_changed_generation_2": sat2["state_transition"]["state_changed"],
            },
            "negative_controls": negatives,
            "gates": gates,
            "scientific_boundary": contract["scientific_boundary"],
        }
        result["integrity_sha256"] = digest(result)
        return result


def load_cnf_document(path: Path) -> tuple[Any, int | None]:
    doc = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(doc, list):
        return doc, None
    if isinstance(doc, Mapping) and isinstance(doc.get("formula"), list):
        return doc["formula"], int(doc["budget"]) if doc.get("budget") is not None else None
    raise ValueError("CNF input must be a JSON list of clauses or {'formula': [...], 'budget': N}")


def main() -> None:
    ap = argparse.ArgumentParser(description="OSIRIS ORIGIN_PRIME spiral compute runtime")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--cnf", type=Path)
    ap.add_argument("--budget", type=int, default=None)
    ap.add_argument("--state", type=Path, default=Path.home() / ".janus" / "osiris_spiral_state.json")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    if args.self_test or args.cnf is None:
        result = self_test()
    else:
        formula, embedded_budget = load_cnf_document(args.cnf)
        budget = args.budget if args.budget is not None else embedded_budget
        if budget is None:
            raise SystemExit("--budget is required when the CNF document does not provide one")
        result = solve_spiral(formula, budget, args.state)

    text = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
    if args.self_test and not str(result.get("status", "")).startswith("PASS_KEEP"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
