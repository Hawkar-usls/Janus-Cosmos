#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import json
import tempfile
from pathlib import Path
from typing import Any, Mapping

import osiris_spiral_runtime as spiral

CONTRACT_PATH = Path(__file__).with_name("OSIRIS_V3_1_RESONANT_LINEAGE_ORION_PYRAMID_FROZEN_CONTRACT.json")
CERT_SCHEMA = "janus.cosmos.osiris_lineage_variable_renaming_certificate.v1"
TRANSFER_CLASS = "VARIABLE_RENAMING_BIJECTION"
RUNTIME_ID = "OSIRIS-V3.1-RESONANT-LINEAGE-2026-08-19-v1"


def _vars(formula: Any) -> list[int]:
    return sorted({abs(int(lit)) for clause in spiral.canonical_cnf(formula) for lit in clause})


def _normalize_map(variable_map: Mapping[Any, Any]) -> dict[int, int]:
    out: dict[int, int] = {}
    for raw_src, raw_dst in variable_map.items():
        src = int(raw_src)
        dst = int(raw_dst)
        if src <= 0 or dst <= 0:
            raise ValueError("VARIABLE_MAP_REQUIRES_POSITIVE_VARIABLE_IDS")
        if src in out:
            raise ValueError("VARIABLE_MAP_DUPLICATE_SOURCE")
        out[src] = dst
    return out


def remap_formula(formula: Any, variable_map: Mapping[Any, Any]) -> Any:
    source = spiral.canonical_cnf(formula)
    mapping = _normalize_map(variable_map)
    source_vars = _vars(source)
    if sorted(mapping) != source_vars:
        raise ValueError("VARIABLE_MAP_NOT_TOTAL_ON_SOURCE")
    if len(set(mapping.values())) != len(mapping):
        raise ValueError("VARIABLE_MAP_NOT_INJECTIVE")
    remapped = []
    for clause in source:
        row = []
        for lit in clause:
            dst = mapping[abs(int(lit))]
            row.append(dst if int(lit) > 0 else -dst)
        remapped.append(row)
    return spiral.canonical_cnf(remapped)


def build_variable_renaming_certificate(source_formula: Any, target_formula: Any, variable_map: Mapping[Any, Any]) -> dict[str, Any]:
    source = spiral.canonical_cnf(source_formula)
    target = spiral.canonical_cnf(target_formula)
    mapping = _normalize_map(variable_map)
    remapped = remap_formula(source, mapping)
    if remapped != target:
        raise ValueError("REMAPPED_SOURCE_DOES_NOT_EQUAL_TARGET")
    if sorted(mapping.values()) != _vars(target):
        raise ValueError("VARIABLE_MAP_NOT_SURJECTIVE_ON_TARGET")
    core = {
        "schema": CERT_SCHEMA,
        "transformation_class": TRANSFER_CLASS,
        "source_formula": [list(clause) for clause in source],
        "source_formula_hash": spiral.cnf_hash(source),
        "target_formula_hash": spiral.cnf_hash(target),
        "variable_map": {str(k): int(v) for k, v in sorted(mapping.items())},
        "literal_sign_preserved": True,
        "authority_delta": 0,
    }
    return {**core, "certificate_sha256": spiral.digest(core)}


def verify_variable_renaming_certificate(certificate: Mapping[str, Any], target_formula: Any) -> tuple[Any, dict[int, int]]:
    if not isinstance(certificate, Mapping) or certificate.get("schema") != CERT_SCHEMA:
        raise ValueError("LINEAGE_CERTIFICATE_SCHEMA_INVALID")
    if certificate.get("transformation_class") != TRANSFER_CLASS:
        raise ValueError("LINEAGE_TRANSFORMATION_CLASS_INVALID")
    if certificate.get("literal_sign_preserved") is not True or certificate.get("authority_delta") != 0:
        raise ValueError("LINEAGE_CERTIFICATE_AUTHORITY_INVALID")
    claimed = certificate.get("certificate_sha256")
    if not isinstance(claimed, str) or len(claimed) != 64:
        raise ValueError("LINEAGE_CERTIFICATE_HASH_INVALID")
    core = dict(certificate)
    core.pop("certificate_sha256", None)
    if spiral.digest(core) != claimed:
        raise ValueError("LINEAGE_CERTIFICATE_TAMPERED")

    source = spiral.canonical_cnf(certificate.get("source_formula") or [])
    target = spiral.canonical_cnf(target_formula)
    mapping = _normalize_map(certificate.get("variable_map") or {})
    if certificate.get("source_formula_hash") != spiral.cnf_hash(source):
        raise ValueError("LINEAGE_SOURCE_FORMULA_HASH_INVALID")
    if certificate.get("target_formula_hash") != spiral.cnf_hash(target):
        raise ValueError("LINEAGE_TARGET_FORMULA_HASH_INVALID")
    if remap_formula(source, mapping) != target:
        raise ValueError("LINEAGE_TARGET_FORMULA_MISMATCH")
    if sorted(mapping.values()) != _vars(target):
        raise ValueError("LINEAGE_VARIABLE_MAP_NOT_BIJECTIVE")
    return source, mapping


def _remap_assignment(assignment: Mapping[str, Any] | None, mapping: Mapping[int, int]) -> dict[str, bool] | None:
    if assignment is None:
        return None
    return {str(mapping[int(k)]): bool(v) for k, v in assignment.items()}


def _remap_var_list(values: Any, mapping: Mapping[int, int]) -> list[int] | None:
    if values is None:
        return None
    return sorted(mapping[int(v)] for v in values)


def _remap_components(values: Any, mapping: Mapping[int, int]) -> list[list[int]]:
    return [sorted(mapping[int(v)] for v in component) for component in (values or [])]


def transform_experience(
    source_record: Mapping[str, Any],
    target_formula: Any,
    budget: int,
    certificate: Mapping[str, Any],
) -> dict[str, Any]:
    target = spiral.canonical_cnf(target_formula)
    source_formula, mapping = verify_variable_renaming_certificate(certificate, target)
    if not spiral.verify_experience_commitment(source_record):
        raise ValueError("LINEAGE_SOURCE_EXPERIENCE_COMMITMENT_INVALID")
    if source_record.get("provider") != spiral.provider_identity():
        raise ValueError("LINEAGE_SOURCE_PROVIDER_MISMATCH")
    if source_record.get("formula_hash") != spiral.cnf_hash(source_formula):
        raise ValueError("LINEAGE_SOURCE_EXPERIENCE_FORMULA_MISMATCH")
    if int(source_record.get("budget", -1)) != int(budget):
        raise ValueError("LINEAGE_BUDGET_MISMATCH")

    target_residual = spiral.guarded.residual_after_parent_preprocessing(target)
    sat_assignment = _remap_assignment(source_record.get("sat_assignment"), mapping)
    if sat_assignment is not None:
        assignment = spiral.assignment_from_json(sat_assignment) or {}
        if not spiral.satisfies(target, assignment):
            raise ValueError("LINEAGE_TRANSFERRED_SAT_WITNESS_REJECTED")

    separator = _remap_var_list(source_record.get("separator"), mapping)
    components = _remap_components(source_record.get("components"), mapping)
    route_reusable = False
    separator_claim = None
    if separator:
        separator_claim = spiral.gate.verify_separator_claim(target_residual, separator, components)
        if not separator_claim.get("passed"):
            raise ValueError("LINEAGE_TRANSFERRED_SEPARATOR_REJECTED")
        route_reusable = True

    core = {
        "schema": spiral.EXPERIENCE_SCHEMA,
        "created_generation": int(source_record.get("created_generation", 0)),
        "formula_hash": spiral.cnf_hash(target),
        "residual_formula_hash": spiral.cnf_hash(target_residual),
        "budget": int(budget),
        "provider": spiral.provider_identity(),
        "prior_status": source_record.get("prior_status"),
        "prior_authorized": bool(source_record.get("prior_authorized")),
        "prior_lane": source_record.get("prior_lane"),
        "found_k": source_record.get("found_k"),
        "separator": separator,
        "components": components,
        "prior_separator_certificate": None,
        "prior_minimality_proof": copy.deepcopy(source_record.get("prior_minimality_proof")),
        "route_reusable_after_revalidation": route_reusable,
        "sat_assignment": sat_assignment,
        "unsat_memory_is_verdict_shortcut": False,
        "lineage_transfer": {
            "class": TRANSFER_CLASS,
            "source_experience_commitment": source_record["experience_commitment"],
            "source_formula_hash": source_record["formula_hash"],
            "transformation_certificate_sha256": certificate["certificate_sha256"],
            "target_sat_witness_revalidated": sat_assignment is not None,
            "target_separator_revalidated": bool(separator_claim and separator_claim.get("passed")),
            "memory_may_propose_not_verdict": True,
            "authority_delta": 0,
        },
    }
    return spiral.committed(core, "experience_commitment")


def propose_lineage_experience(
    state_path: Path,
    target_formula: Any,
    budget: int,
    certificate: Mapping[str, Any],
) -> dict[str, Any]:
    target = spiral.canonical_cnf(target_formula)
    source_formula, _mapping = verify_variable_renaming_certificate(certificate, target)
    store, memory_event = spiral.load_store(state_path)
    if memory_event.get("memory_status") != "VALID":
        raise ValueError("LINEAGE_SOURCE_STORE_NOT_VALID")
    source_key = spiral.experience_key(spiral.cnf_hash(source_formula), int(budget))
    source_record = store.get("experiences", {}).get(source_key)
    if source_record is None:
        raise ValueError("LINEAGE_SOURCE_EXPERIENCE_NOT_FOUND")
    transferred = transform_experience(source_record, target, int(budget), certificate)
    target_key = spiral.experience_key(spiral.cnf_hash(target), int(budget))

    candidate = copy.deepcopy(store)
    candidate["experiences"][target_key] = transferred
    ok, reason = spiral.verify_store(candidate)
    if not ok:
        raise ValueError(f"LINEAGE_PROPOSED_STORE_INVALID:{reason}")
    spiral.save_store(state_path, candidate)
    return {
        "proposed": True,
        "transfer_class": TRANSFER_CLASS,
        "source_formula_hash": spiral.cnf_hash(source_formula),
        "target_formula_hash": spiral.cnf_hash(target),
        "source_experience_commitment": source_record["experience_commitment"],
        "target_experience_commitment": transferred["experience_commitment"],
        "transformation_certificate_sha256": certificate["certificate_sha256"],
        "authority_delta": 0,
    }


def solve_resonant(
    formula: Any,
    budget: int,
    state_path: Path,
    certificate: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    if contract.get("status") != "FROZEN_BEFORE_IMPLEMENTATION":
        raise ValueError("RESONANT_LINEAGE_CONTRACT_NOT_FROZEN")
    lineage = {
        "proposal_attempted": certificate is not None,
        "proposal_accepted": False,
        "fallback_to_fresh_target_solve": False,
        "reason": None,
        "transfer": None,
    }
    if certificate is not None:
        try:
            transfer = propose_lineage_experience(state_path, formula, int(budget), certificate)
            lineage["proposal_accepted"] = True
            lineage["transfer"] = transfer
        except (TypeError, ValueError) as exc:
            lineage["fallback_to_fresh_target_solve"] = True
            lineage["reason"] = str(exc)

    result = spiral.solve_spiral(formula, int(budget), state_path)
    result = dict(result)
    result["artifact_id"] = RUNTIME_ID
    result["resonant_lineage"] = {
        **lineage,
        "memory_may_propose_not_verdict": True,
        "unsat_memory_is_verdict_shortcut": False,
        "authority_delta": 0,
    }
    result.pop("integrity_sha256", None)
    result["integrity_sha256"] = spiral.digest(result)
    return result


def self_test() -> dict[str, Any]:
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    if contract["admitted_transformation_classes"][TRANSFER_CLASS]["status"] != "ADMITTED_WITH_EXACT_VERIFIER":
        raise AssertionError("TRANSFER_CLASS_NOT_ADMITTED")
    _base_contract, _semantic = spiral.load_frozen_inputs()
    specs = {row["id"]: row for row in _base_contract["frozen_fixtures"]}
    source_formula = spiral.fixture_formula(specs["SPIRAL_K4_SAT"])
    mapping = {v: v + 100 for v in _vars(source_formula)}
    target_formula = remap_formula(source_formula, mapping)
    if spiral.cnf_hash(source_formula) == spiral.cnf_hash(target_formula):
        raise AssertionError("SELFTEST_RENAMED_FORMULA_HASH_DID_NOT_CHANGE")
    cert = build_variable_renaming_certificate(source_formula, target_formula, mapping)

    with tempfile.TemporaryDirectory(prefix="osiris-resonant-lineage-") as td:
        state_path = Path(td) / "state.json"
        first = spiral.solve_spiral(source_formula, 50000, state_path)
        if first.get("status") != "SAT" or first.get("authorized") is not True:
            raise AssertionError("SELFTEST_SOURCE_SOLVE_FAILED")
        second = solve_resonant(target_formula, 50000, state_path, cert)
        if second["resonant_lineage"]["proposal_accepted"] is not True:
            raise AssertionError("SELFTEST_LINEAGE_PROPOSAL_NOT_ACCEPTED")
        if second["technical"]["spiral"]["reuse_mode"] != "REVERIFIED_SAT_WITNESS":
            raise AssertionError("SELFTEST_TARGET_WITNESS_NOT_REVERIFIED")
        if second["state_transition"]["generation_advanced"] is not True:
            raise AssertionError("SELFTEST_GENERATION_DID_NOT_ADVANCE")

        negatives: dict[str, bool] = {}
        bad_map = dict(mapping)
        keys = sorted(bad_map)
        bad_map[keys[1]] = bad_map[keys[0]]
        try:
            build_variable_renaming_certificate(source_formula, target_formula, bad_map)
            negatives["non_bijective_map"] = False
        except ValueError:
            negatives["non_bijective_map"] = True

        tampered = copy.deepcopy(cert)
        tampered["target_formula_hash"] = "0" * 64
        try:
            verify_variable_renaming_certificate(tampered, target_formula)
            negatives["certificate_tamper"] = False
        except ValueError:
            negatives["certificate_tamper"] = True

        wrong_target = spiral.canonical_cnf(list(target_formula) + [[999]])
        try:
            verify_variable_renaming_certificate(cert, wrong_target)
            negatives["target_mismatch"] = False
        except ValueError:
            negatives["target_mismatch"] = True

        try:
            transform_experience(
                next(iter(json.loads(state_path.read_text(encoding="utf-8"))["experiences"].values())),
                target_formula,
                49999,
                cert,
            )
            negatives["budget_mismatch"] = False
        except ValueError:
            negatives["budget_mismatch"] = True

        if not all(negatives.values()):
            raise AssertionError(f"RESONANT_LINEAGE_NEGATIVE_FAILED:{negatives}")

        return {
            "status": "PASS_KEEP_OSIRIS_V3_1_RESONANT_LINEAGE",
            "source_formula_hash": spiral.cnf_hash(source_formula),
            "target_formula_hash": spiral.cnf_hash(target_formula),
            "hash_changed": True,
            "transfer_class": TRANSFER_CLASS,
            "target_reuse_mode": second["technical"]["spiral"]["reuse_mode"],
            "lineage_proposal_accepted": True,
            "negative_controls": negatives,
            "authority_delta": 0,
            "P_VS_NP": "OPEN",
        }


def main() -> int:
    parser = argparse.ArgumentParser(description="OSIRIS v3.1 proof-carrying resonant lineage")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not args.self_test:
        parser.error("this entrypoint currently exposes --self-test only; use solve_resonant() as a library API")
    result = self_test()
    text = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
