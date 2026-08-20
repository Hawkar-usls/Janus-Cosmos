#!/usr/bin/env python3
"""Deterministic evidence-preserving research spiral for Janus-Cosmos.

This adapter does not replace the technical OSIRIS solver spiral. It applies the
same non-resetting state law to research hypotheses: every layer is preserved,
content-addressed, parent-bound, and may refine search without gaining truth
authority by memory alone.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

DEFAULT_INPUT = Path("data/giza/JANUS-GIZA-ORION-JET-SAND-REVERSE-SPIRAL-EVIDENCE-v1.json")
DEFAULT_OUTPUT = Path("data/giza/JANUS-GIZA-ORION-JET-SAND-REVERSE-SPIRAL-RUN-001-RECEIPT.json")
CONTRACT = Path("workspace/JANUS_COSMOS_EPISTEMIC_SPIRAL_FROZEN_CONTRACT_v1.json")
EXPECTED_DIRECTIONS = ["BACK", "FORWARD", "LEFT", "RIGHT", "FORWARD_AGAIN", "BACK_AGAIN"]


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def assemble(source: dict[str, Any], contract: dict[str, Any]) -> dict[str, Any]:
    raw_layers = source.get("layers", [])
    if [x.get("direction") for x in raw_layers] != EXPECTED_DIRECTIONS:
        raise AssertionError("Reverse/spiral direction sequence drift")
    if not contract["runtime_law"]["preserve_all_layers"]:
        raise AssertionError("Contract must preserve all layers")
    if contract["runtime_law"]["automatic_winner_selection"]:
        raise AssertionError("Automatic winner selection is forbidden")

    layers = []
    parent_sha = None
    for generation, raw in enumerate(raw_layers):
        core = {
            "generation": generation,
            "layer_id": raw["layer_id"],
            "direction": raw["direction"],
            "question": raw["question"],
            "necessary_predictions": raw.get("necessary_predictions", []),
            "evidence_for": raw.get("evidence_for", []),
            "evidence_against": raw.get("evidence_against", []),
            "alternatives": raw.get("alternatives", []),
            "status": raw["status"],
            "claim_ceiling": raw["claim_ceiling"],
            "parent_layer_sha256": parent_sha,
        }
        layer_sha = sha256(core)
        layers.append({**core, "layer_sha256": layer_sha})
        parent_sha = layer_sha

    source_digest = sha256(source)
    contract_digest = sha256(contract)
    origin_prime_core = {
        "state": source["origin_prime"]["state"],
        "preserved_result": source["origin_prime"]["preserved_result"],
        "next_test": source["origin_prime"]["next_test"],
        "previous_layer_sha256": parent_sha,
        "source_sha256": source_digest,
        "contract_sha256": contract_digest,
    }
    origin_prime_sha = sha256(origin_prime_core)

    receipt_core = {
        "schema": "janus.cosmos.epistemic_spiral.receipt.v1",
        "experiment_id": "JANUS-GIZA-ORION-JET-SAND-REVERSE-SPIRAL-RUN-001",
        "status": "PASS_SPIRAL_ASSEMBLY__PHYSICAL_HYPOTHESIS_FAIL__SYMBOLIC_LAYER_PRESERVED",
        "input_artifact_id": source["artifact_id"],
        "source_sha256": source_digest,
        "contract_sha256": contract_digest,
        "spiral": {
            "order": EXPECTED_DIRECTIONS,
            "preserve_all_layers": True,
            "automatic_winner_selection": False,
            "layer_count": len(layers),
            "layers": layers,
        },
        "origin_prime": {**origin_prime_core, "origin_prime_sha256": origin_prime_sha},
        "hieroglyphic_overlay": source["hieroglyphic_overlay"],
        "gates": {
            "direction_sequence_exact": True,
            "all_prior_layers_preserved": len(layers) == len(raw_layers),
            "sha_parent_chain_valid": all(layers[i]["parent_layer_sha256"] == (None if i == 0 else layers[i-1]["layer_sha256"]) for i in range(len(layers))),
            "negative_geology_layer_preserved": any(x["layer_id"] == "L1_GEOLOGY" and "FAIL" in x["status"] for x in layers),
            "negative_glass_mechanism_layer_preserved": any(x["layer_id"] == "L2_GLASS_TO_SAND" and "FAIL" in x["status"] for x in layers),
            "symbolic_layer_does_not_override_physics": layers[-1]["claim_ceiling"] == "TEXTUAL_SYMBOLIC_PROVENANCE_ONLY",
            "automatic_canonization_disabled": True,
        },
        "scientific_boundary": {
            "orion_jet_vitrified_giza": "NOT_ESTABLISHED__CURRENT_EVIDENCE_GATE_FAIL",
            "giza_sand_from_hydrated_orion_glass": "NOT_ESTABLISHED__CURRENT_EVIDENCE_GATE_FAIL",
            "orion_osiris_star_shining_textual_motif": "SUPPORTED_AS_SYMBOLIC_TEXTUAL_CONTEXT",
            "hieroglyphic_overlay_is_historical_translation": False,
            "memory_is_truth": False,
        },
    }
    if not all(receipt_core["gates"].values()):
        raise AssertionError(f"Gate failure: {receipt_core['gates']}")
    return {**receipt_core, "receipt_sha256": sha256(receipt_core)}


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def self_test(source_path: Path = DEFAULT_INPUT) -> dict[str, Any]:
    source = load(source_path)
    contract = load(CONTRACT)
    a = assemble(source, contract)
    b = assemble(source, contract)
    if a != b:
        raise AssertionError("Spiral assembly must be deterministic")
    if a["receipt_sha256"] != sha256({k: v for k, v in a.items() if k != "receipt_sha256"}):
        raise AssertionError("Receipt SHA binding failed")
    return {
        "status": "PASS",
        "receipt_sha256": a["receipt_sha256"],
        "origin_prime_sha256": a["origin_prime"]["origin_prime_sha256"],
        "layer_count": a["spiral"]["layer_count"],
        "physical_hypothesis": a["scientific_boundary"]["orion_jet_vitrified_giza"],
        "symbolic_layer": a["scientific_boundary"]["orion_osiris_star_shining_textual_motif"],
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        print(json.dumps(self_test(args.input), ensure_ascii=False, indent=2))
        return 0
    result = assemble(load(args.input), load(CONTRACT))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "receipt_sha256": result["receipt_sha256"], "gates": result["gates"]}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
