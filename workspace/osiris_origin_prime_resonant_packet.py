#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import json
import tempfile
from pathlib import Path
from typing import Any, Mapping

import osiris_origin_prime_voice_packet as voice_packet
import osiris_resonant_lineage as resonant
import osiris_spiral_runtime as spiral

EXTENSION_SCHEMA = "janus.cosmos.origin_prime_resonant_representation_extension.v1"
ASTRAL_REPRESENTATION = {
    "anchor_id": "ORION_BELT_SAH_OSIRIS_CONTEXT_v1",
    "star_triplet": ["Mintaka", "Alnilam", "Alnitak"],
    "egyptological_context": "SAH_ORION_OSIRIS_RELIGIOUS_TEXTUAL_CONTEXT",
    "giza_orion_correlation": "HYPOTHESIS_NOT_ASSERTED_AS_ARCHITECTURAL_FACT",
    "janus_rebus_alias": "S𓂸ḥ",
    "janus_rebus_alias_is_historical_transliteration": False,
    "seasonal_visibility": "CONTEXT_ONLY",
    "role": "ASTRAL_CONTEXT_AND_NAVIGATION_REPRESENTATION_ONLY",
    "authority_delta": 0,
    "astral_geometry_is_proof": False,
    "astral_context_changes_solver_correctness": False,
}


def _lineage_summary(experience: Mapping[str, Any] | None) -> dict[str, Any]:
    transfer = experience.get("lineage_transfer") if isinstance(experience, Mapping) else None
    if not isinstance(transfer, Mapping):
        return {
            "transfer_present": False,
            "transfer_class": None,
            "source_experience_commitment": None,
            "transformation_certificate_sha256": None,
            "memory_may_propose_not_verdict": True,
            "authority_delta": 0,
        }
    return {
        "transfer_present": True,
        "transfer_class": transfer.get("class"),
        "source_experience_commitment": transfer.get("source_experience_commitment"),
        "transformation_certificate_sha256": transfer.get("transformation_certificate_sha256"),
        "target_sat_witness_revalidated": bool(transfer.get("target_sat_witness_revalidated")),
        "target_separator_revalidated": bool(transfer.get("target_separator_revalidated")),
        "memory_may_propose_not_verdict": transfer.get("memory_may_propose_not_verdict") is True,
        "authority_delta": 0,
    }


def _representation_binding_core(packet: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "origin_prime_state_commitment": packet["origin_prime"]["state_commitment"],
        "experience_commitment": packet["origin_prime"].get("experience_commitment"),
        "voice_representation": packet["voice_representation"],
        "astral_representation": packet["astral_representation"],
        "lineage_representation": packet["lineage_representation"],
        "authority_delta": 0,
    }


def build_packet(state_path: Path, *, revision: str | None = None) -> dict[str, Any]:
    base = voice_packet.build_packet(state_path, revision=revision)
    core = dict(base)
    core.pop("packet_sha256", None)
    core["representation_extension_schema"] = EXTENSION_SCHEMA
    core["astral_representation"] = copy.deepcopy(ASTRAL_REPRESENTATION)
    core["lineage_representation"] = _lineage_summary(core.get("bound_experience"))
    core["representation_binding_sha256"] = voice_packet.digest(_representation_binding_core(core))
    packet = {**core, "packet_sha256": voice_packet.digest(core)}
    validate_packet(packet)
    return packet


def validate_packet(packet: Mapping[str, Any]) -> None:
    voice_packet.validate_packet(packet)
    if packet.get("representation_extension_schema") != EXTENSION_SCHEMA:
        raise ValueError("RESONANT_PACKET_EXTENSION_SCHEMA_INVALID")
    if packet.get("astral_representation") != ASTRAL_REPRESENTATION:
        raise ValueError("RESONANT_PACKET_ASTRAL_REPRESENTATION_INVALID")
    lineage = packet.get("lineage_representation")
    if not isinstance(lineage, Mapping):
        raise ValueError("RESONANT_PACKET_LINEAGE_REPRESENTATION_INVALID")
    if lineage.get("memory_may_propose_not_verdict") is not True or lineage.get("authority_delta") != 0:
        raise ValueError("RESONANT_PACKET_LINEAGE_AUTHORITY_LEAK")
    if lineage.get("transfer_present"):
        if lineage.get("transfer_class") != resonant.TRANSFER_CLASS:
            raise ValueError("RESONANT_PACKET_TRANSFER_CLASS_INVALID")
        for field in ("source_experience_commitment", "transformation_certificate_sha256"):
            value = lineage.get(field)
            if not isinstance(value, str) or len(value) != 64:
                raise ValueError(f"RESONANT_PACKET_{field.upper()}_INVALID")

    claimed_binding = packet.get("representation_binding_sha256")
    if not isinstance(claimed_binding, str) or len(claimed_binding) != 64:
        raise ValueError("RESONANT_PACKET_REPRESENTATION_BINDING_INVALID")
    if claimed_binding != voice_packet.digest(_representation_binding_core(packet)):
        raise ValueError("RESONANT_PACKET_REPRESENTATION_BINDING_TAMPERED")
    astral = packet["astral_representation"]
    if astral.get("authority_delta") != 0 or astral.get("astral_geometry_is_proof") is not False:
        raise ValueError("RESONANT_PACKET_ASTRAL_AUTHORITY_LEAK")


def verify_packet(packet: Mapping[str, Any]) -> bool:
    try:
        validate_packet(packet)
    except (TypeError, ValueError):
        return False
    return True


def _rehash(candidate: dict[str, Any]) -> None:
    candidate.pop("packet_sha256", None)
    candidate["packet_sha256"] = voice_packet.digest(candidate)


def self_test() -> dict[str, Any]:
    base_contract, _semantic = spiral.load_frozen_inputs()
    specs = {row["id"]: row for row in base_contract["frozen_fixtures"]}
    source_formula = spiral.fixture_formula(specs["SPIRAL_K4_SAT"])
    mapping = {v: v + 200 for v in resonant._vars(source_formula)}
    target_formula = resonant.remap_formula(source_formula, mapping)
    cert = resonant.build_variable_renaming_certificate(source_formula, target_formula, mapping)

    with tempfile.TemporaryDirectory(prefix="osiris-resonant-packet-") as td:
        state_path = Path(td) / "state.json"
        first = spiral.solve_spiral(source_formula, 50000, state_path)
        if first.get("status") != "SAT" or first.get("authorized") is not True:
            raise AssertionError("RESONANT_PACKET_SOURCE_SOLVE_FAILED")
        second = resonant.solve_resonant(target_formula, 50000, state_path, cert)
        if second["resonant_lineage"]["proposal_accepted"] is not True:
            raise AssertionError("RESONANT_PACKET_LINEAGE_TRANSFER_FAILED")
        packet = build_packet(state_path, revision="a" * 40)
        if not verify_packet(packet):
            raise AssertionError("RESONANT_PACKET_VERIFY_FAILED")
        if packet["lineage_representation"]["transfer_present"] is not True:
            raise AssertionError("RESONANT_PACKET_TRANSFER_PROVENANCE_MISSING")

        negatives: dict[str, bool] = {}

        candidate = copy.deepcopy(packet)
        candidate["astral_representation"]["authority_delta"] = 1
        candidate["representation_binding_sha256"] = voice_packet.digest(_representation_binding_core(candidate))
        _rehash(candidate)
        negatives["astral_authority_escalation"] = not verify_packet(candidate)

        candidate = copy.deepcopy(packet)
        candidate["astral_representation"]["star_triplet"][0] = "Tampered"
        candidate["representation_binding_sha256"] = voice_packet.digest(_representation_binding_core(candidate))
        _rehash(candidate)
        negatives["astral_profile_tamper"] = not verify_packet(candidate)

        candidate = copy.deepcopy(packet)
        candidate["representation_binding_sha256"] = "0" * 64
        _rehash(candidate)
        negatives["representation_binding_tamper"] = not verify_packet(candidate)

        candidate = copy.deepcopy(packet)
        candidate["lineage_representation"]["memory_may_propose_not_verdict"] = False
        candidate["representation_binding_sha256"] = voice_packet.digest(_representation_binding_core(candidate))
        _rehash(candidate)
        negatives["lineage_verdict_escalation"] = not verify_packet(candidate)

        if not all(negatives.values()):
            raise AssertionError(f"RESONANT_PACKET_NEGATIVE_FAILED:{negatives}")

        return {
            "status": "PASS_KEEP_OSIRIS_V3_1_RESONANT_ORION_PYRAMID_PACKET",
            "packet_schema": packet["schema"],
            "extension_schema": EXTENSION_SCHEMA,
            "state_commitment": packet["origin_prime"]["state_commitment"],
            "representation_binding_sha256": packet["representation_binding_sha256"],
            "voice_profile": packet["voice_representation"]["profile_id"],
            "orion_anchor": packet["astral_representation"]["anchor_id"],
            "lineage_transfer_present": True,
            "negative_controls": negatives,
            "authority_delta": 0,
            "P_VS_NP": "OPEN",
        }


def main() -> int:
    parser = argparse.ArgumentParser(description="Bind verified ORIGIN_PRIME lineage to Voice/Pyramid and Orion contextual representations")
    parser.add_argument("state", nargs="?", type=Path)
    parser.add_argument("--source-revision")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        result: Any = self_test()
    else:
        if args.state is None:
            parser.error("state is required unless --self-test is used")
        result = build_packet(args.state, revision=args.source_revision)
    text = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
