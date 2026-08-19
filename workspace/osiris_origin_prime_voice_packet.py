#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import re
import sys
import tempfile
from pathlib import Path
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT / "workspace"
if str(WORKSPACE) not in sys.path:
    sys.path.insert(0, str(WORKSPACE))

import osiris_spiral_runtime as spiral  # noqa: E402

CONTRACT_PATH = WORKSPACE / "OSIRIS_V3_1_DEMIHEAD_VOICE_PYRAMID_MEDIATED_LINK_FROZEN_CONTRACT.json"
PACKET_SCHEMA = "janus.cosmos.origin_prime_voice_packet.v1"
COSMOS_REPOSITORY = "Hawkar-usls/Janus-Cosmos"
DEMIHEAD_REPOSITORY = "Hawkar-usls/Demi_Head"
VOICE_REPOSITORY = "Hawkar-usls/The-Voice-of-Janus"
ECHO_REPOSITORY = "Hawkar-usls/Echo-Pyramid"
VOICE_SHA = "e58d65aa46b7e3a64a5131708578a9a3346915c4"
ECHO_SHA = "15712f5b14b123d4e3cb64ddeaa693c5bf6af788"
PROFILE = {
    "profile_id": "PYRAMID_LANGUAGE_117_121_ANCHORED_SPACE_v0.3",
    "anchor_band_hz": [117.0, 121.0],
    "center_hz": 119.0,
    "q": 29.75,
    "gain_db": 11.5,
    "decay_s": 1.65,
    "role": "REPRESENTATION_AND_ACOUSTIC_COLORATION_ONLY",
    "frequencies_create_math_authority": False,
    "audio_output_is_evidence": False,
}
CONTROL = {
    "direct_cosmos_to_echo_route_permitted": False,
    "demihead_mediation_required": True,
    "network_io_required": False,
    "automatic_playback": False,
    "automatic_microphone_start": False,
    "authority_delta": 0,
    "mass_effect_budget_delta": 0,
}
GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
MAX_PACKET_BYTES = 65536


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def source_revision(explicit: str | None = None) -> str:
    candidate = explicit or os.environ.get("GITHUB_SHA") or os.environ.get("JANUS_COSMOS_REVISION")
    if candidate is None:
        try:
            import subprocess
            candidate = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
        except Exception as exc:  # pragma: no cover - environment dependent
            raise ValueError("COSMOS_SOURCE_REVISION_REQUIRED") from exc
    if not isinstance(candidate, str) or GIT_SHA.fullmatch(candidate) is None:
        raise ValueError("COSMOS_SOURCE_REVISION_INVALID")
    return candidate


def load_contract() -> dict[str, Any]:
    value = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    if value.get("status") != "FROZEN_BEFORE_IMPLEMENTATION":
        raise ValueError("VOICE_LINK_CONTRACT_NOT_FROZEN")
    if value.get("parent", {}).get("sha") != "b99949a5529024987396fb8a353e912a21bbfe82":
        raise ValueError("VOICE_LINK_PARENT_MISMATCH")
    return value


def load_verified_store(path: Path) -> dict[str, Any]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    ok, reason = spiral.verify_store(raw)
    if not ok:
        raise ValueError(f"OSIRIS_STATE_STORE_INVALID:{reason}")
    current = raw.get("current_state")
    if not isinstance(current, Mapping) or current.get("state_type") != "ORIGIN_PRIME":
        raise ValueError("ORIGIN_PRIME_STATE_REQUIRED")
    if int(current.get("generation", 0)) < 1:
        raise ValueError("ORIGIN_PRIME_GENERATION_INVALID")
    return raw


def bound_experience(store: Mapping[str, Any], state: Mapping[str, Any]) -> dict[str, Any] | None:
    commitment = state.get("experience_commitment")
    if commitment is None:
        return None
    matches = [
        copy.deepcopy(dict(record))
        for record in store.get("experiences", {}).values()
        if isinstance(record, Mapping) and record.get("experience_commitment") == commitment
    ]
    if len(matches) != 1:
        raise ValueError("STATE_EXPERIENCE_BINDING_NOT_UNIQUE")
    record = matches[0]
    if not spiral.verify_experience_commitment(record):
        raise ValueError("BOUND_EXPERIENCE_COMMITMENT_INVALID")
    return record


def build_packet(state_path: Path, *, revision: str | None = None) -> dict[str, Any]:
    contract = load_contract()
    store = load_verified_store(state_path)
    state = copy.deepcopy(dict(store["current_state"]))
    experience = bound_experience(store, state)
    core: dict[str, Any] = {
        "schema": PACKET_SCHEMA,
        "source": {
            "repository": COSMOS_REPOSITORY,
            "revision": source_revision(revision),
            "canonical_gate": "OSIRIS_V3_ORIGIN_PRIME_SPIRAL_COMPUTE",
            "state_store_schema": store["schema"],
        },
        "origin_prime": state,
        "bound_experience": experience,
        "mediation": {
            "required_mediator": DEMIHEAD_REPOSITORY,
            "voice_repository": VOICE_REPOSITORY,
            "voice_revision": contract["peers"]["voice"]["sha"],
            "physical_body_repository": ECHO_REPOSITORY,
            "physical_body_revision": contract["peers"]["echo_pyramid"]["sha"],
            "route": "COSMOS -> DEMIHEAD -> THE_VOICE_OF_JANUS -> ECHO_PYRAMID",
        },
        "voice_representation": copy.deepcopy(PROFILE),
        "control": copy.deepcopy(CONTROL),
        "scientific_boundary": {
            "P_VS_NP": "OPEN",
            "P_EQUALS_NP": "NOT_ESTABLISHED",
            "P_NOT_EQUALS_NP": "NOT_ESTABLISHED",
            "voice_profile_changes_solver_correctness": False,
            "acoustic_frequencies_are_proof": False,
        },
    }
    packet = {**core, "packet_sha256": digest(core)}
    validate_packet(packet)
    if len(canonical_bytes(packet)) > MAX_PACKET_BYTES:
        raise ValueError("ORIGIN_PRIME_VOICE_PACKET_TOO_LARGE")
    return packet


def validate_packet(packet: Mapping[str, Any]) -> None:
    if not isinstance(packet, Mapping) or packet.get("schema") != PACKET_SCHEMA:
        raise ValueError("VOICE_PACKET_SCHEMA_INVALID")
    source = packet.get("source")
    if not isinstance(source, Mapping) or source.get("repository") != COSMOS_REPOSITORY:
        raise ValueError("VOICE_PACKET_SOURCE_INVALID")
    if not isinstance(source.get("revision"), str) or GIT_SHA.fullmatch(source["revision"]) is None:
        raise ValueError("VOICE_PACKET_SOURCE_REVISION_INVALID")
    if source.get("canonical_gate") != "OSIRIS_V3_ORIGIN_PRIME_SPIRAL_COMPUTE":
        raise ValueError("VOICE_PACKET_GATE_INVALID")

    state = packet.get("origin_prime")
    if not isinstance(state, Mapping) or state.get("state_type") != "ORIGIN_PRIME":
        raise ValueError("VOICE_PACKET_STATE_INVALID")
    if not spiral.verify_state_commitment(state):
        raise ValueError("VOICE_PACKET_STATE_COMMITMENT_INVALID")

    experience = packet.get("bound_experience")
    if experience is None:
        if state.get("experience_commitment") is not None:
            raise ValueError("VOICE_PACKET_EXPERIENCE_MISSING")
    else:
        if not isinstance(experience, Mapping) or not spiral.verify_experience_commitment(experience):
            raise ValueError("VOICE_PACKET_EXPERIENCE_INVALID")
        if experience.get("experience_commitment") != state.get("experience_commitment"):
            raise ValueError("VOICE_PACKET_EXPERIENCE_BINDING_INVALID")

    mediation = packet.get("mediation")
    expected_mediation = {
        "required_mediator": DEMIHEAD_REPOSITORY,
        "voice_repository": VOICE_REPOSITORY,
        "voice_revision": VOICE_SHA,
        "physical_body_repository": ECHO_REPOSITORY,
        "physical_body_revision": ECHO_SHA,
        "route": "COSMOS -> DEMIHEAD -> THE_VOICE_OF_JANUS -> ECHO_PYRAMID",
    }
    if mediation != expected_mediation:
        raise ValueError("VOICE_PACKET_MEDIATION_INVALID")
    if packet.get("voice_representation") != PROFILE:
        raise ValueError("VOICE_PACKET_PROFILE_INVALID")
    if packet.get("control") != CONTROL:
        raise ValueError("VOICE_PACKET_CONTROL_INVALID")
    boundary = packet.get("scientific_boundary")
    if not isinstance(boundary, Mapping) or boundary.get("P_VS_NP") != "OPEN":
        raise ValueError("VOICE_PACKET_SCIENTIFIC_BOUNDARY_INVALID")
    if boundary.get("voice_profile_changes_solver_correctness") is not False or boundary.get("acoustic_frequencies_are_proof") is not False:
        raise ValueError("VOICE_PACKET_AUTHORITY_LEAK")

    claimed = packet.get("packet_sha256")
    if not isinstance(claimed, str) or len(claimed) != 64:
        raise ValueError("VOICE_PACKET_HASH_INVALID")
    body = dict(packet)
    body.pop("packet_sha256", None)
    if claimed != digest(body):
        raise ValueError("VOICE_PACKET_HASH_TAMPERED")


def verify_packet(packet: Mapping[str, Any]) -> bool:
    try:
        validate_packet(packet)
    except (TypeError, ValueError):
        return False
    return True


def self_test() -> dict[str, Any]:
    load_contract()
    with tempfile.TemporaryDirectory(prefix="osiris-voice-link-") as td:
        state_path = Path(td) / "state.json"
        solved = spiral.solve_spiral([[1]], 256, state_path)
        if solved.get("status") != "SAT" or solved.get("authorized") is not True:
            raise AssertionError("SELFTEST_OSIRIS_SOLVE_FAILED")
        packet = build_packet(state_path, revision="a" * 40)
        if not verify_packet(packet):
            raise AssertionError("SELFTEST_PACKET_VERIFY_FAILED")

        negatives: dict[str, bool] = {}
        mutations = {
            "state_commitment_bitflip": lambda p: p["origin_prime"].__setitem__("state_commitment", "0" * 64),
            "experience_commitment_bitflip": lambda p: p["bound_experience"].__setitem__("experience_commitment", "0" * 64),
            "voice_profile_tamper": lambda p: p["voice_representation"].__setitem__("center_hz", 120.0),
            "direct_cosmos_to_echo_route": lambda p: p["control"].__setitem__("direct_cosmos_to_echo_route_permitted", True),
            "authority_escalation": lambda p: p["control"].__setitem__("authority_delta", 1),
            "formula_hash_unbound": lambda p: p["bound_experience"].__setitem__("formula_hash", "0" * 64),
            "packet_hash_tamper": lambda p: p.__setitem__("packet_sha256", "0" * 64),
        }
        for name, mutate in mutations.items():
            candidate = copy.deepcopy(packet)
            mutate(candidate)
            negatives[name] = not verify_packet(candidate)

        raw_store = json.loads(state_path.read_text(encoding="utf-8"))
        raw_store["history"][0]["previous_state_commitment"] = "0" * 64
        bad_store = Path(td) / "bad_state.json"
        bad_store.write_text(json.dumps(raw_store, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        try:
            build_packet(bad_store, revision="a" * 40)
            negatives["broken_state_lineage"] = False
        except ValueError:
            negatives["broken_state_lineage"] = True

        if not all(negatives.values()):
            raise AssertionError(f"VOICE_LINK_NEGATIVE_FAILED:{negatives}")
        return {
            "status": "PASS_KEEP_OSIRIS_V3_1_DEMIHEAD_VOICE_PYRAMID_MEDIATED_LINK",
            "packet_schema": PACKET_SCHEMA,
            "generation": packet["origin_prime"]["generation"],
            "profile": packet["voice_representation"],
            "route": packet["mediation"]["route"],
            "negative_controls": negatives,
            "P_VS_NP": "OPEN",
            "authority_delta": 0,
        }


def main() -> int:
    parser = argparse.ArgumentParser(description="Export a verified OSIRIS ORIGIN_PRIME state for DemiHead-mediated Pyramid Language rendering")
    parser.add_argument("--state", type=Path)
    parser.add_argument("--source-revision")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        value = self_test()
    else:
        if args.state is None:
            parser.error("--state is required unless --self-test is used")
        value = build_packet(args.state, revision=args.source_revision)
    text = json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
