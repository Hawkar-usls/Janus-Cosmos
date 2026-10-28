#!/usr/bin/env python3
"""JANUS Cosmos: independent, no-credit scientific evidence bridge.

Cross-repo evidence transport only; no biological/astronomical inference,
no provider submissions, no changes to OSIRIS or ocean/earth/space gates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from urllib.request import Request, urlopen

BASE = "https://raw.githubusercontent.com/Hawkar-usls/"
SOURCES = {
    "canonical": ("janus-meta-registry", "data/rowan/CANONICAL_STATE.json"),
    "aura": ("aura-oracle-tg", "research/rowan/LATEST_HYPOTHESES.json"),
    "demiurge": ("Janus-Demiurge", "research/rowan/LATEST_COUNCIL.json"),
    "coordination": ("janus-meta-registry", "data/rowan/coordination/LATEST_COORDINATION_STATE.json"),
    "topa": ("janus-meta-registry", "data/rowan/topa/TOPA-SHADOW-86ca72205cc89a03.json"),
    "shadow": ("janus-meta-registry", "data/rowan/learned/SHADOW-BBC1-abac9057873ac3c4.json"),
}
HEX64 = re.compile(r"^[0-9a-f]{64}$")
MAX_BYTES = 120_000

class GateHold(ValueError):
    """Fail-closed: evidence is missing, inconsistent or over-authorized."""

def digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()

def semantic_digest(payload: bytes) -> str:
    """SHA-256 of the source protocol's sorted compact JSON, not file bytes."""
    obj = parse_payload(payload, "semantic")
    return digest(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode())


def require(condition: bool, reason: str) -> None:
    if not condition:
        raise GateHold(reason)

def parse_payload(raw: bytes, name: str) -> dict:
    require(0 < len(raw) <= MAX_BYTES, "SOURCE_SIZE_" + name)
    try:
        obj = json.loads(raw)
    except (ValueError, UnicodeDecodeError) as exc:
        raise GateHold("INVALID_JSON_" + name) from exc
    require(isinstance(obj, dict), "INVALID_OBJECT_" + name)
    return obj

def fetch_public(repo: str, path: str) -> bytes:
    # These names are program-owned allowlisted constants or a validated receipt
    # path; no arbitrary user URL, credentials or provider API access.
    url = BASE + repo + "/main/" + path
    request = Request(url, headers={"User-Agent": "JANUS-Cosmos-Scientific-Witness/1",
                                     "Accept": "application/json"})
    with urlopen(request, timeout=18) as response:
        body = response.read(MAX_BYTES + 1)
    require(0 < len(body) <= MAX_BYTES, "UPSTREAM_TOO_LARGE_OR_EMPTY")
    return body

def acquire() -> dict[str, bytes]:
    raw = {key: fetch_public(repo, path) for key, (repo, path) in SOURCES.items()}
    c = parse_payload(raw["coordination"], "coordination")
    receipt_path = c.get("receipt_path")
    require(isinstance(receipt_path, str) and re.fullmatch(
        r"data/rowan/coordination/receipts/ROWAN-COUNCIL-[0-9]+-[0-9a-f]{16}\.json",
        receipt_path), "UNSAFE_OR_MISSING_COORDINATION_RECEIPT_PATH")
    raw["receipt"] = fetch_public("janus-meta-registry", receipt_path)
    return raw

def reconcile(raw: dict[str, bytes]) -> dict:
    require(set(raw) == set(SOURCES) | {"receipt"}, "INCOMPLETE_SOURCE_SET")
    obj = {key: parse_payload(value, key) for key, value in raw.items()}
    c, a, d, s, t, m, r = (obj[k] for k in
        ("canonical", "aura", "demiurge", "coordination", "topa", "shadow", "receipt"))
    run = c.get("latest_durable_external_run")
    count = c.get("last_completed_mutants")
    require(c.get("schema_id") == "janus.rowan.observer_state.v1" and
            type(run) is int and run > 0 and type(count) is int and 0 <= count <= 8,
            "CANONICAL_IDENTITY_OR_PROGRESS")
    source_sha = semantic_digest(raw["canonical"])
    require(HEX64.fullmatch(str(c.get("underlying_state_is_derived_from", ""))) is not None,
            "CANONICAL_SOURCE_ARTIFACT_UNATTESTED")
    require(a.get("schema_id") == "janus.aura.rowan_hypothesis_packet.v1", "AURA_SCHEMA")
    a_source = a.get("source", {})
    require(a_source.get("state_sha256") == source_sha and
            a_source.get("upstream_run") == run and
            a_source.get("completed_mutants") == count and
            a_source.get("upstream_artifact_sha256") == c["underlying_state_is_derived_from"],
            "AURA_CANONICAL_BINDING_DRIFT")
    require(a.get("source_has_unseen_outcome") is False and
            a.get("epistemic_status") == "HYPOTHESIS_AND_FALSIFIER_PROPOSAL_ONLY" and
            a.get("bounded_authority", {}).get("may_submit_rowan_jobs") is False,
            "AURA_AUTHORITY_LEAK")
    require(d.get("schema_id") == "janus.demiurge.rowan_council_proposal.v1" and
            d.get("scientific_authority") == "NONE__PROPOSAL_NOT_EXPERIMENT",
            "DEMIURGE_AUTHORITY_DRIFT")
    ds = d.get("sources", {})
    require(ds.get("meta_registry", {}).get("canonical_run") == run and
            ds.get("meta_registry", {}).get("sha256") == source_sha and
            ds.get("meta_registry", {}).get("source_artifact_sha256") == c["underlying_state_is_derived_from"] and
            ds.get("aura_oracle", {}).get("packet_sha256") == semantic_digest(raw["aura"]) and
            ds.get("aura_oracle", {}).get("packet_digest") == a.get("packet_digest_sha256"),
            "DEMIURGE_AURA_BINDING_DRIFT")
    require(d.get("gates", {}).get("rowan_job_submission_new") is False and
            d.get("gates", {}).get("production_model_weight_promotion") is False,
            "DEMIURGE_SPEND_OR_MODEL_PROMOTION")
    require(s.get("schema_id") == "janus.meta_registry.rowan_council_state.v1" and
            s.get("meta_canonical_run") == run, "COORDINATION_RUN_DRIFT")
    require(r.get("schema_id") == "janus.meta_registry.rowan_council_provenance.v1" and
            r.get("source_binding_sha256") == s.get("source_binding_sha256") and
            HEX64.fullmatch(str(s.get("source_binding_sha256", ""))) is not None,
            "COORDINATION_RECEIPT_DRIFT")
    rs = r.get("sources", {})
    require(rs.get("meta_canonical_run") == run and
            rs.get("meta_canonical_sha256") == source_sha and
            rs.get("rowan_source_artifact_sha256") == c["underlying_state_is_derived_from"] and
            rs.get("aura_packet_sha256") == digest(raw["aura"]) and
            rs.get("aura_packet_digest") == a.get("packet_digest_sha256") and
            rs.get("demiurge_council_sha256") == semantic_digest(raw["demiurge"]),
            "IMMUTABLE_RECEIPT_SOURCE_DRIFT")
    authority = r.get("authority", {})
    require(all(authority.get(key) is False for key in
                ("new_credit_spend", "new_rowan_jobs", "new_study_authorized",
                 "hypotheses_are_empirical_proof")), "RECEIPT_AUTHORITY_LEAK")
    require(t.get("schema_id") == "janus.topa.independent_shadow_replay.v1" and
            t.get("audit") == "PASS_INDEPENDENT_WITHIN_ASSAY_REPLAY" and
            t.get("promotion_gate") == "SHADOW_ONLY" and
            t.get("model_file_sha256") == digest(raw["shadow"]) and
            t.get("independent_remote_source_sha256") == m.get("source", {}).get("sha256"),
            "TOPA_SHADOW_EVIDENCE_DRIFT")
    require(m.get("model", {}).get("features") == 44 and
            m.get("authority", {}).get("deployment") == "SHADOW_ONLY" and
            m.get("authority", {}).get("rowan_credits_used") == 0,
            "SHADOW_MODEL_AUTHORITY_DRIFT")
    require(t.get("external_lab_replication") == "OPEN" and
            t.get("out_of_assay_generalization") == "NOT_TESTED",
            "UNREVIEWED_SCIENTIFIC_PROMOTION")
    packet = {
        "schema_id": "janus.cosmos.rowan_cross_domain_evidence_witness.v1",
        "status": "PASS_READ_ONLY_SOURCE_BOUND_REPLAY",
        "canonical_rowan_run": run,
        "completed_mutants": count,
        "coordination_receipt": s["receipt_path"],
        "source_sha256": {key: digest(value) for key, value in sorted(raw.items())},
        "coordination_source_binding_sha256": s["source_binding_sha256"],
        "source_artifact_sha256": c["underlying_state_is_derived_from"],
        "topa": {
            "scope": "WITHIN_BBC1_ASSAY_ONLY",
            "audit": t["audit"],
            "test_variant_count": t["heldout_variant_count"],
            "test_mae": t["test_mae_independent"],
            "baseline_mae": t["training_mean_baseline_mae_independent"],
            "model_family": t["model_family"],
            "promotion_gate": "SHADOW_ONLY",
        },
        "aura_next_question": a.get("next_discriminating_question_id"),
        "demiurge_next_control": d.get("next_independent_control"),
        "cosmos_domain_transfer": "PROTOCOL_AND_PROVENANCE_ONLY",
        "astronomical_physical_evidence": False,
        "protein_fitness_or_ddg_established": False,
        "protein_prediction_calibrated_on_other_proteins": False,
        "cosmos_earth_space_ocean_bhb_osiris_gates_modified": False,
        "rowan_provider_jobs_created": 0,
        "rowan_credits_spent": 0,
        "authority": "OBSERVE_AND_ARCHIVE_ONLY",
    }
    packet["observation_key"] = digest(json.dumps(
        packet["source_sha256"], sort_keys=True, separators=(",", ":")).encode())
    return packet

def render(obj: dict) -> bytes:
    return (json.dumps(obj, indent=2, ensure_ascii=False, sort_keys=True,
                       allow_nan=False) + "\n").encode("utf-8")

def persist(root: Path, packet: dict) -> tuple[str, bool]:
    base = root / "data/integrations/rowan"
    receipt = base / "receipts" / ("COSMOS-ROWAN-" +
        str(packet["canonical_rowan_run"]) + "-" + packet["observation_key"][:16] + ".json")
    payload = render(packet)
    if receipt.exists():
        require(receipt.read_bytes() == payload, "IMMUTABLE_RECEIPT_CHANGED")
        changed = False
    else:
        receipt.parent.mkdir(parents=True, exist_ok=True)
        receipt.write_bytes(payload)
        changed = True
    pointer = {"schema_id": "janus.cosmos.rowan_latest_witness.v1",
               "observation_key": packet["observation_key"],
               "receipt_path": str(receipt.relative_to(root)).replace("\\", "/"),
               "canonical_rowan_run": packet["canonical_rowan_run"],
               "status": packet["status"],
               "cosmos_domain_transfer": packet["cosmos_domain_transfer"]}
    latest = base / "LATEST.json"
    latest_bytes = render(pointer)
    if not latest.exists() or latest.read_bytes() != latest_bytes:
        latest.parent.mkdir(parents=True, exist_ok=True)
        latest.write_bytes(latest_bytes)
        changed = True
    return str(receipt), changed

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--persist-root", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    try:
        packet = reconcile(acquire())
        if args.out:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_bytes(render(packet))
        if args.persist_root:
            path, changed = persist(args.persist_root, packet)
            print("COSMOS_ROWAN_RECEIPT=" + path)
            print("COSMOS_ROWAN_PERSIST_CHANGED=" + str(changed).upper())
        print("COSMOS_ROWAN_WITNESS=" + packet["status"])
        print("COSMOS_ROWAN_CANONICAL_RUN=" + str(packet["canonical_rowan_run"]))
        return 0
    except (GateHold, OSError, TimeoutError) as exc:
        print("COSMOS_ROWAN_WITNESS=HOLD_" + type(exc).__name__ +
              (":" + str(exc) if isinstance(exc, GateHold) else ""))
        return 2

if __name__ == "__main__":
    raise SystemExit(main())
