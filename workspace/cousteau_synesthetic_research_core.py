#!/usr/bin/env python3
"""JANUS Cousteau Synesthetic Research Core v2.

Research-grade, fail-closed wrapper around the frozen v1 sensory-passport core.
It adds raw/provenance binding, explicit epistemic state, multiscale memory and
an immutable Cousteau -> DemiHead handshake. Retrieval quality is never truth
confidence; mnemonic similarity is never scientific convergence.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import cousteau_synesthetic_memory_core as base
import cousteau_synesthetic_semantic_overlay as semantic

CORE_ID = "JANUS_COUSTEAU_SYNESTHETIC_RESEARCH_CORE"
CORE_VERSION = "2.0.0"
RESEARCH_SCHEMA = "janus.cosmos.cousteau.synesthetic_research_passport.v2"
HANDSHAKE_SCHEMA = "janus.synesthesia.handshake.packet.v1"
BUNDLE_SCHEMA = "janus.cosmos.cousteau.synesthetic_multiscale_bundle.v2"
PROTOCOL_ID = "JANUS_SYNAESTHETIC_RESEARCH_HANDSHAKE"
PROTOCOL_VERSION = "1.0.0"
PROTOCOL_CONTRACT_SHA256 = "3aec527be027fc280fc9a8ace1255c9a3a7da73fc884d9b4856694a1f1530306"
CONTRACT_PATH = ROOT / ".janus" / "JANUS_SYNAESTHETIC_RESEARCH_HANDSHAKE_V1.json"
EPISTEMIC_STATES = {"OBSERVED", "UNKNOWN", "STALE", "CONTAMINATED", "BLOCKED"}
PROFILES = {"HANNAH_BODC", "SYNTHETIC_TEST", "GENERIC_RESEARCH"}
HEX64 = re.compile(r"^[0-9a-f]{64}$")

AUTHORITY = {
    "memory_equals_truth": False,
    "mnemonic_similarity_is_scientific_similarity": False,
    "association_is_evidence": False,
    "may_change_raw_data": False,
    "may_change_calibration": False,
    "may_change_scientific_verdict": False,
    "may_reorder_review_priority": True,
    "authority_delta": 0,
    "mass_effect_budget_delta": 0,
}


def canonical_json(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)


def digest(obj: Any) -> str:
    return hashlib.sha256(canonical_json(obj).encode("utf-8")).hexdigest()


def _hex64(value: Any) -> bool:
    return isinstance(value, str) and HEX64.fullmatch(value) is not None


def load_and_verify_contract() -> dict[str, Any]:
    obj = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    if digest(obj) != PROTOCOL_CONTRACT_SHA256:
        raise RuntimeError("SYNAESTHETIC_HANDSHAKE_CONTRACT_HASH_MISMATCH")
    if obj.get("protocol_id") != PROTOCOL_ID or obj.get("protocol_version") != PROTOCOL_VERSION:
        raise RuntimeError("SYNAESTHETIC_HANDSHAKE_CONTRACT_ID_MISMATCH")
    return obj


def _raw_binding_valid(raw_binding: Mapping[str, Any] | None) -> bool:
    if not isinstance(raw_binding, Mapping) or not _hex64(raw_binding.get("source_raw_sha256")):
        return False
    parser_sha = raw_binding.get("parser_sha256")
    if parser_sha is not None and not _hex64(parser_sha):
        return False
    window_id = raw_binding.get("window_id")
    return window_id is None or (isinstance(window_id, str) and bool(window_id.strip()))


def _safe_state_key(name: str, feature_names: set[str]) -> bool:
    if not isinstance(name, str) or not name:
        return False
    slug = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
    if any(tok in slug for tok in base.FORBIDDEN_INFLUENCE_TOKENS):
        return False
    if slug in feature_names or slug in base.MEASUREMENT_KEYS:
        return True
    if any(slug.startswith(prefix) for prefix in base.DYNAMIC_MEASUREMENT_PREFIXES):
        return True
    return (slug.endswith("_sin") or slug.endswith("_cos")) and slug[:-4] in base.MEASUREMENT_KEYS


def summarize_epistemic_state(passport: Mapping[str, Any], field_states: Mapping[str, str] | None = None) -> dict[str, Any]:
    fp = passport.get("measurement_fingerprint")
    if not isinstance(fp, Mapping):
        return {
            "overall_state": "BLOCKED",
            "counts": {"OBSERVED": 0, "UNKNOWN": 0, "STALE": 0, "CONTAMINATED": 0, "BLOCKED": 1},
            "coverage_fraction": 0.0,
            "missing_fraction": 1.0,
            "unknown_fraction": 1.0,
            "stale_fraction": 0.0,
            "contaminated_fraction": 0.0,
            "retrieval_quality_score": 0.0,
            "quality_band": "BLOCKED",
            "truth_confidence": None,
            "rule": "RETRIEVAL_QUALITY_NE_TRUTH_CONFIDENCE",
        }

    features = fp.get("features") or {}
    feature_names = set(features)
    states = dict(field_states or {})
    for key, state in states.items():
        if state not in EPISTEMIC_STATES - {"BLOCKED"}:
            raise ValueError(f"invalid epistemic state for {key}: {state}")
        if not _safe_state_key(key, feature_names):
            raise ValueError(f"unsafe or non-measurement epistemic key: {key}")

    counts = {"OBSERVED": len(features), "UNKNOWN": 0, "STALE": 0, "CONTAMINATED": 0, "BLOCKED": 0}
    for key, state in states.items():
        normalized = re.sub(r"[^a-z0-9]+", "_", key.lower()).strip("_")
        if normalized in feature_names and counts["OBSERVED"] > 0:
            counts["OBSERVED"] -= 1
        counts[state] += 1

    total = sum(counts[s] for s in ("OBSERVED", "UNKNOWN", "STALE", "CONTAMINATED"))
    if total <= 0:
        return summarize_epistemic_state({"measurement_fingerprint": None})
    unknown = counts["UNKNOWN"] / total
    stale = counts["STALE"] / total
    contaminated = counts["CONTAMINATED"] / total
    coverage = counts["OBSERVED"] / total
    missing_norm = features.get("missing_fraction")
    missing = max(0.0, min(1.0, (float(missing_norm) + 1.0) / 2.0)) if isinstance(missing_norm, (int, float)) else unknown
    quality = coverage * (1.0 - missing) * max(0.0, 1.0 - 0.5 * stale - 0.9 * contaminated)
    quality = round(max(0.0, min(1.0, quality)), 6)
    overall = "CONTAMINATED" if contaminated else "STALE" if stale else "UNKNOWN" if unknown else "OBSERVED"
    return {
        "overall_state": overall,
        "counts": counts,
        "coverage_fraction": round(coverage, 6),
        "missing_fraction": round(missing, 6),
        "unknown_fraction": round(unknown, 6),
        "stale_fraction": round(stale, 6),
        "contaminated_fraction": round(contaminated, 6),
        "retrieval_quality_score": quality,
        "quality_band": "HIGH" if quality >= 0.85 else "MEDIUM" if quality >= 0.6 else "LOW",
        "truth_confidence": None,
        "rule": "RETRIEVAL_QUALITY_NE_TRUTH_CONFIDENCE",
    }


def _blocked_research(*, event_id: str, reason: str, direction: str, scale: str, provenance: Mapping[str, Any] | None, profile: str) -> dict[str, Any]:
    p = base.build_blocked_passport(
        reason=reason,
        source_receipt={"status": reason, "event_id": event_id},
        provenance=provenance,
    )
    out = {
        "schema": RESEARCH_SCHEMA,
        "core": {"id": CORE_ID, "version": CORE_VERSION, "base_core_version": base.CORE_VERSION},
        "profile": profile,
        "event_id": event_id,
        "context": {"direction": direction, "scale": scale},
        "status": "BLOCKED",
        "source_binding": {"status": "UNBOUND", "source_raw_sha256": None},
        "passport": p,
        "epistemic": summarize_epistemic_state(p),
        "authority": dict(AUTHORITY),
        "scientific_measurement_use_allowed": False,
        "scientific_convergence_claim": False,
    }
    out["research_passport_sha256"] = digest(out)
    return out


def build_research_passport(
    payload: Mapping[str, Any] | None,
    *,
    event_id: str,
    direction: str = "UNKNOWN",
    scale: str = "custom",
    provenance: Mapping[str, Any] | None = None,
    raw_bytes: bytes | None = None,
    raw_binding: Mapping[str, Any] | None = None,
    field_states: Mapping[str, str] | None = None,
    expected_feature_count: int | None = None,
    profile: str = "HANNAH_BODC",
) -> dict[str, Any]:
    if not isinstance(event_id, str) or not event_id.strip():
        raise ValueError("event_id is required")
    if profile not in PROFILES:
        raise ValueError(f"profile must be one of {sorted(PROFILES)}")
    load_and_verify_contract()
    if payload is None:
        return _blocked_research(event_id=event_id, reason="BLOCKED_NO_MEASUREMENT_PAYLOAD", direction=direction, scale=scale, provenance=provenance, profile=profile)

    raw_sha = hashlib.sha256(raw_bytes).hexdigest() if raw_bytes is not None else None
    bound = _raw_binding_valid(raw_binding)
    if profile == "HANNAH_BODC" and raw_sha is None and not bound:
        return _blocked_research(event_id=event_id, reason="BLOCKED_HANNAH_RAW_PROVENANCE_NOT_BOUND", direction=direction, scale=scale, provenance=provenance, profile=profile)

    p = semantic.enrich_passport(base.build_passport(
        payload,
        direction=direction,
        scale=scale,
        provenance=provenance,
        raw_bytes=raw_bytes,
        expected_feature_count=expected_feature_count,
    ))
    if not isinstance(p.get("measurement_fingerprint"), Mapping):
        return _blocked_research(event_id=event_id, reason="BLOCKED_NO_ACCEPTED_MEASUREMENT_FEATURES", direction=direction, scale=scale, provenance=provenance, profile=profile)

    if raw_sha is not None:
        source_binding = {"status": "RAW_BYTES_VERIFIED_THIS_CALL", "source_raw_sha256": raw_sha, "parser_sha256": None, "window_id": None}
    elif bound:
        source_binding = {
            "status": "RAW_SHA_BOUND_DERIVED_WINDOW",
            "source_raw_sha256": raw_binding["source_raw_sha256"],
            "parser_sha256": raw_binding.get("parser_sha256"),
            "window_id": raw_binding.get("window_id"),
        }
    else:
        source_binding = {"status": "UNBOUND_SYNTHETIC_OR_GENERIC", "source_raw_sha256": None, "parser_sha256": None, "window_id": None}

    out = {
        "schema": RESEARCH_SCHEMA,
        "core": {"id": CORE_ID, "version": CORE_VERSION, "base_core_version": base.CORE_VERSION},
        "profile": profile,
        "event_id": event_id,
        "context": {"direction": direction, "scale": scale},
        "status": "READY",
        "source_binding": source_binding,
        "passport": p,
        "epistemic": summarize_epistemic_state(p, field_states=field_states),
        "authority": dict(AUTHORITY),
        "scientific_measurement_use_allowed": source_binding["status"] in {"RAW_BYTES_VERIFIED_THIS_CALL", "RAW_SHA_BOUND_DERIVED_WINDOW"},
        "scientific_convergence_claim": False,
    }
    out["research_passport_sha256"] = digest(out)
    return out


def export_handshake_packet(research_passport: Mapping[str, Any]) -> dict[str, Any]:
    if research_passport.get("schema") != RESEARCH_SCHEMA:
        raise ValueError("invalid research passport schema")
    p = research_passport.get("passport") or {}
    fp = p.get("measurement_fingerprint")
    blocked = research_passport.get("status") == "BLOCKED" or not isinstance(fp, Mapping)
    if blocked:
        fp_export = None
    else:
        embedding = fp.get("embedding") or []
        if len(embedding) != base.EMBED_DIMS:
            raise ValueError("unexpected measurement embedding dimensions")
        fp_export = {
            "sha256": fp.get("sha256"),
            "blake2b_256": p.get("collision_guard_blake2b_256"),
            "embedding": list(embedding),
            "feature_count": int(fp.get("feature_count", len(fp.get("features") or {}))),
            "feature_names_sha256": digest(sorted((fp.get("features") or {}).keys())),
            "units_sha256": digest(fp.get("units") or {}),
        }
        if not _hex64(fp_export["sha256"]) or not _hex64(fp_export["blake2b_256"]):
            raise ValueError("invalid measurement fingerprint digests")

    sensory_material = {
        "sensory_channels": p.get("sensory_channels"),
        "cousteau_semantic_overlay": p.get("cousteau_semantic_overlay"),
    }
    channels = p.get("sensory_channels") or {}
    packet = {
        "schema": HANDSHAKE_SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "protocol_version": PROTOCOL_VERSION,
        "contract_sha256": PROTOCOL_CONTRACT_SHA256,
        "event_id": research_passport.get("event_id"),
        "producer": {
            "repository": "Hawkar-usls/Janus-Cosmos",
            "role": "COUSTEAU_MEASUREMENT_SENSORY_CORE",
            "core_id": CORE_ID,
            "core_version": CORE_VERSION,
        },
        "source_identity": p.get("source_identity") or {},
        "source_binding": research_passport.get("source_binding") or {},
        "context": p.get("context") or research_passport.get("context") or {},
        "measurement_fingerprint": fp_export,
        "sensory_digest": digest(sensory_material),
        "sensory_summary": {
            "color_hex": (channels.get("color") or {}).get("hex"),
            "audio_mode": (channels.get("audio") or {}).get("mode"),
            "audio_frequency_hz": (channels.get("audio") or {}).get("frequency_hz"),
            "texture": (channels.get("texture") or {}).get("label"),
            "semantic_overlay_sha256": (p.get("cousteau_semantic_overlay") or {}).get("semantic_overlay_sha256"),
        },
        "epistemic_state": research_passport.get("epistemic"),
        "authority": dict(AUTHORITY),
        "scientific_measurement_use_allowed": bool(research_passport.get("scientific_measurement_use_allowed")),
        "scientific_convergence_claim": False,
    }
    packet["packet_sha256"] = digest(packet)
    if not verify_handshake_packet(packet):
        raise RuntimeError("COUSTEAU_HANDSHAKE_SELF_VERIFY_FAILED")
    return packet


def verify_handshake_packet(packet: Mapping[str, Any]) -> bool:
    if not isinstance(packet, Mapping) or packet.get("schema") != HANDSHAKE_SCHEMA:
        return False
    if packet.get("protocol_id") != PROTOCOL_ID or packet.get("protocol_version") != PROTOCOL_VERSION or packet.get("contract_sha256") != PROTOCOL_CONTRACT_SHA256:
        return False
    producer = packet.get("producer")
    if not isinstance(producer, Mapping) or producer.get("repository") != "Hawkar-usls/Janus-Cosmos" or producer.get("role") != "COUSTEAU_MEASUREMENT_SENSORY_CORE":
        return False
    if packet.get("authority") != AUTHORITY or packet.get("scientific_convergence_claim") is not False:
        return False
    claimed = packet.get("packet_sha256")
    if not _hex64(claimed):
        return False
    core = dict(packet)
    core.pop("packet_sha256", None)
    if digest(core) != claimed:
        return False
    epistemic = packet.get("epistemic_state")
    if not isinstance(epistemic, Mapping) or epistemic.get("overall_state") not in EPISTEMIC_STATES:
        return False
    fp = packet.get("measurement_fingerprint")
    if fp is None:
        return epistemic.get("overall_state") == "BLOCKED" and packet.get("scientific_measurement_use_allowed") is False
    if not isinstance(fp, Mapping) or not _hex64(fp.get("sha256")) or not _hex64(fp.get("blake2b_256")):
        return False
    emb = fp.get("embedding")
    return isinstance(emb, list) and len(emb) == base.EMBED_DIMS and all(isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(float(x)) for x in emb)


def compare_research_passports(a: Mapping[str, Any], b: Mapping[str, Any]) -> dict[str, Any]:
    cmp = base.compare_passports(a.get("passport") or {}, b.get("passport") or {})
    sim = cmp.get("common_measurement_similarity")
    qa = float((a.get("epistemic") or {}).get("retrieval_quality_score", 0.0))
    qb = float((b.get("epistemic") or {}).get("retrieval_quality_score", 0.0))
    review = None if sim is None else round(float(sim) * math.sqrt(max(0.0, qa) * max(0.0, qb)), 6)
    return {
        "schema": "janus.cosmos.cousteau.synesthetic_research_comparison.v2",
        "event_a": a.get("event_id"),
        "event_b": b.get("event_id"),
        "measurement_similarity": sim,
        "mnemonic_cosine": cmp.get("mnemonic_cosine"),
        "retrieval_quality_a": qa,
        "retrieval_quality_b": qb,
        "quality_adjusted_review_score": review,
        "common_feature_count": cmp.get("common_feature_count"),
        "scientific_convergence_claim": False,
        "authority": "REVIEW_PRIORITY_ONLY",
    }


def rank_cross_front_research_pairs(head: Sequence[Mapping[str, Any]], tail: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = [{"head_event_id": h.get("event_id"), "tail_event_id": t.get("event_id"), **compare_research_passports(h, t)} for h in head for t in tail]
    rows.sort(key=lambda r: -1.0 if r["quality_adjusted_review_score"] is None else r["quality_adjusted_review_score"], reverse=True)
    return rows


def build_multiscale_bundle(passports: Sequence[Mapping[str, Any]], *, bundle_id: str) -> dict[str, Any]:
    if not passports or not isinstance(bundle_id, str) or not bundle_id.strip():
        raise ValueError("non-empty passports and bundle_id required")
    packets = [export_handshake_packet(p) for p in passports]
    entries = [{
        "event_id": p.get("event_id"),
        "scale": ((p.get("passport") or {}).get("context") or p.get("context") or {}).get("scale"),
        "direction": ((p.get("passport") or {}).get("context") or p.get("context") or {}).get("direction"),
        "research_passport_sha256": p.get("research_passport_sha256"),
        "handshake_packet_sha256": hp.get("packet_sha256"),
        "measurement_fingerprint_sha256": (hp.get("measurement_fingerprint") or {}).get("sha256"),
        "epistemic_state": (hp.get("epistemic_state") or {}).get("overall_state"),
    } for p, hp in zip(passports, packets)]
    out = {
        "schema": BUNDLE_SCHEMA,
        "core_version": CORE_VERSION,
        "bundle_id": bundle_id,
        "entries": entries,
        "packet_chain_sha256": digest([p["packet_sha256"] for p in packets]),
        "all_scientific_convergence_claims_false": all(p.get("scientific_convergence_claim") is False for p in packets),
        "authority": "MULTISCALE_MEMORY_INDEX_ONLY",
    }
    out["bundle_sha256"] = digest(out)
    return out


def self_test() -> dict[str, Any]:
    fixture = {
        "em122_depth": 3421.4, "ea600_depth": 3419.1, "em122_minus_ea600": 2.3,
        "latitude": -7.845673, "longitude": -14.48023, "heading": 359.0,
        "cadence_jitter": 0.04, "missing_fraction": 0.02,
        "verdict": "H1_REAL_MORPHOLOGY", "target_label": "0012",
    }
    a = build_research_passport(fixture, event_id="SYNTH-A", direction="HEAD_FORWARD", scale="60s", profile="SYNTHETIC_TEST")
    mutated = dict(fixture); mutated.update(verdict="H0_INSTRUMENT_ONLY", target_label="0037")
    b = build_research_passport(mutated, event_id="SYNTH-B", direction="TAIL_REVERSE", scale="60s", profile="SYNTHETIC_TEST")
    blocked = build_research_passport(fixture, event_id="HANNAH-BLOCKED", profile="HANNAH_BODC")
    ha, hb = export_handshake_packet(a), export_handshake_packet(b)
    checks = {
        "contract_hash": digest(load_and_verify_contract()) == PROTOCOL_CONTRACT_SHA256,
        "label_leakage_blocked": ha["measurement_fingerprint"]["sha256"] == hb["measurement_fingerprint"]["sha256"],
        "blocked_hannah_without_raw": blocked["status"] == "BLOCKED",
        "blocked_handshake_null": export_handshake_packet(blocked)["measurement_fingerprint"] is None,
        "handshake_self_verify": verify_handshake_packet(ha),
        "science_claim_false": ha["scientific_convergence_claim"] is False,
    }
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}


if __name__ == "__main__":
    result = self_test()
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    raise SystemExit(0 if result["status"] == "PASS" else 1)
