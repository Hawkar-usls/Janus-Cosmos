#!/usr/bin/env python3
"""Cousteau Synesthetic Memory Core for the Hannah/BODC lane.

This module creates deterministic cross-modal mnemonic "sensory passports"
from measurement features. It is a read-only sidecar:

    MEASUREMENT -> MNEMONIC PASSPORT -> RETRIEVAL / COMPARISON CUE

It must never become:

    MNEMONIC PASSPORT -> SCIENTIFIC VERDICT

The code intentionally excludes hypothesis/verdict/target/template labels from
the measurement fingerprint. Missing raw bytes produce a BLOCKED/NULL passport,
never synthetic measurements.

Python standard library only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import statistics
from pathlib import Path
from typing import Any, Mapping, Sequence

CORE_ID = "JANUS_COUSTEAU_SYNESTHETIC_MEMORY_CORE"
CORE_VERSION = "1.0.0"
MAPPING_SEED = "JANUS_COUSTEAU_SYNESTHESIA_V1__MEASUREMENT_NOT_STORY"
EMBED_DIMS = 16

AUTHORITY_FIREWALL = {
    "mnemonic_layer_is_verdict_authority": False,
    "memory_equals_truth": False,
    "mnemonic_similarity_is_scientific_similarity": False,
    "blocked_input_may_generate_measurement_values": False,
    "raw_bytes_out_rank_mnemonic_mapping": True,
}

# Frozen from the Hannah/BODC preregistered feature vector. Dynamic status and
# acquisition fields are accepted via prefixes below.
MEASUREMENT_KEYS = {
    "em122_depth",
    "ea600_depth",
    "em122_minus_ea600",
    "latitude",
    "longitude",
    "heading",
    "speed",
    "course",
    "turn_rate",
    "delta_depth_dt",
    "delta2_depth_dt2",
    "rolling_depth_median",
    "rolling_depth_mad",
    "depth_local_range",
    "depth_local_slope",
    "timestamp_cadence",
    "cadence_jitter",
    "missing_fraction",
    "null_run_length",
    "identical_value_run_length",
    "outlier_score",
}
DYNAMIC_MEASUREMENT_PREFIXES = (
    "status_",
    "acquisition_",
    "originator_status_",
    "originator_acquisition_",
)

# These concepts may be retained as provenance/context, but may not affect the
# measurement fingerprint. This is the main anti-confirmation-bias gate.
FORBIDDEN_INFLUENCE_TOKENS = {
    "verdict",
    "hypothesis",
    "interpretation",
    "claim",
    "pyramid",
    "target",
    "candidate",
    "anomaly",
    "h0",
    "h1",
    "h2",
    "artificial",
    "natural",
    "control_label",
    "class_label",
    "expected",
    "prediction",
    "story",
}

DIRECTIONS = {"HEAD_FORWARD", "TAIL_REVERSE", "CENTER", "SPACE_REPLAY", "UNKNOWN"}
SCALES = {
    "1s",
    "10s",
    "60s",
    "300s",
    "1800s",
    "7200s",
    "line",
    "ping",
    "beam",
    "footprint",
    "custom",
}

SHAPES = ("circle", "diamond", "hexagon", "triangle", "square", "wave")
TIMBRES = ("sine", "soft_bell", "reed", "glass", "wood", "bowed")
TEXTURES = ("smooth", "silken", "ridged", "granular", "rough", "fragmented")


def canonical_json(obj: Any) -> str:
    return json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def blake2_text(text: str) -> str:
    return hashlib.blake2b(text.encode("utf-8"), digest_size=32).hexdigest()


def _slug(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")


def _flatten(obj: Any, prefix: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {}
    if isinstance(obj, Mapping):
        for key in sorted(obj, key=lambda x: str(x)):
            path = f"{prefix}.{key}" if prefix else str(key)
            out.update(_flatten(obj[key], path))
    else:
        out[prefix] = obj
    return out


def _leaf(path: str) -> str:
    return _slug(path.rsplit(".", 1)[-1])


def _forbidden_path(path: str) -> bool:
    parts = set(_slug(path).split("_"))
    normalized = _slug(path)
    return any(tok in parts or tok in normalized for tok in FORBIDDEN_INFLUENCE_TOKENS)


def _accepted_measurement_path(path: str) -> bool:
    leaf = _leaf(path)
    if _forbidden_path(path):
        return False
    if leaf in MEASUREMENT_KEYS:
        return True
    return leaf.startswith(DYNAMIC_MEASUREMENT_PREFIXES)


def _extract_value_and_unit(value: Any) -> tuple[Any, str | None]:
    if isinstance(value, Mapping) and "value" in value:
        return value.get("value"), value.get("unit")
    return value, None


def _to_number(value: Any) -> float | None:
    value, _ = _extract_value_and_unit(value)
    if isinstance(value, bool):
        return 1.0 if value else 0.0
    if isinstance(value, (int, float)) and math.isfinite(float(value)):
        return float(value)
    if isinstance(value, str):
        v = value.strip().lower()
        if v in {"true", "ok", "valid", "locked", "on", "pass"}:
            return 1.0
        if v in {"false", "bad", "invalid", "unlocked", "off", "fail"}:
            return 0.0
        try:
            f = float(v)
            if math.isfinite(f):
                return f
        except ValueError:
            pass
    return None


def _signed_log_norm(value: float, scale: float) -> float:
    if value == 0:
        return 0.0
    x = math.log1p(abs(value)) / math.log1p(scale)
    return max(-1.0, min(1.0, math.copysign(x, value)))


def _norm_scalar(name: str, value: float) -> list[tuple[str, float]]:
    """Map a scalar to stable bounded dimensions without data-dependent fitting."""
    name = _slug(name)
    if name in {"heading", "course"}:
        a = math.radians(value % 360.0)
        return [(name + "_sin", math.sin(a)), (name + "_cos", math.cos(a))]
    if name == "latitude":
        return [(name, max(-1.0, min(1.0, value / 90.0)))]
    if name == "longitude":
        return [(name, max(-1.0, min(1.0, value / 180.0)))]
    if name == "missing_fraction":
        return [(name, max(0.0, min(1.0, value)) * 2.0 - 1.0)]
    if name in {"em122_depth", "ea600_depth", "rolling_depth_median"}:
        return [(name, _signed_log_norm(value, 12000.0))]
    if name in {"em122_minus_ea600", "depth_local_range"}:
        return [(name, _signed_log_norm(value, 3000.0))]
    if name in {"delta_depth_dt", "delta2_depth_dt2", "depth_local_slope", "turn_rate"}:
        return [(name, math.tanh(value / 10.0))]
    if name in {
        "timestamp_cadence",
        "cadence_jitter",
        "null_run_length",
        "identical_value_run_length",
    }:
        return [(name, _signed_log_norm(value, 7200.0))]
    if name == "speed":
        return [(name, math.tanh(value / 15.0))]
    if name == "rolling_depth_mad":
        return [(name, _signed_log_norm(value, 1000.0))]
    if name == "outlier_score":
        return [(name, math.tanh(value / 6.0))]
    # Dynamic status/acquisition channels and future predeclared scalars.
    return [(name, math.tanh(value))]


def extract_measurement_features(
    payload: Mapping[str, Any],
) -> tuple[dict[str, float], dict[str, str], list[str]]:
    flat = _flatten(payload)
    features: dict[str, float] = {}
    units: dict[str, str] = {}
    excluded: list[str] = []
    for path, raw in flat.items():
        if not _accepted_measurement_path(path):
            if _forbidden_path(path):
                excluded.append(path)
            continue
        leaf = _leaf(path)
        value, unit = _extract_value_and_unit(raw)
        num = _to_number(value)
        if num is None:
            continue
        for normalized_name, normalized_value in _norm_scalar(leaf, num):
            # If the same leaf appears in nested objects, bind it to its path to
            # prevent silent overwrites while preserving deterministic ordering.
            key = normalized_name
            if key in features:
                key = _slug(path) + (
                    "_" + normalized_name if normalized_name != leaf else ""
                )
            features[key] = round(float(normalized_value), 12)
            if unit:
                units[key] = str(unit)
    return (
        dict(sorted(features.items())),
        dict(sorted(units.items())),
        sorted(set(excluded)),
    )


def _projection_slot(feature_name: str, dim: int) -> tuple[int, float]:
    d = hashlib.sha256(
        f"{MAPPING_SEED}|{feature_name}|{dim}".encode("utf-8")
    ).digest()
    slot = int.from_bytes(d[:4], "big") % EMBED_DIMS
    sign = 1.0 if (d[4] & 1) else -1.0
    return slot, sign


def project_features(features: Mapping[str, float]) -> list[float]:
    vec = [0.0] * EMBED_DIMS
    if not features:
        return vec
    for name, value in sorted(features.items()):
        # Three deterministic signed projections per feature reduce accidental
        # collisions while keeping the embedding compact.
        for j in range(3):
            slot, sign = _projection_slot(name, j)
            vec[slot] += sign * float(value)
    norm = math.sqrt(sum(x * x for x in vec))
    if norm:
        vec = [x / norm for x in vec]
    return [round(x, 12) for x in vec]


def _u01(x: float) -> float:
    return max(0.0, min(1.0, (x + 1.0) / 2.0))


def _rgb_from_embedding(v: Sequence[float]) -> dict[str, Any]:
    # Base mnemonic color depends on measurements only. Epistemic uncertainty is
    # represented as a separate overlay and never changes this base color.
    a, b, c = v[0], v[1], v[2]
    r = int(round(40 + 190 * _u01(a)))
    g = int(round(40 + 190 * _u01(b)))
    bl = int(round(40 + 190 * _u01(c)))
    return {"rgb": [r, g, bl], "hex": f"#{r:02X}{g:02X}{bl:02X}"}


def _sensory_channels(
    vec: list[float], completeness: float, direction: str, scale: str
) -> dict[str, Any]:
    if not any(abs(x) > 1e-15 for x in vec):
        return {
            "color": {"rgb": [128, 128, 128], "hex": "#808080"},
            "audio": {"mode": "SILENCE", "frequency_hz": None, "timbre": None},
            "rhythm": {"bpm": None, "pattern_8": "00000000"},
            "texture": {"label": "fog", "roughness": None},
            "glyph": {"shape": "open_ring", "rotation_deg": 0},
            "direction_overlay": direction,
            "scale_overlay": scale,
            "epistemic_overlay": {"completeness": completeness, "fog": 1.0},
        }

    color = _rgb_from_embedding(vec)
    tone = 110.0 * (2.0 ** (3.0 * _u01(vec[3])))
    timbre = TIMBRES[int(_u01(vec[4]) * (len(TIMBRES) - 1))]
    bpm = int(round(40 + 140 * _u01(vec[5])))
    rough = sum(abs(x) for x in vec[6:10]) / 4.0
    texture = TEXTURES[min(len(TEXTURES) - 1, int(rough * len(TEXTURES)))]
    shape = SHAPES[int(_u01(vec[10]) * (len(SHAPES) - 1))]
    rotation = int(round(359 * _u01(vec[11])))
    rhythm_digest = hashlib.sha256(canonical_json(vec).encode("utf-8")).digest()[0]
    pattern = f"{rhythm_digest:08b}"
    return {
        "color": color,
        "audio": {"mode": "TONE", "frequency_hz": round(tone, 3), "timbre": timbre},
        "rhythm": {"bpm": bpm, "pattern_8": pattern},
        "texture": {"label": texture, "roughness": round(rough, 6)},
        "glyph": {"shape": shape, "rotation_deg": rotation},
        "direction_overlay": direction,
        "scale_overlay": scale,
        "epistemic_overlay": {
            "completeness": round(completeness, 6),
            "fog": round(1.0 - completeness, 6),
            "rule": "OVERLAY_ONLY__DOES_NOT_MODIFY_MEASUREMENT_COLOR_OR_SIMILARITY",
        },
    }


def _source_identity(source: Any, raw_bytes: bytes | None = None) -> dict[str, Any]:
    if raw_bytes is not None:
        return {
            "source_hash_kind": "RAW_BYTES_SHA256",
            "source_sha256": hashlib.sha256(raw_bytes).hexdigest(),
            "source_blake2b_256": hashlib.blake2b(
                raw_bytes, digest_size=32
            ).hexdigest(),
        }
    text = canonical_json(source)
    return {
        "source_hash_kind": "CANONICAL_JSON_SHA256_NOT_RAW_FILE_HASH",
        "source_sha256": sha256_text(text),
        "source_blake2b_256": blake2_text(text),
    }


def build_passport(
    payload: Mapping[str, Any],
    *,
    direction: str = "UNKNOWN",
    scale: str = "custom",
    provenance: Mapping[str, Any] | None = None,
    raw_bytes: bytes | None = None,
    expected_feature_count: int | None = None,
) -> dict[str, Any]:
    if direction not in DIRECTIONS:
        raise ValueError(f"direction must be one of {sorted(DIRECTIONS)}")
    if scale not in SCALES:
        raise ValueError(f"scale must be one of {sorted(SCALES)}")

    features, units, excluded = extract_measurement_features(payload)
    source_id = _source_identity(payload, raw_bytes=raw_bytes)
    measurement_material = {
        "mapping_seed": MAPPING_SEED,
        "features": features,
        "units": units,
    }
    measurement_digest = sha256_text(canonical_json(measurement_material))
    second_digest = blake2_text(canonical_json(measurement_material))
    vec = project_features(features)
    denom = max(
        1,
        expected_feature_count
        if expected_feature_count is not None
        else len(MEASUREMENT_KEYS),
    )
    completeness = min(1.0, len(features) / denom)
    channels = _sensory_channels(vec, completeness, direction, scale)

    passport_material = {
        "core": CORE_ID,
        "version": CORE_VERSION,
        "measurement_digest": measurement_digest,
        "direction": direction,
        "scale": scale,
        "source_sha256": source_id["source_sha256"],
    }
    passport_sha = sha256_text(canonical_json(passport_material))
    return {
        "schema": "janus.cosmos.cousteau.synesthetic_memory_passport.v1",
        "core_id": CORE_ID,
        "core_version": CORE_VERSION,
        "mapping_seed": MAPPING_SEED,
        "status": (
            "MEASUREMENT_MNEMONIC_READY"
            if features
            else "NO_ACCEPTED_MEASUREMENT_FEATURES"
        ),
        "passport_id": f"CSMC1-{passport_sha[:20]}",
        "passport_sha256": passport_sha,
        "collision_guard_blake2b_256": second_digest,
        "source_identity": source_id,
        "measurement_fingerprint": {
            "sha256": measurement_digest,
            "feature_count": len(features),
            "features": features,
            "units": units,
            "embedding": vec,
        },
        "sensory_channels": channels,
        "context": {"direction": direction, "scale": scale},
        "epistemic": {
            "expected_feature_count": expected_feature_count,
            "completeness": round(completeness, 6),
            "excluded_forbidden_influence_paths": excluded,
        },
        "provenance": dict(provenance or {}),
        "authority_firewall": AUTHORITY_FIREWALL,
        "warnings": [
            "MNEMONIC_SIMILARITY_IS_A_RETRIEVAL_CUE_NOT_A_SCIENTIFIC_VERDICT",
            "PASSPORT_IS_NOT_REVERSIBLE_TO_RAW_DATA__USE_SOURCE_POINTER_AND_HASH",
            "HYPOTHESIS_VERDICT_TARGET_AND_TEMPLATE_LABELS_DO_NOT_INFLUENCE_MEASUREMENT_FINGERPRINT",
        ],
    }


def build_blocked_passport(
    *,
    source_receipt: Mapping[str, Any],
    reason: str = "RAW_BYTES_NOT_AVAILABLE",
    provenance: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    source_id = _source_identity(source_receipt)
    material = {
        "core": CORE_ID,
        "version": CORE_VERSION,
        "status": "BLOCKED_NULL",
        "reason": reason,
        "source_sha256": source_id["source_sha256"],
    }
    digest = sha256_text(canonical_json(material))
    channels = _sensory_channels(
        [0.0] * EMBED_DIMS, 0.0, "UNKNOWN", "custom"
    )
    return {
        "schema": "janus.cosmos.cousteau.synesthetic_memory_passport.v1",
        "core_id": CORE_ID,
        "core_version": CORE_VERSION,
        "status": "BLOCKED_NULL",
        "reason": reason,
        "passport_id": f"CSMC1-BLOCKED-{digest[:16]}",
        "passport_sha256": digest,
        "source_identity": source_id,
        "measurement_claims_allowed": False,
        "measurement_fingerprint": None,
        "sensory_channels": channels,
        "provenance": dict(provenance or {}),
        "authority_firewall": AUTHORITY_FIREWALL,
        "hard_rule": "NO_RAW_BYTES_OR_MEASUREMENTS__NO_SYNTHETIC_SENSORY_MEASUREMENT",
    }


def _cosine(a: Sequence[float], b: Sequence[float]) -> float | None:
    if len(a) != len(b) or not a:
        return None
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(x * x for x in b))
    if na == 0 or nb == 0:
        return None
    return sum(x * y for x, y in zip(a, b)) / (na * nb)


def compare_passports(a: Mapping[str, Any], b: Mapping[str, Any]) -> dict[str, Any]:
    fa = (a.get("measurement_fingerprint") or {}).get("features") or {}
    fb = (b.get("measurement_fingerprint") or {}).get("features") or {}
    common = sorted(set(fa) & set(fb))
    union = sorted(set(fa) | set(fb))
    if common:
        # 1 - RMS distance / 2 maps identical -> 1, maximally opposite -> 0.
        rms = math.sqrt(
            sum((float(fa[k]) - float(fb[k])) ** 2 for k in common) / len(common)
        )
        common_similarity = max(0.0, min(1.0, 1.0 - rms / 2.0))
    else:
        common_similarity = None
    va = (a.get("measurement_fingerprint") or {}).get("embedding") or []
    vb = (b.get("measurement_fingerprint") or {}).get("embedding") or []
    mnemonic_cosine = _cosine(va, vb)
    return {
        "schema": "janus.cosmos.cousteau.synesthetic_memory_comparison.v1",
        "status": "COMPARABLE" if common else "INSUFFICIENT_COMMON_FEATURES",
        "common_feature_count": len(common),
        "union_feature_count": len(union),
        "feature_overlap_fraction": (
            round(len(common) / len(union), 6) if union else 0.0
        ),
        "common_measurement_similarity": (
            round(common_similarity, 6) if common_similarity is not None else None
        ),
        "mnemonic_embedding_cosine": (
            round(mnemonic_cosine, 6) if mnemonic_cosine is not None else None
        ),
        "authority": "RETRIEVAL_AND_REVIEW_PRIORITY_ONLY",
        "scientific_convergence_claim": False,
        "hard_rule": "REQUIRE_TIME_SPACE_GROUND_FIXED_REPLICATION_GATES_OUTSIDE_THIS_CORE",
    }


def rank_cross_front_pairs(
    head: Sequence[Mapping[str, Any]],
    tail: Sequence[Mapping[str, Any]],
    *,
    top_k: int = 20,
) -> list[dict[str, Any]]:
    pairs: list[dict[str, Any]] = []
    for ha in head:
        for tb in tail:
            cmp = compare_passports(ha, tb)
            if cmp["common_measurement_similarity"] is None:
                continue
            score = (
                cmp["common_measurement_similarity"]
                * cmp["feature_overlap_fraction"]
            )
            pairs.append(
                {
                    "head_passport_id": ha.get("passport_id"),
                    "tail_passport_id": tb.get("passport_id"),
                    "review_priority_score": round(score, 6),
                    "comparison": cmp,
                    "scientific_convergence_claim": False,
                }
            )
    pairs.sort(
        key=lambda x: (
            -x["review_priority_score"],
            str(x["head_passport_id"]),
            str(x["tail_passport_id"]),
        )
    )
    return pairs[: max(0, top_k)]


def aggregate_numeric_records(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Create a conservative window summary from accepted scalar fields.

    This helper does not invent temporal derivatives. It summarizes only
    accepted numeric fields actually present in records using medians and a
    missingness marker. Derivatives/slope/jitter should be supplied by the
    dedicated Hannah measurement pipeline when available.
    """
    if not records:
        return {"missing_fraction": 1.0}
    series: dict[str, list[float]] = {}
    accepted_leafs = set(MEASUREMENT_KEYS)
    rows_with_measurement = 0
    for rec in records:
        flat = _flatten(rec)
        row_has_measurement = False
        for path, raw in flat.items():
            leaf = _leaf(path)
            if (
                leaf not in accepted_leafs
                and not leaf.startswith(DYNAMIC_MEASUREMENT_PREFIXES)
            ):
                continue
            num = _to_number(raw)
            if num is not None:
                series.setdefault(leaf, []).append(num)
                row_has_measurement = True
        if row_has_measurement:
            rows_with_measurement += 1
    out: dict[str, Any] = {}
    for name, vals in sorted(series.items()):
        if vals:
            out[name] = statistics.median(vals)
    out["missing_fraction"] = 1.0 - (rows_with_measurement / len(records))
    return out


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path | None, obj: Any) -> None:
    text = json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    if path is None:
        print(text, end="")
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")


def self_test() -> dict[str, Any]:
    base = {
        "em122_depth": 3421.4,
        "ea600_depth": 3419.1,
        "latitude": -7.845673,
        "longitude": -14.48023,
        "heading": 359.0,
        "missing_fraction": 0.02,
        "verdict": "H1_REAL_MORPHOLOGY",
        "target_label": "DO_NOT_LEAK",
    }
    shuffled = {
        "target_label": "CHANGED",
        "missing_fraction": 0.02,
        "heading": 359.0,
        "longitude": -14.48023,
        "latitude": -7.845673,
        "ea600_depth": 3419.1,
        "em122_depth": 3421.4,
        "verdict": "H0_INSTRUMENT_ONLY",
    }
    a = build_passport(base, direction="HEAD_FORWARD", scale="60s")
    b = build_passport(shuffled, direction="TAIL_REVERSE", scale="60s")
    c = build_passport(
        {**base, "em122_depth": 4100.0},
        direction="HEAD_FORWARD",
        scale="60s",
    )
    blocked = build_blocked_passport(
        source_receipt={"status": "BLOCKED_RAW_BYTES_NOT_MOUNTED"}
    )
    checks = {
        "verdict_target_leakage_blocked": (
            a["measurement_fingerprint"]["sha256"]
            == b["measurement_fingerprint"]["sha256"]
        ),
        "direction_overlay_does_not_change_measurement": (
            a["measurement_fingerprint"]["sha256"]
            == b["measurement_fingerprint"]["sha256"]
        ),
        "measurement_change_changes_fingerprint": (
            a["measurement_fingerprint"]["sha256"]
            != c["measurement_fingerprint"]["sha256"]
        ),
        "blocked_has_no_measurement_fingerprint": blocked["measurement_fingerprint"]
        is None,
        "blocked_disallows_measurement_claims": blocked["measurement_claims_allowed"]
        is False,
        "self_similarity_is_one": (
            compare_passports(a, a)["common_measurement_similarity"] == 1.0
        ),
        "source_hash_kind_is_explicit": (
            a["source_identity"]["source_hash_kind"]
            == "CANONICAL_JSON_SHA256_NOT_RAW_FILE_HASH"
        ),
        "authority_firewall": (
            a["authority_firewall"]["mnemonic_layer_is_verdict_authority"] is False
        ),
    }
    return {
        "core_id": CORE_ID,
        "core_version": CORE_VERSION,
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="JANUS Cousteau Synesthetic Memory Core"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("self-test")

    p = sub.add_parser("passport")
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path)
    p.add_argument("--direction", choices=sorted(DIRECTIONS), default="UNKNOWN")
    p.add_argument("--scale", choices=sorted(SCALES), default="custom")
    p.add_argument("--expected-features", type=int)

    b = sub.add_parser("blocked")
    b.add_argument("--receipt", type=Path, required=True)
    b.add_argument("--output", type=Path)

    c = sub.add_parser("compare")
    c.add_argument("a", type=Path)
    c.add_argument("b", type=Path)
    c.add_argument("--output", type=Path)

    args = parser.parse_args()
    if args.command == "self-test":
        result = self_test()
        _write_json(None, result)
        return 0 if result["status"] == "PASS" else 1
    if args.command == "passport":
        payload = _read_json(args.input)
        if not isinstance(payload, Mapping):
            raise SystemExit("passport input must be a JSON object")
        result = build_passport(
            payload,
            direction=args.direction,
            scale=args.scale,
            provenance={"input_path": str(args.input)},
            raw_bytes=args.input.read_bytes(),
            expected_feature_count=args.expected_features,
        )
        _write_json(args.output, result)
        return 0
    if args.command == "blocked":
        receipt = _read_json(args.receipt)
        result = build_blocked_passport(
            source_receipt=receipt,
            provenance={"receipt_path": str(args.receipt)},
        )
        _write_json(args.output, result)
        return 0
    if args.command == "compare":
        result = compare_passports(_read_json(args.a), _read_json(args.b))
        _write_json(args.output, result)
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
