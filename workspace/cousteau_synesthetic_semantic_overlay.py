#!/usr/bin/env python3
"""Cousteau-domain semantic overlay for Synesthetic Memory Core v1.

The base core creates a deterministic identity fingerprint. This adapter makes
that fingerprint easier to remember in the Hannah/BODC lane by binding specific
marine-measurement concepts to stable sensory metaphors.

It never modifies measurement_fingerprint and is never used as verdict evidence.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


SCALE_OCTAVE = {
    "7200s": -2,
    "1800s": -1,
    "300s": 0,
    "60s": 1,
    "10s": 2,
    "1s": 3,
    "line": -2,
    "ping": 0,
    "beam": 1,
    "footprint": 2,
    "custom": 0,
}
PAN = {
    "HEAD_FORWARD": -0.75,
    "TAIL_REVERSE": 0.75,
    "CENTER": 0.0,
    "SPACE_REPLAY": 0.25,
    "UNKNOWN": 0.0,
}


def _clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def _f(features: Mapping[str, Any], key: str) -> float | None:
    value = features.get(key)
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _roughness(features: Mapping[str, Any]) -> float | None:
    vals = []
    for key in ("rolling_depth_mad", "depth_local_range", "depth_local_slope", "outlier_score"):
        v = _f(features, key)
        if v is not None:
            vals.append(abs(v))
    if not vals:
        return None
    return _clamp01(sum(vals) / len(vals))


def _texture_label(r: float | None) -> str:
    if r is None:
        return "unknown"
    if r < 0.15:
        return "glass_smooth"
    if r < 0.35:
        return "silken"
    if r < 0.55:
        return "fine_ridged"
    if r < 0.75:
        return "granular"
    return "coarse_rough"


def build_cousteau_semantic_overlay(passport: Mapping[str, Any]) -> dict[str, Any]:
    fp = passport.get("measurement_fingerprint")
    if not isinstance(fp, Mapping):
        return {
            "status": "BLOCKED_OR_NO_MEASUREMENT",
            "depth_register": {"frequency_hz": None, "meaning": "silence"},
            "cross_sensor_beating": {"beat_hz": None, "meaning": "unknown"},
            "texture": {"roughness": None, "label": "fog"},
            "cadence": {"stability": None, "meaning": "no measured pulse"},
            "missingness": {"fraction": None, "fog": 1.0},
            "track_pan": 0.0,
            "mnemonic_sentence": "SILENCE | FOG | NO MEASUREMENT",
            "authority": "NONE",
        }

    features = fp.get("features") or {}
    context = passport.get("context") or {}
    direction = str(context.get("direction", "UNKNOWN"))
    scale = str(context.get("scale", "custom"))

    # Normalized depth is positive for ordinary positive ocean depth. Deeper
    # water moves downward in pitch across roughly three octaves.
    depth = _f(features, "em122_depth")
    depth_source = "EM122"
    if depth is None:
        depth = _f(features, "ea600_depth")
        depth_source = "EA600"
    if depth is None:
        depth_hz = None
        depth_band = "unknown"
    else:
        dn = _clamp01(abs(depth))
        octave = SCALE_OCTAVE.get(scale, 0)
        depth_hz = 880.0 * (2.0 ** (-3.0 * dn)) * (2.0 ** (octave / 12.0))
        depth_hz = round(depth_hz, 3)
        depth_band = "deep_low" if dn > 0.66 else "mid" if dn > 0.33 else "shallow_high"

    # Explicit EM122-EA600 disagreement becomes audible as beat frequency.
    disagreement = _f(features, "em122_minus_ea600")
    if disagreement is None:
        beat_hz = None
        beat_band = "unknown"
    else:
        d = _clamp01(abs(disagreement))
        beat_hz = round(0.5 + 11.5 * d, 3)
        beat_band = "tight" if d < 0.2 else "noticeable" if d < 0.55 else "strong"

    rough = _roughness(features)
    texture = _texture_label(rough)

    jitter = _f(features, "cadence_jitter")
    cadence_stability = None if jitter is None else round(_clamp01(1.0 - abs(jitter)), 6)
    cadence_meaning = (
        "unknown"
        if cadence_stability is None
        else "clocklike" if cadence_stability > 0.85
        else "wavering" if cadence_stability > 0.55
        else "irregular"
    )

    missing_norm = _f(features, "missing_fraction")
    # Base core maps raw missing_fraction [0,1] -> normalized [-1,1].
    missing_fraction = None if missing_norm is None else _clamp01((missing_norm + 1.0) / 2.0)
    fog = 1.0 if missing_fraction is None else missing_fraction

    pan = PAN.get(direction, 0.0)
    bits = [
        f"DEPTH:{depth_band}",
        f"BEAT:{beat_band}",
        f"TEXTURE:{texture}",
        f"PULSE:{cadence_meaning}",
        f"FOG:{fog:.3f}",
        f"PAN:{pan:+.2f}",
        f"SCALE:{scale}",
    ]
    sentence = " | ".join(bits)
    digest = hashlib.sha256(sentence.encode("utf-8")).hexdigest()

    return {
        "status": "READY",
        "depth_register": {
            "source": depth_source if depth is not None else None,
            "frequency_hz": depth_hz,
            "band": depth_band,
            "rule": "DEEPER_MEASURED_DEPTH__LOWER_MNEMONIC_REGISTER",
        },
        "cross_sensor_beating": {
            "source": "EM122_MINUS_EA600",
            "beat_hz": beat_hz,
            "band": beat_band,
            "rule": "LARGER_MEASURED_DISAGREEMENT__STRONGER_MNEMONIC_BEATING",
        },
        "texture": {
            "roughness": None if rough is None else round(rough, 6),
            "label": texture,
            "basis": ["rolling_depth_MAD", "depth_local_range", "depth_local_slope", "outlier_score"],
            "warning": "MNEMONIC_TEXTURE_NE_SEAFLOOR_MORPHOLOGY_CLASS",
        },
        "cadence": {
            "stability": cadence_stability,
            "meaning": cadence_meaning,
            "basis": "cadence_jitter",
        },
        "missingness": {
            "fraction": None if missing_fraction is None else round(missing_fraction, 6),
            "fog": round(fog, 6),
            "rule": "MISSINGNESS_IS_VISIBLE_AS_FOG_NOT_FILLED_WITH_SYNTHETIC_DATA",
        },
        "track_pan": pan,
        "scale_octave_semitones": SCALE_OCTAVE.get(scale, 0),
        "mnemonic_sentence": sentence,
        "semantic_overlay_sha256": digest,
        "authority": "RETRIEVAL_CUE_ONLY",
        "scientific_claim": False,
    }


def enrich_passport(passport: Mapping[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(dict(passport))
    before = json.dumps(out.get("measurement_fingerprint"), sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    out["cousteau_semantic_overlay"] = build_cousteau_semantic_overlay(out)
    after = json.dumps(out.get("measurement_fingerprint"), sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    if before != after:
        raise RuntimeError("semantic overlay mutated measurement_fingerprint")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Add Cousteau semantic mnemonic overlay to a sensory passport")
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    passport = json.loads(args.input.read_text(encoding="utf-8"))
    result = enrich_passport(passport)
    text = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
