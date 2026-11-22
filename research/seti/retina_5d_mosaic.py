"""Research-only JANUS COSMOS 5D mosaic-kaleidoscope advisory.

Radiometry is evaluated ONLY on native numeric time/frequency power arrays.
Mirror views are correlated transforms of the same samples: never extra votes.
No threshold, baseline alert, source classification, or production model is changed.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

import numpy as np

SLOPES = (-0.125, 0.0, 0.125)  # native channel / native time sample; development-only
DIMS = (
    "D1_FORWARD_BACK_TEMPORAL",
    "D2_LEFT_RIGHT_MIRROR",
    "D3_ASCEND_DIVE_NATIVE_MOSAIC",
    "D4_EXPAND_CONTRACT_ALTERNATIVES",
    "D5_CLAIM_SOURCE_PROVENANCE",
)


def _power(value: Any) -> np.ndarray:
    x = np.asarray(value, dtype=np.float64)
    if x.ndim != 2 or x.shape[0] < 16 or x.shape[1] < 16:
        raise ValueError("expected 2D native (time, frequency) with >=16 each")
    if not np.isfinite(x).all():
        raise ValueError("nonfinite native measurement; no imputation")
    return x


def _digest(x: np.ndarray) -> str:
    stable = np.ascontiguousarray(x.astype("<f8", copy=False))
    return hashlib.sha256(json.dumps(list(x.shape)).encode() + stable.tobytes()).hexdigest()


def _normalize_rows(x: np.ndarray) -> np.ndarray:
    median = np.median(x, axis=1, keepdims=True)
    mad = 1.4826 * np.median(np.abs(x - median), axis=1, keepdims=True)
    return (x - median) / np.maximum(mad, 1e-8)


def _path_score(z: np.ndarray, path: np.ndarray) -> tuple[float, list[float]]:
    t = np.arange(z.shape[0])
    local = z[t, path] - (z[t, path - 2] + z[t, path + 2]) * 0.5
    blocks = [float(np.mean(b)) for b in np.array_split(local, 4)]
    # A lone pulse cannot dominate; require temporal support, not a max-pixel vote.
    return float(0.5 * np.mean(local) + 0.5 * min(blocks)), blocks


def _best_path(z: np.ndarray) -> tuple[np.ndarray, float, list[float], float]:
    ntime, nfreq = z.shape
    best: tuple[np.ndarray, float, list[float], float] | None = None
    for slope in SLOPES:
        drift = np.rint(slope * (np.arange(ntime) - (ntime - 1) / 2.0)).astype(int)
        for center in range(2 - int(drift.min()), nfreq - 2 - int(drift.max())):
            path = center + drift
            score, blocks = _path_score(z, path)
            if best is None or score > best[1]:
                best = (path, score, blocks, slope)
    if best is None:
        raise ValueError("no native legal path with two sidebands")
    return best


def _heatmap(z: np.ndarray) -> list[list[float]]:
    rows = np.array_split(np.arange(z.shape[0]), 4)
    cols = np.array_split(np.arange(z.shape[1]), 4)
    return [[round(float(np.mean(np.maximum(z[np.ix_(r, c)], 0))), 5)
             for c in cols] for r in rows]


def analyze_spectrogram(on_power: Any, *, off_power: Any | None = None,
                        source_id: str | None = None, reference_id: str | None = None,
                        native_grid_verified: bool = False,
                        baseline_alert: bool = False,
                        prior_evidence_sha256: str | None = None) -> dict[str, Any]:
    """Advisory 5D pass. Baseline output is returned UNCHANGED, even for RFI.

    OFF must already be matched to the same native time/frequency grid. If not
    independently qualified, its resemblance is only an unresolved hypothesis.
    """
    on = _power(on_power)
    off = _power(off_power) if off_power is not None else None
    if off is not None and off.shape != on.shape:
        raise ValueError("OFF grid mismatch; never resample or shift silently")
    z = _normalize_rows(on)
    path, score, blocks, slope = _best_path(z)

    rev_score, _ = _path_score(z[::-1], path[::-1])
    mirror_score, _ = _path_score(z[:, ::-1], z.shape[1] - 1 - path)
    if not (np.isclose(score, rev_score, atol=1e-10) and
            np.isclose(score, mirror_score, atol=1e-10)):
        raise AssertionError("kaleidoscope inverse-map parity failure")

    off_score = None
    if off is not None:
        off_score, _ = _path_score(_normalize_rows(off), path)

    independent_pair = bool(off is not None and source_id and reference_id and
                            source_id != reference_id and native_grid_verified)
    shared_like = bool(off_score is not None and score > 0 and
                       off_score >= 0.7 * score and off_score > 1.0)
    hypotheses = ["NOISE_OR_NONSTATIONARY_STRUCTURE", "COHERENT_NARROWBAND_STRUCTURE"]
    if off is not None:
        hypotheses.append("SHARED_SOURCE_REFERENCE_STRUCTURE" if shared_like
                          else "SOURCE_REFERENCE_DIFFERENCE_UNCALIBRATED")

    source_hash = _digest(on)
    ref_hash = _digest(off) if off is not None else None
    evidence = {
        "source_sha256": source_hash,
        "reference_sha256": ref_hash,
        "source_id": source_id,
        "reference_id": reference_id,
        "caller_attested_distinct_native_scan_pair": independent_pair,
        "native_shape": list(on.shape),
        "fixed_path_bins": path.tolist(),
        "fixed_path_slope_bins_per_sample": slope,
        "native_research_score": round(score, 10),
        "four_temporal_mosaic_tile_scores": [round(x, 10) for x in blocks],
        "reference_same_path_score": round(off_score, 10) if off_score is not None else None,
        "mirror_parity": True,
        "reverse_parity": True,
        "visual_mosaic_4x4": _heatmap(z),
        "shared_pattern_advisory_only": shared_like,
        "hypotheses_not_verdicts": hypotheses,
    }
    evidence_hash = hashlib.sha256(json.dumps(evidence, sort_keys=True,
                               separators=(",", ":")).encode()).hexdigest()
    return {
        "schema": "janus.cosmos.retina_5d_mosaic.advisory.v1",
        "mode": "RESEARCH_SHADOW_ADVISORY_ONLY",
        "dimensions": DIMS,
        "origin_sha256": source_hash,
        "origin_prime_evidence_sha256": evidence_hash,
        "spiral_return": ("NO_INFORMATION_GAIN__STOP" if
                           prior_evidence_sha256 == evidence_hash else
                           "EVIDENCE_REFINED_NOT_SOURCE_TRUTH"),
        "path_history": list(DIMS),
        "mosaic_views_are_independent_evidence": False,
        "raw_measurement_mutated": False,
        "production_adoption": False,
        "baseline_alert_before": bool(baseline_alert),
        "baseline_alert_after": bool(baseline_alert),
        "science_gate": "HOLD_INDEPENDENT_SOURCES_AND_FRESH_VALIDATION",
        "evidence": evidence,
    }


def inspect_retina_video_tiles(frames: Any, *, tile_rows: int = 4,
                               tile_cols: int = 4) -> dict[str, Any]:
    """Separate visual QC for TOPA Retina frames; never a radio detector.

    Input: same-sized grayscale sampled frames (time, height, width). TOPA's
    original frame hashes/PTS and decoded source video stay authoritative.
    """
    data = np.asarray(frames, dtype=np.float64)
    if data.ndim != 3 or len(data) < 2 or min(data.shape[1:]) < max(tile_rows, tile_cols):
        raise ValueError("video frames must be (time,height,width) and cover tiles")
    if not np.isfinite(data).all() or tile_rows < 2 or tile_cols < 2:
        raise ValueError("invalid visual frames or tile dimensions")
    rows = np.array_split(np.arange(data.shape[1]), tile_rows)
    cols = np.array_split(np.arange(data.shape[2]), tile_cols)
    means = np.array([[np.mean(frame[np.ix_(r, c)]) for c in cols for r in rows]
                      for frame in data])
    # Pixel mosaic order is a visual coordinate, NOT a radio-frequency axis.
    contrast = means - np.median(means, axis=1, keepdims=True)
    consistency = np.mean(contrast > 0, axis=0).reshape(tile_cols, tile_rows).T
    return {
        "schema": "janus.retina.video_mosaic.visual_qc.v1",
        "sampled_frames": int(len(data)),
        "tile_positive_contrast_fraction": consistency.round(6).tolist(),
        "mirror_inverse_map_parity": bool(np.array_equal(means[:, ::-1][:, ::-1], means)),
        "mosaic_votes_independent": False,
        "frame_source_provenance_required": True,
        "class": "VISUAL_PATTERN_QC_ONLY_NOT_SIGNAL_ORIGIN",
    }
