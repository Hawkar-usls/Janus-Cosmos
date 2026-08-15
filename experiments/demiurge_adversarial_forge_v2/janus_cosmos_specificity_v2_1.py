#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np
from scipy import ndimage

import janus_cosmos_core_v2 as core


ROOT = Path(__file__).resolve().parent
PROTOCOL_PATH = ROOT / "SPECIFICITY_PROTOCOL_v2_1.json"
EXPECTED_PATH = ROOT / "EXPECTED_PROTOCOL_v2_1.json"


def canonical_json(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def canonical_sha256(obj) -> str:
    return hashlib.sha256(canonical_json(obj).encode("utf-8")).hexdigest()


def normalized_text_sha256(path: Path) -> str:
    text = path.read_text(encoding="utf-8").replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def angular_separation_deg(a: tuple[float, float], b: tuple[float, float]) -> float:
    ra1, dec1 = map(math.radians, a)
    ra2, dec2 = map(math.radians, b)
    dot = math.sin(dec1) * math.sin(dec2) + math.cos(dec1) * math.cos(dec2) * math.cos(ra1 - ra2)
    return math.degrees(math.acos(max(-1.0, min(1.0, dot))))


def regenerate_sky_controls(protocol: dict) -> list[dict]:
    cfg = protocol["real_sky_controls"]
    gen = cfg["generation"]
    seed = gen["seed_string"]
    dec_cap = math.sin(math.radians(float(gen["declination_cap_deg"])))
    target = protocol["orion_target"]["center_j2000"]
    target_xy = (float(target["ra_deg"]), float(target["dec_deg"]))
    rows: list[dict] = []
    attempt = 0
    while len(rows) < len(cfg["centers"]):
        digest = hashlib.sha256(f"{seed}|{attempt}".encode("utf-8")).digest()
        ra = int.from_bytes(digest[:8], "big") / 2**64 * 360.0
        unit = int.from_bytes(digest[8:16], "big") / 2**64
        dec = math.degrees(math.asin((2.0 * unit - 1.0) * dec_cap))
        point = (ra, dec)
        far_target = angular_separation_deg(point, target_xy) >= float(gen["minimum_orion_separation_deg"])
        separated = all(
            angular_separation_deg(point, (float(row["ra_deg"]), float(row["dec_deg"])))
            >= float(gen["minimum_pairwise_separation_deg"])
            for row in rows
        )
        if far_target and separated:
            rows.append(
                {
                    "id": f"SKYCTRL_{len(rows) + 1:02d}",
                    "ra_deg": round(ra, 12),
                    "dec_deg": round(dec, 12),
                    "generator_attempt": attempt,
                }
            )
        attempt += 1
        if attempt > 100000:
            raise RuntimeError("blind-control regeneration did not converge")
    if attempt != int(gen["accepted_attempt_count"]):
        raise RuntimeError("blind-control accepted-attempt count drift")
    return rows


def verify_protocol_sources(protocol: dict | None = None) -> dict:
    protocol = protocol or json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    expected = json.loads(EXPECTED_PATH.read_text(encoding="utf-8"))
    if protocol.get("status") != "FROZEN":
        raise RuntimeError("v2.1 specificity protocol is not frozen")
    actual_protocol_sha = canonical_sha256(protocol)
    if actual_protocol_sha != expected["protocol_sha256"]:
        raise RuntimeError("v2.1 specificity protocol hash mismatch")
    for name, wanted in expected["source_sha256"].items():
        path = ROOT / name
        if not path.is_file() or normalized_text_sha256(path) != wanted:
            raise RuntimeError(f"v2.1 implementation source hash mismatch: {name}")
    negative = protocol["negative_parent_certificate"]
    negative_path = ROOT / negative["path"]
    if core.sha256_file(negative_path) != negative["sha256"]:
        raise RuntimeError("v2.0.2 negative specificity certificate hash mismatch")
    cert = json.loads(negative_path.read_text(encoding="utf-8"))
    if cert.get("classification") != negative["required_classification"]:
        raise RuntimeError("v2.0.2 negative specificity classification mismatch")
    if regenerate_sky_controls(protocol) != protocol["real_sky_controls"]["centers"]:
        raise RuntimeError("v2.1 real-sky blind-control coordinate drift")
    real_count = len(protocol["real_sky_controls"]["centers"])
    hst_ids = protocol["hst_real_controls"]["field_ids"]
    if real_count < 20 or len(hst_ids) < 20 or len(set(hst_ids)) != len(hst_ids):
        raise RuntimeError("v2.1 requires at least 20 unique controls per target family")
    if protocol["hst_target"]["id"] in hst_ids:
        raise RuntimeError("HST target leaked into real-control cohort")
    return {
        "protocol_sha256": actual_protocol_sha,
        "real_sky_control_count": real_count,
        "hst_control_count": len(hst_ids),
    }


def empirical_test_with_tail(
    x: np.ndarray,
    genome: dict,
    model: str,
    test_nulls: int,
    cal_nulls: int,
    seeds: Sequence[int],
    seed_parts: Sequence,
    progress: Callable[[int, int, int], None] | None = None,
) -> dict:
    center, scale = core.calibrate(x, genome, model, cal_nulls, seed_parts)
    observed = core.weighted_std_dist(core.geometry(x), center, scale, genome)
    base, extra = divmod(test_nulls, len(seeds))
    values: list[float] = []
    chunks: list[dict] = []
    for seed_index, seed in enumerate(seeds):
        count = base + (1 if seed_index < extra else 0)
        rng = np.random.default_rng(core.stable_seed("test-v2.1", *seed_parts, model, seed))
        chunk: list[float] = []
        for index in range(count):
            stat = core.weighted_std_dist(core.geometry(core.surrogate(x, rng, model, genome)), center, scale, genome)
            values.append(stat)
            chunk.append(stat)
            if progress and ((index + 1) % max(1, count // 4) == 0 or index + 1 == count):
                progress(int(seed), index + 1, count)
        chunks.append({"seed": int(seed), "null_count": count, "median_stat": float(np.median(chunk))})
    null = np.asarray(values, dtype=np.float64)
    ge_count = int(np.count_nonzero(null >= observed))
    q95 = float(np.quantile(null, 0.95, method="higher"))
    q99 = float(np.quantile(null, 0.99, method="higher"))
    return {
        "observed_stat": float(observed),
        "null_count": int(null.size),
        "ge_count": ge_count,
        "p_empirical": float((ge_count + 1) / (null.size + 1)),
        "null_min": float(null.min()),
        "null_median": float(np.median(null)),
        "null_q95": q95,
        "null_q99": q99,
        "null_max": float(null.max()),
        "tail_ratio_q99": float(observed / max(q99, 1e-12)),
        "calibration_nulls": int(cal_nulls),
        "seed_chunks": chunks,
    }


def band_tail_effect(models: dict) -> float:
    required = ("phase_iaaft", "block_shuffle")
    return float(min(float(models[name]["tail_ratio_q99"]) for name in required))


def field_tail_effect(band_analyses: Iterable[dict]) -> float:
    values = [band_tail_effect(item["models"]) for item in band_analyses]
    if not values:
        raise ValueError("field tail effect requires at least one band")
    return float(min(values))


def real_field_rank(target_score: float, control_scores: Sequence[float]) -> dict:
    controls = np.asarray(control_scores, dtype=np.float64)
    if controls.size < 20 or not np.all(np.isfinite(controls)) or not np.isfinite(target_score):
        raise RuntimeError("real-field rank requires target plus at least 20 finite controls")
    exceedances = int(np.count_nonzero(controls >= float(target_score)))
    return {
        "target_score": float(target_score),
        "control_count": int(controls.size),
        "control_exceedances": exceedances,
        "p_empirical": float((exceedances + 1) / (controls.size + 1)),
        "control_min": float(controls.min()),
        "control_median": float(np.median(controls)),
        "control_max": float(controls.max()),
        "outperforms_all_controls": bool(exceedances == 0),
    }


def corridor_spec(seed_parts: Sequence, length: float, half_width: float, size: int = core.IMAGE_SIZE) -> dict:
    rng = np.random.default_rng(core.stable_seed(*seed_parts))
    for _ in range(10000):
        angle = float(rng.uniform(0.0, math.pi))
        c, s = abs(math.cos(angle)), abs(math.sin(angle))
        margin_x = c * length / 2.0 + s * half_width + 1.0
        margin_y = s * length / 2.0 + c * half_width + 1.0
        if margin_x >= (size - 1) / 2 or margin_y >= (size - 1) / 2:
            continue
        center_x = float(rng.uniform(margin_x, size - 1 - margin_x))
        center_y = float(rng.uniform(margin_y, size - 1 - margin_y))
        return {"center_xy": [center_x, center_y], "angle_rad": angle, "length": float(length), "half_width": float(half_width)}
    raise RuntimeError("could not sample an in-bounds corridor")


def extract_corridor(x: np.ndarray, spec: dict) -> np.ndarray:
    center_x, center_y = map(float, spec["center_xy"])
    angle = float(spec["angle_rad"])
    length = float(spec["length"])
    half_width = float(spec["half_width"])
    width = max(16, int(math.ceil(length)))
    height = max(8, int(math.ceil(2 * half_width)) + 1)
    along = np.linspace(-length / 2.0, length / 2.0, width)
    across = np.linspace(-half_width, half_width, height)
    a, p = np.meshgrid(along, across)
    axis = np.asarray([math.cos(angle), math.sin(angle)])
    perp = np.asarray([-axis[1], axis[0]])
    xx = center_x + a * axis[0] + p * perp[0]
    yy = center_y + a * axis[1] + p * perp[1]
    crop = ndimage.map_coordinates(np.asarray(x, np.float32), [yy, xx], order=1, mode="reflect")
    side = max(crop.shape)
    pad_y, pad_x = side - crop.shape[0], side - crop.shape[1]
    crop = np.pad(crop, ((pad_y // 2, pad_y - pad_y // 2), (pad_x // 2, pad_x - pad_x // 2)), mode="reflect")
    if side != core.IMAGE_SIZE:
        crop = ndimage.zoom(crop, (core.IMAGE_SIZE / side, core.IMAGE_SIZE / side), order=1)
    return np.clip(crop, 0.0, 1.0).astype(np.float32)


def corridor_local_rank(x: np.ndarray, candidate: np.ndarray, genome: dict, file_sha: str, label: str, cfg: dict) -> dict:
    null_count = int(cfg["null_corridors_per_image"])
    length = float(cfg["length_pixels_normalized"])
    half_width = float(cfg["half_width_pixels_normalized"])
    vectors = []
    for index in range(null_count):
        spec = corridor_spec((cfg["null_seed_domain"], file_sha, label, index), length, half_width)
        vectors.append(core.geometry(extract_corridor(x, spec)))
    matrix = np.asarray(vectors, dtype=np.float64)
    center = np.median(matrix, axis=0)
    mad_scale = 1.4826 * np.median(np.abs(matrix - center), axis=0)
    std_scale = matrix.std(axis=0, ddof=1)
    scale = np.where(mad_scale > 1e-8, mad_scale, np.maximum(std_scale, 1e-8))
    candidate_score = core.weighted_std_dist(core.geometry(candidate), center, scale, genome)
    null_scores = np.asarray([core.weighted_std_dist(row, center, scale, genome) for row in matrix])
    ge_count = int(np.count_nonzero(null_scores >= candidate_score))
    null_q99 = float(np.quantile(null_scores, 0.99, method="higher"))
    return {
        "candidate_score": float(candidate_score),
        "null_count": null_count,
        "ge_count": ge_count,
        "p_empirical": float((ge_count + 1) / (null_count + 1)),
        "null_min": float(null_scores.min()),
        "null_median": float(np.median(null_scores)),
        "null_q99": null_q99,
        "null_max": float(null_scores.max()),
        "tail_ratio_q99": float(candidate_score / max(null_q99, 1e-12)),
        "passes_local_alpha": bool((ge_count + 1) / (null_count + 1) <= float(cfg["local_empirical_alpha"])),
    }


def rank_standardize(x: np.ndarray) -> np.ndarray:
    flat = np.asarray(x, np.float64).ravel()
    order = np.argsort(flat, kind="mergesort")
    ranks = np.empty(flat.size, dtype=np.float64)
    ranks[order] = np.linspace(-1.0, 1.0, flat.size)
    return ranks.reshape(np.asarray(x).shape)


def psf_matched_morphology(x: np.ndarray, cfg: dict) -> np.ndarray:
    ranked = rank_standardize(x)
    smooth = ndimage.gaussian_filter(ranked, float(cfg["gaussian_sigma_pixels"]), mode="reflect")
    size = int(cfg["comparison_size_pixels"])
    if smooth.shape != (size, size):
        smooth = ndimage.zoom(smooth, (size / smooth.shape[0], size / smooth.shape[1]), order=1)
    smooth -= smooth.mean()
    smooth /= max(float(smooth.std()), 1e-12)
    return smooth


def morphology_correlation(a: np.ndarray, b: np.ndarray, cfg: dict) -> float:
    return core.corr(psf_matched_morphology(a, cfg), psf_matched_morphology(b, cfg))


def orion_cross_survey_agreement(images: dict[tuple[str, str], np.ndarray], cfg: dict) -> dict:
    dss = [value for (family, _), value in images.items() if family == "DSS2"]
    tmass = [value for (family, _), value in images.items() if family == "2MASS"]
    if len(dss) != 2 or len(tmass) != 2:
        raise RuntimeError("Orion agreement requires two DSS2 and two 2MASS bands")
    cross = [morphology_correlation(a, b, cfg) for a in dss for b in tmass]
    return {"cross_family_correlations": cross, "score": float(np.median(cross))}


def _longest_true_run(values: np.ndarray) -> tuple[int, int]:
    best = (0, 0)
    start = None
    for index, value in enumerate(np.asarray(values, dtype=bool)):
        if value and start is None:
            start = index
        if start is not None and (not value or index == len(values) - 1):
            end = index + 1 if value and index == len(values) - 1 else index
            if end - start > best[1] - best[0]:
                best = (start, end)
            start = None
    return best


def normalize_masked(raw: np.ndarray, mask: np.ndarray, genome: dict) -> np.ndarray:
    data = np.asarray(raw, dtype=np.float32)
    valid = np.asarray(mask, dtype=bool) & np.isfinite(data)
    if np.count_nonzero(valid) < 64:
        raise RuntimeError("insufficient valid pixels")
    values = data[valid]
    lo = float(np.percentile(values, 1.0))
    hi = float(np.percentile(values, float(genome["clip_high"])))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        raise RuntimeError("invalid masked photometric range")
    fill = float(np.median(values))
    data = np.where(valid, data, fill)
    data = np.clip((data - lo) / (hi - lo), 0.0, 1.0)
    scale = float(genome["asinh_scale"])
    data = np.arcsinh(scale * data) / np.arcsinh(scale)
    if data.shape != (core.IMAGE_SIZE, core.IMAGE_SIZE):
        data = ndimage.zoom(data, (core.IMAGE_SIZE / data.shape[0], core.IMAGE_SIZE / data.shape[1]), order=1)
    return np.clip(data, 0.0, 1.0).astype(np.float32)


def common_valid_support_pair(
    raw_a: np.ndarray,
    raw_b: np.ndarray,
    weight_a: np.ndarray,
    weight_b: np.ndarray,
    genome: dict,
    cfg: dict,
) -> tuple[np.ndarray, np.ndarray, dict]:
    arrays = [np.asarray(item) for item in (raw_a, raw_b, weight_a, weight_b)]
    if len({item.shape for item in arrays}) != 1:
        raise RuntimeError("HST science and weight products have different shapes")
    mask = np.isfinite(arrays[0]) & np.isfinite(arrays[1]) & np.isfinite(arrays[2]) & np.isfinite(arrays[3])
    mask &= arrays[2] > 0
    mask &= arrays[3] > 0
    # WFPC2 chip products can have invalid borders. Locate the contiguous
    # observed chip support first; the stricter 98% coverage gate is applied
    # to the resulting square below.
    row_run = _longest_true_run(np.any(mask, axis=1))
    col_run = _longest_true_run(np.any(mask, axis=0))
    if row_run[1] <= row_run[0] or col_run[1] <= col_run[0]:
        raise RuntimeError("no common HST valid-support rectangle")
    y0, y1 = row_run
    x0, x1 = col_run
    side = min(y1 - y0, x1 - x0)
    cy, cx = (y0 + y1) // 2, (x0 + x1) // 2
    minimum_side = int(cfg["minimum_native_side_pixels"])
    minimum_fraction = float(cfg["minimum_valid_fraction"])
    while side >= minimum_side:
        yy0 = max(0, cy - side // 2)
        xx0 = max(0, cx - side // 2)
        yy1, xx1 = yy0 + side, xx0 + side
        if yy1 <= mask.shape[0] and xx1 <= mask.shape[1]:
            submask = mask[yy0:yy1, xx0:xx1]
            if float(submask.mean()) >= minimum_fraction:
                break
        side = int(side * 0.95)
    else:
        raise RuntimeError("HST common valid support did not reach the frozen coverage gate")
    submask = mask[yy0:yy1, xx0:xx1]
    image_a = normalize_masked(arrays[0][yy0:yy1, xx0:xx1], submask, genome)
    image_b = normalize_masked(arrays[1][yy0:yy1, xx0:xx1], submask, genome)
    return image_a, image_b, {
        "native_shape": list(mask.shape),
        "crop_yx_0based_half_open": [int(yy0), int(yy1), int(xx0), int(xx1)],
        "native_side": int(side),
        "common_valid_fraction": float(submask.mean()),
        "mask_gate_pass": True,
    }
