#!/usr/bin/env python3
"""TOPA label-blind bathymetry benchmark.

Candidate generator only. It never calls a vent active, artificial, or anomalous from
shape alone. Labels/ground truth are intentionally separate from feature extraction.

Supported inputs:
  * ESRI ASCII grid (.asc) without extra dependencies
  * NumPy .npy arrays
  * NetCDF when xarray or netCDF4 is installed

The executable synthetic self-test validates the mechanics and deliberately includes
fault/scarp/pit confounders that should remain morphology candidates rather than be
silently discarded.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy import ndimage

EPS = 1e-12


def robust_z(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    med = np.nanmedian(v)
    mad = np.nanmedian(np.abs(v - med))
    if not np.isfinite(mad) or mad < EPS:
        sd = np.nanstd(v)
        return (v - med) / (sd + EPS)
    return (v - med) / (1.4826 * mad + EPS)


def read_esri_ascii(path: Path) -> Tuple[np.ndarray, Dict[str, float]]:
    meta: Dict[str, float] = {}
    with path.open("r", encoding="utf-8") as f:
        header = [f.readline().strip().split() for _ in range(6)]
        for k, val in header:
            meta[k.lower()] = float(val)
        arr = np.loadtxt(f)
    nodata = meta.get("nodata_value")
    if nodata is not None:
        arr = arr.astype(float)
        arr[np.isclose(arr, nodata)] = np.nan
    return arr, meta


def load_grid(path: Path, variable: str | None = None) -> Tuple[np.ndarray, Dict[str, object]]:
    suffix = path.suffix.lower()
    if suffix == ".npy":
        return np.asarray(np.load(path), dtype=float), {"format": "npy"}
    if suffix in {".asc", ".ascii", ".txt"}:
        arr, meta = read_esri_ascii(path)
        meta["format"] = "esri_ascii"
        return arr, meta
    if suffix in {".nc", ".grd", ".netcdf"}:
        try:
            import xarray as xr  # type: ignore
            ds = xr.open_dataset(path)
            if variable is None:
                candidates = [k for k, da in ds.data_vars.items() if da.ndim >= 2]
                if not candidates:
                    raise ValueError("no 2-D data variable found")
                variable = candidates[0]
            da = ds[variable].squeeze()
            if da.ndim != 2:
                raise ValueError(f"variable {variable!r} is not 2-D after squeeze")
            return np.asarray(da.values, dtype=float), {
                "format": "netcdf",
                "variable": variable,
                "shape": list(da.shape),
            }
        except ImportError as exc:
            raise RuntimeError("NetCDF input requires xarray (or convert to ESRI ASCII/.npy)") from exc
    raise ValueError(f"unsupported input suffix: {suffix}")


def finite_fill(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=float)
    if np.isfinite(a).all():
        return a.copy()
    out = a.copy()
    med = np.nanmedian(out)
    out[~np.isfinite(out)] = med
    return out


def orientation_anisotropy(tile: np.ndarray) -> float:
    gy, gx = np.gradient(tile)
    g = np.column_stack((gx.ravel(), gy.ravel()))
    good = np.isfinite(g).all(axis=1)
    if good.sum() < 8:
        return 0.0
    cov = np.cov(g[good], rowvar=False)
    ev = np.sort(np.linalg.eigvalsh(cov))
    return float((ev[-1] - ev[0]) / (ev[-1] + ev[0] + EPS))


def pinnacle_density(tile: np.ndarray) -> float:
    smooth = ndimage.gaussian_filter(tile, 1.0)
    local_max = smooth == ndimage.maximum_filter(smooth, size=5, mode="nearest")
    hp = smooth - ndimage.gaussian_filter(smooth, 5.0)
    rz = robust_z(hp)
    peaks = local_max & (rz > 2.5)
    return float(peaks.mean())


def multiscale_persistence(tile: np.ndarray) -> float:
    masks = []
    for sigma in (1.0, 2.0, 4.0):
        sm = ndimage.gaussian_filter(tile, sigma)
        hp = sm - ndimage.gaussian_filter(sm, sigma * 4.0)
        masks.append(robust_z(hp) > 2.0)
    intersection = masks[0] & masks[1] & masks[2]
    union = masks[0] | masks[1] | masks[2]
    return float(intersection.sum() / (union.sum() + EPS))


def tile_metrics(tile: np.ndarray) -> Dict[str, float]:
    tile = finite_fill(tile)
    gy, gx = np.gradient(tile)
    slope = np.hypot(gx, gy)
    lap = ndimage.laplace(ndimage.gaussian_filter(tile, 1.0))
    broad = ndimage.gaussian_filter(tile, 8.0)
    hp = tile - broad
    return {
        "local_relief_p95_p05": float(np.percentile(tile, 95) - np.percentile(tile, 5)),
        "slope_p95": float(np.percentile(slope, 95)),
        "abs_curvature_p95": float(np.percentile(np.abs(lap), 95)),
        "rugosity_residual_std": float(np.std(hp)),
        "pinnacle_density": pinnacle_density(tile),
        "multiscale_positive_persistence": multiscale_persistence(tile),
        "orientation_anisotropy": orientation_anisotropy(tile),
    }


def tile_grid(a: np.ndarray, tile_size: int, stride: int) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    h, w = a.shape
    for y0 in range(0, h - tile_size + 1, stride):
        for x0 in range(0, w - tile_size + 1, stride):
            t = a[y0:y0 + tile_size, x0:x0 + tile_size]
            if np.isfinite(t).mean() < 0.8:
                continue
            m = tile_metrics(t)
            rows.append({
                "x0_px": x0,
                "y0_px": y0,
                "x_center_px": x0 + tile_size / 2.0,
                "y_center_px": y0 + tile_size / 2.0,
                "metrics": m,
            })
    return rows


def score_rows(rows: List[Dict[str, object]]) -> None:
    # Candidate score emphasizes multiscale positive relief but deliberately retains
    # rough/linear geological confounders instead of pretending morphology is specific.
    keys = [
        "local_relief_p95_p05",
        "slope_p95",
        "abs_curvature_p95",
        "rugosity_residual_std",
        "pinnacle_density",
        "multiscale_positive_persistence",
    ]
    weights = np.array([0.15, 0.10, 0.15, 0.15, 0.25, 0.20], dtype=float)
    mat = np.array([[float(r["metrics"][k]) for k in keys] for r in rows], dtype=float)
    z = np.column_stack([robust_z(mat[:, i]) for i in range(mat.shape[1])])
    score = z @ weights
    anis = np.array([float(r["metrics"]["orientation_anisotropy"]) for r in rows])
    # Anisotropy is reported as a confounder flag, not used as a positive bonus.
    for i, r in enumerate(rows):
        r["morphology_score"] = float(score[i])
        r["lineation_confounder_flag"] = bool(anis[i] >= 0.55)
        r["classification_ceiling"] = "MORPHOLOGICAL_CANDIDATE_ONLY"
    rows.sort(key=lambda r: float(r["morphology_score"]), reverse=True)
    for rank, r in enumerate(rows, 1):
        r["rank"] = rank


def gaussian_bump(shape: Tuple[int, int], cx: float, cy: float, amp: float, sigma: float) -> np.ndarray:
    yy, xx = np.indices(shape)
    return amp * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2.0 * sigma ** 2))


def synthetic_terrain(seed: int = 230303, n: int = 512) -> Tuple[np.ndarray, Dict[str, object]]:
    rng = np.random.default_rng(seed)
    base = ndimage.gaussian_filter(rng.normal(size=(n, n)), 7.0)
    base *= 5.0 / (np.std(base) + EPS)
    yy, xx = np.indices((n, n))
    base += 35.0 * np.exp(-((xx - n/2) ** 2 + (yy - n/2) ** 2) / (2 * (0.35*n) ** 2))

    # Positive chimney-like cluster.
    vent_center = (112.0, 112.0)
    for dx, dy, amp, sig in [(-12,-8,20,3),(0,0,28,2.5),(9,6,23,3),(16,-11,18,2.5),(-4,15,16,3)]:
        base += gaussian_bump(base.shape, vent_center[0]+dx, vent_center[1]+dy, amp, sig)

    # Strong natural fault/graben confounder.
    base[:, 285:] += 11.0
    base[:, 310:] -= 18.0
    base[:, 340:] += 14.0
    for x in (278, 286, 304, 312, 333, 341):
        base += 7.0 * np.tanh((xx - x) / 1.5)

    # Pit-field confounder.
    pit_center = (120.0, 370.0)
    for _ in range(18):
        cx = pit_center[0] + rng.normal(0, 35)
        cy = pit_center[1] + rng.normal(0, 28)
        base -= gaussian_bump(base.shape, cx, cy, rng.uniform(9, 20), rng.uniform(3, 7))

    truth = {
        "seed": seed,
        "positive_chimney_like_center_px": list(vent_center),
        "fault_graben_confounder_x_range_px": [275, 345],
        "pit_field_center_px": list(pit_center),
        "note": "Synthetic truth is engine QA only; it is not FKt230303 field-data performance."
    }
    return base, truth


def nearest_rank(rows: List[Dict[str, object]], x: float, y: float) -> Dict[str, object]:
    best = min(rows, key=lambda r: (float(r["x_center_px"])-x)**2 + (float(r["y_center_px"])-y)**2)
    return {
        "rank": best["rank"],
        "x_center_px": best["x_center_px"],
        "y_center_px": best["y_center_px"],
        "morphology_score": best["morphology_score"],
        "lineation_confounder_flag": best["lineation_confounder_flag"],
    }


def run_benchmark(a: np.ndarray, tile_size: int, stride: int, top_k: int) -> Dict[str, object]:
    rows = tile_grid(a, tile_size, stride)
    if not rows:
        raise RuntimeError("no valid tiles")
    score_rows(rows)
    return {
        "schema": "JANUS_TOPA_BATHYMETRY_CANDIDATE_OUTPUT",
        "shape": list(a.shape),
        "tile_size_px": tile_size,
        "stride_px": stride,
        "metric_rule": "FROZEN_GENERIC_MORPHOLOGY_V1",
        "classification_ceiling": "MORPHOLOGICAL_CANDIDATE_ONLY",
        "candidate_count": len(rows),
        "top_candidates": rows[:top_k],
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("input", nargs="?", type=Path)
    p.add_argument("--variable")
    p.add_argument("--tile-size", type=int, default=64)
    p.add_argument("--stride", type=int, default=32)
    p.add_argument("--top-k", type=int, default=20)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--synthetic-selftest", action="store_true")
    args = p.parse_args()

    if args.synthetic_selftest:
        arr, truth = synthetic_terrain()
        result = run_benchmark(arr, args.tile_size, args.stride, args.top_k)
        result["mode"] = "SYNTHETIC_ENGINE_QA"
        result["truth_after_scoring"] = truth
        result["posthoc_diagnostics"] = {
            "positive_nearest_tile": nearest_rank(result["top_candidates"] if len(result["top_candidates"]) == result["candidate_count"] else tile_grid_and_score(arr, args.tile_size, args.stride), 112, 112)
        }
    else:
        if args.input is None:
            p.error("input is required unless --synthetic-selftest is used")
        arr, meta = load_grid(args.input, args.variable)
        result = run_benchmark(arr, args.tile_size, args.stride, args.top_k)
        result["mode"] = "FIELD_GRID_CANDIDATE_GENERATION"
        result["input_meta"] = meta

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")


def tile_grid_and_score(a: np.ndarray, tile_size: int, stride: int) -> List[Dict[str, object]]:
    rows = tile_grid(a, tile_size, stride)
    score_rows(rows)
    return rows


if __name__ == "__main__":
    main()
