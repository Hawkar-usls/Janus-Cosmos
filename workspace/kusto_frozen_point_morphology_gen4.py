#!/usr/bin/env python3
"""JANUS Cousteau Gen4 frozen-point morphology gate.

This module deliberately separates:
  historical frozen coordinate -> direct-data provenance -> resolving morphology.
It never recenters on a visually interesting nearby feature and never treats
unmasked/global interpolated bathymetry as direct sonar evidence.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable

import numpy as np

TARGET_LAT = -3.865418
TARGET_LON = 3.854924
CONTRACT_PATH = Path("data/cousteau/JANUS-KUSTO-FROZEN-POINT-MORPHOLOGY-GEN4-CONTRACT-v1.0.json")
GMRT_URL = "https://www.gmrt.org/services/GridServer"
GEBCO_WMS_URL = "https://wms.gebco.net/mapserv"

SCALE_BANDS = [
    {"id": "BASE_0P5_1KM", "min_m": 500.0, "max_m": 1000.0, "width_m": 750.0},
    {"id": "BASE_1_2KM", "min_m": 1000.0, "max_m": 2000.0, "width_m": 1500.0},
    {"id": "BASE_2_4KM", "min_m": 2000.0, "max_m": 4000.0, "width_m": 3000.0},
    {"id": "BASE_4_8KM", "min_m": 4000.0, "max_m": 8000.0, "width_m": 6000.0},
]
MIN_SAMPLES_ACROSS_BASE = 8.0
MIN_FINITE_FRACTION = 0.80
MIN_VALID_CONTROLS = 20
EMPIRICAL_ALPHA = 0.05
PROM_SNR_MIN = 3.0
RADIAL_MONO_MIN = 0.66
FOURFOLD_MIN = 0.20
ORIENT_MAX_DEG = 20.0


def canonical_json(obj: Any) -> str:
    return json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_obj(obj: Any) -> str:
    return sha256_bytes(canonical_json(obj).encode("utf-8"))


def meters_per_degree(lat_deg: float) -> tuple[float, float]:
    lat = math.radians(lat_deg)
    mlat = 111132.92 - 559.82 * math.cos(2 * lat) + 1.175 * math.cos(4 * lat)
    mlon = 111412.84 * math.cos(lat) - 93.5 * math.cos(3 * lat)
    return mlat, mlon


def offset_latlon(lat: float, lon: float, distance_km: float, azimuth_deg: float) -> tuple[float, float]:
    # Local tangent approximation is more than adequate over <= 40 km here.
    mlat, mlon = meters_per_degree(lat)
    a = math.radians(azimuth_deg)
    north_m = distance_km * 1000.0 * math.cos(a)
    east_m = distance_km * 1000.0 * math.sin(a)
    return lat + north_m / mlat, lon + east_m / mlon


def deterministic_controls() -> list[dict[str, float]]:
    controls: list[dict[str, float]] = []
    for ring_km in (10.0, 20.0, 30.0, 40.0):
        for i in range(16):
            az = i * 22.5
            lat, lon = offset_latlon(TARGET_LAT, TARGET_LON, ring_km, az)
            controls.append({"lat": lat, "lon": lon, "ring_km": ring_km, "azimuth_deg": az})
    return controls


def bbox_for_halfwidth_km(lat: float, lon: float, halfwidth_km: float) -> tuple[float, float, float, float]:
    mlat, mlon = meters_per_degree(lat)
    dy = halfwidth_km * 1000.0 / mlat
    dx = halfwidth_km * 1000.0 / mlon
    return lon - dx, lon + dx, lat - dy, lat + dy


def _pick_grid(ds: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Return latitude, longitude, 2D values, source variable name."""
    coord_names = {name.lower(): name for name in ds.coords}
    lat_name = next((coord_names[k] for k in ("lat", "latitude", "y") if k in coord_names), None)
    lon_name = next((coord_names[k] for k in ("lon", "longitude", "x") if k in coord_names), None)

    if lat_name is None or lon_name is None:
        for name in ds.variables:
            low = name.lower()
            if lat_name is None and low in ("lat", "latitude", "y"):
                lat_name = name
            if lon_name is None and low in ("lon", "longitude", "x"):
                lon_name = name
    if lat_name is None or lon_name is None:
        raise RuntimeError(f"Could not identify lat/lon coordinates: {list(ds.variables)}")

    lat = np.asarray(ds[lat_name].values, dtype=float).squeeze()
    lon = np.asarray(ds[lon_name].values, dtype=float).squeeze()
    if lat.ndim != 1 or lon.ndim != 1:
        raise RuntimeError("Expected 1D latitude/longitude coordinates")

    candidate = None
    for name, var in ds.data_vars.items():
        arr = np.asarray(var.values)
        if arr.ndim == 2 and set(var.dims) >= {lat_name, lon_name}:
            candidate = (name, var)
            break
    if candidate is None:
        # Some servers expose dimension names different from coordinate variable names.
        for name, var in ds.data_vars.items():
            arr = np.asarray(var.values)
            if arr.ndim == 2 and sorted(arr.shape) == sorted((lat.size, lon.size)):
                candidate = (name, var)
                break
    if candidate is None:
        raise RuntimeError(f"Could not identify a 2D grid variable: {list(ds.data_vars)}")

    name, var = candidate
    z = np.asarray(var.values, dtype=float)
    dims = list(var.dims)
    if lat_name in dims and lon_name in dims:
        if dims.index(lat_name) > dims.index(lon_name):
            z = z.T
    elif z.shape == (lon.size, lat.size):
        z = z.T
    if z.shape != (lat.size, lon.size):
        raise RuntimeError(f"Unexpected grid shape {z.shape} vs {(lat.size, lon.size)}")

    if lat[0] > lat[-1]:
        lat = lat[::-1]
        z = z[::-1, :]
    if lon[0] > lon[-1]:
        lon = lon[::-1]
        z = z[:, ::-1]
    return lat, lon, z, name


def fetch_gmrt_mask(halfwidth_km: float = 46.0) -> dict[str, Any]:
    import requests
    import xarray as xr

    west, east, south, north = bbox_for_halfwidth_km(TARGET_LAT, TARGET_LON, halfwidth_km)
    params = {
        "west": f"{west:.8f}",
        "east": f"{east:.8f}",
        "south": f"{south:.8f}",
        "north": f"{north:.8f}",
        "layer": "topo-mask",
        "format": "coards",
        "resolution": "max",
    }
    r = requests.get(GMRT_URL, params=params, timeout=180)
    r.raise_for_status()
    payload_sha = sha256_bytes(r.content)
    with xr.open_dataset(io.BytesIO(r.content), engine="netcdf4") as ds:
        lat, lon, z, variable = _pick_grid(ds)
        attrs = {k: str(v) for k, v in ds.attrs.items()}

    mlat, mlon = meters_per_degree(TARGET_LAT)
    dy = abs(float(np.median(np.diff(lat)))) * mlat if lat.size > 1 else float("inf")
    dx = abs(float(np.median(np.diff(lon)))) * mlon if lon.size > 1 else float("inf")
    resolution_m = max(dx, dy)
    return {
        "lat": lat,
        "lon": lon,
        "z": z,
        "resolution_m": resolution_m,
        "grid_dx_m": dx,
        "grid_dy_m": dy,
        "variable": variable,
        "request_url": r.url,
        "http_status": r.status_code,
        "payload_sha256": payload_sha,
        "attrs": attrs,
        "finite_fraction_entire_request": float(np.isfinite(z).mean()),
    }


def fetch_gebco_measured_mask() -> dict[str, Any]:
    """Query the official GEBCO measured-data WMS layer around the frozen point.

    This intentionally returns only a measured-mask indication, never an exact TID code.
    """
    import requests
    from PIL import Image

    delta = 0.0125  # 3 x 15-arcsec cells each side; fixed before inspection.
    params = {
        "service": "WMS",
        "version": "1.3.0",
        "request": "GetMap",
        "layers": "gebco_latest_tid",
        "crs": "EPSG:4326",
        # WMS 1.3.0 axis order for EPSG:4326 is lat,lon.
        "bbox": f"{TARGET_LAT-delta},{TARGET_LON-delta},{TARGET_LAT+delta},{TARGET_LON+delta}",
        "width": "96",
        "height": "96",
        "format": "image/png",
        "transparent": "TRUE",
    }
    r = requests.get(GEBCO_WMS_URL, params=params, timeout=60)
    r.raise_for_status()
    im = Image.open(io.BytesIO(r.content)).convert("RGBA")
    arr = np.asarray(im)
    h, w, _ = arr.shape
    center = arr[h//2-4:h//2+5, w//2-4:w//2+5, :]
    rgb = center[..., :3].astype(float)
    alpha = center[..., 3].astype(float)
    dark = (rgb.mean(axis=2) < 80.0) & (alpha > 0)
    opaque = alpha > 0
    return {
        "request_url": r.url,
        "http_status": r.status_code,
        "payload_sha256": sha256_bytes(r.content),
        "central_dark_fraction": float(dark.mean()),
        "central_opaque_fraction": float(opaque.mean()),
        "measured_mask_present": bool(dark.mean() >= 0.50),
        "interpretation": "MEASURED_MASK_ONLY__NOT_EXACT_TID_CODE",
    }


def grid_local_xy(lat: np.ndarray, lon: np.ndarray, center_lat: float, center_lon: float) -> tuple[np.ndarray, np.ndarray]:
    mlat, mlon = meters_per_degree(center_lat)
    yy = (lat[:, None] - center_lat) * mlat
    xx = (lon[None, :] - center_lon) * mlon
    return xx, yy


def _robust_noise(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan")
    med = float(np.median(values))
    mad = float(np.median(np.abs(values - med)))
    sigma = 1.4826 * mad
    if not math.isfinite(sigma) or sigma < 0.5:
        sigma = max(float(np.std(values)), 0.5)
    return sigma


def _harmonic(values: np.ndarray, theta: np.ndarray, k: int) -> complex:
    good = np.isfinite(values) & np.isfinite(theta)
    if good.sum() < 16:
        return complex(float("nan"), float("nan"))
    v = values[good] - np.median(values[good])
    return complex(np.sum(v * np.exp(-1j * k * theta[good])))


def _orientation_diff_deg(a: float, b: float) -> float:
    # k=4 orientation is periodic every 90 degrees.
    period = math.pi / 2.0
    d = abs((a - b + period / 2.0) % period - period / 2.0)
    return math.degrees(d)


def evaluate_at(lat: np.ndarray, lon: np.ndarray, z: np.ndarray, center_lat: float, center_lon: float, width_m: float) -> dict[str, Any]:
    xx, yy = grid_local_xy(lat, lon, center_lat, center_lon)
    rr = np.hypot(xx, yy)
    theta = np.arctan2(yy, xx)
    test_radius = 0.65 * width_m
    inside = rr <= test_radius
    finite_inside = inside & np.isfinite(z)
    finite_fraction = float(finite_inside.sum() / max(1, inside.sum()))

    iy = int(np.argmin(np.abs(lat - center_lat)))
    ix = int(np.argmin(np.abs(lon - center_lon)))
    center_finite = bool(np.isfinite(z[iy, ix]))
    if not center_finite or finite_fraction < MIN_FINITE_FRACTION:
        return {"admissible": False, "reason": "NO_DIRECT_COVERAGE", "finite_fraction": finite_fraction, "center_finite": center_finite}

    inner = z[(rr <= 0.12 * width_m) & np.isfinite(z)]
    outer = z[(rr >= 0.48 * width_m) & (rr <= 0.65 * width_m) & np.isfinite(z)]
    if inner.size < 4 or outer.size < 12:
        return {"admissible": False, "reason": "INSUFFICIENT_CELLS", "finite_fraction": finite_fraction, "center_finite": center_finite}

    prominence = float(np.median(inner) - np.median(outer))
    noise = _robust_noise(outer)
    prom_snr = float(prominence / noise) if math.isfinite(noise) and noise > 0 else float("nan")

    edges = np.array([0.00, 0.15, 0.30, 0.45, 0.60]) * width_m
    radial_medians: list[float] = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        vals = z[(rr >= lo) & (rr < hi) & np.isfinite(z)]
        radial_medians.append(float(np.median(vals)) if vals.size else float("nan"))
    transitions = []
    for a, b in zip(radial_medians[:-1], radial_medians[1:]):
        if math.isfinite(a) and math.isfinite(b):
            transitions.append(float(b <= a))
    radial_monotonicity = float(np.mean(transitions)) if transitions else float("nan")

    ring = (rr >= 0.20 * width_m) & (rr <= 0.55 * width_m) & np.isfinite(z)
    harmonics = [_harmonic(z[ring], theta[ring], k) for k in range(1, 9)]
    powers = np.array([abs(c) ** 2 if math.isfinite(c.real) and math.isfinite(c.imag) else 0.0 for c in harmonics], dtype=float)
    fourfold_ratio = float(powers[3] / powers.sum()) if powers.sum() > 0 else 0.0

    inner_ring = (rr >= 0.20 * width_m) & (rr <= 0.35 * width_m) & np.isfinite(z)
    outer_ring = (rr > 0.35 * width_m) & (rr <= 0.55 * width_m) & np.isfinite(z)
    c4_inner = _harmonic(z[inner_ring], theta[inner_ring], 4)
    c4_outer = _harmonic(z[outer_ring], theta[outer_ring], 4)
    if all(math.isfinite(x) for x in (c4_inner.real, c4_inner.imag, c4_outer.real, c4_outer.imag)) and abs(c4_inner) > 0 and abs(c4_outer) > 0:
        ori_inner = math.atan2(c4_inner.imag, c4_inner.real) / 4.0
        ori_outer = math.atan2(c4_outer.imag, c4_outer.real) / 4.0
        orientation_diff = _orientation_diff_deg(ori_inner, ori_outer)
    else:
        orientation_diff = 90.0

    prominence_component = 0.0 if not math.isfinite(prom_snr) else float(1.0 / (1.0 + math.exp(-(prom_snr - 2.0))))
    orientation_component = max(0.0, 1.0 - orientation_diff / 30.0)
    rm = radial_monotonicity if math.isfinite(radial_monotonicity) else 0.0
    composite = 0.40 * prominence_component + 0.25 * rm + 0.20 * fourfold_ratio + 0.15 * orientation_component

    hard_floors = {
        "prominence_snr": bool(math.isfinite(prom_snr) and prom_snr >= PROM_SNR_MIN),
        "radial_monotonicity": bool(math.isfinite(radial_monotonicity) and radial_monotonicity >= RADIAL_MONO_MIN),
        "fourfold_power_ratio": bool(fourfold_ratio >= FOURFOLD_MIN),
        "orientation_stability": bool(orientation_diff <= ORIENT_MAX_DEG),
    }
    return {
        "admissible": True,
        "finite_fraction": finite_fraction,
        "center_finite": center_finite,
        "prominence_m": prominence,
        "outer_noise_robust_m": noise,
        "prominence_snr": prom_snr,
        "radial_medians_m": radial_medians,
        "radial_monotonicity": radial_monotonicity,
        "fourfold_power_ratio": fourfold_ratio,
        "orientation_stability_deg": orientation_diff,
        "composite_score": float(composite),
        "hard_metric_floors": hard_floors,
        "all_hard_floors_pass": bool(all(hard_floors.values())),
    }


def evaluate_scale(g: dict[str, Any], band: dict[str, Any]) -> dict[str, Any]:
    resolution_m = float(g["resolution_m"])
    required_max_resolution = band["min_m"] / MIN_SAMPLES_ACROSS_BASE
    if resolution_m > required_max_resolution:
        return {
            "scale_id": band["id"],
            "resolvable": False,
            "reason": "RESOLUTION_INSUFFICIENT",
            "grid_resolution_m": resolution_m,
            "required_max_resolution_m": required_max_resolution,
        }

    target = evaluate_at(g["lat"], g["lon"], g["z"], TARGET_LAT, TARGET_LON, band["width_m"])
    if not target.get("admissible"):
        return {
            "scale_id": band["id"],
            "resolvable": True,
            "target": target,
            "reason": target.get("reason", "NO_DIRECT_COVERAGE"),
        }

    control_scores: list[float] = []
    controls_detail: list[dict[str, Any]] = []
    for c in deterministic_controls():
        ev = evaluate_at(g["lat"], g["lon"], g["z"], c["lat"], c["lon"], band["width_m"])
        item = {"ring_km": c["ring_km"], "azimuth_deg": c["azimuth_deg"], "admissible": bool(ev.get("admissible"))}
        if ev.get("admissible"):
            score = float(ev["composite_score"])
            control_scores.append(score)
            item["composite_score"] = score
        controls_detail.append(item)

    n = len(control_scores)
    if n < MIN_VALID_CONTROLS:
        return {
            "scale_id": band["id"],
            "resolvable": True,
            "target": target,
            "valid_controls": n,
            "reason": "CONTROL_POWER_INSUFFICIENT",
            "controls": controls_detail,
        }

    target_score = float(target["composite_score"])
    ge = sum(1 for s in control_scores if s >= target_score)
    empirical_p = (1.0 + ge) / (1.0 + n)
    signal = bool(empirical_p <= EMPIRICAL_ALPHA and target["all_hard_floors_pass"])
    return {
        "scale_id": band["id"],
        "resolvable": True,
        "target": target,
        "valid_controls": n,
        "empirical_p": empirical_p,
        "target_rank_descending": 1 + ge,
        "signal_at_preregistered_tail": signal,
        "controls": controls_detail,
        "reason": "MORPHOLOGY_SIGNAL_REQUIRES_INDEPENDENT_REPLICATION" if signal else "NO_PREREGISTERED_MORPHOLOGY_SIGNAL",
    }


def run_live() -> dict[str, Any]:
    contract_bytes = CONTRACT_PATH.read_bytes() if CONTRACT_PATH.exists() else b""
    receipt: dict[str, Any] = {
        "schema": "janus.cousteau.kusto_frozen_point_morphology_receipt.v1",
        "artifact_id": "JANUS-KUSTO-FROZEN-POINT-MORPHOLOGY-GEN4-LIVE-RECEIPT",
        "generation": 4,
        "frozen_target": {"latitude_deg": TARGET_LAT, "longitude_deg": TARGET_LON, "no_recentering": True},
        "contract_path": str(CONTRACT_PATH),
        "contract_sha256": sha256_bytes(contract_bytes) if contract_bytes else None,
        "simulation_count_increased": False,
        "claim_ceiling": "MORPHOLOGY_ONLY__NOT_PYRAMID_IDENTITY__NOT_ANOMALY_PROOF",
        "underwater_pyramid_detected": False,
    }

    try:
        tid = fetch_gebco_measured_mask()
    except Exception as exc:  # network/source failure must remain explicit
        tid = {"status": "DATA_SOURCE_UNAVAILABLE", "error": f"{type(exc).__name__}: {exc}", "measured_mask_present": None}
    receipt["gebco_2026_tid_wms"] = tid

    try:
        g = fetch_gmrt_mask()
    except Exception as exc:
        receipt["gmrt_topo_mask"] = {"status": "DATA_SOURCE_UNAVAILABLE", "error": f"{type(exc).__name__}: {exc}"}
        receipt["decision"] = "DATA_SOURCE_UNAVAILABLE"
        receipt["result_sha256"] = sha256_obj({k: v for k, v in receipt.items() if k != "result_sha256"})
        return receipt

    receipt["gmrt_topo_mask"] = {
        "http_status": g["http_status"],
        "request_url": g["request_url"],
        "payload_sha256": g["payload_sha256"],
        "variable": g["variable"],
        "grid_resolution_m": g["resolution_m"],
        "grid_dx_m": g["grid_dx_m"],
        "grid_dy_m": g["grid_dy_m"],
        "finite_fraction_entire_request": g["finite_fraction_entire_request"],
        "attrs": g["attrs"],
    }

    scales = [evaluate_scale(g, b) for b in SCALE_BANDS]
    receipt["scale_results"] = scales

    reasons = [s.get("reason") for s in scales]
    if any(r == "MORPHOLOGY_SIGNAL_REQUIRES_INDEPENDENT_REPLICATION" for r in reasons):
        decision = "MORPHOLOGY_SIGNAL_REQUIRES_INDEPENDENT_REPLICATION"
    elif all(r == "RESOLUTION_INSUFFICIENT" for r in reasons):
        decision = "RESOLUTION_INSUFFICIENT"
    elif any(r == "NO_PREREGISTERED_MORPHOLOGY_SIGNAL" for r in reasons):
        decision = "NO_PREREGISTERED_MORPHOLOGY_SIGNAL"
    elif any(r == "CONTROL_POWER_INSUFFICIENT" for r in reasons):
        decision = "CONTROL_POWER_INSUFFICIENT"
    elif any(r == "NO_DIRECT_COVERAGE" for r in reasons):
        decision = "NO_DIRECT_COVERAGE"
    else:
        decision = "RESOLUTION_INSUFFICIENT"

    receipt["decision"] = decision
    receipt["interpretation"] = {
        "positive_signal_means": "target morphology is unusual versus preregistered local blind controls at at least one resolvable scale; independent resolving replication still required",
        "negative_signal_means": "no preregistered morphology signal in admitted direct GMRT data at resolvable scales",
        "insufficient_means": "data coverage, resolution, or control power does not permit the morphology question to be answered",
        "never_means": "proof of a pyramid, artificial structure, anomaly, or LOVE/EDEM identity",
    }
    receipt["result_sha256"] = sha256_obj({k: v for k, v in receipt.items() if k != "result_sha256"})
    return receipt


def self_test() -> dict[str, Any]:
    # 100 m grid, deterministic square-pyramid positive control and smooth negative control.
    axis = np.arange(-5000.0, 5000.1, 100.0)
    lat0, lon0 = TARGET_LAT, TARGET_LON
    mlat, mlon = meters_per_degree(lat0)
    lat = lat0 + axis / mlat
    lon = lon0 + axis / mlon
    xx, yy = np.meshgrid(axis, axis)
    width = 3000.0
    pyramid = 150.0 * np.clip(1.0 - np.maximum(np.abs(xx), np.abs(yy)) / (width / 2.0), 0.0, 1.0)
    ripple = 1.5 * np.sin(xx / 700.0) + 1.0 * np.cos(yy / 900.0)
    positive = pyramid + ripple
    negative = 3.0 * np.sin(xx / 1200.0) + 2.0 * np.cos(yy / 1700.0)
    p = evaluate_at(lat, lon, positive, lat0, lon0, width)
    n = evaluate_at(lat, lon, negative, lat0, lon0, width)
    if not p.get("admissible") or not n.get("admissible"):
        raise AssertionError("synthetic grids unexpectedly inadmissible")
    if not (p["composite_score"] > n["composite_score"] and p["prominence_snr"] > n["prominence_snr"]):
        raise AssertionError(f"positive control did not dominate negative control: p={p}, n={n}")
    if deterministic_controls() != deterministic_controls():
        raise AssertionError("controls are not deterministic")
    return {"status": "PASS", "positive_score": p["composite_score"], "negative_score": n["composite_score"], "positive_prominence_snr": p["prominence_snr"], "negative_prominence_snr": n["prominence_snr"], "controls": len(deterministic_controls())}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--live", action="store_true")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    if args.self_test:
        result = self_test()
        print("JANUS_GEN4_SELF_TEST=" + canonical_json(result))
    if args.live:
        result = run_live()
        text = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
        if args.out:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(text, encoding="utf-8")
        print("JANUS_GEN4_DECISION=" + str(result.get("decision")))
        print("JANUS_GEN4_RESULT=" + canonical_json(result))
    if not args.self_test and not args.live:
        ap.error("choose --self-test and/or --live")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
