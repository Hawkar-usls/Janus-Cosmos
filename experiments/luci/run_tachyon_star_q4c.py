#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits

from experiments.luci.run_tachyon_star_q4b import (
    FRAME_SHA256,
    FWHM,
    download_exact,
    forced_gaussian_plane,
    sha256_file,
    temporal_contrast,
)

FRAME_ORDER = ("A_BEFORE", "B_CANDIDATE", "C_AFTER")
ROLE_TO_FILE = {
    "A_BEFORE": "luci1.20220319.0114.fits.gz",
    "B_CANDIDATE": "luci1.20220319.0116.fits.gz",
    "C_AFTER": "luci1.20220319.0120.fits.gz",
}
TARGET_X = 1437.216755721343
TARGET_Y = 2038.7227543099102
GRID_X = tuple(float(x) for x in range(32, 2017, 32) if abs(float(x) - TARGET_X) >= 128.0)
EXPECTED_CONTROLS = 55
RESIDUAL_Z_MIN = 3.0
RANK_MAX = 0.05
INJECTION_CONTROL_SIGMA = 5.0
Q4B_EXPECTED_TARGET_AMPLITUDE = {
    "A_BEFORE": 2392.0767796086225,
    "B_CANDIDATE": 2325.455455012585,
    "C_AFTER": 2243.341603642889,
}
Q4B_REPLAY_ABS_TOL = 1.0e-3
CLAIM = "DENSE_SAME_Y_LOCAL_RESIDUAL_ADJUDICATION_ONLY__POST_Q4B__NO_SPECIFIC_ARTIFACT_MECHANISM__NO_TACHYON_OR_FTL_OR_UAP_OR_ARTIFICIAL_ORIGIN_CLAIM"


def robust_center_scale(values: list[float]) -> tuple[float, float, str]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size != EXPECTED_CONTROLS:
        raise RuntimeError(f"expected {EXPECTED_CONTROLS} finite control contrasts, got {x.size}")
    center = float(np.median(x))
    mad = float(np.median(np.abs(x - center)))
    scale = 1.4826 * mad
    mode = "MAD"
    if not math.isfinite(scale) or scale <= 0:
        scale = float(np.std(x, ddof=1))
        mode = "STD_FALLBACK"
    if not math.isfinite(scale) or scale <= 0:
        raise RuntimeError("nonpositive dense-control scale")
    return center, scale, mode


def empirical_rank_fraction(target: float, controls: list[float]) -> dict:
    vals = [float(v) for v in controls if math.isfinite(float(v))]
    if len(vals) != EXPECTED_CONTROLS:
        raise RuntimeError(f"expected {EXPECTED_CONTROLS} finite controls, got {len(vals)}")
    ge = sum(v >= float(target) for v in vals)
    gt = sum(v > float(target) for v in vals)
    return {
        "controls_ge_target": ge,
        "rank_descending_1_is_largest": 1 + gt,
        "empirical_one_sided_rank_fraction": float((1 + ge) / (EXPECTED_CONTROLS + 1)),
        "minimum_possible_rank_fraction": 1.0 / (EXPECTED_CONTROLS + 1),
    }


def unique_gate(target_contrast: float, center: float, scale: float, controls: list[float]) -> dict:
    residual_z = float((float(target_contrast) - center) / scale)
    rank = empirical_rank_fraction(target_contrast, controls)
    passed = bool(
        float(target_contrast) > 0
        and residual_z >= RESIDUAL_Z_MIN
        and rank["empirical_one_sided_rank_fraction"] <= RANK_MAX
    )
    return {
        "passed": passed,
        "target_temporal_contrast": float(target_contrast),
        "control_center": center,
        "control_scale": scale,
        "target_residual_z": residual_z,
        "rank": rank,
        "thresholds": {
            "positive_temporal_contrast_required": True,
            "residual_z_min": RESIDUAL_Z_MIN,
            "rank_fraction_max": RANK_MAX,
        },
    }


def add_fixed_gaussian(image: np.ndarray, x0: float, y0: float, fwhm_px: float, peak_amp: float) -> np.ndarray:
    out = np.array(image, dtype=float, copy=True)
    h, w = out.shape
    cx, cy = int(round(float(x0))), int(round(float(y0)))
    half = 8
    x1, x2 = cx - half, cx + half + 1
    y1, y2 = cy - half, cy + half + 1
    if x1 < 0 or y1 < 0 or x2 > w or y2 > h:
        raise RuntimeError("injection patch out of bounds")
    yy, xx = np.indices((17, 17), dtype=float)
    xx += x1
    yy += y1
    sigma = float(fwhm_px) / 2.354820045
    g = np.exp(-0.5 * ((xx - float(x0)) ** 2 + (yy - float(y0)) ** 2) / (sigma * sigma))
    out[y1:y2, x1:x2] += float(peak_amp) * g
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", default="results/tachyon_star_q4c")
    ap.add_argument("--cache-dir", default=".cache/tachyon_star_q4c")
    args = ap.parse_args()
    out = Path(args.output_dir)
    cache = Path(args.cache_dir)
    out.mkdir(parents=True, exist_ok=True)

    if len(GRID_X) != EXPECTED_CONTROLS:
        raise RuntimeError(f"dense control grid changed: {len(GRID_X)}")
    if any(abs(x - TARGET_X) < 128.0 for x in GRID_X):
        raise RuntimeError("target-neighborhood exclusion violated")

    images = {}
    provenance = []
    for role in FRAME_ORDER:
        filename = ROLE_TO_FILE[role]
        path = download_exact(filename, cache)
        got = sha256_file(path)
        if got != FRAME_SHA256[filename]:
            raise RuntimeError(f"frame SHA mismatch: {filename}")
        image = np.asarray(fits.getdata(path, ext=0), dtype=float)
        if image.shape != (2048, 2048):
            raise RuntimeError(f"unexpected shape for {filename}: {image.shape}")
        images[role] = image
        provenance.append({"role": role, "filename": filename, "sha256": got, "shape": list(image.shape)})

    target = {}
    for role in FRAME_ORDER:
        filename = ROLE_TO_FILE[role]
        target[role] = forced_gaussian_plane(images[role], TARGET_X, TARGET_Y, FWHM[filename])
        if abs(float(target[role]["amplitude"]) - Q4B_EXPECTED_TARGET_AMPLITUDE[role]) > Q4B_REPLAY_ABS_TOL:
            raise RuntimeError(
                f"Q4B target amplitude replay mismatch {role}: {target[role]['amplitude']} vs {Q4B_EXPECTED_TARGET_AMPLITUDE[role]}"
            )
    target_contrast = temporal_contrast(target["A_BEFORE"], target["B_CANDIDATE"], target["C_AFTER"])

    controls = []
    for i, x in enumerate(GRID_X, start=1):
        fits_by_role = {}
        for role in FRAME_ORDER:
            filename = ROLE_TO_FILE[role]
            fits_by_role[role] = forced_gaussian_plane(images[role], x, TARGET_Y, FWHM[filename])
        contrast = temporal_contrast(fits_by_role["A_BEFORE"], fits_by_role["B_CANDIDATE"], fits_by_role["C_AFTER"])
        controls.append({
            "control_id": f"Q4C-DENSE-{i:02d}",
            "x_zero_based": x,
            "y_zero_based": TARGET_Y,
            "fits": fits_by_role,
            "temporal_contrast": contrast,
        })

    control_contrasts = [float(q["temporal_contrast"]) for q in controls]
    center, scale, scale_mode = robust_center_scale(control_contrasts)
    raw_gate = unique_gate(target_contrast, center, scale, control_contrasts)

    injected_peak = INJECTION_CONTROL_SIGMA * scale
    injected_b = add_fixed_gaussian(
        images["B_CANDIDATE"], TARGET_X, TARGET_Y, FWHM[ROLE_TO_FILE["B_CANDIDATE"]], injected_peak
    )
    injected_b_fit = forced_gaussian_plane(
        injected_b, TARGET_X, TARGET_Y, FWHM[ROLE_TO_FILE["B_CANDIDATE"]]
    )
    injected_contrast = temporal_contrast(target["A_BEFORE"], injected_b_fit, target["C_AFTER"])
    injection_gate = unique_gate(injected_contrast, center, scale, control_contrasts)

    if raw_gate["passed"]:
        status = "UNIQUE_LOCAL_B_RESIDUAL_CANDIDATE_REQUIRES_ADJUDICATION"
    elif injection_gate["passed"]:
        status = "QUALIFIED_NO_UNIQUE_LOCAL_B_RESIDUAL_AT_5SIGMA_CONTROL_SCALE"
    else:
        status = "UNRESOLVED_DENSE_CONTROL_SENSITIVITY"

    measurements = {
        "schema": "janus.cosmos.tachyon_star.q4c.measurements.v1",
        "experiment_id": "JANUS-TACHYON-STAR-Q4C-LUCI-DENSE-SAME-Y-EDGE-RESIDUAL",
        "scope": "POST_Q4B_ADJUDICATION__NOT_BLIND_DISCOVERY",
        "frame_provenance": provenance,
        "target": target,
        "target_temporal_contrast": target_contrast,
        "dense_controls": controls,
        "dense_control_count": len(controls),
        "dense_control_center": center,
        "dense_control_scale": scale,
        "dense_control_scale_mode": scale_mode,
        "raw_unique_gate": raw_gate,
        "sensitivity_injection": {
            "injected_peak_amplitude": injected_peak,
            "injected_control_sigma": INJECTION_CONTROL_SIGMA,
            "injected_b_fit": injected_b_fit,
            "injected_temporal_contrast": injected_contrast,
            "control_distribution_recomputed": False,
            "gate": injection_gate,
        },
    }
    mpath = out / "measurements.json"
    mpath.write_text(json.dumps(measurements, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    rec = {
        "schema": "janus.cosmos.tachyon_star.q4c.receipt.v1",
        "experiment_id": "JANUS-TACHYON-STAR-Q4C-LUCI-DENSE-SAME-Y-EDGE-RESIDUAL",
        "status": status,
        "scope": "POST_Q4B_ADJUDICATION__NOT_BLIND_DISCOVERY",
        "control_count": len(controls),
        "control_center": center,
        "control_scale": scale,
        "control_scale_mode": scale_mode,
        "target_temporal_contrast": target_contrast,
        "target_residual_z": raw_gate["target_residual_z"],
        "target_rank": raw_gate["rank"],
        "raw_unique_gate_passed": raw_gate["passed"],
        "sensitivity_injected_control_sigma": INJECTION_CONTROL_SIGMA,
        "sensitivity_injected_peak_amplitude": injected_peak,
        "sensitivity_injection_gate_passed": injection_gate["passed"],
        "sensitivity_injection_residual_z": injection_gate["target_residual_z"],
        "sensitivity_injection_rank": injection_gate["rank"],
        "measurements_sha256": sha256_file(mpath),
        "interpretation_limit": "A qualified no-unique-residual result strengthens a broad detector/background explanation at the preregistered residual scale but does not exclude a smaller one-frame sky transient or identify a specific detector mechanism.",
        "claim_ceiling": CLAIM,
    }
    (out / "receipt.json").write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(rec, indent=2, sort_keys=True))
    return 0 if status != "UNRESOLVED_DENSE_CONTROL_SENSITIVITY" else 3


if __name__ == "__main__":
    raise SystemExit(main())
