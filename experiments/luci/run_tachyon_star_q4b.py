#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits

BASE_URL = "https://archive.lbto.org/files/lbt/"
FRAME_ORDER = ("A_BEFORE", "B_CANDIDATE", "C_AFTER")
FRAME_SHA256 = {
    "luci1.20220319.0114.fits.gz": "217ecc39c5d7753a1647536c91006f04b7c4abc73e3de86a8e029bd0d8e40f4a",
    "luci1.20220319.0116.fits.gz": "251a1221646347e7415eeaa9fead6a40e5fda8c9b04a06a3bdb855924783637a",
    "luci1.20220319.0120.fits.gz": "9e98b0cc1e08b0b1e8e16fb10d9c58609411f6deeaeb2df86dced16dde6ceeb6",
}
FWHM = {
    "luci1.20220319.0114.fits.gz": 6.026402,
    "luci1.20220319.0116.fits.gz": 6.0313418764,
    "luci1.20220319.0120.fits.gz": 5.945551,
}
PATCH_HALF = 8
MIN_FINITE = 250
COND_MAX = 1.0e8
B_Z_MIN = 3.0
COMPAT_Z_MAX = 3.0
CLAIM = "THREE_FRAME_FIXED_COORDINATE_ORIGIN_ADJUDICATION_ONLY__NO_TACHYON_OR_FTL_OR_UAP_OR_ARTIFICIAL_ORIGIN_CLAIM"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download_exact(filename: str, cache: Path) -> Path:
    cache.mkdir(parents=True, exist_ok=True)
    out = cache / filename
    if out.exists() and sha256_file(out) == FRAME_SHA256[filename]:
        return out
    req = urllib.request.Request(BASE_URL + filename, headers={"User-Agent": "Janus-Cosmos-Q4B/1.0"})
    with urllib.request.urlopen(req, timeout=300) as r:
        out.write_bytes(r.read())
    got = sha256_file(out)
    if got != FRAME_SHA256[filename]:
        raise RuntimeError(f"exact frame SHA mismatch for {filename}: {got}")
    return out


def robust_sigma(values: np.ndarray) -> tuple[float, str]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 10:
        raise RuntimeError("too few residuals for robust sigma")
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    sig = 1.4826 * mad
    mode = "MAD"
    if not math.isfinite(sig) or sig <= 0:
        sig = float(np.std(x, ddof=1))
        mode = "STD_FALLBACK"
    if not math.isfinite(sig) or sig <= 0:
        raise RuntimeError("nonpositive residual sigma")
    return sig, mode


def forced_gaussian_plane(image: np.ndarray, x0: float, y0: float, fwhm_px: float) -> dict:
    a = np.asarray(image, dtype=float)
    if a.ndim != 2:
        raise RuntimeError(f"expected 2-D image, got {a.shape}")
    h, w = a.shape
    cx, cy = int(round(float(x0))), int(round(float(y0)))
    x1, x2 = cx - PATCH_HALF, cx + PATCH_HALF + 1
    y1, y2 = cy - PATCH_HALF, cy + PATCH_HALF + 1
    if x1 < 0 or y1 < 0 or x2 > w or y2 > h:
        raise RuntimeError(f"17x17 patch out of bounds at x={x0},y={y0}")
    patch = a[y1:y2, x1:x2]
    yy, xx = np.indices(patch.shape, dtype=float)
    xx += x1
    yy += y1
    finite = np.isfinite(patch)
    n = int(finite.sum())
    if n < MIN_FINITE:
        raise RuntimeError(f"insufficient finite patch pixels: {n}/{MIN_FINITE}")
    dx = (xx - float(x0))[finite]
    dy = (yy - float(y0))[finite]
    vals = patch[finite]
    sig_psf = float(fwhm_px) / 2.354820045
    if not math.isfinite(sig_psf) or sig_psf <= 0:
        raise RuntimeError("invalid fixed FWHM")
    g = np.exp(-0.5 * (dx * dx + dy * dy) / (sig_psf * sig_psf))
    X = np.column_stack([g, np.ones_like(g), dx, dy])
    cond = float(np.linalg.cond(X))
    if not math.isfinite(cond) or cond > COND_MAX:
        raise RuntimeError(f"fit condition number too high: {cond}")
    beta, _, _, _ = np.linalg.lstsq(X, vals, rcond=None)
    pred = X @ beta
    resid = vals - pred
    resid_sigma, sigma_mode = robust_sigma(resid)
    xtx_inv = np.linalg.pinv(X.T @ X)
    amp_sigma = float(resid_sigma * math.sqrt(max(0.0, float(xtx_inv[0, 0]))))
    if not math.isfinite(amp_sigma) or amp_sigma <= 0:
        raise RuntimeError("invalid amplitude uncertainty")
    amp = float(beta[0])
    z = float(amp / amp_sigma)

    X0 = X[:, 1:]
    beta0, _, _, _ = np.linalg.lstsq(X0, vals, rcond=None)
    resid0 = vals - X0 @ beta0
    sse = float(np.sum(resid * resid))
    sse0 = float(np.sum(resid0 * resid0))
    return {
        "x_zero_based": float(x0),
        "y_zero_based": float(y0),
        "patch_integer_center": [cx, cy],
        "patch_bounds_inclusive": [x1, y1, x2 - 1, y2 - 1],
        "finite_pixels": n,
        "fixed_fwhm_px": float(fwhm_px),
        "fit_condition_number": cond,
        "amplitude": amp,
        "amplitude_sigma": amp_sigma,
        "forced_z": z,
        "background_b0": float(beta[1]),
        "background_bx": float(beta[2]),
        "background_by": float(beta[3]),
        "residual_sigma": resid_sigma,
        "residual_sigma_mode": sigma_mode,
        "gaussian_plus_plane_sse": sse,
        "plane_only_sse": sse0,
        "delta_sse_plane_minus_gaussian": sse0 - sse,
    }


def compatibility_z(b: dict, x: dict) -> float:
    denom = math.hypot(float(b["amplitude_sigma"]), float(x["amplitude_sigma"]))
    if not math.isfinite(denom) or denom <= 0:
        raise RuntimeError("invalid compatibility denominator")
    return float((float(b["amplitude"]) - float(x["amplitude"])) / denom)


def temporal_contrast(a: dict, b: dict, c: dict) -> float:
    return float(float(b["amplitude"]) - 0.5 * (float(a["amplitude"]) + float(c["amplitude"])))


def empirical_rank_p(target: float, controls: list[float]) -> dict:
    vals = [float(v) for v in controls if math.isfinite(float(v))]
    if len(vals) != 13:
        raise RuntimeError(f"all 13 frozen controls required, got {len(vals)}")
    ge = sum(v >= float(target) for v in vals)
    return {
        "target": float(target),
        "control_count": len(vals),
        "controls_ge_target": ge,
        "one_sided_rank_p": float((1 + ge) / (len(vals) + 1)),
        "rank_descending_1_is_largest": 1 + sum(v > float(target) for v in vals),
        "minimum_possible_p": 1.0 / 14.0,
    }


def classify(sky: dict, detector: dict) -> dict:
    b = sky["B_CANDIDATE"]
    if float(b["amplitude"]) <= 0 or float(b["forced_z"]) < B_Z_MIN:
        return {
            "classification": "B_FIXED_EXCESS_NOT_CONFIRMED",
            "b_fixed_excess_confirmed": False,
            "sky_persistence_compatible": None,
            "detector_persistence_compatible": None,
        }
    sky_za = compatibility_z(b, sky["A_BEFORE"])
    sky_zc = compatibility_z(b, sky["C_AFTER"])
    det_za = compatibility_z(b, detector["A_BEFORE"])
    det_zc = compatibility_z(b, detector["C_AFTER"])
    sky_ok = abs(sky_za) <= COMPAT_Z_MAX and abs(sky_zc) <= COMPAT_Z_MAX
    det_ok = abs(det_za) <= COMPAT_Z_MAX and abs(det_zc) <= COMPAT_Z_MAX
    if sky_ok and not det_ok:
        cls = "EVIDENCE_FAVORS_SKY_FIXED_PERSISTENCE"
    elif det_ok and not sky_ok:
        cls = "EVIDENCE_FAVORS_DETECTOR_FIXED_PERSISTENCE"
    elif sky_ok and det_ok:
        cls = "TRACKS_NOT_DISCRIMINATED"
    else:
        cls = "ONE_FRAME_EXCESS_ORIGIN_INDETERMINATE"
    return {
        "classification": cls,
        "b_fixed_excess_confirmed": True,
        "sky_persistence_compatible": sky_ok,
        "detector_persistence_compatible": det_ok,
        "sky_constant_flux_difference_z": {"A": sky_za, "C": sky_zc},
        "detector_constant_flux_difference_z": {"A": det_za, "C": det_zc},
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry-manifest", required=True)
    ap.add_argument("--output-dir", default="results/tachyon_star_q4b")
    ap.add_argument("--cache-dir", default=".cache/tachyon_star_q4b")
    args = ap.parse_args()
    geometry_path = Path(args.geometry_manifest)
    out = Path(args.output_dir)
    cache = Path(args.cache_dir)
    out.mkdir(parents=True, exist_ok=True)

    geometry = json.loads(geometry_path.read_text(encoding="utf-8"))
    if geometry.get("schema") != "janus.cosmos.tachyon_star.q4a.geometry_manifest.v1":
        raise RuntimeError("unexpected Q4A geometry schema")
    if geometry.get("image_pixel_arrays_dereferenced") is not False or geometry.get("new_forced_photometry_performed") is not False:
        raise RuntimeError("Q4A geometry lineage indicates prior pixel measurement")
    if len(geometry.get("detector_controls", [])) != 13:
        raise RuntimeError("Q4A detector-control cardinality changed")

    by_role = {r["role"]: r for r in geometry["frame_headers"]}
    sky_pos = {r["role"]: r for r in geometry["sky_fixed_track"]}
    det_pos = {r["role"]: r for r in geometry["detector_fixed_track"]}

    images = {}
    frame_provenance = []
    for role in FRAME_ORDER:
        filename = by_role[role]["filename"]
        path = download_exact(filename, cache)
        got = sha256_file(path)
        if got != FRAME_SHA256[filename]:
            raise RuntimeError(f"frame SHA changed: {filename}")
        image = np.asarray(fits.getdata(path, ext=0), dtype=float)
        if image.shape != (2048, 2048):
            raise RuntimeError(f"unexpected image shape {image.shape} for {filename}")
        images[role] = image
        frame_provenance.append({"role": role, "filename": filename, "sha256": got, "shape": list(image.shape)})

    sky = {}
    detector = {}
    for role in FRAME_ORDER:
        filename = by_role[role]["filename"]
        sky[role] = forced_gaussian_plane(images[role], sky_pos[role]["x_zero_based"], sky_pos[role]["y_zero_based"], FWHM[filename])
        detector[role] = forced_gaussian_plane(images[role], det_pos[role]["x_zero_based"], det_pos[role]["y_zero_based"], FWHM[filename])

    controls = []
    for ci, pos in enumerate(geometry["detector_controls"]):
        row = {"control_id": f"Q4A-DETCTRL-{ci+1:02d}", "x_zero_based": pos["x_zero_based"], "y_zero_based": pos["y_zero_based"], "fits": {}}
        for role in FRAME_ORDER:
            filename = by_role[role]["filename"]
            row["fits"][role] = forced_gaussian_plane(images[role], pos["x_zero_based"], pos["y_zero_based"], FWHM[filename])
        row["temporal_contrast"] = temporal_contrast(row["fits"]["A_BEFORE"], row["fits"]["B_CANDIDATE"], row["fits"]["C_AFTER"])
        controls.append(row)

    cls = classify(sky, detector)
    b = sky["B_CANDIDATE"]
    sky_contrast = temporal_contrast(sky["A_BEFORE"], b, sky["C_AFTER"])
    det_contrast = temporal_contrast(detector["A_BEFORE"], b, detector["C_AFTER"])
    control_contrasts = [q["temporal_contrast"] for q in controls]
    sky_rank = empirical_rank_p(sky_contrast, control_contrasts)
    det_rank = empirical_rank_p(det_contrast, control_contrasts)

    detectability = {
        "A_SKY_B_EQUIVALENT_SNR": float(b["amplitude"] / sky["A_BEFORE"]["amplitude_sigma"]),
        "C_SKY_B_EQUIVALENT_SNR": float(b["amplitude"] / sky["C_AFTER"]["amplitude_sigma"]),
        "A_DETECTOR_B_EQUIVALENT_SNR": float(b["amplitude"] / detector["A_BEFORE"]["amplitude_sigma"]),
        "C_DETECTOR_B_EQUIVALENT_SNR": float(b["amplitude"] / detector["C_AFTER"]["amplitude_sigma"]),
    }

    measurements = {
        "schema": "janus.cosmos.tachyon_star.q4b.measurements.v1",
        "experiment_id": "JANUS-TACHYON-STAR-Q4B-LUCI-FIXED-COORDINATE-FORCED-PHOTOMETRY",
        "frame_provenance": frame_provenance,
        "sky_fixed": sky,
        "detector_fixed": detector,
        "detector_controls": controls,
        "sky_temporal_contrast": sky_contrast,
        "detector_temporal_contrast": det_contrast,
        "sky_temporal_contrast_rank": sky_rank,
        "detector_temporal_contrast_rank": det_rank,
        "b_equivalent_detectability": detectability,
        "classification_metrics": cls,
    }
    mpath = out / "measurements.json"
    mpath.write_text(json.dumps(measurements, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    rec = {
        "schema": "janus.cosmos.tachyon_star.q4b.receipt.v1",
        "experiment_id": "JANUS-TACHYON-STAR-Q4B-LUCI-FIXED-COORDINATE-FORCED-PHOTOMETRY",
        "status": cls["classification"],
        "b_fixed_amplitude": b["amplitude"],
        "b_fixed_amplitude_sigma": b["amplitude_sigma"],
        "b_fixed_forced_z": b["forced_z"],
        "sky_persistence_compatible": cls.get("sky_persistence_compatible"),
        "detector_persistence_compatible": cls.get("detector_persistence_compatible"),
        "sky_constant_flux_difference_z": cls.get("sky_constant_flux_difference_z"),
        "detector_constant_flux_difference_z": cls.get("detector_constant_flux_difference_z"),
        "sky_temporal_contrast": sky_contrast,
        "detector_temporal_contrast": det_contrast,
        "sky_temporal_rank": sky_rank,
        "detector_temporal_rank": det_rank,
        "b_equivalent_detectability": detectability,
        "control_count": len(controls),
        "measurements_sha256": sha256_file(mpath),
        "origin_conclusion_limit": "If neither persistence track is compatible, one-frame sky transient and one-frame detector/cosmic-ray/background explanations remain unresolved.",
        "claim_ceiling": CLAIM,
    }
    (out / "receipt.json").write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(rec, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
