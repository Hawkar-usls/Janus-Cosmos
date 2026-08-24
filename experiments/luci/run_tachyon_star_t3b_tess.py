#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

from experiments.luci.run_tachyon_star_t3a_tess import _extract_single_fits, _fetch, sha256_file, tesscut_url

PARENT_MANIFEST_SHA256 = "be8bc23021e9b9a6fa6efa65815e79bc39a1bb011192acab2825c013c1ac4202"
EXPECTED_SOURCES = 42
TIME_TOL_DAY = 1e-9
CONTROL_RADIUS = 120
CONTROL_EXCLUDE_B_DISTANCE = 3
MIN_CONTROLS = 30
Z_MIN = 8.0
NEIGHBOR_DELTA_MIN_SIGMA = 4.0
INJECTION_SNRS = (8.0, 10.0, 12.0)
CLAIM = "PREPOINTED_TESS_CADENCE_NATIVE_APERTURE_EVENT_REPLICATION_ONLY__NO_SOURCE_IDENTITY__NO_FTL__NO_RETROCAUSALITY__NO_TACHYON_PARTICLE_IDENTIFICATION__NO_UAP_OR_ARTIFICIAL_ORIGIN_CLAIM"


def _rows(path: Path) -> list[dict]:
    return list(csv.DictReader(path.open("r", encoding="utf-8", newline="")))


def read_flux_table(path: Path) -> dict:
    with fits.open(path, mode="readonly", memmap=True, lazy_load_hdus=True) as hdul:
        chosen = None
        for i, hdu in enumerate(hdul):
            cols = list(getattr(getattr(hdu, "columns", None), "names", []) or [])
            upper = {str(c).upper(): str(c) for c in cols}
            if "TIME" in upper and "QUALITY" in upper and "FLUX" in upper:
                chosen = (i, hdu, cols, upper)
                break
        if chosen is None:
            raise RuntimeError("no TESS table containing TIME, QUALITY and FLUX")
        hdu_index, hdu, cols, upper = chosen
        data = hdu.data
        t = np.asarray(data[upper["TIME"]], dtype=float).copy()
        q = np.asarray(data[upper["QUALITY"]], dtype=np.int64).copy()
        if "CADENCENO" in upper:
            c = np.asarray(data[upper["CADENCENO"]], dtype=np.int64).copy()
        else:
            c = np.arange(len(t), dtype=np.int64)
        flux = np.asarray(data[upper["FLUX"]], dtype=float).copy()
        if "FLUX_BKG" in upper:
            bkg = np.asarray(data[upper["FLUX_BKG"]], dtype=float).copy()
        else:
            bkg = None
        return {
            "hdu": int(hdu_index),
            "columns": cols,
            "time": t,
            "quality": q,
            "cadenceno": c,
            "flux": flux,
            "flux_bkg": bkg,
        }


def aperture_series(flux: np.ndarray, flux_bkg: np.ndarray | None = None) -> tuple[np.ndarray, dict]:
    a = np.asarray(flux, dtype=float)
    if a.ndim != 3:
        raise RuntimeError(f"expected 3-D TESS FLUX cube, got shape {a.shape}")
    n, h, w = a.shape
    if h < 9 or w < 9 or h % 2 == 0 or w % 2 == 0:
        raise RuntimeError(f"expected odd >=9 pixel cutout, got {h}x{w}")
    if flux_bkg is not None:
        b = np.asarray(flux_bkg, dtype=float)
        if b.shape != a.shape:
            raise RuntimeError("FLUX_BKG shape mismatch")
        net = a - b
        bkg_mode = "FLUX_MINUS_FLUX_BKG"
    else:
        net = a
        bkg_mode = "FLUX_ONLY"
    cy, cx = h // 2, w // 2
    target = np.zeros((h, w), dtype=bool)
    target[cy - 1:cy + 2, cx - 1:cx + 2] = True
    yy, xx = np.indices((h, w))
    outer = np.maximum(np.abs(yy - cy), np.abs(xx - cx)) >= 4
    out = np.full(n, np.nan, dtype=float)
    for i in range(n):
        frame = net[i]
        tv = frame[target]
        ov = frame[outer]
        if np.all(np.isfinite(tv)) and np.count_nonzero(np.isfinite(ov)) >= 20:
            out[i] = float(np.sum(tv) - 9.0 * np.nanmedian(ov))
    return out, {
        "shape": [int(h), int(w)],
        "center_xy": [int(cx), int(cy)],
        "target_pixels": int(target.sum()),
        "background_pixels": int(outer.sum()),
        "background_mode": bkg_mode,
    }


def isolated_contrast(series: np.ndarray, i: int) -> float:
    s = np.asarray(series, dtype=float)
    return float(s[i] - 0.5 * (s[i - 1] + s[i + 1]))


def robust_center_sigma(values: list[float]) -> tuple[float, float, str]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < MIN_CONTROLS:
        raise RuntimeError(f"insufficient control contrasts: {x.size}/{MIN_CONTROLS}")
    center = float(np.median(x))
    mad = float(np.median(np.abs(x - center)))
    sigma = 1.4826 * mad
    mode = "MAD"
    if not math.isfinite(sigma) or sigma <= 0:
        sigma = float(np.std(x, ddof=1))
        mode = "STD_FALLBACK"
    if not math.isfinite(sigma) or sigma <= 0:
        raise RuntimeError("nonpositive robust control sigma")
    return center, sigma, mode


def control_contrasts(series: np.ndarray, quality: np.ndarray, b_row: int) -> list[float]:
    s = np.asarray(series, dtype=float)
    q = np.asarray(quality, dtype=np.int64)
    lo = max(1, int(b_row) - CONTROL_RADIUS)
    hi = min(len(s) - 1, int(b_row) + CONTROL_RADIUS + 1)
    out = []
    for i in range(lo, hi):
        if abs(i - int(b_row)) <= CONTROL_EXCLUDE_B_DISTANCE:
            continue
        if not (int(q[i - 1]) == 0 and int(q[i]) == 0 and int(q[i + 1]) == 0):
            continue
        if not (np.isfinite(s[i - 1]) and np.isfinite(s[i]) and np.isfinite(s[i + 1])):
            continue
        out.append(isolated_contrast(s, i))
    return out


def candidate_metrics(series: np.ndarray, b_row: int, center: float, sigma: float) -> dict:
    s = np.asarray(series, dtype=float)
    a, b, c = int(b_row) - 1, int(b_row), int(b_row) + 1
    if not (0 <= a < b < c < len(s)):
        return {"passed": False, "reason": "ABC_OUT_OF_RANGE"}
    if not (np.isfinite(s[a]) and np.isfinite(s[b]) and np.isfinite(s[c])):
        return {"passed": False, "reason": "NONFINITE_ABC_APERTURE"}
    contrast = isolated_contrast(s, b)
    z = float((contrast - center) / sigma)
    delta_a = float((s[b] - s[a]) / sigma)
    delta_c = float((s[b] - s[c]) / sigma)
    passed = bool(z >= Z_MIN and delta_a >= NEIGHBOR_DELTA_MIN_SIGMA and delta_c >= NEIGHBOR_DELTA_MIN_SIGMA)
    return {
        "passed": passed,
        "reason": "PASS" if passed else "THRESHOLD_NOT_MET",
        "a_flux": float(s[a]),
        "b_flux": float(s[b]),
        "c_flux": float(s[c]),
        "isolated_contrast": contrast,
        "contrast_z": z,
        "b_minus_a_sigma": delta_a,
        "b_minus_c_sigma": delta_c,
        "thresholds": {
            "contrast_z_min": Z_MIN,
            "b_minus_a_sigma_min": NEIGHBOR_DELTA_MIN_SIGMA,
            "b_minus_c_sigma_min": NEIGHBOR_DELTA_MIN_SIGMA,
        },
    }


def verify_frozen_cadence(row: dict, table: dict) -> dict:
    idx = {k: int(row[k]) for k in ("a_row", "b_row", "c_row")}
    n = len(table["time"])
    if not (0 <= idx["a_row"] < idx["b_row"] < idx["c_row"] < n):
        raise RuntimeError("frozen A/B/C rows outside regenerated table")
    if not (idx["b_row"] == idx["a_row"] + 1 and idx["c_row"] == idx["b_row"] + 1):
        raise RuntimeError("frozen A/B/C rows are not consecutive")
    for prefix in ("a", "b", "c"):
        i = idx[f"{prefix}_row"]
        if int(table["quality"][i]) != 0:
            raise RuntimeError(f"frozen {prefix.upper()} QUALITY changed from zero")
        if abs(float(table["time"][i]) - float(row[f"{prefix}_time"])) > TIME_TOL_DAY:
            raise RuntimeError(f"frozen {prefix.upper()} TIME changed")
        if int(table["cadenceno"][i]) != int(row[f"{prefix}_cadenceno"]):
            raise RuntimeError(f"frozen {prefix.upper()} CADENCENO changed")
    return idx


def classify_trial(series: np.ndarray, quality: np.ndarray, b_row: int) -> dict:
    controls = control_contrasts(series, quality, b_row)
    center, sigma, sigma_mode = robust_center_sigma(controls)
    raw = candidate_metrics(series, b_row, center, sigma)
    injections = []
    for snr in INJECTION_SNRS:
        z = np.array(series, copy=True, dtype=float)
        z[int(b_row)] += float(snr) * sigma
        m = candidate_metrics(z, b_row, center, sigma)
        injections.append({"injected_snr": snr, "candidate_recovered": bool(m.get("passed")), "metrics": m})
    recovered = {float(q["injected_snr"]): bool(q["candidate_recovered"]) for q in injections}
    if raw.get("passed"):
        classification = "CANDIDATE_REQUIRES_ADJUDICATION"
    elif raw.get("contrast_z", -math.inf) >= Z_MIN:
        classification = "NON_ISOLATED_OR_VARIABLE_PATTERN"
    elif recovered.get(10.0, False) and recovered.get(12.0, False):
        classification = "QUALIFIED_NO_ISOLATED_B_EVENT"
    else:
        classification = "UNRESOLVED_TRIAL"
    return {
        "classification": classification,
        "control_count": len(controls),
        "control_center": center,
        "control_sigma": sigma,
        "control_sigma_mode": sigma_mode,
        "raw": raw,
        "injections": injections,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepoint-manifest", required=True)
    ap.add_argument("--output-dir", default="results/tachyon_star_t3b_tess")
    ap.add_argument("--cache-dir", default=".cache/tachyon_star_t3b_tess")
    args = ap.parse_args()

    manifest = Path(args.prepoint_manifest)
    out = Path(args.output_dir)
    cache = Path(args.cache_dir)
    out.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)

    if sha256_file(manifest) != PARENT_MANIFEST_SHA256:
        raise RuntimeError("T3A prepointed cadence manifest SHA mismatch")
    rows = _rows(manifest)
    if len(rows) != EXPECTED_SOURCES or len({r["src_id"] for r in rows}) != EXPECTED_SOURCES:
        raise RuntimeError("T3A prepoint manifest cardinality changed")

    results = []
    for i, row in enumerate(rows):
        sid = row["src_id"]
        safe = f"{i:02d}_{hashlib.sha256(sid.encode()).hexdigest()[:12]}"
        url = tesscut_url(row)
        try:
            payload = _fetch(url)
            zip_path = cache / f"{safe}.zip"
            zip_path.write_bytes(payload)
            fits_path, source_name = _extract_single_fits(payload, cache, safe)
            if source_name != row["source_fits_name"]:
                raise RuntimeError(f"source FITS identity changed: {source_name} != {row['source_fits_name']}")
            table = read_flux_table(fits_path)
            idx = verify_frozen_cadence(row, table)
            series, aperture = aperture_series(table["flux"], table["flux_bkg"])
            q = classify_trial(series, table["quality"], idx["b_row"])
            q.update({
                "trial_id": f"TS-T3B-TESS-{i+1:02d}",
                "src_id": sid,
                "sector": int(row["sector"]),
                "ffi_cadence_s": int(row["ffi_cadence_s"]),
                "a_row": idx["a_row"],
                "b_row": idx["b_row"],
                "c_row": idx["c_row"],
                "a_cadenceno": int(row["a_cadenceno"]),
                "b_cadenceno": int(row["b_cadenceno"]),
                "c_cadenceno": int(row["c_cadenceno"]),
                "a_time": float(row["a_time"]),
                "b_time": float(row["b_time"]),
                "c_time": float(row["c_time"]),
                "tesscut_url": url,
                "source_fits_name": source_name,
                "parent_fits_sha256": row["tesscut_fits_sha256"],
                "replayed_fits_sha256": sha256_file(fits_path),
                "exact_fits_sha256_match": sha256_file(fits_path) == row["tesscut_fits_sha256"],
                "parent_zip_sha256": row["tesscut_zip_sha256"],
                "replayed_zip_sha256": sha256_file(zip_path),
                "exact_zip_sha256_match": sha256_file(zip_path) == row["tesscut_zip_sha256"],
                "table_hdu": table["hdu"],
                "aperture": aperture,
            })
        except Exception as exc:
            q = {
                "trial_id": f"TS-T3B-TESS-{i+1:02d}",
                "src_id": sid,
                "sector": int(row["sector"]),
                "classification": "UNRESOLVED_TRIAL",
                "error": f"{type(exc).__name__}: {exc}",
            }
        results.append(q)
        time.sleep(0.25)

    candidates = [r for r in results if r["classification"] == "CANDIDATE_REQUIRES_ADJUDICATION"]
    nulls = [r for r in results if r["classification"] == "QUALIFIED_NO_ISOLATED_B_EVENT"]
    nonisolated = [r for r in results if r["classification"] == "NON_ISOLATED_OR_VARIABLE_PATTERN"]
    unresolved = [r for r in results if r["classification"] == "UNRESOLVED_TRIAL"]
    if candidates:
        status = "CANDIDATE_REQUIRES_ADJUDICATION"
    elif len(nulls) == EXPECTED_SOURCES and not nonisolated and not unresolved:
        status = "PASS_ALL_42_PREPOINTED_NULL"
    else:
        status = "BLOCKED_PREPOINTED_REPLICATION_INCOMPLETE"

    (out / "trial_results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    rec = {
        "schema": "janus.cosmos.tachyon_star.t3b.tess_prepointed_flux.receipt.v1",
        "experiment_id": "JANUS-TACHYON-STAR-T3B-TESS-PREPOINTED-FLUX-REPLICATION",
        "status": status,
        "parent_manifest_sha256": PARENT_MANIFEST_SHA256,
        "trials_total": len(results),
        "qualified_null_trials": len(nulls),
        "candidate_trials": len(candidates),
        "nonisolated_trials": len(nonisolated),
        "unresolved_trials": len(unresolved),
        "candidate_trial_ids": [r["trial_id"] for r in candidates],
        "nonisolated_trial_ids": [r["trial_id"] for r in nonisolated],
        "unresolved_trial_ids": [r["trial_id"] for r in unresolved],
        "exact_fits_sha_matches": sum(bool(r.get("exact_fits_sha256_match")) for r in results),
        "exact_zip_sha_matches": sum(bool(r.get("exact_zip_sha256_match")) for r in results),
        "trial_results_sha256": sha256_file(out / "trial_results.json"),
        "claim_ceiling": CLAIM,
    }
    (out / "receipt.json").write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(rec, indent=2, sort_keys=True))
    return 0 if status == "PASS_ALL_42_PREPOINTED_NULL" else 3


if __name__ == "__main__":
    raise SystemExit(main())
