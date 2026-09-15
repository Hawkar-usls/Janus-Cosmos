#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
from urllib.parse import urlencode

import numpy as np

from janus_cosmos.pipeline import EventWriter, download_source, sha256_file
from janus_cosmos.luci_psf import _inject_gaussian, robust_background
from janus_cosmos.luci_psf_r1 import estimate_native_psf_fwhm, measure_psf_at
from experiments.luci.run_tachyon_star_t2b_ztf import (
    MANIFEST_GIT_BLOB_SHA1,
    MANIFEST_SHA256_OBSERVED_PRE_PIXEL,
    _git_blob_sha1,
    _mask_clean,
    _quality_ok,
    _read_image_wcs,
    _rows,
    _xy,
)

PARENT_RESULT_COMMIT = "7971f48275dddcd4ae1b2d13962ac0ae7a40efb5"
EXPECTED_PARENT_NULLS = 20
UNRESOLVED_IDS = (
    "TS-T2B-ZTF-03",
    "TS-T2B-ZTF-09",
    "TS-T2B-ZTF-10",
    "TS-T2B-ZTF-16",
    "TS-T2B-ZTF-17",
    "TS-T2B-ZTF-27",
    "TS-T2B-ZTF-39",
    "TS-T2B-ZTF-40",
)
CUTOUT_ARCSEC = 360
MASK_RADIUS_PX = 3.0
RETAINED_FRACTION_MIN = 0.995
MEASURE_RADIUS_PX = 6
FWHM_RATIOS = (0.8, 1.0, 1.2)
SNR_GRID = (8.0, 12.0)
HOT_PIXEL_SNR = 12.0
CLAIM = "TARGETED_ZTF_NATIVE_SENSITIVITY_RECOVERY_ONLY__NO_FTL__NO_RETROCAUSALITY__NO_TACHYON_PARTICLE_IDENTIFICATION__NO_NUCLEAR_CAUSALITY__NO_UAP_OR_ARTIFICIAL_ORIGIN_CLAIM"


def _product_base(row: dict, epoch: str, suffix: str) -> str:
    ff = row[f"{epoch}_filefracday"]
    year, mmdd, frac = ff[:4], ff[4:8], ff[8:14]
    field = str(row["field"]).zfill(6)
    ccd = str(row["ccdid"]).zfill(2)
    qid = str(row["qid"])
    filt = row["filtercode"]
    itype = row[f"{epoch}_imgtypecode"]
    return f"https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/{year}/{mmdd}/{frac}/ztf_{ff}_{field}_{filt}_c{ccd}_{itype}_q{qid}_{suffix}"


def _product_urls(row: dict, epoch: str, suffix: str) -> list[tuple[str, str]]:
    base = _product_base(row, epoch, suffix)
    cut = base + "?" + urlencode({
        "center": f"{row['ra_deg']},{row['dec_deg']}",
        "size": f"{CUTOUT_ARCSEC}arcsec",
        "gzip": "false",
    })
    return [("CUTOUT_360_ARCSEC", cut), ("SAME_PRODUCT_FULL_FRAME", base)]


def _download_same_product(row, epoch, suffix, cache, events, trial_id, label):
    failures = []
    for mode, url in _product_urls(row, epoch, suffix):
        try:
            path, meta = download_source(
                url, cache, events,
                target=trial_id,
                filter_name=f"T2C_{epoch}_{label}_{mode}",
            )
            return path, meta, mode, url, failures
        except Exception as exc:
            failures.append({
                "mode": mode,
                "url": url,
                "error": f"{type(exc).__name__}: {exc}",
            })
    raise RuntimeError("same frozen ZTF product unavailable in preregistered retrieval modes: " + json.dumps(failures, sort_keys=True))


def gaussian_retained_fraction(shape: tuple[int, int], x: float, y: float, fwhm_px: float) -> float:
    h, w = shape
    sig = float(fwhm_px) / 2.354820045
    if not math.isfinite(sig) or sig <= 0:
        return 0.0
    s2 = math.sqrt(2.0) * sig
    fx = 0.5 * (math.erf((w - 0.5 - x) / s2) - math.erf((-0.5 - x) / s2))
    fy = 0.5 * (math.erf((h - 0.5 - y) / s2) - math.erf((-0.5 - y) / s2))
    return float(max(0.0, min(1.0, fx * fy)))


def local_native_sensitivity_gate(image: np.ndarray, x: float, y: float, native_fwhm_px: float) -> dict:
    a = np.asarray(image, dtype=float)
    h, w = a.shape
    base = max(1.5, min(8.0, float(native_fwhm_px)))
    retained = gaussian_retained_fraction(a.shape, x, y, max(base * max(FWHM_RATIOS), base))
    window_ok = (
        MEASURE_RADIUS_PX <= x < w - MEASURE_RADIUS_PX
        and MEASURE_RADIUS_PX <= y < h - MEASURE_RADIUS_PX
    )
    if retained < RETAINED_FRACTION_MIN or not window_ok:
        return {
            "passed": False,
            "reason": "TARGET_EDGE_NOT_ADMISSIBLE",
            "retained_fraction": retained,
            "threshold": RETAINED_FRACTION_MIN,
            "measurement_window_ok": bool(window_ok),
            "native_fwhm_px": base,
        }

    med, sigma = robust_background(a)
    hot = np.array(a, copy=True)
    iy, ix = int(round(y)), int(round(x))
    if not (0 <= iy < h and 0 <= ix < w):
        return {"passed": False, "reason": "TARGET_PIXEL_OUTSIDE_IMAGE"}
    hot[iy, ix] += HOT_PIXEL_SNR * sigma
    hot_accepted = measure_psf_at(hot, y, x) is not None

    trials = []
    for ratio in FWHM_RATIOS:
        fwhm = max(1.5, min(10.0, base * ratio))
        for snr in SNR_GRID:
            z = np.array(a, copy=True)
            _inject_gaussian(z, y, x, fwhm, snr * sigma)
            q = measure_psf_at(z, y, x)
            trials.append({
                "fwhm_ratio": ratio,
                "fwhm_px": fwhm,
                "snr": snr,
                "recovered": q is not None,
                "measured": None if q is None else q.to_dict(),
            })

    passed = bool(all(t["recovered"] for t in trials) and not hot_accepted)
    return {
        "schema": "janus.cosmos.tachyon_star.t2c.local_native_sensitivity.v1",
        "passed": passed,
        "reason": "PASS" if passed else ("HOT_PIXEL_ACCEPTED" if hot_accepted else "LOCAL_INJECTION_RECOVERY_FAILED"),
        "background_median": med,
        "background_sigma": sigma,
        "native_fwhm_px": base,
        "fwhm_ratios": list(FWHM_RATIOS),
        "snr_grid": list(SNR_GRID),
        "retained_fraction": retained,
        "retained_fraction_threshold": RETAINED_FRACTION_MIN,
        "hot_pixel_snr": HOT_PIXEL_SNR,
        "hot_pixel_accepted": bool(hot_accepted),
        "trials": trials,
    }


def classify_epoch_native(image, wcs, mask, mask_wcs, ra, dec):
    x, y = _xy(wcs, ra, dec)
    h, w = image.shape
    if not (math.isfinite(x) and math.isfinite(y) and 0 <= x < w and 0 <= y < h):
        return {"status": "BLOCKED_TARGET_OUTSIDE_WCS", "x": x, "y": y, "shape": [h, w]}

    mg = _mask_clean(mask, mask_wcs, ra, dec)
    target = measure_psf_at(image, y, x)
    if target is not None:
        return {
            "status": "SOURCE_PRESENT" if mg.get("passed") else "SOURCE_PRESENT_MASK_BLOCKED",
            "x": x, "y": y, "mask_gate": mg, "target": target.to_dict(),
        }

    native_fwhm, native_n = estimate_native_psf_fwhm(image)
    if native_fwhm is None:
        return {
            "status": "BLOCKED_NATIVE_PSF_REFERENCE",
            "x": x, "y": y, "mask_gate": mg, "native_psf_source_count": native_n,
        }

    local = local_native_sensitivity_gate(image, x, y, native_fwhm)
    if not mg.get("passed"):
        return {
            "status": "BLOCKED_TARGET_MASK",
            "x": x, "y": y, "mask_gate": mg,
            "native_psf_source_count": native_n,
            "native_psf_fwhm_px": native_fwhm,
            "local_gate": local,
        }
    if not local.get("passed"):
        return {
            "status": "BLOCKED_LOCAL_NATIVE_SENSITIVITY",
            "x": x, "y": y, "mask_gate": mg,
            "native_psf_source_count": native_n,
            "native_psf_fwhm_px": native_fwhm,
            "local_gate": local,
        }
    return {
        "status": "QUALIFIED_ABSENCE",
        "x": x, "y": y, "mask_gate": mg,
        "native_psf_source_count": native_n,
        "native_psf_fwhm_px": native_fwhm,
        "local_gate": local,
    }


def classify_trial(a: dict, b: dict, c: dict) -> str:
    if a["status"] == "QUALIFIED_ABSENCE" and b["status"] == "SOURCE_PRESENT" and c["status"] == "QUALIFIED_ABSENCE":
        return "ISOLATED_B_L0"
    if all(q["status"] == "QUALIFIED_ABSENCE" for q in (a, b, c)):
        return "NO_ISOLATED_B_EVENT"
    if "SOURCE_PRESENT" in a["status"] or "SOURCE_PRESENT" in c["status"]:
        return "NON_ISOLATED_SOURCE_PATTERN"
    return "UNRESOLVED_TRIAL"


def _blocked_epoch(exc: Exception) -> dict:
    return {"status": "BLOCKED_PIXEL_REPLAY_ERROR", "error": f"{type(exc).__name__}: {exc}"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="data/tachyon_star/JANUS-TACHYON-STAR-T2B-ZTF-PREPOINTED-MANIFEST.csv")
    ap.add_argument("--output-dir", default="results/tachyon_star_t2c_ztf")
    ap.add_argument("--cache-dir", default=".cache/tachyon_star_t2c_ztf")
    args = ap.parse_args()

    manifest = Path(args.manifest)
    out = Path(args.output_dir)
    cache = Path(args.cache_dir)
    out.mkdir(parents=True, exist_ok=True)

    if _git_blob_sha1(manifest) != MANIFEST_GIT_BLOB_SHA1:
        raise RuntimeError("parent T2B manifest Git blob changed")
    if sha256_file(manifest) != MANIFEST_SHA256_OBSERVED_PRE_PIXEL:
        raise RuntimeError("parent T2B manifest SHA256 changed")

    rows = _rows(manifest)
    by_id = {r["trial_id"]: r for r in rows}
    if set(UNRESOLVED_IDS) - set(by_id):
        raise RuntimeError("frozen unresolved trial missing from manifest")
    selected = [by_id[k] for k in UNRESOLVED_IDS]
    if len(selected) != 8 or len({r["trial_id"] for r in selected}) != 8:
        raise RuntimeError("T2C unresolved selection cardinality changed")
    if not all(_quality_ok(r) for r in selected):
        raise RuntimeError("T2C frozen unresolved set unexpectedly includes metadata-quality-blocked trial")

    events = EventWriter(out / "events.jsonl")
    results = []
    for ti, row in enumerate(selected):
        epochs = {}
        for epoch in "abc":
            try:
                sci, smeta, sci_mode, sci_url, sci_fail = _download_same_product(
                    row, epoch, "sciimg.fits", cache, events, row["trial_id"], "SCI"
                )
                msk, mmeta, msk_mode, msk_url, msk_fail = _download_same_product(
                    row, epoch, "mskimg.fits", cache, events, row["trial_id"], "MASK"
                )
                image, wcs, sci_hdu = _read_image_wcs(sci)
                mask, mwcs, mask_hdu = _read_image_wcs(msk)
                er = classify_epoch_native(image, wcs, mask, mwcs, row["ra_deg"], row["dec_deg"])
                er.update({
                    "science_retrieval_mode": sci_mode,
                    "mask_retrieval_mode": msk_mode,
                    "science_url": sci_url,
                    "mask_url": msk_url,
                    "science_prior_failures": sci_fail,
                    "mask_prior_failures": msk_fail,
                    "science_hdu": sci_hdu,
                    "mask_hdu": mask_hdu,
                    "science_sha256": smeta["sha256"],
                    "mask_sha256": mmeta["sha256"],
                })
            except Exception as exc:
                er = _blocked_epoch(exc)
            epochs[epoch] = er
        cls = classify_trial(epochs["a"], epochs["b"], epochs["c"])
        results.append({
            "trial_id": row["trial_id"],
            "src_id": row["src_id"],
            "classification": cls,
            "metadata": row,
            "epochs": epochs,
        })

    candidates = [r for r in results if r["classification"] == "ISOLATED_B_L0"]
    nulls = [r for r in results if r["classification"] == "NO_ISOLATED_B_EVENT"]
    nonisolated = [r for r in results if r["classification"] == "NON_ISOLATED_SOURCE_PATTERN"]
    unresolved = [r for r in results if r["classification"] == "UNRESOLVED_TRIAL"]
    combined_nulls = EXPECTED_PARENT_NULLS + len(nulls)

    if candidates:
        status = "CANDIDATE_REQUIRES_ADJUDICATION"
    elif len(nulls) == len(UNRESOLVED_IDS) and not unresolved and not nonisolated:
        status = "PASS_ALL_28_NULL"
    else:
        status = "BLOCKED_RECOVERY_INCOMPLETE"

    rec = {
        "schema": "janus.cosmos.tachyon_star.t2c.ztf_native_sensitivity_recovery.receipt.v1",
        "experiment_id": "JANUS-TACHYON-STAR-T2C-ZTF-NATIVE-SENSITIVITY-RECOVERY",
        "status": status,
        "parent_result_commit": PARENT_RESULT_COMMIT,
        "parent_qualified_null_trials_immutable": EXPECTED_PARENT_NULLS,
        "parent_qualified_null_trials_reopened": 0,
        "targeted_unresolved_trials": len(UNRESOLVED_IDS),
        "recovered_no_isolated_b_event": len(nulls),
        "recovery_candidates": len(candidates),
        "recovery_nonisolated": len(nonisolated),
        "recovery_unresolved": len(unresolved),
        "combined_qualified_null_trials_out_of_28": combined_nulls,
        "candidate_trial_ids": [r["trial_id"] for r in candidates],
        "remaining_unresolved_trial_ids": [r["trial_id"] for r in unresolved],
        "nonisolated_trial_ids": [r["trial_id"] for r in nonisolated],
        "results": results,
        "claim_ceiling": CLAIM,
    }
    (out / "receipt.json").write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: rec[k] for k in (
        "status",
        "parent_qualified_null_trials_immutable",
        "targeted_unresolved_trials",
        "recovered_no_isolated_b_event",
        "recovery_candidates",
        "recovery_nonisolated",
        "recovery_unresolved",
        "combined_qualified_null_trials_out_of_28",
        "candidate_trial_ids",
        "remaining_unresolved_trial_ids",
        "nonisolated_trial_ids",
    )}, indent=2))
    return 0 if status == "PASS_ALL_28_NULL" else 3


if __name__ == "__main__":
    raise SystemExit(main())
