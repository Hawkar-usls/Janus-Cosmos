#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, math, time, urllib.parse, urllib.request
from pathlib import Path

import numpy as np

from janus_cosmos.pipeline import EventWriter, download_source, sha256_file
from janus_cosmos.luci_psf import _inject_gaussian, robust_background
from janus_cosmos.luci_psf_r1 import estimate_native_psf_fwhm, measure_psf_at
from experiments.luci.run_tachyon_star_t2b_ztf import (
    MANIFEST_GIT_BLOB_SHA1,
    MANIFEST_SHA256_OBSERVED_PRE_PIXEL,
    _git_blob_sha1,
    _mask_clean,
    _read_image_wcs,
    _rows,
    _xy,
)
from experiments.luci.run_tachyon_star_t2c_ztf import _product_urls, gaussian_retained_fraction

PARENT_RESULT_COMMIT = "40d391df27a928cea8db38226c30581aec9a5599"
EXPECTED_PARENT_QUALIFIED_NULLS = 23
SENSITIVITY_TARGETS = {
    "TS-T2B-ZTF-03": ("b", "c"),
    "TS-T2B-ZTF-09": ("b",),
    "TS-T2B-ZTF-10": ("a", "b", "c"),
}
PROVENANCE_TARGETS = {
    "TS-T2B-ZTF-16": ("b",),
    "TS-T2B-ZTF-39": ("a", "b"),
}
FWHM_RATIOS = (0.8, 1.0, 1.2)
SNR_GRID = (8.0, 10.0, 12.0, 15.0, 20.0)
HOT_PIXEL_SNR = 12.0
MEASURE_RADIUS_PX = 6
RETAINED_FRACTION_MIN = 0.995
CLAIM = "ZTF_RESIDUAL_FORENSIC_CHARACTERIZATION_ONLY__NO_NEW_NEGATIVE_UPGRADE__NO_FTL__NO_RETROCAUSALITY__NO_TACHYON_PARTICLE_IDENTIFICATION__NO_NUCLEAR_CAUSALITY"


def _get(url: str, timeout: int = 120, retries: int = 2) -> bytes:
    err = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "Janus-Cosmos-TachyonStar-T2D/1.0", "Accept": "text/csv,text/plain,*/*"},
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as exc:
            err = f"{type(exc).__name__}: {exc}"
            if attempt < retries:
                time.sleep(2 * attempt)
    raise RuntimeError(f"GET failed after {retries} attempts: {url}: {err}")


def _download_cutout(row: dict, epoch: str, suffix: str, cache: Path, events: EventWriter, label: str):
    mode, url = _product_urls(row, epoch, suffix)[0]
    path, meta = download_source(url, cache, events, target=row["trial_id"], filter_name=f"T2D_{label}_{epoch}_{mode}")
    return path, meta, url


def sensitivity_floor(all_by_snr: dict[float, bool]) -> float | None:
    grid = sorted(float(s) for s in all_by_snr)
    for s in grid:
        if all(bool(all_by_snr[q]) for q in grid if q >= s):
            return s
    return None


def sensitivity_curve(image: np.ndarray, wcs, mask: np.ndarray, mask_wcs, ra, dec) -> dict:
    x, y = _xy(wcs, ra, dec)
    h, w = image.shape
    if not (math.isfinite(x) and math.isfinite(y) and 0 <= x < w and 0 <= y < h):
        return {"status": "CURVE_BLOCKED_TARGET_OUTSIDE_WCS", "x": x, "y": y, "shape": [h, w]}

    mg = _mask_clean(mask, mask_wcs, ra, dec)
    if not mg.get("passed"):
        return {"status": "CURVE_BLOCKED_TARGET_MASK", "x": x, "y": y, "mask_gate": mg}

    target = measure_psf_at(image, y, x)
    if target is not None:
        return {"status": "CURVE_BLOCKED_SOURCE_PRESENT", "x": x, "y": y, "mask_gate": mg, "target": target.to_dict()}

    native_fwhm, native_n = estimate_native_psf_fwhm(image)
    if native_fwhm is None:
        return {"status": "CURVE_BLOCKED_NATIVE_PSF_REFERENCE", "x": x, "y": y, "mask_gate": mg, "native_psf_source_count": native_n}

    base = max(1.5, min(8.0, float(native_fwhm)))
    retained = gaussian_retained_fraction(image.shape, x, y, base * max(FWHM_RATIOS))
    window_ok = (
        MEASURE_RADIUS_PX <= x < w - MEASURE_RADIUS_PX
        and MEASURE_RADIUS_PX <= y < h - MEASURE_RADIUS_PX
    )
    if retained < RETAINED_FRACTION_MIN or not window_ok:
        return {
            "status": "CURVE_BLOCKED_TARGET_EDGE",
            "x": x, "y": y,
            "native_psf_source_count": native_n,
            "native_psf_fwhm_px": base,
            "retained_fraction": retained,
            "measurement_window_ok": bool(window_ok),
        }

    med, sigma = robust_background(image)
    hot = np.array(image, copy=True)
    iy, ix = int(round(y)), int(round(x))
    hot[iy, ix] += HOT_PIXEL_SNR * sigma
    hot_accepted = measure_psf_at(hot, y, x) is not None

    trials = []
    for snr in SNR_GRID:
        for ratio in FWHM_RATIOS:
            fwhm = max(1.5, min(10.0, base * ratio))
            z = np.array(image, copy=True)
            _inject_gaussian(z, y, x, fwhm, snr * sigma)
            q = measure_psf_at(z, y, x)
            trials.append({
                "snr": snr,
                "fwhm_ratio": ratio,
                "fwhm_px": fwhm,
                "recovered": q is not None,
                "measured": None if q is None else q.to_dict(),
            })

    all_by_snr = {
        float(snr): all(t["recovered"] for t in trials if float(t["snr"]) == float(snr))
        for snr in SNR_GRID
    }
    floor = sensitivity_floor(all_by_snr)
    status = "CURVE_COMPLETE" if not hot_accepted else "CURVE_BLOCKED_HOT_PIXEL_ACCEPTED"
    return {
        "status": status,
        "x": x, "y": y,
        "mask_gate": mg,
        "native_psf_source_count": native_n,
        "native_psf_fwhm_px": base,
        "background_median": med,
        "background_sigma": sigma,
        "retained_fraction": retained,
        "hot_pixel_snr": HOT_PIXEL_SNR,
        "hot_pixel_accepted": bool(hot_accepted),
        "snr_grid": list(SNR_GRID),
        "fwhm_ratios": list(FWHM_RATIOS),
        "all_ratios_recovered_by_snr": {str(k): bool(v) for k, v in all_by_snr.items()},
        "sensitivity_floor_snr": floor,
        "sensitivity_floor_label": f"{floor:g}" if floor is not None else "GREATER_THAN_20",
        "trials": trials,
        "interpretation": "Sensitivity characterization only; this does not upgrade T2C to QUALIFIED_ABSENCE.",
    }


def _same_int(a, b) -> bool:
    try:
        return int(str(a).strip()) == int(str(b).strip())
    except Exception:
        return str(a).strip() == str(b).strip()


def exact_metadata_match(row: dict, target: dict) -> bool:
    if str(row.get("filefracday", "")).strip() != str(target["filefracday"]).strip():
        return False
    for key in ("field", "ccdid", "qid"):
        if not _same_int(row.get(key, ""), target[key]):
            return False
    if str(row.get("filtercode", "")).strip() != str(target["filtercode"]).strip():
        return False
    if "pid" in row and str(row.get("pid", "")).strip():
        if str(row.get("pid", "")).strip() != str(target["pid"]).strip():
            return False
    return True


def _path_fields(row: dict) -> dict:
    out = {}
    for k, v in row.items():
        kl = str(k).lower()
        s = str(v or "")
        if any(token in kl for token in ("file", "path", "url", "uri")) or ".fits" in s.lower():
            if s:
                out[k] = s
    return out


def provenance_audit(row: dict, epoch: str, out: Path) -> dict:
    ra, dec = row["ra_deg"], row["dec_deg"]
    url = "https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?" + urllib.parse.urlencode({
        "POS": f"{ra},{dec}",
        "ct": "csv",
    })
    body = _get(url)
    raw_path = out / f"provenance_{row['trial_id']}_{epoch}.csv"
    raw_path.write_bytes(body)
    parsed = list(csv.DictReader(body.decode("utf-8-sig").splitlines()))
    target = {
        "filefracday": row[f"{epoch}_filefracday"],
        "pid": row[f"{epoch}_pid"],
        "field": row["field"],
        "ccdid": row["ccdid"],
        "qid": row["qid"],
        "filtercode": row["filtercode"],
    }
    exact = [q for q in parsed if exact_metadata_match(q, target)]
    return {
        "status": "EXACT_METADATA_MATCH_FOUND" if exact else "NO_EXACT_METADATA_MATCH",
        "trial_id": row["trial_id"],
        "epoch": epoch,
        "query_url": url,
        "query_sha256": sha256_file(raw_path),
        "query_rows": len(parsed),
        "query_columns": list(parsed[0].keys()) if parsed else [],
        "target_identity": target,
        "exact_match_count": len(exact),
        "exact_rows": exact,
        "path_like_fields": [_path_fields(q) for q in exact],
        "pixels_opened": False,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="data/tachyon_star/JANUS-TACHYON-STAR-T2B-ZTF-PREPOINTED-MANIFEST.csv")
    ap.add_argument("--output-dir", default="results/tachyon_star_t2d_ztf")
    ap.add_argument("--cache-dir", default=".cache/tachyon_star_t2d_ztf")
    args = ap.parse_args()

    manifest = Path(args.manifest)
    out = Path(args.output_dir)
    cache = Path(args.cache_dir)
    out.mkdir(parents=True, exist_ok=True)

    if _git_blob_sha1(manifest) != MANIFEST_GIT_BLOB_SHA1:
        raise RuntimeError("parent T2B manifest Git blob changed")
    if sha256_file(manifest) != MANIFEST_SHA256_OBSERVED_PRE_PIXEL:
        raise RuntimeError("parent T2B manifest SHA256 changed")
    by_id = {q["trial_id"]: q for q in _rows(manifest)}

    events = EventWriter(out / "events.jsonl")
    curve_receipts = []
    for trial_id, epochs in SENSITIVITY_TARGETS.items():
        row = by_id[trial_id]
        for epoch in epochs:
            try:
                sci, smeta, sci_url = _download_cutout(row, epoch, "sciimg.fits", cache, events, "SCI")
                msk, mmeta, msk_url = _download_cutout(row, epoch, "mskimg.fits", cache, events, "MASK")
                image, wcs, sci_hdu = _read_image_wcs(sci)
                mask, mwcs, mask_hdu = _read_image_wcs(msk)
                curve = sensitivity_curve(image, wcs, mask, mwcs, row["ra_deg"], row["dec_deg"])
                curve.update({
                    "trial_id": trial_id,
                    "src_id": row["src_id"],
                    "epoch": epoch,
                    "science_url": sci_url,
                    "mask_url": msk_url,
                    "science_sha256": smeta["sha256"],
                    "mask_sha256": mmeta["sha256"],
                    "science_hdu": sci_hdu,
                    "mask_hdu": mask_hdu,
                })
            except Exception as exc:
                curve = {
                    "status": "CURVE_TRANSPORT_OR_REPLAY_ERROR",
                    "trial_id": trial_id,
                    "src_id": row["src_id"],
                    "epoch": epoch,
                    "error": f"{type(exc).__name__}: {exc}",
                }
            curve_receipts.append(curve)

    provenance_receipts = []
    for trial_id, epochs in PROVENANCE_TARGETS.items():
        row = by_id[trial_id]
        for epoch in epochs:
            try:
                q = provenance_audit(row, epoch, out)
            except Exception as exc:
                q = {
                    "status": "PROVENANCE_QUERY_ERROR",
                    "trial_id": trial_id,
                    "src_id": row["src_id"],
                    "epoch": epoch,
                    "error": f"{type(exc).__name__}: {exc}",
                    "pixels_opened": False,
                }
            provenance_receipts.append(q)

    curves_complete = sum(q["status"] == "CURVE_COMPLETE" for q in curve_receipts)
    provenance_complete = sum(q["status"] in ("EXACT_METADATA_MATCH_FOUND", "NO_EXACT_METADATA_MATCH") for q in provenance_receipts)
    status = (
        "PASS_FORENSIC_CHARACTERIZATION_COMPLETE"
        if curves_complete == 6 and provenance_complete == 3
        else "BLOCKED_FORENSIC_CHARACTERIZATION_INCOMPLETE"
    )
    rec = {
        "schema": "janus.cosmos.tachyon_star.t2d.ztf_residual_forensics.receipt.v1",
        "experiment_id": "JANUS-TACHYON-STAR-T2D-ZTF-RESIDUAL-FORENSICS",
        "status": status,
        "parent_result_commit": PARENT_RESULT_COMMIT,
        "parent_qualified_null_trials_immutable": EXPECTED_PARENT_QUALIFIED_NULLS,
        "parent_qualified_null_trials_reopened": 0,
        "sensitivity_curve_targets": 6,
        "sensitivity_curves_complete": curves_complete,
        "product_provenance_targets": 3,
        "product_provenance_audits_complete": provenance_complete,
        "sensitivity_curves": curve_receipts,
        "product_provenance_audits": provenance_receipts,
        "claim_ceiling": CLAIM,
    }
    (out / "receipt.json").write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": status,
        "sensitivity_curves_complete": curves_complete,
        "product_provenance_audits_complete": provenance_complete,
        "sensitivity_floors": [
            {
                "trial_id": q.get("trial_id"),
                "epoch": q.get("epoch"),
                "status": q.get("status"),
                "floor": q.get("sensitivity_floor_label"),
            } for q in curve_receipts
        ],
        "provenance": [
            {
                "trial_id": q.get("trial_id"),
                "epoch": q.get("epoch"),
                "status": q.get("status"),
                "exact_match_count": q.get("exact_match_count"),
            } for q in provenance_receipts
        ],
    }, indent=2))
    return 0 if status == "PASS_FORENSIC_CHARACTERIZATION_COMPLETE" else 3


if __name__ == "__main__":
    raise SystemExit(main())
