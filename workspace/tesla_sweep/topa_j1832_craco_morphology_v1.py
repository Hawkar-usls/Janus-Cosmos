#!/usr/bin/env python3
"""Frozen TOPA CRACO morphology/DM replay for ASKAP J1832-0911.

Input is the author-exported Zenodo NPY only after checksum/header validation.
No pickle is executed. This script characterizes one discovery filterbank; it does
not classify origin, intelligence, or physical object type.
"""
from __future__ import annotations

import hashlib
import json
import tempfile
from pathlib import Path

import numpy as np
import requests
from scipy import signal

URL = "https://zenodo.org/records/15228816/files/SB55237_CRACO_filterbank.npy?download=1"
EXPECTED_MD5 = "41fca6b1dccc464e65e948ed6c82695e"
EXPECTED_SHAPE = (49152, 288)
DT = 0.013824
DF_MHZ = 1.0
FC_MHZ = 887.49074074
KDM_SEC_MHZ2 = 4.148808e3
NSUB = 24
COARSE_DMS = np.arange(0.0, 800.0 + 0.1, 10.0)
OUT = Path("data/tesla-sweep/results/TOPA-HUNT-005C3-J1832-CRACO-NUMERIC-MORPHOLOGY-RUN-001.json")
EPS = 1e-12


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def robust_center_scale(x: np.ndarray, axis=None):
    med = np.nanmedian(x, axis=axis, keepdims=True)
    mad = np.nanmedian(np.abs(x - med), axis=axis, keepdims=True)
    scale = 1.4826 * mad
    if np.ndim(scale) == 0:
        if not np.isfinite(scale) or scale < EPS:
            scale = np.nanstd(x) + EPS
    else:
        bad = (~np.isfinite(scale)) | (scale < EPS)
        if np.any(bad):
            sd = np.nanstd(x, axis=axis, keepdims=True)
            scale = np.where(bad, sd + EPS, scale)
    return med, scale


def robust_z_1d(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale < EPS:
        scale = np.std(x) + EPS
    return (x - med) / scale


def frequency_vector(nfreq: int, orientation: str) -> np.ndarray:
    idx = np.arange(nfreq, dtype=float) - (nfreq - 1) / 2.0
    f = FC_MHZ + idx * DF_MHZ
    if orientation == "descending":
        f = f[::-1]
    return f


def shift_earlier_accumulate(profiles: np.ndarray, freqs: np.ndarray, dm: float):
    """Dedisperse subband profiles to highest-frequency arrival time."""
    nt, nf = profiles.shape
    ref = float(np.max(freqs))
    delays = KDM_SEC_MHZ2 * dm * (freqs ** -2 - ref ** -2)
    shifts = np.rint(delays / DT).astype(int)
    acc = np.zeros(nt, dtype=np.float64)
    count = np.zeros(nt, dtype=np.float64)
    for j, sh in enumerate(shifts):
        p = profiles[:, j]
        if sh > 0:
            acc[: nt - sh] += p[sh:]
            count[: nt - sh] += 1.0
        elif sh < 0:
            s = -sh
            acc[s:] += p[: nt - s]
            count[s:] += 1.0
        else:
            acc += p
            count += 1.0
    good = count > 0
    out = np.zeros(nt, dtype=np.float64)
    out[good] = acc[good] / count[good]
    return out, shifts, delays


def score_dm(sub: np.ndarray, freqs: np.ndarray, dm: float) -> float:
    p, _, _ = shift_earlier_accumulate(sub, freqs, dm)
    z = robust_z_1d(p)
    edge = int(np.ceil((KDM_SEC_MHZ2 * 800 * ((np.min(freqs) ** -2) - (np.max(freqs) ** -2))) / DT)) + 5
    if edge * 2 < len(z):
        z = z[edge:-edge]
    return float(np.max(z))


def dm_search(sub: np.ndarray, freqs: np.ndarray):
    coarse = [(float(dm), score_dm(sub, freqs, float(dm))) for dm in COARSE_DMS]
    best_coarse = max(coarse, key=lambda t: t[1])
    lo = max(0.0, best_coarse[0] - 30.0)
    hi = min(800.0, best_coarse[0] + 30.0)
    fine_dms = np.arange(lo, hi + 0.1, 1.0)
    fine = [(float(dm), score_dm(sub, freqs, float(dm))) for dm in fine_dms]
    best = max(fine, key=lambda t: t[1])
    return {
        "best_coarse_dm": best_coarse[0],
        "best_coarse_score": best_coarse[1],
        "best_dm": best[0],
        "best_score": best[1],
        "score_dm_430": score_dm(sub, freqs, 430.0),
        "score_dm_458": score_dm(sub, freqs, 458.0),
        "coarse_top5": sorted(coarse, key=lambda t: t[1], reverse=True)[:5],
        "fine_top10": sorted(fine, key=lambda t: t[1], reverse=True)[:10],
    }


def contiguous_halfmax_width(z: np.ndarray, peak: int):
    target = 0.5 * z[peak]
    left = peak
    while left > 0 and z[left - 1] >= target:
        left -= 1
    right = peak
    while right + 1 < len(z) and z[right + 1] >= target:
        right += 1
    return left, right, (right - left + 1) * DT


def aligned_channel_onpulse(xnorm: np.ndarray, freqs: np.ndarray, dm: float, left: int, right: int):
    ref = float(np.max(freqs))
    delays = KDM_SEC_MHZ2 * dm * (freqs ** -2 - ref ** -2)
    shifts = np.rint(delays / DT).astype(int)
    means = np.full(xnorm.shape[1], np.nan, dtype=float)
    nwin = max(1, right - left + 1)
    for j, sh in enumerate(shifts):
        # Aligned index i corresponds to original i+shift for positive DM delay.
        a = left + sh
        b = right + sh + 1
        a = max(0, a)
        b = min(xnorm.shape[0], b)
        if b > a:
            means[j] = float(np.mean(xnorm[a:b, j]))
    significance = means * np.sqrt(nwin)
    return means, significance, shifts


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    receipt = {
        "schema": "JANUS_TOPA_J1832_CRACO_NUMERIC_MORPHOLOGY_RUN",
        "version": "1.0",
        "hunt_id": "TOPA-TESLA-HUNT-005C3",
        "parent_prereg": "data/tesla-sweep/results/TOPA-HUNT-005C3-J1832-CRACO-NUMERIC-MORPHOLOGY-PREREG-2026-08-23-v1.0.json",
        "status": "STARTED",
        "security": {"allow_pickle": False, "external_pickle_execution": False},
    }
    try:
        with tempfile.TemporaryDirectory(prefix="topa_j1832_morph_") as td:
            path = Path(td) / "SB55237_CRACO_filterbank.npy"
            with requests.get(URL, stream=True, timeout=(20, 120)) as r:
                r.raise_for_status()
                with path.open("wb") as f:
                    for chunk in r.iter_content(1024 * 1024):
                        if chunk:
                            f.write(chunk)
            got = md5_file(path)
            receipt["computed_md5"] = got
            receipt["md5_match"] = got == EXPECTED_MD5
            if got != EXPECTED_MD5:
                raise RuntimeError("checksum mismatch")
            a = np.load(path, mmap_mode="r", allow_pickle=False)
            if a.shape != EXPECTED_SHAPE:
                raise RuntimeError(f"shape mismatch {a.shape} != {EXPECTED_SHAPE}")
            if a.dtype.hasobject:
                raise RuntimeError("object dtype rejected")
            x = np.asarray(a, dtype=np.float32)

        raw_std = np.std(x, axis=0, dtype=np.float64)
        med_std = float(np.median(raw_std))
        lowvar_fraction = float(np.mean(raw_std < 0.05 * med_std)) if med_std > 0 else 0.0

        med, scale = robust_center_scale(x, axis=0)
        xnorm = ((x - med) / scale).astype(np.float32)
        sub = xnorm.reshape(xnorm.shape[0], NSUB, xnorm.shape[1] // NSUB).mean(axis=2)

        orientation_results = {}
        for orientation in ("ascending", "descending"):
            fchan = frequency_vector(x.shape[1], orientation)
            fsub = fchan.reshape(NSUB, x.shape[1] // NSUB).mean(axis=1)
            orientation_results[orientation] = dm_search(sub, fsub)

        best_orientation = max(orientation_results, key=lambda k: orientation_results[k]["best_score"])
        best_dm = float(orientation_results[best_orientation]["best_dm"])
        fchan = frequency_vector(x.shape[1], best_orientation)
        fsub = fchan.reshape(NSUB, x.shape[1] // NSUB).mean(axis=1)
        profile, sub_shifts, sub_delays = shift_earlier_accumulate(sub, fsub, best_dm)
        z = robust_z_1d(profile)
        smooth = signal.savgol_filter(z, 101, 3, mode="interp")
        peak = int(np.argmax(smooth))
        left, right, fwhm_sec = contiguous_halfmax_width(smooth, peak)
        means, chan_sig, chan_shifts = aligned_channel_onpulse(xnorm, fchan, best_dm, left, right)
        valid = np.isfinite(means)
        occ = float(np.mean(chan_sig[valid] > 3.0)) if valid.any() else None
        m = means[valid]
        patch = float(np.std(m) / (abs(np.mean(m)) + EPS)) if m.size else None

        peaks, props = signal.find_peaks(smooth, height=5.0, distance=max(1, int(round(1.0 / DT))))

        receipt.update({
            "status": "PASS_SAFE_NUMERIC_CRACO_MORPHOLOGY_CHARACTERIZATION__CONTROL_COMPARISON_PENDING",
            "array": {
                "shape": list(x.shape),
                "dtype": str(x.dtype),
                "time_resolution_ms": DT * 1000.0,
                "frequency_resolution_MHz": DF_MHZ,
                "reference_center_frequency_MHz": FC_MHZ,
                "window_duration_seconds": x.shape[0] * DT,
                "nominal_bandwidth_MHz": x.shape[1] * DF_MHZ,
                "low_variance_channel_fraction": lowvar_fraction,
            },
            "DM_replay": {
                "orientation_results": orientation_results,
                "best_frequency_orientation": best_orientation,
                "best_DM_pc_cm3": best_dm,
                "published_discovery_dynamic_spectrum_DM_pc_cm3": "458 +/- 14",
                "public_timing_model_DM_pc_cm3": 430.0,
                "important_boundary": "Grid-peak DM is a replay diagnostic; no uncertainty interval is claimed from score width alone."
            },
            "morphology": {
                "dedispersed_broadband_peak_robust_z": float(smooth[peak]),
                "peak_time_seconds_from_array_start": peak * DT,
                "FWHM_seconds": float(fwhm_sec),
                "FWHM_fraction_of_filterbank_window": float(fwhm_sec / (x.shape[0] * DT)),
                "onpulse_frequency_channel_occupancy_fraction_z_gt_3": occ,
                "spectral_patchiness_coefficient": patch,
                "prominent_time_subpeaks_z_gt_5_separated_by_at_least_1s": int(len(peaks)),
                "prominent_subpeak_times_seconds": [float(i * DT) for i in peaks[:50]],
                "classification_ceiling": "NUMERIC_ASTROPHYSICAL_MORPHOLOGY_CHARACTERIZED__CONTROL_COMPARISON_PENDING"
            },
            "raw_telescope_voltage_processing_completed": False,
            "author_exported_filterbank_processing_completed": True,
            "cross_telescope_morphology_replication_completed": False,
            "next_gate": "TOPA-HUNT-005C4-LPT-FRB-PULSAR-RFI-CONTROL-MATRIX",
        })
    except Exception as exc:
        receipt.update({
            "status": "BLOCKED_OR_FAILED_NUMERIC_MORPHOLOGY",
            "error_type": type(exc).__name__,
            "error": str(exc),
            "rule": "FAIL_CLOSED_AND_KEEP_RECEIPT",
        })

    OUT.write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print(json.dumps(receipt, indent=2))
    return 0 if receipt["status"].startswith("PASS_") else 2


if __name__ == "__main__":
    raise SystemExit(main())
