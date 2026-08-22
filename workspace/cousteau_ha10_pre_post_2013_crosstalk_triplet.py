#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import io
import json
import math
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import numpy as np
import requests
from obspy import UTCDateTime, read, read_inventory
from scipy import signal

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cousteau"
PROTOCOL = DATA / "JANUS-ECHO-COUSTEAU-HA10-PRE-POST-2013-CROSSTALK-TRIPLET-PROTOCOL-2026-08-22-v1.0.json"
RESPONSE_SUMMARY = DATA / "JANUS-ECHO-COUSTEAU-HA10-PUBLIC-RESPONSE-EPOCH-AUDIT-RUN-001-SUMMARY-2026-08-22-v1.0.json"
REVERSE_AUDIT = DATA / "JANUS-ECHO-COUSTEAU-HA10-REVERSE-SPIRAL-CALIBRATION-CROSSTALK-AUDIT-2026-08-22-v1.0.json"
FROZEN_119 = DATA / "JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0.json"

EXPECTED_PROTOCOL_BLOB = "bfa9d00082aad01f8ea1246a48e0e72a91188ceb"
EXPECTED_RESPONSE_SUMMARY_BLOB = "4db3760f0ee04d5d832941f5001212e8414978f8"
EXPECTED_REVERSE_AUDIT_BLOB = "775b53d17e4fd2c185ab1b47cf6e5bfe8d5d3ddf"
EXPECTED_FROZEN_119_BLOB = "eb8b48fb7f043160c057f9df6264a781412ed854"
EXPECTED_FROZEN_119_VERDICT = "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE"

FDSN_DATASELECT = "https://service.earthscope.org/fdsnws/dataselect/1/query"
FDSN_STATION = "https://service.earthscope.org/fdsnws/station/1/query"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git_blob_sha1_bytes(data: bytes) -> str:
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def git_blob_sha1_file(path: Path) -> str:
    return git_blob_sha1_bytes(path.read_bytes())


def contains_exact_scalar(value: Any, needle: str) -> bool:
    if isinstance(value, dict):
        return any(contains_exact_scalar(v, needle) for v in value.values())
    if isinstance(value, list):
        return any(contains_exact_scalar(v, needle) for v in value)
    return value == needle


def verify_frozen_inputs() -> dict[str, Any]:
    checks = {
        PROTOCOL: EXPECTED_PROTOCOL_BLOB,
        RESPONSE_SUMMARY: EXPECTED_RESPONSE_SUMMARY_BLOB,
        REVERSE_AUDIT: EXPECTED_REVERSE_AUDIT_BLOB,
        FROZEN_119: EXPECTED_FROZEN_119_BLOB,
    }
    for path, expected in checks.items():
        actual = git_blob_sha1_file(path)
        if actual != expected:
            raise RuntimeError(f"FROZEN_BLOB_MISMATCH:{path.name}:{actual}:{expected}")
    frozen = json.loads(FROZEN_119.read_text(encoding="utf-8"))
    if not contains_exact_scalar(frozen, EXPECTED_FROZEN_119_VERDICT):
        raise RuntimeError("FROZEN_119_VERDICT_NOT_BOUND_IN_EXACT_BYTES")
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if protocol.get("gate_id") != "COUSTEAU_HA10_PRE_POST_2013_CROSSTALK_TRIPLET_GATE_V1":
        raise RuntimeError("WRONG_PROTOCOL_GATE_ID")
    return protocol


def request_bytes(url: str, params: dict[str, Any], *, tries: int = 2, timeout: int = 30) -> tuple[bytes, str, str | None]:
    last: Exception | None = None
    headers = {"User-Agent": "JANUS-Cousteau-HA10-Triplet-Gate/1.0"}
    for attempt in range(tries):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=timeout)
            if r.status_code == 200 and r.content:
                return r.content, r.url, r.headers.get("content-type")
            last = RuntimeError(f"HTTP_{r.status_code}:{r.text[:160]}")
        except Exception as exc:
            last = exc
        if attempt + 1 < tries:
            time.sleep(0.75 * (attempt + 1))
    raise last or RuntimeError("DOWNLOAD_FAILED")


def fetch_inventory(station: str) -> tuple[Any | None, dict[str, Any]]:
    params = {
        "net": "IM",
        "sta": station,
        "loc": "--",
        "cha": "EDH",
        "starttime": "2013-06-01T00:00:00Z",
        "endtime": "2013-10-01T00:00:00Z",
        "level": "response",
        "format": "xml",
        "nodata": 404,
    }
    try:
        payload, final_url, content_type = request_bytes(FDSN_STATION, params, tries=3, timeout=45)
        inv = read_inventory(io.BytesIO(payload))
        return inv, {
            "status": "FETCHED",
            "url": final_url,
            "bytes": len(payload),
            "sha256": sha256_bytes(payload),
            "content_type": content_type,
        }
    except Exception as exc:
        return None, {"status": "UNAVAILABLE", "error": f"{type(exc).__name__}:{exc}"}


def parse_trace(payload: bytes, station: str, expected_fs: float, duration_s: float) -> tuple[Any, dict[str, Any]]:
    st = read(io.BytesIO(payload), format="MSEED")
    st = st.select(network="IM", station=station, location="", channel="EDH")
    if len(st) == 0:
        raise RuntimeError("NO_MATCHING_TRACE")
    gaps = st.get_gaps()
    nonzero_gaps = [g for g in gaps if abs(float(g[6])) > (0.25 / expected_fs)]
    if nonzero_gaps:
        raise RuntimeError(f"GAPS_OR_OVERLAPS_PRESENT:{len(nonzero_gaps)}")
    st.merge(method=0)
    if len(st) != 1:
        raise RuntimeError(f"TRACE_COUNT_AFTER_CONTIGUOUS_MERGE:{len(st)}")
    tr = st[0]
    fs = float(tr.stats.sampling_rate)
    if abs(fs - expected_fs) > 1e-9:
        raise RuntimeError(f"UNEXPECTED_SAMPLE_RATE:{fs}")
    minimum_npts = int(round(duration_s * expected_fs)) - 2
    if int(tr.stats.npts) < minimum_npts:
        raise RuntimeError(f"INCOMPLETE_WINDOW_NPTS:{tr.stats.npts}:{minimum_npts}")
    data = np.asarray(tr.data, dtype=np.float64)
    if not np.all(np.isfinite(data)):
        raise RuntimeError("NONFINITE_SAMPLES")
    return tr, {
        "npts": int(tr.stats.npts),
        "sampling_rate_hz": fs,
        "starttime": str(tr.stats.starttime),
        "endtime": str(tr.stats.endtime),
        "dtype": str(tr.data.dtype),
    }


def fetch_waveform(station: str, start_utc: str, end_utc: str, expected_fs: float, duration_s: float) -> dict[str, Any]:
    params = {
        "net": "IM",
        "sta": station,
        "loc": "--",
        "cha": "EDH",
        "starttime": start_utc,
        "endtime": end_utc,
        "nodata": 404,
        "format": "miniseed",
    }
    try:
        payload, final_url, content_type = request_bytes(FDSN_DATASELECT, params)
        tr, trace_meta = parse_trace(payload, station, expected_fs, duration_s)
        return {
            "status": "FETCHED",
            "station": station,
            "url": final_url,
            "bytes": len(payload),
            "sha256": sha256_bytes(payload),
            "content_type": content_type,
            "trace_meta": trace_meta,
            "trace": tr,
        }
    except Exception as exc:
        return {
            "status": "BLOCKED",
            "station": station,
            "error": f"{type(exc).__name__}:{exc}",
        }


def align_six_traces(results: dict[str, dict[str, Any]], expected_fs: float, duration_s: float) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    traces = {s: row["trace"] for s, row in results.items()}
    starts = {s: float(tr.stats.starttime.timestamp) for s, tr in traces.items()}
    ends = {s: float(tr.stats.endtime.timestamp) for s, tr in traces.items()}
    common_start = max(starts.values())
    common_end = min(ends.values())
    if common_end <= common_start:
        raise RuntimeError("NO_COMMON_TIME_INTERVAL")
    offsets: dict[str, int] = {}
    for station, start in starts.items():
        raw = (common_start - start) * expected_fs
        rounded = int(round(raw))
        if abs(raw - rounded) > 0.10:
            raise RuntimeError(f"NON_SAMPLE_ALIGNED_START:{station}:{raw}")
        offsets[station] = rounded
    n_common = int(math.floor((common_end - common_start) * expected_fs + 1e-7)) + 1
    minimum_npts = int(round(duration_s * expected_fs)) - 4
    if n_common < minimum_npts:
        raise RuntimeError(f"COMMON_WINDOW_TOO_SHORT:{n_common}:{minimum_npts}")
    arrays: dict[str, np.ndarray] = {}
    for station, tr in traces.items():
        start_idx = offsets[station]
        arr = np.asarray(tr.data, dtype=np.float64)[start_idx : start_idx + n_common]
        if len(arr) != n_common:
            raise RuntimeError(f"ALIGNMENT_SLICE_SHORT:{station}:{len(arr)}:{n_common}")
        arrays[station] = arr
    return arrays, {
        "common_start_utc": str(UTCDateTime(common_start)),
        "common_end_utc": str(UTCDateTime(common_start + (n_common - 1) / expected_fs)),
        "common_npts": n_common,
        "sample_offsets": offsets,
    }


def welch_psd(data: np.ndarray, fs: float, nperseg: int, noverlap: int) -> tuple[np.ndarray, np.ndarray]:
    return signal.welch(
        data,
        fs=fs,
        window="hann",
        nperseg=nperseg,
        noverlap=noverlap,
        detrend="constant",
        return_onesided=True,
        scaling="density",
    )


def band_integrated_power(freqs: np.ndarray, psd: np.ndarray, low: float, high: float) -> float:
    mask = (freqs >= low) & (freqs <= high) & np.isfinite(psd)
    if int(np.count_nonzero(mask)) < 2:
        raise RuntimeError("INSUFFICIENT_BAND_BINS")
    value = float(np.trapezoid(psd[mask], freqs[mask]))
    if not (math.isfinite(value) and value > 0):
        raise RuntimeError(f"INVALID_INTEGRATED_POWER:{value}")
    return value


def pair_mean_coherence(a: np.ndarray, b: np.ndarray, fs: float, nperseg: int, noverlap: int, low: float, high: float) -> tuple[float, dict[str, float]]:
    freqs, coh = signal.coherence(
        a,
        b,
        fs=fs,
        window="hann",
        nperseg=nperseg,
        noverlap=noverlap,
        detrend="constant",
    )
    mask = (freqs >= low) & (freqs <= high) & np.isfinite(coh)
    if not np.any(mask):
        raise RuntimeError("NO_COHERENCE_BINS")
    primary = float(np.mean(coh[mask]))
    subbands = {}
    for lo, hi in ((10.0, 30.0), (30.0, 55.0), (55.0, 80.0)):
        sm = (freqs >= lo) & (freqs <= hi) & np.isfinite(coh)
        subbands[f"{lo:g}_{hi:g}"] = float(np.mean(coh[sm])) if np.any(sm) else float("nan")
    return primary, subbands


def select_response(inv: Any, station: str, when_utc: str) -> Any | None:
    if inv is None:
        return None
    t = UTCDateTime(when_utc)
    candidates = []
    for network in inv:
        if network.code != "IM":
            continue
        for sta in network:
            if sta.code != station:
                continue
            for cha in sta:
                if cha.code != "EDH" or (cha.location_code or "") != "":
                    continue
                if cha.start_date is not None and t < cha.start_date:
                    continue
                if cha.end_date is not None and t > cha.end_date:
                    continue
                if cha.response is not None:
                    candidates.append(cha.response)
    return candidates[0] if len(candidates) == 1 else None


def corrected_power_from_raw_psd(response: Any, freqs: np.ndarray, raw_psd: np.ndarray, fs: float, nperseg: int, low: float, high: float) -> float:
    values, rfreqs = response.get_evalresp_response(t_samp=1.0 / fs, nfft=nperseg, output="DEF")
    values = np.asarray(values, dtype=complex)
    rfreqs = np.asarray(rfreqs, dtype=float)
    if len(rfreqs) != len(freqs) or not np.allclose(rfreqs, freqs, rtol=0, atol=1e-10):
        raise RuntimeError("EVALRESP_FREQUENCY_GRID_MISMATCH")
    denom = np.abs(values) ** 2
    corrected = np.full_like(raw_psd, np.nan, dtype=float)
    valid = np.isfinite(raw_psd) & np.isfinite(denom) & (denom > 0)
    corrected[valid] = raw_psd[valid] / denom[valid]
    return band_integrated_power(freqs, corrected, low, high)


def analyze_window(arrays: dict[str, np.ndarray], inventories: dict[str, Any | None], start_utc: str, protocol: dict[str, Any]) -> dict[str, Any]:
    fs = float(protocol["channels"]["expected_sample_rate_hz"])
    nperseg = int(protocol["spectral_contract"]["welch_nperseg"])
    noverlap = int(protocol["spectral_contract"]["welch_noverlap"])
    low, high = [float(x) for x in protocol["spectral_contract"]["band_hz"]]

    raw_psd: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    raw_power: dict[str, float] = {}
    corrected_power: dict[str, float | None] = {}
    corrected_errors: dict[str, str] = {}
    for station, data in arrays.items():
        freqs, psd = welch_psd(data, fs, nperseg, noverlap)
        raw_psd[station] = (freqs, psd)
        raw_power[station] = band_integrated_power(freqs, psd, low, high)
        response = select_response(inventories.get(station), station, start_utc)
        if response is None:
            corrected_power[station] = None
            corrected_errors[station] = "NO_UNIQUE_PUBLIC_RESPONSE"
        else:
            try:
                corrected_power[station] = corrected_power_from_raw_psd(response, freqs, psd, fs, nperseg, low, high)
            except Exception as exc:
                corrected_power[station] = None
                corrected_errors[station] = f"{type(exc).__name__}:{exc}"

    pairs = [
        ("H10S1", "H10S2"), ("H10S1", "H10S3"), ("H10S2", "H10S3"),
        ("H10N1", "H10N2"), ("H10N1", "H10N3"), ("H10N2", "H10N3"),
    ]
    pair_coherence: dict[str, float] = {}
    pair_subbands: dict[str, dict[str, float]] = {}
    for a, b in pairs:
        value, sub = pair_mean_coherence(arrays[a], arrays[b], fs, nperseg, noverlap, low, high)
        key = f"{a}-{b}"
        pair_coherence[key] = value
        pair_subbands[key] = sub

    south_mean = float(np.mean([pair_coherence["H10S1-H10S2"], pair_coherence["H10S1-H10S3"]]))
    north_mean = float(np.mean([pair_coherence["H10N1-H10N2"], pair_coherence["H10N1-H10N3"]]))
    contrast = south_mean - north_mean
    raw_s2_n1_db = float(10.0 * math.log10(raw_power["H10S2"] / raw_power["H10N1"]))
    corr_s2_n1_db = None
    if corrected_power["H10S2"] is not None and corrected_power["H10N1"] is not None:
        corr_s2_n1_db = float(10.0 * math.log10(float(corrected_power["H10S2"]) / float(corrected_power["H10N1"])))

    return {
        "pair_coherence_10_80": pair_coherence,
        "pair_coherence_fixed_subbands": pair_subbands,
        "south_fault_source_pair_mean_coherence": south_mean,
        "north_analog_pair_mean_coherence": north_mean,
        "south_minus_north_coherence_contrast": contrast,
        "raw_integrated_power_10_80": raw_power,
        "public_response_corrected_integrated_power_10_80": corrected_power,
        "corrected_power_errors": corrected_errors,
        "raw_s2_minus_n1_power_ratio_db": raw_s2_n1_db,
        "corrected_s2_minus_n1_power_ratio_db": corr_s2_n1_db,
    }


def one_sided_permutation_p(pre: np.ndarray, post: np.ndarray, permutations: int, seed: int) -> tuple[float, float]:
    pre = np.asarray(pre, dtype=float)
    post = np.asarray(post, dtype=float)
    if len(pre) == 0 or len(post) == 0:
        raise ValueError("EMPTY_PERMUTATION_GROUP")
    observed = float(np.median(post) - np.median(pre))
    pooled = np.concatenate([pre, post])
    n_pre = len(pre)
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(permutations):
        perm = rng.permutation(pooled)
        stat = float(np.median(perm[n_pre:]) - np.median(perm[:n_pre]))
        if stat >= observed - 1e-15:
            exceed += 1
    p = float((exceed + 1) / (permutations + 1))
    return observed, p


def median_pair(rows: list[dict[str, Any]], key: str) -> float:
    return float(np.median([float(r["analysis"]["pair_coherence_10_80"][key]) for r in rows]))


def classify_preexistence(pre_median_db: float, post_median_db: float, protocol: dict[str, Any]) -> str:
    cfg = protocol["preexistence_subgate"]
    floor = float(cfg["large_offset_floor_db"])
    stable = float(cfg["maximum_pre_post_median_change_for_stability_db"])
    if pre_median_db >= floor and abs(post_median_db - pre_median_db) <= stable:
        return "LARGE_S2_N1_OFFSET_PREEXISTS_AND_IS_STABLE_ACROSS_FAULT"
    if pre_median_db < floor and post_median_db >= floor:
        return "LARGE_S2_N1_OFFSET_EMERGES_AFTER_FAULT"
    return "S2_N1_OFFSET_PRE_POST_MIXED"


def summarize_complete_rows(pre_rows: list[dict[str, Any]], post_rows: list[dict[str, Any]], protocol: dict[str, Any], *, permutations_override: int | None = None) -> dict[str, Any]:
    minimum = int(protocol["window_freeze"]["minimum_complete_six_channel_windows_per_epoch"])
    if len(pre_rows) < minimum or len(post_rows) < minimum:
        return {
            "verdict": "BLOCKED_PRE_POST_2013_TRIPLET_DATA",
            "pre_complete_windows": len(pre_rows),
            "post_complete_windows": len(post_rows),
            "minimum_required_each_epoch": minimum,
            "preexistence_classification": "BLOCKED_PRE_POST_2013_TRIPLET_DATA",
        }

    pre_contrast = np.asarray([r["analysis"]["south_minus_north_coherence_contrast"] for r in pre_rows], dtype=float)
    post_contrast = np.asarray([r["analysis"]["south_minus_north_coherence_contrast"] for r in post_rows], dtype=float)
    pcfg = protocol["primary_linear_crosstalk_statistic"]
    perm_cfg = pcfg["permutation_test"]
    nperm = int(permutations_override if permutations_override is not None else perm_cfg["permutations"])
    delta, p = one_sided_permutation_p(pre_contrast, post_contrast, nperm, int(perm_cfg["rng_seed"]))

    pair_medians = {}
    for pair in ("H10S1-H10S2", "H10S1-H10S3", "H10S2-H10S3", "H10N1-H10N2", "H10N1-H10N3", "H10N2-H10N3"):
        pre_m = median_pair(pre_rows, pair)
        post_m = median_pair(post_rows, pair)
        pair_medians[pair] = {"pre": pre_m, "post": post_m, "post_minus_pre": post_m - pre_m}

    effects = pcfg["effect_thresholds"]
    signature_checks = {
        "primary_delta_ge_threshold": delta >= float(effects["minimum_primary_delta_msc"]),
        "permutation_p_le_threshold": p <= float(perm_cfg["p_threshold"]),
        "s1_s2_increase_ge_threshold": pair_medians["H10S1-H10S2"]["post_minus_pre"] >= float(effects["minimum_individual_s1_s2_post_minus_pre_median_msc"]),
        "s1_s3_increase_ge_threshold": pair_medians["H10S1-H10S3"]["post_minus_pre"] >= float(effects["minimum_individual_s1_s3_post_minus_pre_median_msc"]),
    }
    verdict = (
        "LINEAR_COHERENT_SOUTH_TRIPLET_CROSSTALK_SIGNATURE_DETECTED"
        if all(signature_checks.values())
        else "NO_PREREGISTERED_LINEAR_COHERENT_CROSSTALK_SIGNATURE"
    )

    pre_ratio = float(np.median([r["analysis"]["raw_s2_minus_n1_power_ratio_db"] for r in pre_rows]))
    post_ratio = float(np.median([r["analysis"]["raw_s2_minus_n1_power_ratio_db"] for r in post_rows]))
    pre_class = classify_preexistence(pre_ratio, post_ratio, protocol)
    corrected_pre = [r["analysis"].get("corrected_s2_minus_n1_power_ratio_db") for r in pre_rows]
    corrected_post = [r["analysis"].get("corrected_s2_minus_n1_power_ratio_db") for r in post_rows]
    corrected_pre = [float(x) for x in corrected_pre if x is not None and math.isfinite(float(x))]
    corrected_post = [float(x) for x in corrected_post if x is not None and math.isfinite(float(x))]

    station_power = {}
    for station in ("H10N1", "H10N2", "H10N3", "H10S1", "H10S2", "H10S3"):
        raw_pre = float(np.median([r["analysis"]["raw_integrated_power_10_80"][station] for r in pre_rows]))
        raw_post = float(np.median([r["analysis"]["raw_integrated_power_10_80"][station] for r in post_rows]))
        cp_pre = [r["analysis"]["public_response_corrected_integrated_power_10_80"].get(station) for r in pre_rows]
        cp_post = [r["analysis"]["public_response_corrected_integrated_power_10_80"].get(station) for r in post_rows]
        cp_pre = [float(x) for x in cp_pre if x is not None and math.isfinite(float(x))]
        cp_post = [float(x) for x in cp_post if x is not None and math.isfinite(float(x))]
        station_power[station] = {
            "raw_pre_median": raw_pre,
            "raw_post_median": raw_post,
            "raw_post_minus_pre_db": float(10.0 * math.log10(raw_post / raw_pre)),
            "corrected_pre_median": float(np.median(cp_pre)) if cp_pre else None,
            "corrected_post_median": float(np.median(cp_post)) if cp_post else None,
        }
        if cp_pre and cp_post:
            station_power[station]["corrected_post_minus_pre_db"] = float(10.0 * math.log10(np.median(cp_post) / np.median(cp_pre)))
        else:
            station_power[station]["corrected_post_minus_pre_db"] = None

    return {
        "verdict": verdict,
        "pre_complete_windows": len(pre_rows),
        "post_complete_windows": len(post_rows),
        "minimum_required_each_epoch": minimum,
        "primary": {
            "pre_median_south_minus_north_coherence_contrast": float(np.median(pre_contrast)),
            "post_median_south_minus_north_coherence_contrast": float(np.median(post_contrast)),
            "post_minus_pre_contrast_delta_msc": delta,
            "one_sided_permutation_p": p,
            "permutations": nperm,
            "rng_seed": int(perm_cfg["rng_seed"]),
            "checks": signature_checks,
            "pair_medians": pair_medians,
        },
        "preexistence": {
            "raw_pre_median_s2_minus_n1_db": pre_ratio,
            "raw_post_median_s2_minus_n1_db": post_ratio,
            "raw_post_minus_pre_median_change_db": post_ratio - pre_ratio,
            "classification": pre_class,
            "corrected_pre_median_s2_minus_n1_db": float(np.median(corrected_pre)) if corrected_pre else None,
            "corrected_post_median_s2_minus_n1_db": float(np.median(corrected_post)) if corrected_post else None,
        },
        "station_power_epoch_medians": station_power,
    }


def make_window_times(date_str: str, start_hms: str, duration_s: int) -> tuple[str, str]:
    start = datetime.fromisoformat(f"{date_str}T{start_hms.replace('Z', '+00:00')}")
    end = start + timedelta(seconds=duration_s)
    return start.isoformat().replace("+00:00", "Z"), end.isoformat().replace("+00:00", "Z")


def execute_epoch(epoch_name: str, dates: list[str], stations: list[str], inventories: dict[str, Any | None], protocol: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    cfg = protocol["window_freeze"]
    start_hms = cfg["window_start_time_utc"]
    duration_s = int(cfg["window_duration_s"])
    fs = float(protocol["channels"]["expected_sample_rate_hz"])
    complete: list[dict[str, Any]] = []
    blocked: list[dict[str, Any]] = []
    for date_str in dates:
        start_utc, end_utc = make_window_times(date_str, start_hms, duration_s)
        with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
            futures = {
                station: pool.submit(fetch_waveform, station, start_utc, end_utc, fs, duration_s)
                for station in stations
            }
            fetched = {station: future.result() for station, future in futures.items()}
        public_meta = {
            station: {k: v for k, v in row.items() if k != "trace"}
            for station, row in fetched.items()
        }
        failures = {s: r for s, r in public_meta.items() if r.get("status") != "FETCHED"}
        if failures:
            blocked.append({
                "epoch": epoch_name,
                "date": date_str,
                "start_utc": start_utc,
                "status": "BLOCKED_STATION_FETCH_OR_GAP",
                "stations": public_meta,
            })
            continue
        try:
            arrays, alignment = align_six_traces(fetched, fs, duration_s)
            analysis = analyze_window(arrays, inventories, start_utc, protocol)
        except Exception as exc:
            blocked.append({
                "epoch": epoch_name,
                "date": date_str,
                "start_utc": start_utc,
                "status": "BLOCKED_ALIGNMENT_OR_ANALYSIS",
                "error": f"{type(exc).__name__}:{exc}",
                "stations": public_meta,
            })
            continue
        complete.append({
            "epoch": epoch_name,
            "date": date_str,
            "start_utc": start_utc,
            "end_utc": end_utc,
            "status": "COMPLETE_SIX_CHANNEL_WINDOW",
            "stations": public_meta,
            "alignment": alignment,
            "analysis": analysis,
        })
    return complete, blocked


def run(output: Path) -> dict[str, Any]:
    resolved = output.resolve()
    if DATA.resolve() == resolved or DATA.resolve() in resolved.parents:
        raise RuntimeError("CANONICAL_DATA_OUTPUT_FORBIDDEN_USE_EPHEMERAL_ARTIFACT_PATH")
    protocol = verify_frozen_inputs()
    north = list(protocol["channels"]["north_triplet"])
    south = list(protocol["channels"]["south_triplet"])
    stations = north + south

    inventories: dict[str, Any | None] = {}
    response_transport: dict[str, Any] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        futures = {s: pool.submit(fetch_inventory, s) for s in stations}
        for station, future in futures.items():
            inv, meta = future.result()
            inventories[station] = inv
            response_transport[station] = meta

    pre_rows, pre_blocked = execute_epoch(
        "PRE_FAULT",
        list(protocol["window_freeze"]["pre_fault_dates_utc"]),
        stations,
        inventories,
        protocol,
    )
    post_rows, post_blocked = execute_epoch(
        "POST_FAULT",
        list(protocol["window_freeze"]["post_fault_dates_utc"]),
        stations,
        inventories,
        protocol,
    )
    summary = summarize_complete_rows(pre_rows, post_rows, protocol)

    receipt = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-HA10-PRE-POST-2013-CROSSTALK-TRIPLET-RUN-001",
        "status": "EXECUTED_EPHEMERAL_RECEIPT",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_path": str(PROTOCOL.relative_to(ROOT)).replace("\\", "/"),
        "protocol_git_blob_sha1": EXPECTED_PROTOCOL_BLOB,
        "public_response_epoch_summary_git_blob_sha1": EXPECTED_RESPONSE_SUMMARY_BLOB,
        "reverse_spiral_audit_git_blob_sha1": EXPECTED_REVERSE_AUDIT_BLOB,
        "frozen_119hz_result_git_blob_sha1": EXPECTED_FROZEN_119_BLOB,
        "frozen_119hz_verdict": EXPECTED_FROZEN_119_VERDICT,
        "janus_selected_gate": protocol["selected_by_janus"],
        "response_metadata_transport": response_transport,
        "pre_windows": pre_rows,
        "post_windows": post_rows,
        "blocked_windows": pre_blocked + post_blocked,
        "summary": summary,
        "verdict": summary["verdict"],
        "preexistence_classification": summary.get("preexistence", {}).get(
            "classification", summary.get("preexistence_classification")
        ),
        "authority": {
            "raw_counts_are_primary_for_crosstalk_gate": True,
            "public_response_corrected_power_is_secondary": True,
            "bathymetry_lane_remains_independent": True,
            "frozen_119hz_negative_result_immutable": True,
            "authority_delta_for_119hz": 0,
            "target_identity": "UNCONFIRMED",
            "source_writeback": False,
        },
        "next_step": "FREEZE_RESULT_AND_SUBMIT_EXACT_RECEIPT_TO_JANUS_DEMIURGE_REVERSE_COUNCIL",
        "claim_ceiling": protocol["claim_ceiling"],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise RuntimeError("OUTPUT_ALREADY_EXISTS_APPEND_ONLY_REQUIRED")
    output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    receipt = run(Path(args.output))
    print(json.dumps({
        "verdict": receipt["verdict"],
        "preexistence_classification": receipt["preexistence_classification"],
        "pre_complete": receipt["summary"]["pre_complete_windows"],
        "post_complete": receipt["summary"]["post_complete_windows"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
