#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import time
from collections import defaultdict
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import requests
from obspy import UTCDateTime, read, read_inventory
from scipy import signal

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cousteau"
PROTOCOL = DATA / "JANUS-ECHO-COUSTEAU-HA10-FIRST-AVAILABLE-SIX-CHANNEL-BASELINE-PROTOCOL-2026-08-22-v1.0.json"
REFERENCE = DATA / "JANUS-ECHO-COUSTEAU-HA10-N1-S2-TEMPORAL-SCALE-STABILITY-DIAGNOSTIC-RUN-001-SUMMARY-2026-08-22-v1.0.json"
FROZEN_119 = DATA / "JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0.json"

EXPECTED_PROTOCOL_BLOB = "f040c5df715fb1cc3dca3d478be9334c5c09d271"
EXPECTED_REFERENCE_BLOB = "45f2901adde866df610a563f27c50177d32be414"
EXPECTED_FROZEN_119_BLOB = "eb8b48fb7f043160c057f9df6264a781412ed854"
EXPECTED_FROZEN_119_VERDICT = "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE"

DATASELECT = "https://service.earthscope.org/fdsnws/dataselect/1/query"
STATION = "https://service.earthscope.org/fdsnws/station/1/query"
UA = "JANUS-Cousteau-First-Six-Channel-Baseline/1.0"
HTTP_TIMEOUT_S = 60
HTTP_ATTEMPTS = 3


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git_blob_sha1_file(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def verify_frozen_inputs() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    for path, expected in {
        PROTOCOL: EXPECTED_PROTOCOL_BLOB,
        REFERENCE: EXPECTED_REFERENCE_BLOB,
        FROZEN_119: EXPECTED_FROZEN_119_BLOB,
    }.items():
        actual = git_blob_sha1_file(path)
        if actual != expected:
            raise RuntimeError(f"FROZEN_BLOB_DRIFT:{path.name}:{actual}:{expected}")
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    frozen = json.loads(FROZEN_119.read_text(encoding="utf-8"))
    if protocol.get("status") != "PREREGISTERED_BEFORE_AVAILABILITY_PROBES_OR_NEW_WAVEFORM_VALUE_INSPECTION":
        raise RuntimeError("PROTOCOL_NOT_PREREGISTERED")
    if protocol.get("gate_id") != "COUSTEAU_HA10_FIRST_AVAILABLE_SIX_CHANNEL_BASELINE_V1":
        raise RuntimeError("GATE_ID_DRIFT")
    if protocol.get("authority", {}).get("authority_delta_for_119hz") != 0:
        raise RuntimeError("119_AUTHORITY_DRIFT")
    if protocol.get("spectral_contract", {}).get("119hz_excluded") is not True:
        raise RuntimeError("119_EXCLUSION_DRIFT")
    if frozen.get("summary", {}).get("verdict") != EXPECTED_FROZEN_119_VERDICT:
        raise RuntimeError("FROZEN_119_VERDICT_DRIFT")
    return protocol, reference, frozen


def ensure_noncanonical_output(path: Path) -> None:
    try:
        path.resolve().relative_to(DATA.resolve())
    except ValueError:
        return
    raise RuntimeError("CANONICAL_DATA_WRITE_FORBIDDEN_BY_FIRST_AVAILABLE_BASELINE_GATE")


def parse_channel_id(cid: str) -> tuple[str, str, str, str]:
    parts = cid.split(".")
    if len(parts) != 4:
        raise ValueError(f"BAD_CHANNEL_ID:{cid}")
    return parts[0], parts[1], parts[2], parts[3]


def request_post(session: requests.Session, body: str) -> tuple[bytes | None, dict[str, Any]]:
    last: Exception | None = None
    for attempt in range(HTTP_ATTEMPTS):
        try:
            response = session.post(
                DATASELECT,
                data=body.encode("utf-8"),
                headers={"Content-Type": "text/plain", "User-Agent": UA},
                timeout=HTTP_TIMEOUT_S,
            )
            meta = {
                "url": response.url,
                "http_status": response.status_code,
                "content_type": response.headers.get("content-type"),
                "attempt": attempt + 1,
            }
            if response.status_code == 200 and response.content:
                payload = response.content
                meta.update({"bytes": len(payload), "sha256": sha256_bytes(payload)})
                return payload, meta
            if response.status_code in {204, 404}:
                meta.update({"bytes": len(response.content), "sha256": sha256_bytes(response.content)})
                return None, meta
            last = RuntimeError(f"HTTP_{response.status_code}:{response.text[:240]}")
        except Exception as exc:
            last = exc
        if attempt + 1 < HTTP_ATTEMPTS:
            time.sleep(1.0 * (attempt + 1))
    raise last or RuntimeError("DATASELECT_POST_FAILED")


def selection_body(stations: Iterable[str], windows: Iterable[tuple[datetime, datetime]]) -> str:
    lines = ["quality=M", "format=miniseed"]
    for start, end in windows:
        st = start.strftime("%Y-%m-%dT%H:%M:%S")
        en = end.strftime("%Y-%m-%dT%H:%M:%S")
        for station in stations:
            lines.append(f"IM {station} -- EDH {st} {en}")
    return "\n".join(lines) + "\n"


def intervals_from_headonly(payload: bytes) -> dict[tuple[str, str], list[tuple[float, float, float, int]]]:
    """Return only transport/header information. Sample amplitudes are never accessed."""
    stream = read(io.BytesIO(payload), format="MSEED", headonly=True)
    out: dict[tuple[str, str], list[tuple[float, float, float, int]]] = defaultdict(list)
    for tr in stream:
        station = str(tr.stats.station)
        day = str(tr.stats.starttime.date)
        fs = float(tr.stats.sampling_rate)
        start = float(tr.stats.starttime.timestamp)
        # endtime is timestamp of final sample; add one sample interval for half-open coverage.
        end_exclusive = float(tr.stats.endtime.timestamp) + (1.0 / fs if fs > 0 else 0.0)
        out[(day, station)].append((start, end_exclusive, fs, int(tr.stats.npts)))
    return out


def merge_intervals(intervals: list[tuple[float, float, float, int]], tolerance_s: float) -> list[tuple[float, float]]:
    spans = sorted((float(a), float(b)) for a, b, _fs, _n in intervals if b > a)
    if not spans:
        return []
    merged: list[list[float]] = [[spans[0][0], spans[0][1]]]
    for start, end in spans[1:]:
        if start <= merged[-1][1] + tolerance_s:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(a, b) for a, b in merged]


def header_covers_window(
    intervals: list[tuple[float, float, float, int]],
    start: datetime,
    end: datetime,
    expected_fs: float = 250.0,
) -> bool:
    if not intervals:
        return False
    fs_values = [row[2] for row in intervals]
    if any(abs(fs - expected_fs) > 1e-6 for fs in fs_values):
        return False
    tolerance = 1.5 / expected_fs
    merged = merge_intervals(intervals, tolerance)
    target_start = start.timestamp()
    target_end = end.timestamp()
    for a, b in merged:
        if a <= target_start + tolerance and b >= target_end - tolerance:
            return True
    return False


def daily_probe_scan(
    session: requests.Session,
    stations: list[str],
    scan_start: datetime,
    scan_end: datetime,
    clock: tuple[int, int, int],
    duration_s: int,
) -> tuple[date | None, list[dict[str, Any]]]:
    current_year = scan_start.year
    audit: list[dict[str, Any]] = []
    while current_year <= (scan_end - timedelta(seconds=1)).year:
        year_start = max(scan_start, datetime(current_year, 1, 1, tzinfo=timezone.utc))
        year_end = min(scan_end, datetime(current_year + 1, 1, 1, tzinfo=timezone.utc))
        windows: list[tuple[datetime, datetime]] = []
        d = year_start.date()
        last_date = (year_end - timedelta(seconds=1)).date()
        while d <= last_date:
            start = datetime(d.year, d.month, d.day, clock[0], clock[1], clock[2], tzinfo=timezone.utc)
            if scan_start <= start < scan_end:
                windows.append((start, start + timedelta(seconds=duration_s)))
            d += timedelta(days=1)
        payload, meta = request_post(session, selection_body(stations, windows))
        row: dict[str, Any] = {
            "year": current_year,
            "date_count": len(windows),
            "request": meta,
            "all_six_probe_dates": [],
        }
        if payload:
            headers = intervals_from_headonly(payload)
            station_sets: dict[str, set[str]] = defaultdict(set)
            for (day, station), intervals in headers.items():
                # One-second probe only requires returned header at that date/station.
                if intervals:
                    station_sets[day].add(station)
            all_six = sorted(day for day, ss in station_sets.items() if set(stations).issubset(ss))
            row["all_six_probe_dates"] = all_six
            row["header_trace_keys"] = len(headers)
            if all_six:
                audit.append(row)
                return date.fromisoformat(all_six[0]), audit
        audit.append(row)
        current_year += 1
    return None, audit


def discover_full_windows(
    session: requests.Session,
    stations: list[str],
    earliest: date,
    clock: tuple[int, int, int],
    duration_s: int,
    target_count: int,
    max_days: int,
) -> tuple[list[str], list[dict[str, Any]]]:
    selected: list[str] = []
    audit: list[dict[str, Any]] = []
    for offset in range(max_days):
        d = earliest + timedelta(days=offset)
        start = datetime(d.year, d.month, d.day, clock[0], clock[1], clock[2], tzinfo=timezone.utc)
        end = start + timedelta(seconds=duration_s)
        payload, meta = request_post(session, selection_body(stations, [(start, end)]))
        complete: dict[str, bool] = {station: False for station in stations}
        if payload:
            headers = intervals_from_headonly(payload)
            for station in stations:
                complete[station] = header_covers_window(headers.get((d.isoformat(), station), []), start, end)
        all_complete = all(complete.values())
        audit.append({
            "date": d.isoformat(),
            "start_utc": start.isoformat().replace("+00:00", "Z"),
            "end_utc": end.isoformat().replace("+00:00", "Z"),
            "request": meta,
            "complete_by_station": complete,
            "complete_all_six": all_complete,
        })
        if all_complete:
            selected.append(d.isoformat())
            if len(selected) >= target_count:
                break
    return selected, audit


def fetch_inventory(session: requests.Session, station: str, start: str, end: str) -> tuple[Any | None, dict[str, Any]]:
    params = {
        "net": "IM", "sta": station, "loc": "--", "cha": "EDH",
        "starttime": start, "endtime": end, "level": "response",
        "format": "xml", "nodata": 404,
    }
    last: Exception | None = None
    for attempt in range(HTTP_ATTEMPTS):
        try:
            r = session.get(STATION, params=params, timeout=HTTP_TIMEOUT_S, headers={"User-Agent": UA})
            if r.status_code == 200 and r.content:
                payload = r.content
                return read_inventory(io.BytesIO(payload)), {
                    "url": r.url, "http_status": r.status_code, "bytes": len(payload),
                    "sha256": sha256_bytes(payload), "content_type": r.headers.get("content-type"),
                    "attempt": attempt + 1,
                }
            last = RuntimeError(f"HTTP_{r.status_code}:{r.text[:240]}")
        except Exception as exc:
            last = exc
        if attempt + 1 < HTTP_ATTEMPTS:
            time.sleep(1.0 * (attempt + 1))
    return None, {"status": "UNAVAILABLE", "error": f"{type(last).__name__}:{last}"}


def parse_full_six(payload: bytes, stations: list[str], start: datetime, end: datetime, expected_fs: float = 250.0) -> tuple[dict[str, Any], dict[str, Any]]:
    stream = read(io.BytesIO(payload), format="MSEED")
    traces: dict[str, Any] = {}
    meta: dict[str, Any] = {}
    for station in stations:
        ss = stream.select(network="IM", station=station, location="", channel="EDH").copy()
        if len(ss) == 0:
            raise RuntimeError(f"NO_TRACE:{station}")
        gaps = [g for g in ss.get_gaps() if abs(float(g[6])) > (0.25 / expected_fs)]
        if gaps:
            raise RuntimeError(f"GAPS_OR_OVERLAPS:{station}:{len(gaps)}")
        ss.merge(method=0)
        if len(ss) != 1:
            raise RuntimeError(f"TRACE_COUNT:{station}:{len(ss)}")
        tr = ss[0]
        fs = float(tr.stats.sampling_rate)
        if abs(fs - expected_fs) > 1e-6:
            raise RuntimeError(f"SAMPLE_RATE:{station}:{fs}")
        # Trim to common requested bounds. No interpolation.
        tr.trim(UTCDateTime(start), UTCDateTime(end), nearest_sample=True, pad=False)
        if int(tr.stats.npts) < int((end - start).total_seconds() * expected_fs) - 4:
            raise RuntimeError(f"INCOMPLETE_FULL_WINDOW:{station}:{tr.stats.npts}")
        traces[station] = tr
        meta[station] = {
            "npts": int(tr.stats.npts), "sampling_rate_hz": fs,
            "starttime": str(tr.stats.starttime), "endtime": str(tr.stats.endtime),
        }
    # Intersect time range across all six and truncate arrays to identical length.
    common_start = max(tr.stats.starttime for tr in traces.values())
    common_end = min(tr.stats.endtime for tr in traces.values())
    arrays: dict[str, np.ndarray] = {}
    for station, tr in traces.items():
        tt = tr.copy().trim(common_start, common_end, nearest_sample=True, pad=False)
        arrays[station] = np.asarray(tt.data, dtype=np.float64)
    n = min(len(a) for a in arrays.values())
    if n < int((end - start).total_seconds() * expected_fs) - 8:
        raise RuntimeError(f"COMMON_WINDOW_TOO_SHORT:{n}")
    arrays = {k: v[:n] for k, v in arrays.items()}
    meta["common"] = {"starttime": str(common_start), "endtime": str(common_end), "npts": n}
    return arrays, meta


def welch_psd(data: np.ndarray, fs: float, nperseg: int, noverlap: int) -> tuple[np.ndarray, np.ndarray]:
    return signal.welch(
        data, fs=fs, window="hann", nperseg=nperseg, noverlap=noverlap,
        detrend="constant", return_onesided=True, scaling="density",
    )


def band_power(freq: np.ndarray, psd: np.ndarray, low: float, high: float) -> float:
    mask = (freq >= low) & (freq <= high) & np.isfinite(psd)
    if int(np.count_nonzero(mask)) < 2:
        raise RuntimeError("INSUFFICIENT_BAND_BINS")
    value = float(np.trapezoid(psd[mask], freq[mask]))
    if not (math.isfinite(value) and value > 0):
        raise RuntimeError(f"INVALID_BAND_POWER:{value}")
    return value


def corrected_band_power(
    freq: np.ndarray,
    raw_psd: np.ndarray,
    inventory: Any | None,
    station: str,
    when: datetime,
    low: float,
    high: float,
) -> float | None:
    if inventory is None:
        return None
    mask = (freq >= low) & (freq <= high) & np.isfinite(raw_psd)
    ff = freq[mask]
    if len(ff) < 2:
        return None
    try:
        response = inventory.get_response(f"IM.{station}..EDH", UTCDateTime(when))
        complex_resp = response.get_evalresp_response_for_frequencies(ff, output="DEF")
        amp2 = np.abs(complex_resp) ** 2
        valid = np.isfinite(amp2) & (amp2 > 0)
        if int(np.count_nonzero(valid)) < 2:
            return None
        value = float(np.trapezoid(raw_psd[mask][valid] / amp2[valid], ff[valid]))
        return value if math.isfinite(value) and value > 0 else None
    except Exception:
        return None


def mean_band_coherence(a: np.ndarray, b: np.ndarray, fs: float, nperseg: int, noverlap: int, low: float, high: float) -> float:
    freq, coh = signal.coherence(
        a, b, fs=fs, window="hann", nperseg=nperseg, noverlap=noverlap,
        detrend="constant",
    )
    mask = (freq >= low) & (freq <= high) & np.isfinite(coh)
    if not np.any(mask):
        raise RuntimeError("NO_COHERENCE_BINS")
    return float(np.mean(coh[mask]))


def ratio_db(num: float | None, den: float | None) -> float | None:
    if num is None or den is None or num <= 0 or den <= 0:
        return None
    value = 10.0 * math.log10(num / den)
    return float(value) if math.isfinite(value) else None


def median_iqr(values: list[float]) -> dict[str, float] | None:
    if not values:
        return None
    arr = np.asarray(values, dtype=float)
    q25, q75 = np.percentile(arr, [25, 75], method="linear")
    return {"median": float(np.median(arr)), "q25": float(q25), "q75": float(q75), "iqr": float(q75 - q25)}


def classify_baseline(early_raw_median: float, reference_raw_median: float, floor_db: float, max_diff_db: float) -> str:
    if early_raw_median >= floor_db and abs(early_raw_median - reference_raw_median) <= max_diff_db:
        return "EARLY_LARGE_S2_N1_BASELINE_SIMILAR_TO_2014"
    if early_raw_median < floor_db and reference_raw_median >= floor_db:
        return "EARLY_S2_N1_BASELINE_MATERIALLY_LOWER_THAN_2014"
    return "EARLY_S2_N1_BASELINE_DIFFERENT_OR_MIXED"


def analyze_selected(
    session: requests.Session,
    protocol: dict[str, Any],
    selected_dates: list[str],
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    stations = [parse_channel_id(cid)[1] for cid in protocol["channels"]]
    clock = tuple(int(x) for x in protocol["discovery_contract"]["probe_clock_utc_each_day"].replace("Z", "").split(":"))
    duration_s = 600
    sc = protocol["spectral_contract"]
    low, high = map(float, sc["primary_band_hz"])
    nperseg = int(sc["welch_nperseg"])
    noverlap = int(sc["welch_noverlap"])

    inv_start = selected_dates[0] + "T00:00:00Z"
    inv_end = (date.fromisoformat(selected_dates[-1]) + timedelta(days=2)).isoformat() + "T00:00:00Z"
    inventories: dict[str, Any | None] = {}
    inventory_meta: dict[str, Any] = {}
    for station in stations:
        inv, meta = fetch_inventory(session, station, inv_start, inv_end)
        inventories[station] = inv
        inventory_meta[station] = meta

    rows: list[dict[str, Any]] = []
    pair_names = {
        "S1_S2": ("H10S1", "H10S2"), "S1_S3": ("H10S1", "H10S3"), "S2_S3": ("H10S2", "H10S3"),
        "N1_N2": ("H10N1", "H10N2"), "N1_N3": ("H10N1", "H10N3"), "N2_N3": ("H10N2", "H10N3"),
    }
    for ds in selected_dates:
        d = date.fromisoformat(ds)
        start = datetime(d.year, d.month, d.day, *clock, tzinfo=timezone.utc)
        end = start + timedelta(seconds=duration_s)
        payload, request_meta = request_post(session, selection_body(stations, [(start, end)]))
        row: dict[str, Any] = {
            "date": ds, "start_utc": start.isoformat().replace("+00:00", "Z"),
            "end_utc": end.isoformat().replace("+00:00", "Z"), "request": request_meta,
            "status": "BLOCKED", "stations": {}, "coherence": {},
        }
        if not payload:
            row["error"] = "SELECTED_WINDOW_DISAPPEARED_AT_ANALYSIS_FETCH"
            rows.append(row)
            continue
        try:
            arrays, trace_meta = parse_full_six(payload, stations, start, end)
            for station in stations:
                freq, raw_psd = welch_psd(arrays[station], 250.0, nperseg, noverlap)
                raw = band_power(freq, raw_psd, low, high)
                corrected = corrected_band_power(freq, raw_psd, inventories[station], station, start, low, high)
                row["stations"][station] = {
                    "trace": trace_meta[station],
                    "raw_integrated_power": raw,
                    "corrected_integrated_power_pa2": corrected,
                }
            for name, (a, b) in pair_names.items():
                row["coherence"][name] = mean_band_coherence(arrays[a], arrays[b], 250.0, nperseg, noverlap, low, high)
            row["raw_s2_minus_n1_db"] = ratio_db(
                row["stations"]["H10S2"]["raw_integrated_power"],
                row["stations"]["H10N1"]["raw_integrated_power"],
            )
            row["corrected_s2_minus_n1_db"] = ratio_db(
                row["stations"]["H10S2"]["corrected_integrated_power_pa2"],
                row["stations"]["H10N1"]["corrected_integrated_power_pa2"],
            )
            row["common_trace"] = trace_meta["common"]
            row["status"] = "ANALYZED"
        except Exception as exc:
            row["error"] = f"{type(exc).__name__}:{exc}"
        rows.append(row)

    analyzed = [r for r in rows if r["status"] == "ANALYZED"]
    summary: dict[str, Any] = {"analyzed_windows": len(analyzed)}
    if analyzed:
        raw_ratios = [float(r["raw_s2_minus_n1_db"]) for r in analyzed if r.get("raw_s2_minus_n1_db") is not None]
        corr_ratios = [float(r["corrected_s2_minus_n1_db"]) for r in analyzed if r.get("corrected_s2_minus_n1_db") is not None]
        summary["raw_s2_minus_n1_db"] = median_iqr(raw_ratios)
        summary["corrected_s2_minus_n1_db"] = median_iqr(corr_ratios)
        summary["coherence_epoch_medians"] = {
            name: float(np.median([r["coherence"][name] for r in analyzed])) for name in pair_names
        }
        summary["per_channel_median_power"] = {}
        for station in stations:
            raw = [float(r["stations"][station]["raw_integrated_power"]) for r in analyzed]
            corrected = [r["stations"][station]["corrected_integrated_power_pa2"] for r in analyzed]
            corrected = [float(x) for x in corrected if x is not None]
            summary["per_channel_median_power"][station] = {
                "raw": float(np.median(raw)),
                "corrected_pa2": float(np.median(corrected)) if corrected else None,
            }
    return rows, summary, inventory_meta


def run(output: Path) -> dict[str, Any]:
    ensure_noncanonical_output(output)
    protocol, reference, frozen = verify_frozen_inputs()
    channels = list(protocol["channels"])
    expected = [
        "IM.H10N1..EDH", "IM.H10N2..EDH", "IM.H10N3..EDH",
        "IM.H10S1..EDH", "IM.H10S2..EDH", "IM.H10S3..EDH",
    ]
    if channels != expected:
        raise RuntimeError("CHANNEL_SET_DRIFT")
    stations = [parse_channel_id(cid)[1] for cid in channels]
    dc = protocol["discovery_contract"]
    scan_start = datetime.fromisoformat(dc["scan_start_utc"].replace("Z", "+00:00"))
    scan_end = datetime.fromisoformat(dc["scan_end_utc_exclusive"].replace("Z", "+00:00"))
    clock = tuple(int(x) for x in dc["probe_clock_utc_each_day"].replace("Z", "").split(":"))
    probe_duration = int(dc["daily_probe_duration_s"])
    min_windows = int(dc["minimum_complete_windows_for_analysis"])
    target_windows = int(dc["target_window_count"])

    session = requests.Session()
    session.headers["User-Agent"] = UA

    earliest, annual_probe_audit = daily_probe_scan(
        session, stations, scan_start, scan_end, clock, probe_duration
    )
    selected_dates: list[str] = []
    full_window_audit: list[dict[str, Any]] = []
    if earliest is not None:
        selected_dates, full_window_audit = discover_full_windows(
            session, stations, earliest, clock, 600, target_windows, 120
        )

    rows: list[dict[str, Any]] = []
    analysis_summary: dict[str, Any] | None = None
    inventory_meta: dict[str, Any] = {}
    classification: str
    if earliest is None:
        classification = "NO_SIMULTANEOUS_SIX_CHANNEL_PROBE_FOUND_BEFORE_REFERENCE_HORIZON"
    elif len(selected_dates) < min_windows:
        classification = "BLOCKED_FIRST_AVAILABLE_SIX_CHANNEL_BASELINE_COMPLETENESS"
    else:
        rows, analysis_summary, inventory_meta = analyze_selected(session, protocol, selected_dates)
        analyzed = int(analysis_summary.get("analyzed_windows", 0)) if analysis_summary else 0
        if analyzed < min_windows or not analysis_summary or not analysis_summary.get("raw_s2_minus_n1_db"):
            classification = "BLOCKED_FIRST_AVAILABLE_SIX_CHANNEL_BASELINE_COMPLETENESS"
        else:
            early_raw = float(analysis_summary["raw_s2_minus_n1_db"]["median"])
            ref_raw = float(protocol["frozen_later_reference"]["reference_raw_median_s2_minus_n1_db"])
            pc = protocol["predeclared_baseline_classification"]
            classification = classify_baseline(
                early_raw, ref_raw, float(pc["large_offset_floor_db"]),
                float(pc["maximum_abs_early_vs_2014_raw_median_difference_for_similar_baseline_db"]),
            )

    fault_date = date(2013, 7, 19)
    relative = None
    if earliest is not None:
        relative = "BEFORE_2013_07_19_FAULT_ANCHOR" if earliest < fault_date else "ON_OR_AFTER_2013_07_19_FAULT_ANCHOR"

    receipt = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-HA10-FIRST-AVAILABLE-SIX-CHANNEL-BASELINE-RUN",
        "created_utc": utc_now(),
        "gate_id": protocol["gate_id"],
        "protocol_git_blob_sha1": EXPECTED_PROTOCOL_BLOB,
        "reference_git_blob_sha1": EXPECTED_REFERENCE_BLOB,
        "frozen_119hz_result_git_blob_sha1": EXPECTED_FROZEN_119_BLOB,
        "frozen_119hz_verdict": frozen["summary"]["verdict"],
        "authority_delta_for_119hz": 0,
        "target_identity": "UNCONFIRMED",
        "source_writeback": False,
        "selection_used_waveform_values": False,
        "discovery": {
            "scan_start_utc": dc["scan_start_utc"],
            "scan_end_utc_exclusive": dc["scan_end_utc_exclusive"],
            "probe_clock_utc": dc["probe_clock_utc_each_day"],
            "probe_duration_s": probe_duration,
            "earliest_simultaneous_six_channel_probe_date": earliest.isoformat() if earliest else None,
            "earliest_relative_to_2013_fault": relative,
            "annual_probe_audit": annual_probe_audit,
            "full_window_audit": full_window_audit,
            "selected_dates": selected_dates,
            "selected_count": len(selected_dates),
        },
        "analysis": {
            "classification": classification,
            "summary": analysis_summary,
            "windows": rows,
            "inventory_metadata": inventory_meta,
            "frozen_2014_reference": {
                "verdict": reference.get("result", {}).get("verdict"),
                "raw_median_s2_minus_n1_db": protocol["frozen_later_reference"]["reference_raw_median_s2_minus_n1_db"],
                "corrected_median_s2_minus_n1_db": protocol["frozen_later_reference"]["reference_corrected_median_s2_minus_n1_db"],
            },
        },
        "interpretation_firewall": protocol["interpretation_firewall"],
        "canonical_answer": (
            f"{classification}; earliest simultaneous six-channel probe="
            f"{earliest.isoformat() if earliest else 'NONE'}; this is a temporal/archive baseline, not an instrument-clean label."
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise RuntimeError("OUTPUT_ALREADY_EXISTS_APPEND_ONLY_REQUIRED")
    output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "classification": classification,
        "earliest_probe": earliest.isoformat() if earliest else None,
        "selected_windows": len(selected_dates),
        "analyzed_windows": (analysis_summary or {}).get("analyzed_windows", 0),
        "authority_delta_for_119hz": 0,
    }, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    run(Path(args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
