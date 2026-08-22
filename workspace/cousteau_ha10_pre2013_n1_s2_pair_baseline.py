#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import json
import math
from collections import defaultdict
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import numpy as np
import requests
from obspy import UTCDateTime, read

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cousteau"
PROTOCOL = DATA / "JANUS-ECHO-COUSTEAU-HA10-PRE2013-N1-S2-PAIR-BASELINE-PROTOCOL-2026-08-22-v1.0.json"
REFERENCE = DATA / "JANUS-ECHO-COUSTEAU-HA10-N1-S2-TEMPORAL-SCALE-STABILITY-DIAGNOSTIC-RUN-001-SUMMARY-2026-08-22-v1.0.json"
FROZEN_119 = DATA / "JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0.json"
V1_PATH = Path(__file__).with_name("cousteau_ha10_first_available_six_channel_baseline.py")

EXPECTED_PROTOCOL_BLOB = "b4cb0298672bdb34147737e1efcde1b61becb9b7"
EXPECTED_REFERENCE_BLOB = "45f2901adde866df610a563f27c50177d32be414"
EXPECTED_FROZEN_119_BLOB = "eb8b48fb7f043160c057f9df6264a781412ed854"
EXPECTED_V1_BLOB = "89eb62a693cf5dd0fb1eada5ccb4477ff494e776"

spec = importlib.util.spec_from_file_location("cousteau_first_six_v1", V1_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError("CANNOT_LOAD_V1_HELPERS")
v1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v1)


def git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def verify_inputs() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    for path, expected in {
        PROTOCOL: EXPECTED_PROTOCOL_BLOB,
        REFERENCE: EXPECTED_REFERENCE_BLOB,
        FROZEN_119: EXPECTED_FROZEN_119_BLOB,
        V1_PATH: EXPECTED_V1_BLOB,
    }.items():
        actual = git_blob_sha1(path)
        if actual != expected:
            raise RuntimeError(f"FROZEN_BLOB_DRIFT:{path.name}:{actual}:{expected}")
    p = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    ref = json.loads(REFERENCE.read_text(encoding="utf-8"))
    frozen = json.loads(FROZEN_119.read_text(encoding="utf-8"))
    if p.get("status") != "PREREGISTERED_BEFORE_PRE2013_PAIR_WAVEFORM_VALUE_INSPECTION":
        raise RuntimeError("PROTOCOL_STATUS_DRIFT")
    if p.get("channels") != ["IM.H10N1..EDH", "IM.H10S2..EDH"]:
        raise RuntimeError("PAIR_CHANNEL_DRIFT")
    if p.get("authority", {}).get("authority_delta_for_119hz") != 0:
        raise RuntimeError("119_AUTHORITY_DRIFT")
    return p, ref, frozen


def month_chunks(start: datetime, end: datetime):
    cur = datetime(start.year, start.month, 1, tzinfo=timezone.utc)
    while cur < end:
        nxt = datetime(cur.year + 1, 1, 1, tzinfo=timezone.utc) if cur.month == 12 else datetime(cur.year, cur.month + 1, 1, tzinfo=timezone.utc)
        a, b = max(cur, start), min(nxt, end)
        if a < b:
            yield a, b
        cur = nxt


def find_earliest_pair(session: requests.Session, start: datetime, end: datetime, stations: list[str]) -> tuple[date | None, list[dict[str, Any]]]:
    audit: list[dict[str, Any]] = []
    for idx, (a, b) in enumerate(month_chunks(start, end), 1):
        windows = []
        d = a.date()
        last = (b - timedelta(microseconds=1)).date()
        while d <= last:
            s = datetime(d.year, d.month, d.day, 12, 0, 0, tzinfo=timezone.utc)
            if a <= s < b:
                windows.append((s, s + timedelta(seconds=1)))
            d += timedelta(days=1)
        payload, meta = v1.request_post(session, v1.selection_body(stations, windows))
        row: dict[str, Any] = {"chunk_index": idx, "start": a.isoformat(), "end": b.isoformat(), "request": meta, "pair_dates": []}
        if payload:
            headers = v1.intervals_from_headonly(payload)
            sets: dict[str, set[str]] = defaultdict(set)
            for (day, station), intervals in headers.items():
                if intervals:
                    sets[day].add(station)
            pair_dates = sorted(day for day, ss in sets.items() if set(stations).issubset(ss))
            row["pair_dates"] = pair_dates
            row["header_trace_keys"] = len(headers)
            audit.append(row)
            if pair_dates:
                return date.fromisoformat(pair_dates[0]), audit
        else:
            audit.append(row)
    return None, audit


def freeze_pair_windows(session: requests.Session, earliest: date, end_exclusive: datetime, stations: list[str], target: int) -> tuple[list[str], list[dict[str, Any]]]:
    selected: list[str] = []
    audit: list[dict[str, Any]] = []
    d = earliest
    while len(selected) < target:
        s = datetime(d.year, d.month, d.day, 12, 0, 0, tzinfo=timezone.utc)
        if s >= end_exclusive:
            break
        e = s + timedelta(seconds=600)
        payload, meta = v1.request_post(session, v1.selection_body(stations, [(s, e)]))
        complete = {st: False for st in stations}
        if payload:
            headers = v1.intervals_from_headonly(payload)
            for st in stations:
                complete[st] = v1.header_covers_window(headers.get((d.isoformat(), st), []), s, e)
        ok = all(complete.values())
        audit.append({"date": d.isoformat(), "request": meta, "complete": complete, "pair_complete": ok})
        if ok:
            selected.append(d.isoformat())
        d += timedelta(days=1)
    return selected, audit


def parse_pair(payload: bytes, stations: list[str], start: datetime, end: datetime) -> dict[str, np.ndarray]:
    stream = read(io.BytesIO(payload), format="MSEED")
    out: dict[str, np.ndarray] = {}
    traces = {}
    for st in stations:
        ss = stream.select(network="IM", station=st, location="", channel="EDH").copy()
        if not ss:
            raise RuntimeError(f"NO_TRACE:{st}")
        gaps = [g for g in ss.get_gaps() if abs(float(g[6])) > (0.25 / 250.0)]
        if gaps:
            raise RuntimeError(f"GAPS:{st}:{len(gaps)}")
        ss.merge(method=0)
        if len(ss) != 1:
            raise RuntimeError(f"TRACE_COUNT:{st}:{len(ss)}")
        tr = ss[0]
        if abs(float(tr.stats.sampling_rate) - 250.0) > 1e-6:
            raise RuntimeError(f"SAMPLE_RATE:{st}:{tr.stats.sampling_rate}")
        tr.trim(UTCDateTime(start), UTCDateTime(end), nearest_sample=True, pad=False)
        traces[st] = tr
    common_start = max(t.stats.starttime for t in traces.values())
    common_end = min(t.stats.endtime for t in traces.values())
    n = None
    for st, tr in traces.items():
        tt = tr.copy().trim(common_start, common_end, nearest_sample=True, pad=False)
        arr = np.asarray(tt.data, dtype=np.float64)
        n = len(arr) if n is None else min(n, len(arr))
        out[st] = arr
    if n is None or n < 149990:
        raise RuntimeError(f"COMMON_WINDOW_TOO_SHORT:{n}")
    return {st: arr[:n] for st, arr in out.items()}


def stats(values: list[float]) -> dict[str, float]:
    a = np.asarray(values, dtype=float)
    q25, q75 = np.percentile(a, [25, 75])
    return {"median": float(np.median(a)), "q25": float(q25), "q75": float(q75), "iqr": float(q75-q25), "min": float(np.min(a)), "max": float(np.max(a))}


def run(output: Path) -> dict[str, Any]:
    v1.ensure_noncanonical_output(output)
    p, ref, frozen = verify_inputs()
    dc = p["discovery_contract"]
    start = datetime.fromisoformat(dc["scan_start_utc"].replace("Z", "+00:00"))
    end = datetime.fromisoformat(dc["scan_end_utc_exclusive"].replace("Z", "+00:00"))
    stations = ["H10N1", "H10S2"]
    session = requests.Session()
    session.headers["User-Agent"] = "JANUS-Cousteau-Pre2013-N1S2-Solo/1.0"

    earliest, discovery_audit = find_earliest_pair(session, start, end, stations)
    selected: list[str] = []
    freeze_audit: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    summary = None
    classification: str

    if earliest is None:
        classification = dc["if_no_pair"]
    else:
        selected, freeze_audit = freeze_pair_windows(session, earliest, end, stations, int(dc["target_window_count"]))
        if len(selected) < int(dc["minimum_complete_windows_for_analysis"]):
            classification = dc["if_insufficient_windows"]
        else:
            invs = {}
            inv_meta = {}
            inv_start = selected[0] + "T00:00:00Z"
            inv_end = (date.fromisoformat(selected[-1]) + timedelta(days=1)).isoformat() + "T23:59:59Z"
            for st in stations:
                invs[st], inv_meta[st] = v1.fetch_inventory(session, st, inv_start, inv_end)
            raw_ratios, corr_ratios, coherences = [], [], []
            for ds in selected:
                dd = date.fromisoformat(ds)
                s = datetime(dd.year, dd.month, dd.day, 12, 0, 0, tzinfo=timezone.utc)
                e = s + timedelta(seconds=600)
                payload, meta = v1.request_post(session, v1.selection_body(stations, [(s, e)]))
                row: dict[str, Any] = {"date": ds, "request": meta, "status": "BLOCKED"}
                if not payload:
                    row["error"] = "WINDOW_DISAPPEARED"
                    rows.append(row)
                    continue
                try:
                    arr = parse_pair(payload, stations, s, e)
                    powers = {}
                    corrected = {}
                    for st in stations:
                        f, psd = v1.welch_psd(arr[st], 250.0, 8192, 4096)
                        powers[st] = v1.band_power(f, psd, 10.0, 80.0)
                        corrected[st] = v1.corrected_band_power(f, psd, invs[st], st, s, 10.0, 80.0)
                    rr = v1.ratio_db(powers["H10S2"], powers["H10N1"])
                    cr = v1.ratio_db(corrected["H10S2"], corrected["H10N1"])
                    coh = v1.mean_band_coherence(arr["H10N1"], arr["H10S2"], 250.0, 8192, 4096, 10.0, 80.0)
                    row.update({"status":"ANALYZED", "raw_power":powers, "corrected_power_pa2":corrected, "raw_s2_minus_n1_db":rr, "corrected_s2_minus_n1_db":cr, "n1_s2_coherence":coh})
                    if rr is not None: raw_ratios.append(float(rr))
                    if cr is not None: corr_ratios.append(float(cr))
                    coherences.append(float(coh))
                except Exception as exc:
                    row["error"] = f"{type(exc).__name__}:{exc}"
                rows.append(row)
            analyzed = sum(r.get("status") == "ANALYZED" for r in rows)
            if analyzed < int(dc["minimum_complete_windows_for_analysis"]) or len(raw_ratios) < int(dc["minimum_complete_windows_for_analysis"]):
                classification = dc["if_insufficient_windows"]
            else:
                raw_stat = stats(raw_ratios)
                corr_stat = stats(corr_ratios) if corr_ratios else None
                summary = {"analyzed_windows": analyzed, "raw_s2_minus_n1_db": raw_stat, "corrected_s2_minus_n1_db": corr_stat, "median_n1_s2_coherence": float(np.median(coherences)), "inventory_metadata": inv_meta}
                ref_raw = float(p["frozen_later_reference"]["reference_raw_median_s2_minus_n1_db"])
                floor = float(p["predeclared_classification"]["large_offset_floor_db"])
                maxdiff = float(p["predeclared_classification"]["maximum_abs_difference_from_2014_reference_for_similar_db"])
                med = raw_stat["median"]
                if med >= floor and abs(med-ref_raw) <= maxdiff:
                    classification = "PRE2013_LARGE_OFFSET_SIMILAR_TO_2014"
                elif med < floor and ref_raw >= floor:
                    classification = "PRE2013_OFFSET_MATERIALLY_LOWER_THAN_2014"
                else:
                    classification = "PRE2013_OFFSET_DIFFERENT_OR_MIXED"

    receipt = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-HA10-PRE2013-N1-S2-PAIR-BASELINE-RUN",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "gate_id": p["gate_id"],
        "protocol_git_blob_sha1": EXPECTED_PROTOCOL_BLOB,
        "reference_git_blob_sha1": EXPECTED_REFERENCE_BLOB,
        "frozen_119hz_result_git_blob_sha1": EXPECTED_FROZEN_119_BLOB,
        "selection_used_waveform_values": False,
        "source_writeback": False,
        "authority_delta_for_119hz": 0,
        "target_identity": "UNCONFIRMED",
        "discovery": {"earliest_pre2013_pair_date": earliest.isoformat() if earliest else None, "selected_dates": selected, "selected_count": len(selected), "monthly_audit": discovery_audit, "freeze_audit": freeze_audit},
        "analysis": {"classification": classification, "summary": summary, "windows": rows, "frozen_2014_reference_raw_db": p["frozen_later_reference"]["reference_raw_median_s2_minus_n1_db"]},
        "interpretation_firewall": p["interpretation_firewall"],
        "frozen_119hz_verdict": frozen.get("summary", {}).get("verdict"),
        "canonical_answer": f"{classification}; earliest pre-2013 simultaneous public N1/S2 pair={earliest.isoformat() if earliest else 'NONE'}."
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise RuntimeError("OUTPUT_ALREADY_EXISTS_APPEND_ONLY_REQUIRED")
    output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps({"classification":classification, "earliest":earliest.isoformat() if earliest else None, "selected":len(selected), "raw_median_db":(summary or {}).get("raw_s2_minus_n1_db",{}).get("median"), "authority_delta_for_119hz":0}, sort_keys=True))
    return receipt


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    run(Path(args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
