#!/usr/bin/env python3
"""Transport-repaired execution wrapper for the frozen earliest-six-channel gate.

V1 remains byte-preserved. This wrapper changes only how the value-blind daily
probe selectors are transported to EarthScope: one logical year is partitioned
into deterministic chronological calendar-month POSTs after run 32578915494
showed that a whole-year POST can hit a 60 s read timeout. Scientific selection,
analysis, thresholds and authority are inherited unchanged from V1.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from collections import defaultdict
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Iterator

ROOT = Path(__file__).resolve().parents[1]
V1_PATH = Path(__file__).with_name("cousteau_ha10_first_available_six_channel_baseline.py")
ADDENDUM = ROOT / "data" / "cousteau" / "JANUS-ECHO-COUSTEAU-HA10-FIRST-AVAILABLE-SIX-CHANNEL-BASELINE-TRANSPORT-ADDENDUM-2026-08-22-v1.0.json"
EXPECTED_V1_BLOB = "89eb62a693cf5dd0fb1eada5ccb4477ff494e776"
EXPECTED_ADDENDUM_BLOB = "106ee138ff66f5d3204ad37bbee6a0ff225fc577"

_SPEC = importlib.util.spec_from_file_location("cousteau_first_six_v1_frozen", V1_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError("CANNOT_LOAD_FROZEN_V1")
v1 = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(v1)


def git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def verify_transport_addendum() -> dict[str, Any]:
    if git_blob_sha1(V1_PATH) != EXPECTED_V1_BLOB:
        raise RuntimeError("V1_RUNNER_BLOB_DRIFT")
    if git_blob_sha1(ADDENDUM) != EXPECTED_ADDENDUM_BLOB:
        raise RuntimeError("TRANSPORT_ADDENDUM_BLOB_DRIFT")
    payload = json.loads(ADDENDUM.read_text(encoding="utf-8"))
    if payload.get("status") != "APPEND_ONLY_TRANSPORT_REPAIR_AFTER_INFRASTRUCTURE_FAILURE_BEFORE_SUCCESSFUL_ACQUISITION":
        raise RuntimeError("TRANSPORT_ADDENDUM_STATUS_DRIFT")
    repair = payload.get("transport_only_repair", {})
    required_true = [
        "logical_scan_order_unchanged",
        "scan_start_unchanged",
        "scan_end_unchanged",
        "probe_clock_unchanged",
        "probe_duration_unchanged",
        "six_channel_set_unchanged",
        "earliest_date_selection_rule_unchanged",
        "waveform_value_blinding_unchanged",
        "scientific_thresholds_unchanged",
        "spectral_contract_unchanged",
    ]
    if not all(repair.get(k) is True for k in required_true):
        raise RuntimeError("TRANSPORT_REPAIR_CHANGED_SCIENTIFIC_CONTRACT")
    if payload.get("epistemic_firewall", {}).get("authority_delta_for_119hz") != 0:
        raise RuntimeError("TRANSPORT_ADDENDUM_119_AUTHORITY_DRIFT")
    return payload


def month_start(dt: datetime) -> datetime:
    return datetime(dt.year, dt.month, 1, tzinfo=timezone.utc)


def next_month(dt: datetime) -> datetime:
    if dt.month == 12:
        return datetime(dt.year + 1, 1, 1, tzinfo=timezone.utc)
    return datetime(dt.year, dt.month + 1, 1, tzinfo=timezone.utc)


def iter_month_chunks(scan_start: datetime, scan_end: datetime) -> Iterator[tuple[datetime, datetime]]:
    if scan_start.tzinfo is None or scan_end.tzinfo is None:
        raise ValueError("MONTH_CHUNKS_REQUIRE_AWARE_UTC_DATETIMES")
    cursor = month_start(scan_start)
    while cursor < scan_end:
        chunk_start = max(cursor, scan_start)
        chunk_end = min(next_month(cursor), scan_end)
        if chunk_start < chunk_end:
            yield chunk_start, chunk_end
        cursor = next_month(cursor)


def daily_probe_scan_chunked(
    session: Any,
    stations: list[str],
    scan_start: datetime,
    scan_end: datetime,
    clock: tuple[int, int, int],
    duration_s: int,
) -> tuple[date | None, list[dict[str, Any]]]:
    """Value-blind earliest-date discovery using chronological month transport chunks.

    No chunk is skipped. Any exhausted transport failure raises and therefore cannot
    be mistaken for a scientific no-data result or allow a later date to win.
    """
    audit: list[dict[str, Any]] = []
    for chunk_index, (chunk_start, chunk_end) in enumerate(iter_month_chunks(scan_start, scan_end), start=1):
        windows: list[tuple[datetime, datetime]] = []
        d = chunk_start.date()
        last_date = (chunk_end - timedelta(microseconds=1)).date()
        while d <= last_date:
            start = datetime(d.year, d.month, d.day, clock[0], clock[1], clock[2], tzinfo=timezone.utc)
            if scan_start <= start < scan_end and chunk_start <= start < chunk_end:
                windows.append((start, start + timedelta(seconds=duration_s)))
            d += timedelta(days=1)
        if not windows:
            continue

        payload, meta = v1.request_post(session, v1.selection_body(stations, windows))
        row: dict[str, Any] = {
            "transport_mode": "DETERMINISTIC_CALENDAR_MONTH_CHUNK",
            "chunk_index": chunk_index,
            "logical_year": chunk_start.year,
            "month": chunk_start.month,
            "chunk_start_utc": chunk_start.isoformat().replace("+00:00", "Z"),
            "chunk_end_utc_exclusive": chunk_end.isoformat().replace("+00:00", "Z"),
            "date_count": len(windows),
            "request": meta,
            "all_six_probe_dates": [],
        }
        if payload:
            headers = v1.intervals_from_headonly(payload)
            station_sets: dict[str, set[str]] = defaultdict(set)
            for (day, station), intervals in headers.items():
                if intervals:
                    station_sets[day].add(station)
            all_six = sorted(day for day, ss in station_sets.items() if set(stations).issubset(ss))
            row["all_six_probe_dates"] = all_six
            row["header_trace_keys"] = len(headers)
            audit.append(row)
            if all_six:
                return date.fromisoformat(all_six[0]), audit
        else:
            audit.append(row)
    return None, audit


def run(output: Path) -> dict[str, Any]:
    addendum = verify_transport_addendum()
    original = v1.daily_probe_scan
    try:
        v1.daily_probe_scan = daily_probe_scan_chunked
        receipt = v1.run(output)
    finally:
        v1.daily_probe_scan = original
    receipt["transport_repair"] = {
        "wrapper": "cousteau_ha10_first_available_six_channel_baseline_v2.py",
        "frozen_v1_git_blob_sha1": EXPECTED_V1_BLOB,
        "transport_addendum_git_blob_sha1": EXPECTED_ADDENDUM_BLOB,
        "failed_predecessor_run_id": addendum["failed_execution"]["github_actions_run_id"],
        "mode": "DETERMINISTIC_CALENDAR_MONTH_CHUNKS",
        "scientific_contract_unchanged": True,
        "selection_used_waveform_values": False,
        "authority_delta_for_119hz": 0,
    }
    # V1 already wrote the receipt before this wrapper annotation. Rewrite only the
    # ephemeral run-scoped artifact path, never canonical data/cousteau.
    output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    run(Path(args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
