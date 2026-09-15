#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import numpy as np
import requests

from workspace import cousteau_ha10_tphase_inband_positive_control as spectral_v1
from workspace import cousteau_ha10_tphase_inband_positive_control_v2 as transport_v2

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cousteau"
PROTOCOL = DATA / "JANUS-ECHO-COUSTEAU-HA10-N1-S2-TEMPORAL-SCALE-STABILITY-DIAGNOSTIC-PROTOCOL-2026-08-22-v1.0.json"
SHAPE_SUMMARY = DATA / "JANUS-ECHO-COUSTEAU-HA10-N1-S2-SPECTRAL-SHAPE-VS-SCALE-DIAGNOSTIC-RUN-001-SUMMARY-2026-08-22-v1.0.json"
FROZEN_119 = DATA / "JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0.json"
SPECTRAL_V1 = ROOT / "workspace" / "cousteau_ha10_tphase_inband_positive_control.py"
TRANSPORT_V2 = ROOT / "workspace" / "cousteau_ha10_tphase_inband_positive_control_v2.py"

EXPECTED_PROTOCOL_BLOB = "bd7df63f62b16428e30b3cb37ee76cc556e9917b"
EXPECTED_SHAPE_SUMMARY_BLOB = "26c7b3cf8a5e0707b5ca254ee9ae9d35c064c18b"
EXPECTED_FROZEN_119_BLOB = "eb8b48fb7f043160c057f9df6264a781412ed854"
EXPECTED_SPECTRAL_V1_BLOB = "244169a1331a16529e2a963585a5659dae109c66"
EXPECTED_TRANSPORT_V2_BLOB = "c3d37fe4f5a75514390862be4c1d2870c4c01fe9"
EXPECTED_SHAPE_VERDICT = "BROADBAND_SCALE_LIKE_H10S2_ASYMMETRY"
EXPECTED_FROZEN_119_VERDICT = "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE"


def git_blob_sha1_file(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def ratio_db(numerator: float, denominator: float) -> float:
    if not (
        numerator > 0
        and denominator > 0
        and math.isfinite(numerator)
        and math.isfinite(denominator)
    ):
        return float("nan")
    return 10.0 * math.log10(numerator / denominator)


def parse_start(date_text: str, hhmmss: str) -> datetime:
    return datetime.fromisoformat(f"{date_text}T{hhmmss.replace('Z', '+00:00')}")


def ensure_noncanonical_output(path: Path) -> None:
    try:
        path.resolve().relative_to(DATA.resolve())
    except ValueError:
        return
    raise RuntimeError("CANONICAL_DATA_WRITE_FORBIDDEN_BY_TEMPORAL_STABILITY_DIAGNOSTIC")


def verify_frozen_contracts() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    expected = {
        PROTOCOL: EXPECTED_PROTOCOL_BLOB,
        SHAPE_SUMMARY: EXPECTED_SHAPE_SUMMARY_BLOB,
        FROZEN_119: EXPECTED_FROZEN_119_BLOB,
        SPECTRAL_V1: EXPECTED_SPECTRAL_V1_BLOB,
        TRANSPORT_V2: EXPECTED_TRANSPORT_V2_BLOB,
    }
    for path, blob in expected.items():
        if git_blob_sha1_file(path) != blob:
            raise RuntimeError(f"FROZEN_BLOB_DRIFT:{path.name}")

    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    shape = json.loads(SHAPE_SUMMARY.read_text(encoding="utf-8"))
    frozen = json.loads(FROZEN_119.read_text(encoding="utf-8"))

    if protocol.get("status") != "PREREGISTERED_BEFORE_NEW_0600UTC_WINDOWS_ARE_DOWNLOADED":
        raise RuntimeError("TEMPORAL_STABILITY_PROTOCOL_NOT_PREREGISTERED")
    if shape.get("result", {}).get("verdict") != EXPECTED_SHAPE_VERDICT:
        raise RuntimeError("SHAPE_VERDICT_BINDING_DRIFT")
    if frozen.get("summary", {}).get("verdict") != EXPECTED_FROZEN_119_VERDICT:
        raise RuntimeError("FROZEN_119_VERDICT_DRIFT")
    if protocol.get("epistemic_position", {}).get("authority_delta_for_119hz") != 0:
        raise RuntimeError("AUTHORITY_DELTA_FOR_119HZ_NOT_ZERO")
    if protocol.get("frequency_contract", {}).get("119hz_or_117_121hz_bins_may_be_used") is not False:
        raise RuntimeError("TARGET_BAND_EXCLUSION_NOT_FROZEN")
    if len(protocol["new_window_contract"]["dates"]) != 20:
        raise RuntimeError("NEW_WINDOW_COUNT_DRIFT")
    if protocol["new_window_contract"]["window_start_time_utc_each_day"] != "06:00:00Z":
        raise RuntimeError("NEW_WINDOW_CLOCK_DRIFT")
    return protocol, shape, frozen


def integrated_power(trace, inventory, *, band: list[float], nperseg: int, noverlap: int) -> dict[str, float]:
    frequencies, raw_psd, corrected_psd = spectral_v1.psd_with_response(
        trace,
        inventory,
        nperseg=nperseg,
        noverlap=noverlap,
    )
    raw_power, _ = spectral_v1.integrated_band_power(frequencies, raw_psd, band)
    corrected_power, _ = spectral_v1.integrated_band_power(
        frequencies, corrected_psd, band
    )
    return {
        "raw_integrated_power": raw_power,
        "corrected_integrated_power_pa2": corrected_power,
    }


def stability_statistics(rows: list[dict[str, Any]], floor_db: float) -> dict[str, Any]:
    raw = np.asarray([float(row["raw_s2_minus_n1_db"]) for row in rows], dtype=float)
    corrected = np.asarray(
        [float(row["corrected_s2_minus_n1_db"]) for row in rows], dtype=float
    )
    raw_q25, raw_q75 = np.percentile(raw, [25, 75], method="linear")
    corrected_q25, corrected_q75 = np.percentile(
        corrected, [25, 75], method="linear"
    )
    raw_above = raw >= floor_db
    corrected_above = corrected >= floor_db
    both_above = raw_above & corrected_above
    return {
        "percentile_method": "numpy_linear",
        "median_raw_s2_minus_n1_db": float(np.median(raw)),
        "median_corrected_s2_minus_n1_db": float(np.median(corrected)),
        "raw_q25_db": float(raw_q25),
        "raw_q75_db": float(raw_q75),
        "corrected_q25_db": float(corrected_q25),
        "corrected_q75_db": float(corrected_q75),
        "raw_iqr_db": float(raw_q75 - raw_q25),
        "corrected_iqr_db": float(corrected_q75 - corrected_q25),
        "fraction_raw_windows_ge_floor": float(np.mean(raw_above)),
        "fraction_corrected_windows_ge_floor": float(np.mean(corrected_above)),
        "fraction_both_raw_and_corrected_windows_ge_floor": float(
            np.mean(both_above)
        ),
        "raw_min_db": float(np.min(raw)),
        "raw_max_db": float(np.max(raw)),
        "corrected_min_db": float(np.min(corrected)),
        "corrected_max_db": float(np.max(corrected)),
    }


def decide(
    rows: list[dict[str, Any]],
    *,
    minimum_pairs: int,
    floor_db: float,
    max_iqr_db: float,
    minimum_fraction: float,
) -> tuple[str, dict[str, Any]]:
    if len(rows) < minimum_pairs:
        return "BLOCKED_TEMPORAL_SCALE_STABILITY_DATA_ACCESS", {
            "complete_paired_windows": len(rows),
            "statistics": None,
        }

    stats = stability_statistics(rows, floor_db)
    medians_ok = (
        stats["median_raw_s2_minus_n1_db"] >= floor_db
        and stats["median_corrected_s2_minus_n1_db"] >= floor_db
    )
    iqr_ok = (
        stats["raw_iqr_db"] <= max_iqr_db
        and stats["corrected_iqr_db"] <= max_iqr_db
    )
    fractions_ok = (
        stats["fraction_raw_windows_ge_floor"] >= minimum_fraction
        and stats["fraction_corrected_windows_ge_floor"] >= minimum_fraction
        and stats["fraction_both_raw_and_corrected_windows_ge_floor"]
        >= minimum_fraction
    )
    stats.update(
        {
            "complete_paired_windows": len(rows),
            "medians_ge_floor": medians_ok,
            "iqrs_le_maximum": iqr_ok,
            "fractions_ge_minimum": fractions_ok,
        }
    )
    verdict = (
        "TEMPORALLY_STABLE_BROADBAND_H10S2_SCALE_OFFSET"
        if medians_ok and iqr_ok and fractions_ok
        else "TEMPORALLY_VARIABLE_H10S2_SCALE_OFFSET"
    )
    return verdict, stats


def run(output: Path) -> dict[str, Any]:
    ensure_noncanonical_output(output)
    protocol, shape, frozen = verify_frozen_contracts()

    wc = protocol["new_window_contract"]
    estimator = protocol["spectral_estimator"]
    thresholds = protocol["classification_thresholds"]
    band = [float(x) for x in protocol["frequency_contract"]["band_hz"]]
    nperseg = int(estimator["welch_nperseg_samples"])
    noverlap = int(estimator["welch_noverlap_samples"])
    duration_s = int(wc["window_duration_s"])
    minimum_pairs = int(wc["minimum_complete_paired_windows"])
    floor_db = float(thresholds["large_asymmetry_floor_db"])
    max_iqr_db = float(thresholds["maximum_stable_iqr_db"])
    minimum_fraction = float(thresholds["minimum_fraction_windows_above_floor"])
    channels = list(protocol["channels"])
    if channels != ["IM.H10N1..EDH", "IM.H10S2..EDH"]:
        raise RuntimeError("CHANNEL_SET_DRIFT")

    session = requests.Session()
    session.headers["User-Agent"] = (
        "Janus-Echo-Cousteau/1.0 preregistered H10 temporal scale stability diagnostic"
    )

    first = parse_start(wc["dates"][0], wc["window_start_time_utc_each_day"])
    last = parse_start(wc["dates"][-1], wc["window_start_time_utc_each_day"])
    inv_start = (first - timedelta(days=1)).isoformat().replace("+00:00", "Z")
    inv_end = (last + timedelta(days=2)).isoformat().replace("+00:00", "Z")

    inventories: dict[str, Any] = {}
    inventory_meta: dict[str, Any] = {}
    inventory_errors: dict[str, str] = {}
    for cid in channels:
        try:
            inventory, metadata = transport_v2.fetch_inventory_bounded(
                session, cid, inv_start, inv_end
            )
            inventories[cid] = inventory
            inventory_meta[cid] = metadata
        except Exception as exc:
            inventory_errors[cid] = f"{type(exc).__name__}:{exc}"

    windows: list[dict[str, Any]] = []
    complete: list[dict[str, Any]] = []
    data_errors: list[dict[str, Any]] = []

    for date_text in wc["dates"]:
        start_dt = parse_start(date_text, wc["window_start_time_utc_each_day"])
        end_dt = start_dt + timedelta(seconds=duration_s)
        start = start_dt.isoformat().replace("+00:00", "Z")
        end = end_dt.isoformat().replace("+00:00", "Z")
        row: dict[str, Any] = {
            "date": date_text,
            "start_utc": start,
            "end_utc": end,
            "stations": {},
            "pair_complete": False,
        }
        for cid in channels:
            try:
                if cid not in inventories:
                    raise RuntimeError(
                        "RESPONSE_UNAVAILABLE:" + inventory_errors.get(cid, "UNKNOWN")
                    )
                trace, waveform = transport_v2.fetch_trace_bounded(
                    session, cid, start, end
                )
                power = integrated_power(
                    trace,
                    inventories[cid],
                    band=band,
                    nperseg=nperseg,
                    noverlap=noverlap,
                )
                row["stations"][cid] = {
                    "data_status": "ANALYZED",
                    "waveform": waveform,
                    **power,
                }
            except Exception as exc:
                row["stations"][cid] = {
                    "data_status": "BLOCKED",
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                }
                data_errors.append(
                    {
                        "date": date_text,
                        "station": cid,
                        "error_type": type(exc).__name__,
                        "error": str(exc),
                    }
                )

        n1 = row["stations"].get(channels[0], {})
        s2 = row["stations"].get(channels[1], {})
        if n1.get("data_status") == "ANALYZED" and s2.get("data_status") == "ANALYZED":
            row["pair_complete"] = True
            row["raw_s2_minus_n1_db"] = ratio_db(
                float(s2["raw_integrated_power"]),
                float(n1["raw_integrated_power"]),
            )
            row["corrected_s2_minus_n1_db"] = ratio_db(
                float(s2["corrected_integrated_power_pa2"]),
                float(n1["corrected_integrated_power_pa2"]),
            )
            complete.append(row)
        windows.append(row)

    verdict, statistics = decide(
        complete,
        minimum_pairs=minimum_pairs,
        floor_db=floor_db,
        max_iqr_db=max_iqr_db,
        minimum_fraction=minimum_fraction,
    )

    receipt = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-HA10-N1-S2-TEMPORAL-SCALE-STABILITY-DIAGNOSTIC-RUN",
        "created_utc": utc_now(),
        "gate_id": protocol["gate_id"],
        "protocol_git_blob_sha1": EXPECTED_PROTOCOL_BLOB,
        "shape_summary_git_blob_sha1": EXPECTED_SHAPE_SUMMARY_BLOB,
        "shape_verdict": shape["result"]["verdict"],
        "frozen_119hz_result_git_blob_sha1": EXPECTED_FROZEN_119_BLOB,
        "frozen_119hz_verdict": frozen["summary"]["verdict"],
        "spectral_v1_helper_git_blob_sha1": EXPECTED_SPECTRAL_V1_BLOB,
        "bounded_transport_v2_git_blob_sha1": EXPECTED_TRANSPORT_V2_BLOB,
        "authority_delta_for_119hz": 0,
        "target_identity": "UNCONFIRMED",
        "target_evidence_delta": "NONE_FROM_TEMPORAL_STABILITY_DIAGNOSTIC",
        "source_writeback": False,
        "frequency_band_hz": band,
        "network_budget": {
            "attempts_per_request": transport_v2.HTTP_ATTEMPTS,
            "timeout_per_attempt_s": transport_v2.HTTP_TIMEOUT_S,
            "backoff_s": transport_v2.HTTP_BACKOFF_S,
        },
        "classification_thresholds": {
            "large_asymmetry_floor_db": floor_db,
            "maximum_stable_iqr_db": max_iqr_db,
            "minimum_fraction_windows_above_floor": minimum_fraction,
            "minimum_complete_paired_windows": minimum_pairs,
        },
        "inventory_metadata": inventory_meta,
        "inventory_errors": inventory_errors,
        "windows": windows,
        "statistics": statistics,
        "summary": {
            "frozen_window_count": len(wc["dates"]),
            "complete_paired_windows": len(complete),
            "blocked_station_windows": len(data_errors),
            "verdict": verdict,
        },
        "data_errors": data_errors,
        "claim_ceiling": protocol["claim_ceiling"],
        "hard_rules": protocol["hard_rules"],
        "status": "DIAGNOSTIC_RUN_COMPLETE",
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"NO_OVERWRITE:{output}")
    output.write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return receipt


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    receipt = run(args.output)
    print(json.dumps(receipt["summary"] | {"statistics": receipt["statistics"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
