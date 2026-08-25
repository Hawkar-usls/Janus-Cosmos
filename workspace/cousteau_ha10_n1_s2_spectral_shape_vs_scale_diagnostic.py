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
PROTOCOL = DATA / "JANUS-ECHO-COUSTEAU-HA10-N1-S2-SPECTRAL-SHAPE-VS-SCALE-DIAGNOSTIC-PROTOCOL-2026-08-22-v1.0.json"
REPLICATION = DATA / "JANUS-ECHO-COUSTEAU-HA10-N1-S2-INDEPENDENT-NOISE-ASYMMETRY-REPLICATION-RUN-001-SUMMARY-2026-08-22-v1.0.json"
FROZEN_119 = DATA / "JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0.json"
SPECTRAL_V1 = ROOT / "workspace" / "cousteau_ha10_tphase_inband_positive_control.py"
TRANSPORT_V2 = ROOT / "workspace" / "cousteau_ha10_tphase_inband_positive_control_v2.py"

EXPECTED_PROTOCOL_BLOB = "1dd7a9135d74f7b2ed262b761543f40565db1dea"
EXPECTED_REPLICATION_BLOB = "8add6bdb329f4758a67c0a1525dfd281ec46de61"
EXPECTED_FROZEN_119_BLOB = "eb8b48fb7f043160c057f9df6264a781412ed854"
EXPECTED_SPECTRAL_V1_BLOB = "244169a1331a16529e2a963585a5659dae109c66"
EXPECTED_TRANSPORT_V2_BLOB = "c3d37fe4f5a75514390862be4c1d2870c4c01fe9"
EXPECTED_REPLICATION_VERDICT = "REPLICATED_LARGE_H10S2_BASELINE_POWER_ASYMMETRY"
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
    resolved = path.resolve()
    try:
        resolved.relative_to(DATA.resolve())
    except ValueError:
        return
    raise RuntimeError("CANONICAL_DATA_WRITE_FORBIDDEN_BY_SPECTRAL_SHAPE_DIAGNOSTIC")


def verify_frozen_contracts() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    expected = {
        PROTOCOL: EXPECTED_PROTOCOL_BLOB,
        REPLICATION: EXPECTED_REPLICATION_BLOB,
        FROZEN_119: EXPECTED_FROZEN_119_BLOB,
        SPECTRAL_V1: EXPECTED_SPECTRAL_V1_BLOB,
        TRANSPORT_V2: EXPECTED_TRANSPORT_V2_BLOB,
    }
    for path, blob in expected.items():
        if git_blob_sha1_file(path) != blob:
            raise RuntimeError(f"FROZEN_BLOB_DRIFT:{path.name}")

    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    replication = json.loads(REPLICATION.read_text(encoding="utf-8"))
    frozen_119 = json.loads(FROZEN_119.read_text(encoding="utf-8"))

    if protocol.get("status") != "PREREGISTERED_BEFORE_NEW_DIAGNOSTIC_WINDOWS_ARE_DOWNLOADED":
        raise RuntimeError("SPECTRAL_SHAPE_PROTOCOL_NOT_PREREGISTERED")
    if replication.get("result", {}).get("verdict") != EXPECTED_REPLICATION_VERDICT:
        raise RuntimeError("REPLICATION_VERDICT_BINDING_DRIFT")
    if frozen_119.get("summary", {}).get("verdict") != EXPECTED_FROZEN_119_VERDICT:
        raise RuntimeError("FROZEN_119_VERDICT_DRIFT")
    if protocol.get("epistemic_position", {}).get("authority_delta_for_119hz") != 0:
        raise RuntimeError("AUTHORITY_DELTA_FOR_119HZ_NOT_ZERO")
    if protocol.get("frequency_contract", {}).get("119hz_or_117_121hz_bins_may_be_used") is not False:
        raise RuntimeError("TARGET_BAND_EXCLUSION_NOT_FROZEN")
    if len(protocol["new_window_contract"]["dates"]) != 20:
        raise RuntimeError("NEW_WINDOW_COUNT_DRIFT")
    if protocol["new_window_contract"]["window_start_time_utc_each_day"] != "00:00:00Z":
        raise RuntimeError("NEW_WINDOW_CLOCK_DRIFT")
    return protocol, replication, frozen_119


def analyze_trace(
    trace,
    inventory,
    *,
    bands: dict[str, list[float]],
    nperseg: int,
    noverlap: int,
) -> dict[str, Any]:
    frequencies, raw_psd, corrected_psd = spectral_v1.psd_with_response(
        trace,
        inventory,
        nperseg=nperseg,
        noverlap=noverlap,
    )
    result: dict[str, Any] = {}
    for name, band in bands.items():
        raw_power, raw_bins = spectral_v1.integrated_band_power(
            frequencies, raw_psd, band
        )
        corrected_power, corrected_bins = spectral_v1.integrated_band_power(
            frequencies, corrected_psd, band
        )
        result[name] = {
            "band_hz": band,
            "raw_integrated_power": raw_power,
            "corrected_integrated_power_pa2": corrected_power,
            "raw_bin_count": raw_bins,
            "corrected_bin_count": corrected_bins,
        }
    return result


def summarize_complete_pairs(
    complete_pairs: list[dict[str, Any]],
    *,
    subband_names: list[str],
) -> dict[str, Any]:
    medians: dict[str, Any] = {}
    raw_values: list[float] = []
    corrected_values: list[float] = []
    for name in subband_names:
        raw = [float(row["subband_ratios_db"][name]["raw_s2_minus_n1_db"]) for row in complete_pairs]
        corrected = [
            float(row["subband_ratios_db"][name]["corrected_s2_minus_n1_db"])
            for row in complete_pairs
        ]
        raw_median = float(np.median(raw))
        corrected_median = float(np.median(corrected))
        raw_values.append(raw_median)
        corrected_values.append(corrected_median)
        medians[name] = {
            "median_raw_s2_minus_n1_db": raw_median,
            "median_corrected_s2_minus_n1_db": corrected_median,
            "median_corrected_minus_raw_ratio_delta_db": float(
                np.median(np.asarray(corrected) - np.asarray(raw))
            ),
        }
    return {
        "subband_medians": medians,
        "raw_median_ratio_spread_db": float(max(raw_values) - min(raw_values)),
        "corrected_median_ratio_spread_db": float(
            max(corrected_values) - min(corrected_values)
        ),
    }


def decide(
    complete_pairs: list[dict[str, Any]],
    *,
    subband_names: list[str],
    minimum_pairs: int,
    floor_db: float,
    max_spread_db: float,
) -> tuple[str, dict[str, Any]]:
    if len(complete_pairs) < minimum_pairs:
        return "BLOCKED_SPECTRAL_SHAPE_DIAGNOSTIC_DATA_ACCESS", {
            "complete_paired_windows": len(complete_pairs),
            "subband_medians": None,
            "raw_median_ratio_spread_db": None,
            "corrected_median_ratio_spread_db": None,
        }

    aggregate = summarize_complete_pairs(
        complete_pairs,
        subband_names=subband_names,
    )
    medians = aggregate["subband_medians"]
    all_above_floor = all(
        medians[name]["median_raw_s2_minus_n1_db"] >= floor_db
        and medians[name]["median_corrected_s2_minus_n1_db"] >= floor_db
        for name in subband_names
    )
    spread_ok = (
        aggregate["raw_median_ratio_spread_db"] <= max_spread_db
        and aggregate["corrected_median_ratio_spread_db"] <= max_spread_db
    )
    aggregate.update(
        {
            "complete_paired_windows": len(complete_pairs),
            "all_raw_and_corrected_subband_medians_ge_floor": all_above_floor,
            "raw_and_corrected_spreads_le_maximum": spread_ok,
        }
    )
    verdict = (
        "BROADBAND_SCALE_LIKE_H10S2_ASYMMETRY"
        if all_above_floor and spread_ok
        else "FREQUENCY_SELECTIVE_OR_NONUNIFORM_H10S2_ASYMMETRY"
    )
    return verdict, aggregate


def extract_sensitivity(inventory, cid: str, when: str) -> dict[str, Any] | None:
    try:
        response = inventory.get_response(cid, spectral_v1.UTCDateTime(when))
        sensitivity = response.instrument_sensitivity
        if sensitivity is None:
            return None
        return {
            "value": float(sensitivity.value),
            "frequency_hz": float(sensitivity.frequency),
            "input_units": sensitivity.input_units,
            "output_units": sensitivity.output_units,
        }
    except Exception:
        return None


def run(output: Path) -> dict[str, Any]:
    ensure_noncanonical_output(output)
    protocol, replication, frozen_119 = verify_frozen_contracts()

    windows_contract = protocol["new_window_contract"]
    estimator = protocol["spectral_estimator"]
    frequency = protocol["frequency_contract"]
    subbands = {
        str(name): [float(x) for x in band]
        for name, band in frequency["fixed_subbands_hz"].items()
    }
    bands = {"broadband": [float(x) for x in frequency["broadband_hz"]], **subbands}
    subband_names = list(subbands)
    duration_s = int(windows_contract["window_duration_s"])
    nperseg = int(estimator["welch_nperseg_samples"])
    noverlap = int(estimator["welch_noverlap_samples"])
    minimum_pairs = int(windows_contract["minimum_complete_paired_windows"])
    floor_db = float(
        protocol["classification_thresholds"]["large_asymmetry_floor_db_each_subband"]
    )
    max_spread_db = float(
        protocol["classification_thresholds"]["maximum_scale_like_spread_db"]
    )
    channels = list(protocol["channels"])
    if channels != ["IM.H10N1..EDH", "IM.H10S2..EDH"]:
        raise RuntimeError("CHANNEL_SET_DRIFT")

    session = requests.Session()
    session.headers["User-Agent"] = (
        "Janus-Echo-Cousteau/1.0 preregistered H10 spectral shape-vs-scale diagnostic"
    )

    first_start = parse_start(
        windows_contract["dates"][0],
        windows_contract["window_start_time_utc_each_day"],
    )
    last_start = parse_start(
        windows_contract["dates"][-1],
        windows_contract["window_start_time_utc_each_day"],
    )
    inv_start = (first_start - timedelta(days=1)).isoformat().replace("+00:00", "Z")
    inv_end = (last_start + timedelta(days=2)).isoformat().replace("+00:00", "Z")

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

    rows: list[dict[str, Any]] = []
    complete_pairs: list[dict[str, Any]] = []
    data_errors: list[dict[str, Any]] = []

    for date_text in windows_contract["dates"]:
        start_dt = parse_start(
            date_text, windows_contract["window_start_time_utc_each_day"]
        )
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
                trace, waveform_meta = transport_v2.fetch_trace_bounded(
                    session, cid, start, end
                )
                band_power = analyze_trace(
                    trace,
                    inventories[cid],
                    bands=bands,
                    nperseg=nperseg,
                    noverlap=noverlap,
                )
                row["stations"][cid] = {
                    "data_status": "ANALYZED",
                    "waveform": waveform_meta,
                    "band_power": band_power,
                    "instrument_sensitivity_diagnostic": extract_sensitivity(
                        inventories[cid], cid, start
                    ),
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
            row["subband_ratios_db"] = {}
            for name in subband_names:
                n1_band = n1["band_power"][name]
                s2_band = s2["band_power"][name]
                row["subband_ratios_db"][name] = {
                    "raw_s2_minus_n1_db": ratio_db(
                        float(s2_band["raw_integrated_power"]),
                        float(n1_band["raw_integrated_power"]),
                    ),
                    "corrected_s2_minus_n1_db": ratio_db(
                        float(s2_band["corrected_integrated_power_pa2"]),
                        float(n1_band["corrected_integrated_power_pa2"]),
                    ),
                }
            complete_pairs.append(row)
        rows.append(row)

    verdict, aggregate = decide(
        complete_pairs,
        subband_names=subband_names,
        minimum_pairs=minimum_pairs,
        floor_db=floor_db,
        max_spread_db=max_spread_db,
    )

    receipt = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-HA10-N1-S2-SPECTRAL-SHAPE-VS-SCALE-DIAGNOSTIC-RUN",
        "created_utc": utc_now(),
        "gate_id": protocol["gate_id"],
        "protocol_git_blob_sha1": EXPECTED_PROTOCOL_BLOB,
        "replication_summary_git_blob_sha1": EXPECTED_REPLICATION_BLOB,
        "replication_verdict": replication["result"]["verdict"],
        "frozen_119hz_result_git_blob_sha1": EXPECTED_FROZEN_119_BLOB,
        "frozen_119hz_verdict": frozen_119["summary"]["verdict"],
        "spectral_v1_helper_git_blob_sha1": EXPECTED_SPECTRAL_V1_BLOB,
        "bounded_transport_v2_git_blob_sha1": EXPECTED_TRANSPORT_V2_BLOB,
        "authority_delta_for_119hz": 0,
        "target_identity": "UNCONFIRMED",
        "target_evidence_delta": "NONE_FROM_SPECTRAL_SHAPE_DIAGNOSTIC",
        "source_writeback": False,
        "network_budget": {
            "attempts_per_request": transport_v2.HTTP_ATTEMPTS,
            "timeout_per_attempt_s": transport_v2.HTTP_TIMEOUT_S,
            "backoff_s": transport_v2.HTTP_BACKOFF_S,
        },
        "frequency_contract": {
            "broadband_hz": bands["broadband"],
            "fixed_subbands_hz": subbands,
            "best_subband_selection_performed": False,
            "whole_spectrum_peak_search_performed": False,
            "119hz_or_117_121hz_bins_used": False,
        },
        "classification_thresholds": {
            "large_asymmetry_floor_db_each_subband": floor_db,
            "maximum_scale_like_spread_db": max_spread_db,
            "minimum_complete_paired_windows": minimum_pairs,
        },
        "inventory_metadata": inventory_meta,
        "inventory_errors": inventory_errors,
        "windows": rows,
        "aggregate": aggregate,
        "summary": {
            "frozen_window_count": len(windows_contract["dates"]),
            "complete_paired_windows": len(complete_pairs),
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
    print(
        json.dumps(
            receipt["summary"] | receipt["aggregate"],
            indent=2,
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
