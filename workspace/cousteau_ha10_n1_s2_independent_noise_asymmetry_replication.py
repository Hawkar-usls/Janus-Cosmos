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

from workspace import cousteau_ha10_tphase_inband_positive_control as spectral

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cousteau"
PROTOCOL = DATA / "JANUS-ECHO-COUSTEAU-HA10-N1-S2-INDEPENDENT-NOISE-ASYMMETRY-REPLICATION-PROTOCOL-2026-08-22-v1.0.json"
SUMMARY = DATA / "JANUS-ECHO-COUSTEAU-HA10-TPHASE-INBAND-POSITIVE-CONTROL-RUN-001-SUMMARY-2026-08-22-v1.0.json"
FROZEN_119 = DATA / "JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0.json"
SPECTRAL_HELPER = ROOT / "workspace" / "cousteau_ha10_tphase_inband_positive_control.py"

EXPECTED_PROTOCOL_BLOB = "70e58ca56e202639c7e2d4d647e3ffb70bc53495"
EXPECTED_SUMMARY_BLOB = "848c53633fe662e37df615626600953267119eff"
EXPECTED_FROZEN_119_BLOB = "eb8b48fb7f043160c057f9df6264a781412ed854"
EXPECTED_SPECTRAL_HELPER_BLOB = "244169a1331a16529e2a963585a5659dae109c66"
EXPECTED_POSITIVE_CONTROL_VERDICT = "FAIL_HA10_INBAND_TPHASE_PIPELINE_CONTROL"
EXPECTED_FROZEN_119_VERDICT = "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE"


def git_blob_sha1_file(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def parse_start(date_text: str, hhmmss: str) -> datetime:
    return datetime.fromisoformat(f"{date_text}T{hhmmss.replace('Z', '+00:00')}")


def ratio_db(numerator: float, denominator: float) -> float:
    if not (
        numerator > 0
        and denominator > 0
        and math.isfinite(numerator)
        and math.isfinite(denominator)
    ):
        return float("nan")
    return 10.0 * math.log10(numerator / denominator)


def ensure_noncanonical_output(path: Path) -> None:
    resolved = path.resolve()
    try:
        resolved.relative_to(DATA.resolve())
    except ValueError:
        return
    raise RuntimeError("CANONICAL_DATA_WRITE_FORBIDDEN_BY_NOISE_REPLICATION_RUNNER")


def verify_frozen_contracts() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    expected = {
        PROTOCOL: EXPECTED_PROTOCOL_BLOB,
        SUMMARY: EXPECTED_SUMMARY_BLOB,
        FROZEN_119: EXPECTED_FROZEN_119_BLOB,
        SPECTRAL_HELPER: EXPECTED_SPECTRAL_HELPER_BLOB,
    }
    for path, blob in expected.items():
        if git_blob_sha1_file(path) != blob:
            raise RuntimeError(f"FROZEN_BLOB_DRIFT:{path.name}")

    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    frozen_119 = json.loads(FROZEN_119.read_text(encoding="utf-8"))

    if protocol.get("status") != "PREREGISTERED_BEFORE_INDEPENDENT_NOISE_WINDOWS_ARE_DOWNLOADED":
        raise RuntimeError("NOISE_REPLICATION_PROTOCOL_NOT_PREREGISTERED")
    if summary.get("control_summary", {}).get("verdict") != EXPECTED_POSITIVE_CONTROL_VERDICT:
        raise RuntimeError("POSITIVE_CONTROL_FAILURE_BINDING_DRIFT")
    if frozen_119.get("summary", {}).get("verdict") != EXPECTED_FROZEN_119_VERDICT:
        raise RuntimeError("FROZEN_119_VERDICT_DRIFT")
    if protocol.get("epistemic_position", {}).get("authority_delta_for_119hz") != 0:
        raise RuntimeError("AUTHORITY_DELTA_FOR_119HZ_NOT_ZERO")
    if protocol.get("frequency_contract", {}).get("119hz_or_117_121hz_bins_may_be_used") is not False:
        raise RuntimeError("TARGET_BAND_EXCLUSION_NOT_FROZEN")
    if len(protocol["independent_window_contract"]["dates"]) != 20:
        raise RuntimeError("INDEPENDENT_WINDOW_COUNT_DRIFT")
    return protocol, summary, frozen_119


def analyze_trace(trace, inventory, *, band: list[float], nperseg: int, noverlap: int):
    frequencies, raw_psd, corrected_psd = spectral.psd_with_response(
        trace,
        inventory,
        nperseg=nperseg,
        noverlap=noverlap,
    )
    raw_power, raw_bins = spectral.integrated_band_power(frequencies, raw_psd, band)
    corrected_power, corrected_bins = spectral.integrated_band_power(
        frequencies, corrected_psd, band
    )
    return {
        "raw_integrated_power": raw_power,
        "corrected_integrated_power_pa2": corrected_power,
        "raw_bin_count": raw_bins,
        "corrected_bin_count": corrected_bins,
    }


def decide(
    complete_pairs: list[dict[str, Any]],
    *,
    minimum_pairs: int,
    threshold_db: float,
) -> tuple[str, dict[str, Any]]:
    if len(complete_pairs) < minimum_pairs:
        return "BLOCKED_INDEPENDENT_NOISE_REPLICATION_DATA_ACCESS", {
            "complete_paired_windows": len(complete_pairs),
            "median_per_window_raw_s2_minus_n1_db": None,
            "median_per_window_corrected_s2_minus_n1_db": None,
        }

    raw_ratios = [float(row["s2_minus_n1_raw_power_db"]) for row in complete_pairs]
    corrected_ratios = [
        float(row["s2_minus_n1_corrected_power_db"]) for row in complete_pairs
    ]
    raw_median = float(np.median(raw_ratios))
    corrected_median = float(np.median(corrected_ratios))
    verdict = (
        "REPLICATED_LARGE_H10S2_BASELINE_POWER_ASYMMETRY"
        if raw_median >= threshold_db and corrected_median >= threshold_db
        else "NOT_REPLICATED_LARGE_H10S2_BASELINE_POWER_ASYMMETRY"
    )
    return verdict, {
        "complete_paired_windows": len(complete_pairs),
        "median_per_window_raw_s2_minus_n1_db": raw_median,
        "median_per_window_corrected_s2_minus_n1_db": corrected_median,
    }


def run(output: Path) -> dict[str, Any]:
    ensure_noncanonical_output(output)
    protocol, summary, frozen_119 = verify_frozen_contracts()

    window_contract = protocol["independent_window_contract"]
    estimator = protocol["spectral_estimator"]
    band = [float(x) for x in protocol["frequency_contract"]["band_hz"]]
    duration_s = int(window_contract["window_duration_s"])
    nperseg = int(estimator["welch_nperseg_samples"])
    noverlap = int(estimator["welch_noverlap_samples"])
    minimum_pairs = int(window_contract["minimum_complete_paired_windows"])
    threshold_db = float(
        protocol["replication_threshold"]["large_asymmetry_effect_size_db"]
    )
    channels = list(protocol["channels"])
    if channels != ["IM.H10N1..EDH", "IM.H10S2..EDH"]:
        raise RuntimeError("CHANNEL_SET_DRIFT")

    session = requests.Session()
    session.headers["User-Agent"] = (
        "Janus-Echo-Cousteau/1.0 independent H10 N1-S2 noise replication"
    )

    first_start = parse_start(
        window_contract["dates"][0], window_contract["window_start_time_utc_each_day"]
    )
    last_start = parse_start(
        window_contract["dates"][-1], window_contract["window_start_time_utc_each_day"]
    )
    inv_start = (first_start - timedelta(days=1)).isoformat().replace("+00:00", "Z")
    inv_end = (last_start + timedelta(days=2)).isoformat().replace("+00:00", "Z")

    inventories: dict[str, Any] = {}
    inventory_meta: dict[str, Any] = {}
    inventory_errors: dict[str, str] = {}
    for cid in channels:
        try:
            inventory, metadata = spectral.fetch_inventory(
                session, cid, inv_start, inv_end
            )
            inventories[cid] = inventory
            inventory_meta[cid] = metadata
        except Exception as exc:
            inventory_errors[cid] = f"{type(exc).__name__}:{exc}"

    windows: list[dict[str, Any]] = []
    complete_pairs: list[dict[str, Any]] = []
    data_errors: list[dict[str, Any]] = []

    for date_text in window_contract["dates"]:
        start_dt = parse_start(date_text, window_contract["window_start_time_utc_each_day"])
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
            station: dict[str, Any] = {}
            try:
                if cid not in inventories:
                    raise RuntimeError(
                        "RESPONSE_UNAVAILABLE:" + inventory_errors.get(cid, "UNKNOWN")
                    )
                trace, waveform_meta = spectral.fetch_trace(
                    session, cid, start, end
                )
                power = analyze_trace(
                    trace,
                    inventories[cid],
                    band=band,
                    nperseg=nperseg,
                    noverlap=noverlap,
                )
                station = {
                    "data_status": "ANALYZED",
                    "waveform": waveform_meta,
                    **power,
                }
            except Exception as exc:
                station = {
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
            row["stations"][cid] = station

        n1 = row["stations"][channels[0]]
        s2 = row["stations"][channels[1]]
        if n1.get("data_status") == "ANALYZED" and s2.get("data_status") == "ANALYZED":
            row["pair_complete"] = True
            row["s2_minus_n1_raw_power_db"] = ratio_db(
                float(s2["raw_integrated_power"]),
                float(n1["raw_integrated_power"]),
            )
            row["s2_minus_n1_corrected_power_db"] = ratio_db(
                float(s2["corrected_integrated_power_pa2"]),
                float(n1["corrected_integrated_power_pa2"]),
            )
            complete_pairs.append(row)
        windows.append(row)

    verdict, aggregate = decide(
        complete_pairs,
        minimum_pairs=minimum_pairs,
        threshold_db=threshold_db,
    )

    station_aggregate: dict[str, Any] = {}
    if complete_pairs:
        for cid in channels:
            raw = [
                float(row["stations"][cid]["raw_integrated_power"])
                for row in complete_pairs
            ]
            corrected = [
                float(row["stations"][cid]["corrected_integrated_power_pa2"])
                for row in complete_pairs
            ]
            station_aggregate[cid] = {
                "median_raw_integrated_power": float(np.median(raw)),
                "median_corrected_integrated_power_pa2": float(np.median(corrected)),
            }

    receipt = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-HA10-N1-S2-INDEPENDENT-NOISE-ASYMMETRY-REPLICATION-RUN",
        "created_utc": utc_now(),
        "gate_id": protocol["gate_id"],
        "protocol_path": str(PROTOCOL.relative_to(ROOT)),
        "protocol_git_blob_sha1": EXPECTED_PROTOCOL_BLOB,
        "positive_control_summary_git_blob_sha1": EXPECTED_SUMMARY_BLOB,
        "positive_control_verdict": summary["control_summary"]["verdict"],
        "frozen_119hz_result_git_blob_sha1": EXPECTED_FROZEN_119_BLOB,
        "frozen_119hz_verdict": frozen_119["summary"]["verdict"],
        "spectral_helper_git_blob_sha1": EXPECTED_SPECTRAL_HELPER_BLOB,
        "authority_delta_for_119hz": 0,
        "target_identity": "UNCONFIRMED",
        "target_evidence_delta": "NONE_FROM_NOISE_REPLICATION",
        "frequency_band_hz": band,
        "independent_window_count_frozen": len(window_contract["dates"]),
        "minimum_complete_paired_windows": minimum_pairs,
        "large_asymmetry_effect_size_db": threshold_db,
        "inventory_metadata": inventory_meta,
        "inventory_errors": inventory_errors,
        "windows": windows,
        "station_aggregate": station_aggregate,
        "aggregate": aggregate,
        "summary": {
            "frozen_window_count": len(window_contract["dates"]),
            "complete_paired_windows": len(complete_pairs),
            "blocked_station_windows": len(data_errors),
            "verdict": verdict,
        },
        "data_errors": data_errors,
        "claim_ceiling": protocol["claim_ceiling"],
        "hard_rules": protocol["hard_rules"],
        "source_writeback": False,
        "status": "REPLICATION_RUN_COMPLETE",
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
    print(json.dumps(receipt["summary"] | receipt["aggregate"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
