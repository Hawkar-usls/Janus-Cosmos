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
PROTOCOL = DATA / "JANUS-ECHO-COUSTEAU-HA10-N1-S2-INDEPENDENT-SCALAR-COLLAPSE-HOLDOUT-PROTOCOL-2026-08-22-v1.0.json"
TEMPORAL_SUMMARY = DATA / "JANUS-ECHO-COUSTEAU-HA10-N1-S2-TEMPORAL-SCALE-STABILITY-DIAGNOSTIC-RUN-001-SUMMARY-2026-08-22-v1.0.json"
FROZEN_119 = DATA / "JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0.json"
SPECTRAL_V1 = ROOT / "workspace" / "cousteau_ha10_tphase_inband_positive_control.py"
TRANSPORT_V2 = ROOT / "workspace" / "cousteau_ha10_tphase_inband_positive_control_v2.py"

EXPECTED_PROTOCOL_BLOB = "c642837c2f1efa77699c76260fe6f158f8fb3d1c"
EXPECTED_TEMPORAL_SUMMARY_BLOB = "45f2901adde866df610a563f27c50177d32be414"
EXPECTED_FROZEN_119_BLOB = "eb8b48fb7f043160c057f9df6264a781412ed854"
EXPECTED_SPECTRAL_V1_BLOB = "244169a1331a16529e2a963585a5659dae109c66"
EXPECTED_TRANSPORT_V2_BLOB = "c3d37fe4f5a75514390862be4c1d2870c4c01fe9"
EXPECTED_TEMPORAL_VERDICT = "TEMPORALLY_STABLE_BROADBAND_H10S2_SCALE_OFFSET"
EXPECTED_FROZEN_119_VERDICT = "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE"
EXPECTED_SCALAR_DB = 12.616188132659623


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
    raise RuntimeError("CANONICAL_DATA_WRITE_FORBIDDEN_BY_SCALAR_COLLAPSE_HOLDOUT")


def verify_frozen_contracts() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    expected = {
        PROTOCOL: EXPECTED_PROTOCOL_BLOB,
        TEMPORAL_SUMMARY: EXPECTED_TEMPORAL_SUMMARY_BLOB,
        FROZEN_119: EXPECTED_FROZEN_119_BLOB,
        SPECTRAL_V1: EXPECTED_SPECTRAL_V1_BLOB,
        TRANSPORT_V2: EXPECTED_TRANSPORT_V2_BLOB,
    }
    for path, blob in expected.items():
        if git_blob_sha1_file(path) != blob:
            raise RuntimeError(f"FROZEN_BLOB_DRIFT:{path.name}")

    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    temporal = json.loads(TEMPORAL_SUMMARY.read_text(encoding="utf-8"))
    frozen = json.loads(FROZEN_119.read_text(encoding="utf-8"))

    if protocol.get("status") != "PREREGISTERED_BEFORE_NEW_1800UTC_HOLDOUT_WINDOWS_ARE_DOWNLOADED":
        raise RuntimeError("SCALAR_COLLAPSE_PROTOCOL_NOT_PREREGISTERED")
    if temporal.get("result", {}).get("verdict") != EXPECTED_TEMPORAL_VERDICT:
        raise RuntimeError("TEMPORAL_STABILITY_VERDICT_BINDING_DRIFT")
    if frozen.get("summary", {}).get("verdict") != EXPECTED_FROZEN_119_VERDICT:
        raise RuntimeError("FROZEN_119_VERDICT_DRIFT")
    if protocol.get("epistemic_position", {}).get("authority_delta_for_119hz") != 0:
        raise RuntimeError("AUTHORITY_DELTA_FOR_119HZ_NOT_ZERO")
    if protocol.get("frequency_contract", {}).get("119hz_or_117_121hz_bins_may_be_used") is not False:
        raise RuntimeError("TARGET_BAND_EXCLUSION_NOT_FROZEN")
    scalar = float(protocol["frozen_scalar"]["source_value_db"])
    upstream_scalar = float(temporal["result"]["median_corrected_s2_minus_n1_db"])
    if scalar != EXPECTED_SCALAR_DB or scalar != upstream_scalar:
        raise RuntimeError("FROZEN_SCALAR_BINDING_DRIFT")
    if protocol["frozen_scalar"]["holdout_refit_allowed"] is not False:
        raise RuntimeError("HOLDOUT_REFIT_MUST_BE_FORBIDDEN")
    if protocol["frozen_scalar"]["alternate_scalar_search_allowed"] is not False:
        raise RuntimeError("ALTERNATE_SCALAR_SEARCH_MUST_BE_FORBIDDEN")
    wc = protocol["holdout_window_contract"]
    if len(wc["dates"]) != 20:
        raise RuntimeError("HOLDOUT_WINDOW_COUNT_DRIFT")
    if wc["window_start_time_utc_each_day"] != "18:00:00Z":
        raise RuntimeError("HOLDOUT_WINDOW_CLOCK_DRIFT")
    return protocol, temporal, frozen


def powers_for_bands(
    trace,
    inventory,
    *,
    bands: dict[str, list[float]],
    nperseg: int,
    noverlap: int,
) -> dict[str, dict[str, float]]:
    frequencies, raw_psd, corrected_psd = spectral_v1.psd_with_response(
        trace,
        inventory,
        nperseg=nperseg,
        noverlap=noverlap,
    )
    out: dict[str, dict[str, float]] = {}
    for name, band in bands.items():
        raw_power, raw_bins = spectral_v1.integrated_band_power(frequencies, raw_psd, band)
        corrected_power, corrected_bins = spectral_v1.integrated_band_power(
            frequencies, corrected_psd, band
        )
        out[name] = {
            "band_hz": [float(band[0]), float(band[1])],
            "raw_integrated_power": raw_power,
            "corrected_integrated_power_pa2": corrected_power,
            "raw_bin_count": int(raw_bins),
            "corrected_bin_count": int(corrected_bins),
        }
    return out


def residual_statistics(
    rows: list[dict[str, Any]],
    *,
    tolerance_db: float,
) -> dict[str, Any]:
    raw = np.asarray([float(row["broadband"]["raw_residual_db"]) for row in rows], dtype=float)
    corrected = np.asarray(
        [float(row["broadband"]["corrected_residual_db"]) for row in rows],
        dtype=float,
    )
    raw_q25, raw_q75 = np.percentile(raw, [25, 75], method="linear")
    corrected_q25, corrected_q75 = np.percentile(corrected, [25, 75], method="linear")
    raw_within = np.abs(raw) <= tolerance_db
    corrected_within = np.abs(corrected) <= tolerance_db
    both_within = raw_within & corrected_within
    return {
        "percentile_method": "numpy_linear",
        "median_raw_residual_db": float(np.median(raw)),
        "median_corrected_residual_db": float(np.median(corrected)),
        "raw_q25_residual_db": float(raw_q25),
        "raw_q75_residual_db": float(raw_q75),
        "corrected_q25_residual_db": float(corrected_q25),
        "corrected_q75_residual_db": float(corrected_q75),
        "raw_residual_iqr_db": float(raw_q75 - raw_q25),
        "corrected_residual_iqr_db": float(corrected_q75 - corrected_q25),
        "fraction_raw_windows_within_tolerance": float(np.mean(raw_within)),
        "fraction_corrected_windows_within_tolerance": float(np.mean(corrected_within)),
        "fraction_both_windows_within_tolerance": float(np.mean(both_within)),
        "raw_residual_min_db": float(np.min(raw)),
        "raw_residual_max_db": float(np.max(raw)),
        "corrected_residual_min_db": float(np.min(corrected)),
        "corrected_residual_max_db": float(np.max(corrected)),
    }


def fixed_subband_statistics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    names = ["low", "mid", "high"]
    out: dict[str, Any] = {}
    for name in names:
        raw = np.asarray(
            [float(row["subbands"][name]["raw_residual_db"]) for row in rows],
            dtype=float,
        )
        corrected = np.asarray(
            [float(row["subbands"][name]["corrected_residual_db"]) for row in rows],
            dtype=float,
        )
        out[name] = {
            "band_hz": rows[0]["subbands"][name]["band_hz"],
            "median_raw_residual_db": float(np.median(raw)),
            "median_corrected_residual_db": float(np.median(corrected)),
        }
    return out


def decide(
    rows: list[dict[str, Any]],
    *,
    minimum_pairs: int,
    max_abs_median_db: float,
    max_iqr_db: float,
    tolerance_db: float,
    minimum_fraction: float,
    max_abs_subband_median_db: float,
) -> tuple[str, dict[str, Any] | None, dict[str, Any] | None, dict[str, bool]]:
    if len(rows) < minimum_pairs:
        return (
            "BLOCKED_SCALAR_COLLAPSE_HOLDOUT_DATA_ACCESS",
            None,
            None,
            {
                "minimum_complete_pairs": False,
                "broadband_medians": False,
                "broadband_iqrs": False,
                "window_fractions": False,
                "fixed_subbands": False,
            },
        )

    broadband = residual_statistics(rows, tolerance_db=tolerance_db)
    subbands = fixed_subband_statistics(rows)
    median_ok = (
        abs(broadband["median_raw_residual_db"]) <= max_abs_median_db
        and abs(broadband["median_corrected_residual_db"]) <= max_abs_median_db
    )
    iqr_ok = (
        broadband["raw_residual_iqr_db"] <= max_iqr_db
        and broadband["corrected_residual_iqr_db"] <= max_iqr_db
    )
    fractions_ok = (
        broadband["fraction_raw_windows_within_tolerance"] >= minimum_fraction
        and broadband["fraction_corrected_windows_within_tolerance"] >= minimum_fraction
        and broadband["fraction_both_windows_within_tolerance"] >= minimum_fraction
    )
    subbands_ok = all(
        abs(float(row["median_raw_residual_db"])) <= max_abs_subband_median_db
        and abs(float(row["median_corrected_residual_db"])) <= max_abs_subband_median_db
        for row in subbands.values()
    )
    gates = {
        "minimum_complete_pairs": True,
        "broadband_medians": bool(median_ok),
        "broadband_iqrs": bool(iqr_ok),
        "window_fractions": bool(fractions_ok),
        "fixed_subbands": bool(subbands_ok),
    }
    verdict = (
        "SCALAR_COLLAPSE_HOLDOUT_PASS"
        if all(gates.values())
        else "SCALAR_COLLAPSE_HOLDOUT_FAIL"
    )
    return verdict, broadband, subbands, gates


def run(output: Path) -> dict[str, Any]:
    ensure_noncanonical_output(output)
    protocol, temporal, frozen = verify_frozen_contracts()

    wc = protocol["holdout_window_contract"]
    estimator = protocol["spectral_estimator"]
    thresholds = protocol["classification_thresholds"]
    scalar_db = float(protocol["frozen_scalar"]["source_value_db"])
    nperseg = int(estimator["welch_nperseg_samples"])
    noverlap = int(estimator["welch_noverlap_samples"])
    duration_s = int(wc["window_duration_s"])
    minimum_pairs = int(wc["minimum_complete_paired_windows"])
    max_abs_median_db = float(thresholds["maximum_absolute_broadband_residual_median_db"])
    max_iqr_db = float(thresholds["maximum_broadband_residual_iqr_db"])
    tolerance_db = float(thresholds["window_absolute_residual_tolerance_db"])
    minimum_fraction = float(thresholds["minimum_fraction_windows_within_tolerance"])
    max_abs_subband_median_db = float(
        thresholds["maximum_absolute_fixed_subband_residual_median_db"]
    )
    channels = list(protocol["channels"])
    if channels != ["IM.H10N1..EDH", "IM.H10S2..EDH"]:
        raise RuntimeError("CHANNEL_SET_DRIFT")

    bands: dict[str, list[float]] = {
        "broadband": [float(x) for x in protocol["frequency_contract"]["broadband_hz"]],
        **{
            name: [float(x) for x in band]
            for name, band in protocol["frequency_contract"]["fixed_subbands_hz"].items()
        },
    }

    session = requests.Session()
    session.headers["User-Agent"] = (
        "Janus-Echo-Cousteau/1.0 preregistered H10 independent scalar-collapse holdout"
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
                powers = powers_for_bands(
                    trace,
                    inventories[cid],
                    bands=bands,
                    nperseg=nperseg,
                    noverlap=noverlap,
                )
                row["stations"][cid] = {
                    "data_status": "ANALYZED",
                    "waveform": waveform,
                    "powers": powers,
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
            row["broadband"] = {}
            row["subbands"] = {}
            for name in ["broadband", "low", "mid", "high"]:
                n1p = n1["powers"][name]
                s2p = s2["powers"][name]
                raw_ratio = ratio_db(
                    float(s2p["raw_integrated_power"]),
                    float(n1p["raw_integrated_power"]),
                )
                corrected_ratio = ratio_db(
                    float(s2p["corrected_integrated_power_pa2"]),
                    float(n1p["corrected_integrated_power_pa2"]),
                )
                target = row["broadband"] if name == "broadband" else row["subbands"].setdefault(name, {})
                target.update(
                    {
                        "band_hz": s2p["band_hz"],
                        "raw_uncollapsed_s2_minus_n1_db": raw_ratio,
                        "corrected_uncollapsed_s2_minus_n1_db": corrected_ratio,
                        "raw_residual_db": raw_ratio - scalar_db,
                        "corrected_residual_db": corrected_ratio - scalar_db,
                    }
                )
            complete.append(row)
        windows.append(row)

    verdict, broadband_stats, subband_stats, gates = decide(
        complete,
        minimum_pairs=minimum_pairs,
        max_abs_median_db=max_abs_median_db,
        max_iqr_db=max_iqr_db,
        tolerance_db=tolerance_db,
        minimum_fraction=minimum_fraction,
        max_abs_subband_median_db=max_abs_subband_median_db,
    )

    receipt = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-HA10-N1-S2-INDEPENDENT-SCALAR-COLLAPSE-HOLDOUT-RUN",
        "created_utc": utc_now(),
        "status": "DIAGNOSTIC_RUN_COMPLETE",
        "gate_id": protocol["gate_id"],
        "protocol_git_blob_sha1": EXPECTED_PROTOCOL_BLOB,
        "temporal_stability_summary_git_blob_sha1": EXPECTED_TEMPORAL_SUMMARY_BLOB,
        "temporal_stability_verdict": temporal["result"]["verdict"],
        "frozen_119hz_result_git_blob_sha1": EXPECTED_FROZEN_119_BLOB,
        "frozen_119hz_verdict": frozen["summary"]["verdict"],
        "spectral_v1_helper_git_blob_sha1": EXPECTED_SPECTRAL_V1_BLOB,
        "bounded_transport_v2_git_blob_sha1": EXPECTED_TRANSPORT_V2_BLOB,
        "authority_delta_for_119hz": 0,
        "target_identity": "UNCONFIRMED",
        "target_evidence_delta": "NONE_FROM_SCALAR_COLLAPSE_HOLDOUT",
        "source_writeback": False,
        "frozen_scalar": protocol["frozen_scalar"],
        "frequency_contract": protocol["frequency_contract"],
        "classification_thresholds": thresholds,
        "network_budget": {
            "attempts_per_request": transport_v2.HTTP_ATTEMPTS,
            "timeout_per_attempt_s": transport_v2.HTTP_TIMEOUT_S,
            "backoff_s": transport_v2.HTTP_BACKOFF_S,
        },
        "inventory": {
            "request_start_utc": inv_start,
            "request_end_utc": inv_end,
            "metadata": inventory_meta,
            "errors": inventory_errors,
        },
        "holdout_window_contract": wc,
        "windows": windows,
        "data_errors": data_errors,
        "summary": {
            "frozen_window_count": len(wc["dates"]),
            "complete_paired_windows": len(complete),
            "blocked_station_windows": len(data_errors),
            "verdict": verdict,
            "broadband_residual_statistics": broadband_stats,
            "fixed_subband_residual_statistics": subband_stats,
            "gate_checks": gates,
        },
        "claim_ceiling": protocol["claim_ceiling"],
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = run(args.output)
    summary = receipt["summary"]
    print(
        json.dumps(
            {
                "frozen_window_count": summary["frozen_window_count"],
                "complete_paired_windows": summary["complete_paired_windows"],
                "blocked_station_windows": summary["blocked_station_windows"],
                "verdict": summary["verdict"],
                "frozen_scalar_db": receipt["frozen_scalar"]["source_value_db"],
                "broadband_residual_statistics": summary["broadband_residual_statistics"],
                "fixed_subband_residual_statistics": summary["fixed_subband_residual_statistics"],
                "gate_checks": summary["gate_checks"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
