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
from scipy import signal

from workspace import cousteau_ha10_tphase_inband_positive_control as spectral_v1
from workspace import cousteau_ha10_tphase_inband_positive_control_v2 as transport_v2

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cousteau"
PROTOCOL = DATA / "JANUS-ECHO-COUSTEAU-HA10-PRE-POST-2013-CROSSTALK-TRIPLET-GATE-PROTOCOL-2026-08-22-v1.0.json"
FROZEN_119 = DATA / "JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0.json"
SPECTRAL_V1 = ROOT / "workspace" / "cousteau_ha10_tphase_inband_positive_control.py"
TRANSPORT_V2 = ROOT / "workspace" / "cousteau_ha10_tphase_inband_positive_control_v2.py"

EXPECTED_PROTOCOL_BLOB = "c6006a450b7aae4da86e830c4a25c59713e828cb"
EXPECTED_FROZEN_119_BLOB = "eb8b48fb7f043160c057f9df6264a781412ed854"
EXPECTED_SPECTRAL_V1_BLOB = "244169a1331a16529e2a963585a5659dae109c66"
EXPECTED_TRANSPORT_V2_BLOB = "c3d37fe4f5a75514390862be4c1d2870c4c01fe9"
EXPECTED_FROZEN_119_VERDICT = "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE"


def git_blob_sha1_file(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def ratio_db(num: float, den: float) -> float:
    if not (num > 0 and den > 0 and math.isfinite(num) and math.isfinite(den)):
        return float("nan")
    return 10.0 * math.log10(num / den)


def parse_start(date_text: str, clock: str) -> datetime:
    return datetime.fromisoformat(f"{date_text}T{clock.replace('Z', '+00:00')}")


def ensure_noncanonical_output(path: Path) -> None:
    try:
        path.resolve().relative_to(DATA.resolve())
    except ValueError:
        return
    raise RuntimeError("CANONICAL_DATA_WRITE_FORBIDDEN_BY_PRE_POST_TRIPLET_GATE")


def verify_frozen_contracts() -> tuple[dict[str, Any], dict[str, Any]]:
    expected = {
        PROTOCOL: EXPECTED_PROTOCOL_BLOB,
        FROZEN_119: EXPECTED_FROZEN_119_BLOB,
        SPECTRAL_V1: EXPECTED_SPECTRAL_V1_BLOB,
        TRANSPORT_V2: EXPECTED_TRANSPORT_V2_BLOB,
    }
    for path, blob in expected.items():
        if git_blob_sha1_file(path) != blob:
            raise RuntimeError(f"FROZEN_BLOB_DRIFT:{path.name}")
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    frozen = json.loads(FROZEN_119.read_text(encoding="utf-8"))
    if protocol.get("status") != "PREREGISTERED_BEFORE_NEW_PRE_POST_2013_TRIPLET_WINDOWS_ARE_DOWNLOADED":
        raise RuntimeError("PROTOCOL_NOT_PREREGISTERED")
    if protocol.get("gate_id") != "COUSTEAU_HA10_PRE_POST_2013_CROSSTALK_TRIPLET_GATE_V1":
        raise RuntimeError("GATE_ID_DRIFT")
    if protocol.get("authority", {}).get("authority_delta_for_119hz") != 0:
        raise RuntimeError("AUTHORITY_DELTA_FOR_119HZ_NOT_ZERO")
    if protocol.get("frequency_contract", {}).get("119hz_or_117_121hz_bins_may_be_used") is not False:
        raise RuntimeError("TARGET_BAND_EXCLUSION_DRIFT")
    if frozen.get("summary", {}).get("verdict") != EXPECTED_FROZEN_119_VERDICT:
        raise RuntimeError("FROZEN_119_VERDICT_DRIFT")
    return protocol, frozen


def integrated_powers(trace, inventory, band: list[float], nperseg: int, noverlap: int) -> dict[str, float]:
    f, raw_psd, corrected_psd = spectral_v1.psd_with_response(
        trace, inventory, nperseg=nperseg, noverlap=noverlap
    )
    mask = (f >= band[0]) & (f <= band[1])
    if int(np.count_nonzero(mask)) < 2:
        raise RuntimeError("INSUFFICIENT_BAND_BINS")
    return {
        "raw_integrated_power": float(np.trapezoid(raw_psd[mask], f[mask])),
        "corrected_integrated_power_pa2": float(np.trapezoid(corrected_psd[mask], f[mask])),
    }


def _bandpassed(trace, band: list[float]) -> np.ndarray:
    fs = float(trace.stats.sampling_rate)
    data = np.asarray(trace.data, dtype=float)
    sos = signal.butter(4, band, btype="bandpass", fs=fs, output="sos")
    return signal.sosfiltfilt(sos, data)


def pair_coupling(trace_a, trace_b, band: list[float], nperseg: int, noverlap: int) -> dict[str, float]:
    fs_a = float(trace_a.stats.sampling_rate)
    fs_b = float(trace_b.stats.sampling_rate)
    if abs(fs_a - fs_b) > 1e-9:
        raise RuntimeError("PAIR_SAMPLE_RATE_MISMATCH")
    n = min(len(trace_a.data), len(trace_b.data))
    if n < nperseg * 2:
        raise RuntimeError("PAIR_TOO_SHORT")
    xa = np.asarray(trace_a.data[:n], dtype=float)
    xb = np.asarray(trace_b.data[:n], dtype=float)
    f, coh = signal.coherence(
        xa,
        xb,
        fs=fs_a,
        window="hann",
        nperseg=nperseg,
        noverlap=noverlap,
        detrend="constant",
    )
    mask = (f >= band[0]) & (f <= band[1])
    if int(np.count_nonzero(mask)) < 2:
        raise RuntimeError("INSUFFICIENT_COHERENCE_BINS")
    ba = _bandpassed(trace_a, band)[:n]
    bb = _bandpassed(trace_b, band)[:n]
    corr = float(np.corrcoef(ba, bb)[0, 1])
    if not math.isfinite(corr):
        raise RuntimeError("NONFINITE_CORRELATION")
    return {
        "median_magnitude_squared_coherence": float(np.median(coh[mask])),
        "abs_zero_lag_pearson_correlation": abs(corr),
    }


def median(values: list[float]) -> float:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return float("nan")
    return float(np.median(arr))


def epoch_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    def vals(key: str) -> list[float]:
        return [float(r[key]) for r in rows]

    power_keys = [
        "raw_s1_minus_n1_db", "raw_s2_minus_n1_db", "raw_s3_minus_n1_db",
        "corrected_s1_minus_n1_db", "corrected_s2_minus_n1_db", "corrected_s3_minus_n1_db",
    ]
    out: dict[str, Any] = {
        "complete_all_six_windows": len(rows),
        "south_driver_coherence_median": median(vals("south_driver_coherence")),
        "north_driver_coherence_median": median(vals("north_driver_coherence")),
        "south_driver_abs_correlation_median": median(vals("south_driver_abs_correlation")),
        "north_driver_abs_correlation_median": median(vals("north_driver_abs_correlation")),
    }
    for key in power_keys:
        out[key + "_median"] = median(vals(key))
    return out


def decide(pre: list[dict[str, Any]], post: list[dict[str, Any]], protocol: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    th = protocol["predeclared_classification_thresholds"]
    minimum = int(th["minimum_epoch_complete_windows"])
    if len(pre) < minimum or len(post) < minimum:
        return "BLOCKED_PRE_POST_2013_TRIPLET_DATA_ACCESS", {
            "pre_complete": len(pre),
            "post_complete": len(post),
            "minimum_per_epoch": minimum,
        }

    pre_s = epoch_summary(pre)
    post_s = epoch_summary(post)
    coh_south_inc = post_s["south_driver_coherence_median"] - pre_s["south_driver_coherence_median"]
    coh_north_inc = post_s["north_driver_coherence_median"] - pre_s["north_driver_coherence_median"]
    coh_did = coh_south_inc - coh_north_inc
    corr_south_inc = post_s["south_driver_abs_correlation_median"] - pre_s["south_driver_abs_correlation_median"]
    corr_north_inc = post_s["north_driver_abs_correlation_median"] - pre_s["north_driver_abs_correlation_median"]
    corr_did = corr_south_inc - corr_north_inc

    checks = {
        "south_coherence_increase": coh_south_inc >= float(th["minimum_post_minus_pre_south_driver_coherence_increase"]),
        "coherence_difference_in_difference": coh_did >= float(th["minimum_coherence_difference_in_difference_vs_north"]),
        "south_abs_correlation_increase": corr_south_inc >= float(th["minimum_post_minus_pre_south_driver_abs_correlation_increase"]),
        "abs_correlation_difference_in_difference": corr_did >= float(th["minimum_abs_correlation_difference_in_difference_vs_north"]),
    }
    verdict = (
        "SUPPORTED_POST_2013_SOUTH_TRIPLET_COUPLING_INCREASE"
        if all(checks.values())
        else "NO_STRONG_POST_2013_SOUTH_TRIPLET_COUPLING_INCREASE"
    )
    diagnostics = {
        "pre": pre_s,
        "post": post_s,
        "post_minus_pre": {
            "south_driver_coherence": coh_south_inc,
            "north_driver_coherence": coh_north_inc,
            "coherence_difference_in_difference_vs_north": coh_did,
            "south_driver_abs_correlation": corr_south_inc,
            "north_driver_abs_correlation": corr_north_inc,
            "abs_correlation_difference_in_difference_vs_north": corr_did,
            "raw_s2_minus_n1_db": post_s["raw_s2_minus_n1_db_median"] - pre_s["raw_s2_minus_n1_db_median"],
            "corrected_s2_minus_n1_db": post_s["corrected_s2_minus_n1_db_median"] - pre_s["corrected_s2_minus_n1_db_median"],
            "raw_s3_minus_n1_db": post_s["raw_s3_minus_n1_db_median"] - pre_s["raw_s3_minus_n1_db_median"],
            "corrected_s3_minus_n1_db": post_s["corrected_s3_minus_n1_db_median"] - pre_s["corrected_s3_minus_n1_db_median"],
        },
        "gate_checks": checks,
    }
    return verdict, diagnostics


def run(output: Path) -> dict[str, Any]:
    ensure_noncanonical_output(output)
    protocol, frozen = verify_frozen_contracts()
    wc = protocol["window_contract"]
    band = [float(x) for x in protocol["frequency_contract"]["primary_band_hz"]]
    nperseg = int(protocol["spectral_estimator"]["welch_nperseg_samples"])
    noverlap = int(protocol["spectral_estimator"]["welch_noverlap_samples"])
    duration = int(wc["window_duration_s"])
    channels = list(protocol["channels"])
    expected_channels = [
        "IM.H10N1..EDH", "IM.H10N2..EDH", "IM.H10N3..EDH",
        "IM.H10S1..EDH", "IM.H10S2..EDH", "IM.H10S3..EDH",
    ]
    if channels != expected_channels:
        raise RuntimeError("CHANNEL_SET_DRIFT")

    all_dates = list(wc["pre_fault_dates"]) + list(wc["post_fault_dates"])
    first = parse_start(all_dates[0], wc["window_start_time_utc_each_day"])
    last = parse_start(all_dates[-1], wc["window_start_time_utc_each_day"])
    inv_start = (first - timedelta(days=2)).isoformat().replace("+00:00", "Z")
    inv_end = (last + timedelta(days=2)).isoformat().replace("+00:00", "Z")

    session = requests.Session()
    session.headers["User-Agent"] = "Janus-Echo-Cousteau/3.0 pre-post-2013 triplet gate"

    inventories: dict[str, Any] = {}
    inventory_meta: dict[str, Any] = {}
    inventory_errors: dict[str, str] = {}
    for cid in channels:
        try:
            inv, meta = transport_v2.fetch_inventory_bounded(session, cid, inv_start, inv_end)
            inventories[cid] = inv
            inventory_meta[cid] = meta
        except Exception as exc:
            inventory_errors[cid] = f"{type(exc).__name__}:{exc}"

    epochs = {
        "PRE_MARCH_2013": list(wc["pre_fault_dates"]),
        "POST_SEPTEMBER_2013": list(wc["post_fault_dates"]),
    }
    rows_by_epoch: dict[str, list[dict[str, Any]]] = {key: [] for key in epochs}
    all_windows: list[dict[str, Any]] = []
    data_errors: list[dict[str, Any]] = []

    for epoch, dates in epochs.items():
        for date_text in dates:
            start_dt = parse_start(date_text, wc["window_start_time_utc_each_day"])
            end_dt = start_dt + timedelta(seconds=duration)
            start = start_dt.isoformat().replace("+00:00", "Z")
            end = end_dt.isoformat().replace("+00:00", "Z")
            row: dict[str, Any] = {
                "epoch": epoch,
                "date": date_text,
                "start_utc": start,
                "end_utc": end,
                "stations": {},
                "complete_all_six": False,
            }
            traces: dict[str, Any] = {}
            for cid in channels:
                try:
                    if cid not in inventories:
                        raise RuntimeError("RESPONSE_UNAVAILABLE:" + inventory_errors.get(cid, "UNKNOWN"))
                    tr, waveform = transport_v2.fetch_trace_bounded(session, cid, start, end)
                    powers = integrated_powers(tr, inventories[cid], band, nperseg, noverlap)
                    traces[cid] = tr
                    row["stations"][cid] = {"data_status": "ANALYZED", "waveform": waveform, **powers}
                except Exception as exc:
                    row["stations"][cid] = {
                        "data_status": "BLOCKED",
                        "error_type": type(exc).__name__,
                        "error": str(exc),
                    }
                    data_errors.append({
                        "epoch": epoch,
                        "date": date_text,
                        "station": cid,
                        "error_type": type(exc).__name__,
                        "error": str(exc),
                    })

            if len(traces) == 6:
                try:
                    pair_specs = {
                        "S1_S2": ("IM.H10S1..EDH", "IM.H10S2..EDH"),
                        "S1_S3": ("IM.H10S1..EDH", "IM.H10S3..EDH"),
                        "N1_N2": ("IM.H10N1..EDH", "IM.H10N2..EDH"),
                        "N1_N3": ("IM.H10N1..EDH", "IM.H10N3..EDH"),
                    }
                    pairs: dict[str, Any] = {}
                    for name, (a, b) in pair_specs.items():
                        pairs[name] = pair_coupling(traces[a], traces[b], band, nperseg, noverlap)
                    row["pair_coupling"] = pairs
                    row["south_driver_coherence"] = float(np.mean([
                        pairs["S1_S2"]["median_magnitude_squared_coherence"],
                        pairs["S1_S3"]["median_magnitude_squared_coherence"],
                    ]))
                    row["north_driver_coherence"] = float(np.mean([
                        pairs["N1_N2"]["median_magnitude_squared_coherence"],
                        pairs["N1_N3"]["median_magnitude_squared_coherence"],
                    ]))
                    row["south_driver_abs_correlation"] = float(np.mean([
                        pairs["S1_S2"]["abs_zero_lag_pearson_correlation"],
                        pairs["S1_S3"]["abs_zero_lag_pearson_correlation"],
                    ]))
                    row["north_driver_abs_correlation"] = float(np.mean([
                        pairs["N1_N2"]["abs_zero_lag_pearson_correlation"],
                        pairs["N1_N3"]["abs_zero_lag_pearson_correlation"],
                    ]))
                    n1 = row["stations"]["IM.H10N1..EDH"]
                    for sid, short in [
                        ("IM.H10S1..EDH", "s1"),
                        ("IM.H10S2..EDH", "s2"),
                        ("IM.H10S3..EDH", "s3"),
                    ]:
                        ss = row["stations"][sid]
                        row[f"raw_{short}_minus_n1_db"] = ratio_db(
                            float(ss["raw_integrated_power"]), float(n1["raw_integrated_power"])
                        )
                        row[f"corrected_{short}_minus_n1_db"] = ratio_db(
                            float(ss["corrected_integrated_power_pa2"]),
                            float(n1["corrected_integrated_power_pa2"]),
                        )
                    row["complete_all_six"] = True
                    rows_by_epoch[epoch].append(row)
                except Exception as exc:
                    row["coupling_status"] = "BLOCKED"
                    row["coupling_error_type"] = type(exc).__name__
                    row["coupling_error"] = str(exc)
                    data_errors.append({
                        "epoch": epoch,
                        "date": date_text,
                        "station": "PAIR_COUPLING",
                        "error_type": type(exc).__name__,
                        "error": str(exc),
                    })
            all_windows.append(row)

    pre = rows_by_epoch["PRE_MARCH_2013"]
    post = rows_by_epoch["POST_SEPTEMBER_2013"]
    verdict, diagnostic = decide(pre, post, protocol)

    receipt = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-HA10-PRE-POST-2013-CROSSTALK-TRIPLET-GATE-RUN",
        "created_utc": utc_now(),
        "gate_id": protocol["gate_id"],
        "protocol_git_blob_sha1": EXPECTED_PROTOCOL_BLOB,
        "frozen_119hz_result_git_blob_sha1": EXPECTED_FROZEN_119_BLOB,
        "frozen_119hz_verdict": frozen["summary"]["verdict"],
        "spectral_helper_git_blob_sha1": EXPECTED_SPECTRAL_V1_BLOB,
        "transport_helper_git_blob_sha1": EXPECTED_TRANSPORT_V2_BLOB,
        "historical_fault_boundary": protocol["historical_anchor"],
        "analysis_order": protocol["analysis_order"],
        "frequency_band_hz": band,
        "window_contract": wc,
        "response_metadata": inventory_meta,
        "response_metadata_errors": inventory_errors,
        "epochs": {
            "PRE_MARCH_2013": {"complete_windows": len(pre)},
            "POST_SEPTEMBER_2013": {"complete_windows": len(post)},
        },
        "verdict": verdict,
        "diagnostic": diagnostic,
        "windows": all_windows,
        "data_errors": data_errors,
        "authority_delta_for_119hz": 0,
        "target_identity": "UNCONFIRMED",
        "target_evidence_delta": "NONE_FROM_PRE_POST_CROSSTALK_GATE",
        "source_writeback": False,
        "interpretation_ceiling": protocol["interpretation_firewall"],
        "janus_next_step_required": True,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True), encoding="utf-8")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    receipt = run(Path(args.output))
    print(json.dumps({
        "gate_id": receipt["gate_id"],
        "verdict": receipt["verdict"],
        "pre_complete": receipt["epochs"]["PRE_MARCH_2013"]["complete_windows"],
        "post_complete": receipt["epochs"]["POST_SEPTEMBER_2013"]["complete_windows"],
        "diagnostic": receipt["diagnostic"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
