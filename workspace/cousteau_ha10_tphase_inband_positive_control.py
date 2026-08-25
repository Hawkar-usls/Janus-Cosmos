#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import requests
from obspy import UTCDateTime, read, read_inventory
from scipy.signal import welch

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cousteau"
PROTOCOL = DATA / "JANUS-ECHO-COUSTEAU-HA10-TPHASE-INBAND-ARRIVAL-POSITIVE-CONTROL-PROTOCOL-2026-08-22-v1.0.json"
WINDOWS = DATA / "JANUS-ECHO-COUSTEAU-HA10-CONFIRMATORY-WINDOW-FREEZE-2026-08-22-v1.0.json"
FROZEN_119 = DATA / "JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0.json"
EXPECTED_PROTOCOL_BLOB = "d1800856076fc46cb04b5461bb1f2bf786af532e"
EXPECTED_WINDOWS_BLOB = "f5fbc4b155d514094867b88bb33af71d3b458f76"
EXPECTED_FROZEN_119_BLOB = "eb8b48fb7f043160c057f9df6264a781412ed854"
EXPECTED_WINDOW_FREEZE_SHA256 = "d5f3c29c1dc4f7d7862724d1225688f9ee88460266d32fb0ec99fabe52cf2671"
EXPECTED_FROZEN_119_VERDICT = "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE"
FDSN = "https://service.earthscope.org/fdsnws"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git_blob_sha1_bytes(data: bytes) -> str:
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def git_blob_sha1_file(path: Path) -> str:
    return git_blob_sha1_bytes(path.read_bytes())


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def db_ratio(a: float, b: float) -> float:
    if not (a > 0 and b > 0 and math.isfinite(a) and math.isfinite(b)):
        return float("nan")
    return 10.0 * math.log10(a / b)


def parse_channel_id(cid: str) -> tuple[str, str, str, str]:
    parts = cid.split(".")
    if len(parts) != 4:
        raise ValueError(f"INVALID_CHANNEL_ID:{cid}")
    return parts[0], parts[1], parts[2], parts[3]


def get_bytes(
    session: requests.Session,
    url: str,
    params: dict[str, Any],
    *,
    tries: int = 4,
    timeout: int = 120,
) -> requests.Response:
    last: Exception | None = None
    for attempt in range(tries):
        try:
            response = session.get(url, params=params, timeout=timeout)
            if response.status_code == 200 and response.content:
                return response
            last = RuntimeError(
                f"HTTP_{response.status_code}:{response.text[:300]}"
            )
        except Exception as exc:  # network boundary: preserve failure as data-access state
            last = exc
        if attempt + 1 < tries:
            time.sleep(1.5 * (attempt + 1))
    raise last or RuntimeError("DOWNLOAD_FAILED")


def fetch_inventory(
    session: requests.Session, cid: str, start: str, end: str
):
    net, sta, loc, cha = parse_channel_id(cid)
    params = {
        "net": net,
        "sta": sta,
        "loc": "--" if loc == "" else loc,
        "cha": cha,
        "starttime": start,
        "endtime": end,
        "level": "response",
        "format": "xml",
        "nodata": 404,
    }
    response = get_bytes(session, FDSN + "/station/1/query", params, timeout=90)
    inventory = read_inventory(io.BytesIO(response.content))
    return inventory, {
        "url": response.url,
        "bytes": len(response.content),
        "sha256": sha256_bytes(response.content),
        "content_type": response.headers.get("content-type"),
    }


def fetch_trace(
    session: requests.Session, cid: str, start: str, end: str
):
    net, sta, loc, cha = parse_channel_id(cid)
    params = {
        "net": net,
        "sta": sta,
        "loc": "--" if loc == "" else loc,
        "cha": cha,
        "starttime": start,
        "endtime": end,
        "nodata": 404,
        "format": "miniseed",
    }
    response = get_bytes(session, FDSN + "/dataselect/1/query", params, timeout=180)
    stream = read(io.BytesIO(response.content))
    stream = stream.select(network=net, station=sta, channel=cha)
    if len(stream) == 0:
        raise RuntimeError(f"NO_TRACE:{cid}")
    stream.merge(method=0, fill_value=None)
    if len(stream) != 1:
        raise RuntimeError(f"MULTIPLE_UNMERGED_TRACES:{cid}:{len(stream)}")
    trace = stream[0]
    if np.ma.isMaskedArray(trace.data) and np.any(np.ma.getmaskarray(trace.data)):
        raise RuntimeError(f"GAPS_PRESENT_INTERPOLATION_FORBIDDEN:{cid}")
    trace = trace.copy()
    trace.detrend("demean")
    trace.detrend("linear")
    sampling_rate = float(trace.stats.sampling_rate)
    expected_samples = (UTCDateTime(end) - UTCDateTime(start)) * sampling_rate
    if len(trace.data) < expected_samples * 0.98:
        raise RuntimeError(
            f"INCOMPLETE_WINDOW:{cid}:{len(trace.data)}<{0.98 * expected_samples}"
        )
    return trace, {
        "url": response.url,
        "bytes": len(response.content),
        "sha256": sha256_bytes(response.content),
        "npts": int(trace.stats.npts),
        "sampling_rate_hz": sampling_rate,
        "starttime": str(trace.stats.starttime),
        "endtime": str(trace.stats.endtime),
    }


def psd_with_response(trace, inventory, *, nperseg: int, noverlap: int):
    values = np.asarray(trace.data, dtype=np.float64)
    sampling_rate = float(trace.stats.sampling_rate)
    if len(values) < nperseg:
        raise RuntimeError(f"NPTS_BELOW_NPERSEG:{len(values)}<{nperseg}")
    frequencies, raw_psd = welch(
        values,
        fs=sampling_rate,
        window="hann",
        nperseg=nperseg,
        noverlap=noverlap,
        detrend=False,
        scaling="density",
        return_onesided=True,
    )
    response = inventory.get_response(trace.id, trace.stats.starttime)
    transfer, response_frequencies = response.get_evalresp_response(
        t_samp=1.0 / sampling_rate,
        nfft=nperseg,
        output="DEF",
    )
    if len(response_frequencies) != len(frequencies):
        raise RuntimeError("RESPONSE_FREQUENCY_GRID_LENGTH_MISMATCH")
    if np.max(np.abs(response_frequencies - frequencies)) > 1e-9:
        raise RuntimeError("RESPONSE_FREQUENCY_GRID_MISMATCH")
    amplitude = np.abs(transfer)
    corrected_psd = np.full_like(raw_psd, np.nan, dtype=float)
    valid = np.isfinite(amplitude) & (amplitude > 0) & np.isfinite(raw_psd)
    corrected_psd[valid] = raw_psd[valid] / (amplitude[valid] ** 2)
    return frequencies, raw_psd, corrected_psd


def integrated_band_power(
    frequencies: np.ndarray, psd: np.ndarray, band: list[float]
) -> tuple[float, int]:
    lo, hi = map(float, band)
    mask = (
        (frequencies >= lo)
        & (frequencies <= hi)
        & np.isfinite(psd)
    )
    count = int(np.sum(mask))
    if count < 2:
        raise RuntimeError(f"INSUFFICIENT_FINITE_BINS:{band}:{count}")
    x = frequencies[mask]
    y = psd[mask]
    power = float(np.trapz(y, x))
    if not (power > 0 and math.isfinite(power)):
        raise RuntimeError(f"INVALID_INTEGRATED_POWER:{band}:{power}")
    return power, count


def summarize_control_spectrum(
    frequencies: np.ndarray,
    raw_psd: np.ndarray,
    corrected_psd: np.ndarray,
    protocol: dict[str, Any],
) -> dict[str, Any]:
    contract = protocol["control_frequency_contract"]
    broadband = contract["broadband_hz"]
    raw_power, raw_bins = integrated_band_power(frequencies, raw_psd, broadband)
    corrected_power, corrected_bins = integrated_band_power(
        frequencies, corrected_psd, broadband
    )
    subbands: dict[str, Any] = {}
    for name, band in contract["diagnostic_subbands_hz"].items():
        raw_value, n_raw = integrated_band_power(frequencies, raw_psd, band)
        corrected_value, n_corrected = integrated_band_power(
            frequencies, corrected_psd, band
        )
        subbands[name] = {
            "band_hz": band,
            "raw_integrated_power": raw_value,
            "corrected_integrated_power_pa2": corrected_value,
            "raw_bin_count": n_raw,
            "corrected_bin_count": n_corrected,
        }
    return {
        "broadband_hz": broadband,
        "raw_integrated_power": raw_power,
        "corrected_integrated_power_pa2": corrected_power,
        "raw_bin_count": raw_bins,
        "corrected_bin_count": corrected_bins,
        "diagnostic_subbands": subbands,
        "peak_search_performed": False,
    }


def decide_control(
    events: list[dict[str, Any]], minimum_events: int
) -> tuple[str, int, int]:
    complete_analyzed = 0
    replicated_pass = 0
    for event in events:
        station_rows = list(event.get("stations", {}).values())
        if len(station_rows) == 2 and all(
            row.get("data_status") == "ANALYZED" for row in station_rows
        ):
            complete_analyzed += 1
        if event.get("positive_control_replicated_both_stations") is True:
            replicated_pass += 1
    if complete_analyzed < minimum_events:
        return "BLOCKED_POSITIVE_CONTROL_DATA_ACCESS_OR_RESPONSE", complete_analyzed, replicated_pass
    if replicated_pass >= minimum_events:
        return "PASS_HA10_INBAND_TPHASE_PIPELINE_CONTROL", complete_analyzed, replicated_pass
    return "FAIL_HA10_INBAND_TPHASE_PIPELINE_CONTROL", complete_analyzed, replicated_pass


def verify_frozen_contracts() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    if git_blob_sha1_file(PROTOCOL) != EXPECTED_PROTOCOL_BLOB:
        raise RuntimeError("POSITIVE_CONTROL_PROTOCOL_BLOB_DRIFT")
    if git_blob_sha1_file(WINDOWS) != EXPECTED_WINDOWS_BLOB:
        raise RuntimeError("WINDOW_FREEZE_BLOB_DRIFT")
    if git_blob_sha1_file(FROZEN_119) != EXPECTED_FROZEN_119_BLOB:
        raise RuntimeError("FROZEN_119_RESULT_BLOB_DRIFT")

    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    windows = json.loads(WINDOWS.read_text(encoding="utf-8"))
    frozen_119 = json.loads(FROZEN_119.read_text(encoding="utf-8"))

    if protocol.get("status") != "PREREGISTERED_POST_NEGATIVE_PIPELINE_CONTROL_BEFORE_CONTROL_RUN":
        raise RuntimeError("POSITIVE_CONTROL_PROTOCOL_STATUS_MISMATCH")
    if windows.get("status") != "WINDOWS_FROZEN_READY_FOR_FFT":
        raise RuntimeError("WINDOW_FREEZE_NOT_READY")
    if windows.get("window_freeze_sha256") != EXPECTED_WINDOW_FREEZE_SHA256:
        raise RuntimeError("WINDOW_FREEZE_SEMANTIC_HASH_DRIFT")
    if frozen_119.get("summary", {}).get("verdict") != EXPECTED_FROZEN_119_VERDICT:
        raise RuntimeError("FROZEN_119_NEGATIVE_VERDICT_DRIFT")
    if protocol["epistemic_position"].get("authority_delta_for_119hz") != 0:
        raise RuntimeError("POSITIVE_CONTROL_AUTHORITY_DELTA_NOT_ZERO")
    if protocol["control_frequency_contract"].get("119hz_or_117_121hz_bins_may_be_used") is not False:
        raise RuntimeError("POSITIVE_CONTROL_TARGET_BAND_EXCLUSION_NOT_FROZEN")
    return protocol, windows, frozen_119


def ensure_noncanonical_output(path: Path) -> None:
    resolved = path.resolve()
    data_root = DATA.resolve()
    try:
        resolved.relative_to(data_root)
    except ValueError:
        return
    raise RuntimeError("CANONICAL_DATA_WRITE_FORBIDDEN_BY_POSITIVE_CONTROL_RUNNER")


def run(output: Path) -> dict[str, Any]:
    ensure_noncanonical_output(output)
    protocol, windows, frozen_119 = verify_frozen_contracts()
    estimator = protocol["spectral_estimator"]
    thresholds = protocol["admission_thresholds"]
    nperseg = int(estimator["welch_nperseg_samples"])
    noverlap = int(estimator["welch_noverlap_samples"])

    session = requests.Session()
    session.headers["User-Agent"] = "Janus-Echo-Cousteau/1.0 in-band T-phase positive control"

    inventories: dict[str, Any] = {}
    response_metadata: dict[str, Any] = {}
    inventory_errors: dict[str, str] = {}
    for channel in protocol["channels"]:
        cid = channel["id"]
        try:
            inventory, metadata = fetch_inventory(
                session,
                cid,
                "2014-12-11T00:00:00Z",
                "2015-01-13T00:00:00Z",
            )
            inventories[cid] = inventory
            response_metadata[cid] = metadata
        except Exception as exc:
            inventory_errors[cid] = f"{type(exc).__name__}:{exc}"

    results: list[dict[str, Any]] = []
    data_errors: list[dict[str, Any]] = []
    selected_events = [
        event
        for event in windows["selected_events"]
        if event.get("complete_on_both_stations")
    ]

    for event in selected_events:
        event_result: dict[str, Any] = {
            "selection_rank": event["selection_rank"],
            "source_time_code": event["source_time_code"],
            "origin_utc": event["origin_utc"],
            "stations": {},
        }
        for cid, station_window in event["stations"].items():
            station_result: dict[str, Any] = {
                "predicted_arrival_utc": station_window["predicted_arrival_utc"],
                "event_window": station_window["event_window"],
                "noise_windows": station_window["noise_windows"],
            }
            try:
                if cid not in inventories:
                    raise RuntimeError(
                        "RESPONSE_UNAVAILABLE:" + inventory_errors.get(cid, "UNKNOWN")
                    )
                event_trace, event_metadata = fetch_trace(
                    session,
                    cid,
                    station_window["event_window"]["start_utc"],
                    station_window["event_window"]["end_utc"],
                )
                f, raw, corrected = psd_with_response(
                    event_trace,
                    inventories[cid],
                    nperseg=nperseg,
                    noverlap=noverlap,
                )
                event_spectrum = summarize_control_spectrum(
                    f, raw, corrected, protocol
                )

                noise_spectra: list[dict[str, Any]] = []
                noise_metadata: list[dict[str, Any]] = []
                for noise_window in station_window["noise_windows"]:
                    noise_trace, metadata = fetch_trace(
                        session,
                        cid,
                        noise_window["start_utc"],
                        noise_window["end_utc"],
                    )
                    nf, nraw, ncorrected = psd_with_response(
                        noise_trace,
                        inventories[cid],
                        nperseg=nperseg,
                        noverlap=noverlap,
                    )
                    if np.max(np.abs(nf - f)) > 1e-9:
                        raise RuntimeError("NOISE_FREQUENCY_GRID_MISMATCH")
                    noise_spectra.append(
                        summarize_control_spectrum(nf, nraw, ncorrected, protocol)
                    )
                    noise_metadata.append(metadata)

                raw_noise = float(
                    np.median(
                        [row["raw_integrated_power"] for row in noise_spectra]
                    )
                )
                corrected_noise = float(
                    np.median(
                        [
                            row["corrected_integrated_power_pa2"]
                            for row in noise_spectra
                        ]
                    )
                )
                raw_snr = db_ratio(event_spectrum["raw_integrated_power"], raw_noise)
                corrected_snr = db_ratio(
                    event_spectrum["corrected_integrated_power_pa2"], corrected_noise
                )

                subband_scores: dict[str, Any] = {}
                for name in protocol["control_frequency_contract"][
                    "diagnostic_subbands_hz"
                ]:
                    event_sub = event_spectrum["diagnostic_subbands"][name]
                    raw_sub_noise = float(
                        np.median(
                            [
                                row["diagnostic_subbands"][name]["raw_integrated_power"]
                                for row in noise_spectra
                            ]
                        )
                    )
                    corrected_sub_noise = float(
                        np.median(
                            [
                                row["diagnostic_subbands"][name][
                                    "corrected_integrated_power_pa2"
                                ]
                                for row in noise_spectra
                            ]
                        )
                    )
                    subband_scores[name] = {
                        "raw_event_vs_noise_snr_db": db_ratio(
                            event_sub["raw_integrated_power"], raw_sub_noise
                        ),
                        "corrected_event_vs_noise_snr_db": db_ratio(
                            event_sub["corrected_integrated_power_pa2"],
                            corrected_sub_noise,
                        ),
                    }

                per_station = thresholds["per_station_event"]
                checks = {
                    "raw_broadband_snr_ge_3db": raw_snr
                    >= float(
                        per_station["minimum_raw_broadband_event_vs_noise_snr_db"]
                    ),
                    "corrected_broadband_snr_ge_3db": corrected_snr
                    >= float(
                        per_station[
                            "minimum_corrected_broadband_event_vs_noise_snr_db"
                        ]
                    ),
                }
                station_result.update(
                    {
                        "data_status": "ANALYZED",
                        "event_waveform": event_metadata,
                        "noise_waveforms": noise_metadata,
                        "event_spectrum": event_spectrum,
                        "noise_spectra": noise_spectra,
                        "raw_broadband_event_vs_noise_snr_db": raw_snr,
                        "corrected_broadband_event_vs_noise_snr_db": corrected_snr,
                        "diagnostic_subband_scores": subband_scores,
                        "per_station_event_checks": checks,
                        "per_station_event_pass": all(checks.values()),
                    }
                )
            except Exception as exc:
                station_result.update(
                    {
                        "data_status": "BLOCKED_PAIR",
                        "error_type": type(exc).__name__,
                        "error": str(exc),
                        "per_station_event_pass": False,
                    }
                )
                data_errors.append(
                    {
                        "source_time_code": event["source_time_code"],
                        "station": cid,
                        "error_type": type(exc).__name__,
                        "error": str(exc),
                    }
                )
            event_result["stations"][cid] = station_result

        station_ids = [row["id"] for row in protocol["channels"]]
        event_result["positive_control_replicated_both_stations"] = bool(
            event_result["stations"].get(station_ids[0], {}).get(
                "per_station_event_pass"
            )
            and event_result["stations"].get(station_ids[1], {}).get(
                "per_station_event_pass"
            )
        )
        results.append(event_result)

    minimum_events = int(
        thresholds["cross_station"]["minimum_independent_events_passing_both_stations"]
    )
    verdict, complete_analyzed, replicated_pass = decide_control(
        results, minimum_events
    )

    receipt = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-HA10-TPHASE-INBAND-ARRIVAL-POSITIVE-CONTROL-RUN",
        "created_utc": utc_now(),
        "gate_id": protocol["gate_id"],
        "protocol_path": str(PROTOCOL.relative_to(ROOT)),
        "protocol_git_blob_sha1": EXPECTED_PROTOCOL_BLOB,
        "window_freeze_path": str(WINDOWS.relative_to(ROOT)),
        "window_freeze_git_blob_sha1": EXPECTED_WINDOWS_BLOB,
        "window_freeze_semantic_sha256": EXPECTED_WINDOW_FREEZE_SHA256,
        "frozen_119hz_result_path": str(FROZEN_119.relative_to(ROOT)),
        "frozen_119hz_result_git_blob_sha1": EXPECTED_FROZEN_119_BLOB,
        "frozen_119hz_verdict": frozen_119["summary"]["verdict"],
        "frozen_119hz_negative_result_immutable": True,
        "authority_delta_for_119hz": 0,
        "target_identity": "UNCONFIRMED",
        "target_evidence_delta": "NONE_FROM_POSITIVE_CONTROL",
        "response_metadata": response_metadata,
        "inventory_errors": inventory_errors,
        "processing_contract": {
            "broadband_hz": protocol["control_frequency_contract"]["broadband_hz"],
            "diagnostic_subbands_hz": protocol["control_frequency_contract"][
                "diagnostic_subbands_hz"
            ],
            "whole_spectrum_peak_search_performed": False,
            "best_subband_selection_performed": False,
            "119hz_or_117_121hz_bins_used": False,
            "welch_nperseg": nperseg,
            "welch_noverlap": noverlap,
            "response_correction": estimator["response_correction"],
        },
        "events": results,
        "summary": {
            "selected_complete_events": len(selected_events),
            "complete_events_analyzed_on_both_stations": complete_analyzed,
            "events_passing_both_stations": replicated_pass,
            "minimum_events_required": minimum_events,
            "blocked_station_event_pairs": len(data_errors),
            "verdict": verdict,
        },
        "data_errors": data_errors,
        "claim_ceiling": protocol["claim_ceiling"],
        "hard_rules": protocol["hard_rules"],
        "status": "CONTROL_RUN_COMPLETE",
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
    print(json.dumps(receipt["summary"], indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
