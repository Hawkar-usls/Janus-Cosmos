#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import io
import itertools
import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import requests
from obspy import UTCDateTime, read, read_inventory

from workspace import cousteau_ha10_tphase_inband_positive_control as v1

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cousteau"
PROTOCOL = DATA / "JANUS-ECHO-COUSTEAU-HA10-TPHASE-INBAND-ARRIVAL-POSITIVE-CONTROL-PROTOCOL-2026-08-22-v1.0.json"
ADDENDUM = DATA / "JANUS-ECHO-COUSTEAU-HA10-TPHASE-INBAND-POSITIVE-CONTROL-INDEPENDENCE-ADDENDUM-2026-08-22-v1.0.json"
WINDOWS = DATA / "JANUS-ECHO-COUSTEAU-HA10-CONFIRMATORY-WINDOW-FREEZE-2026-08-22-v1.0.json"
FROZEN_119 = DATA / "JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0.json"
V1_HELPER = ROOT / "workspace" / "cousteau_ha10_tphase_inband_positive_control.py"
FDSN = "https://service.earthscope.org/fdsnws"

EXPECTED_PROTOCOL_BLOB = "d1800856076fc46cb04b5461bb1f2bf786af532e"
EXPECTED_ADDENDUM_BLOB = "f517c170948ed16e6097fc6fca38af7b7fd6edda"
EXPECTED_WINDOWS_BLOB = "f5fbc4b155d514094867b88bb33af71d3b458f76"
EXPECTED_FROZEN_119_BLOB = "eb8b48fb7f043160c057f9df6264a781412ed854"
EXPECTED_V1_HELPER_BLOB = "244169a1331a16529e2a963585a5659dae109c66"
EXPECTED_WINDOW_FREEZE_SHA256 = "d5f3c29c1dc4f7d7862724d1225688f9ee88460266d32fb0ec99fabe52cf2671"
EXPECTED_FROZEN_119_VERDICT = "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE"
HTTP_ATTEMPTS = 2
HTTP_TIMEOUT_S = 15
HTTP_BACKOFF_S = 1.0


def git_blob_sha1_file(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def bounded_get(
    session: requests.Session,
    url: str,
    params: dict[str, Any],
    *,
    attempts: int = HTTP_ATTEMPTS,
    timeout_s: int = HTTP_TIMEOUT_S,
) -> requests.Response:
    last: Exception | None = None
    for attempt in range(attempts):
        try:
            response = session.get(url, params=params, timeout=timeout_s)
            if response.status_code == 200 and response.content:
                return response
            last = RuntimeError(
                f"HTTP_{response.status_code}:{response.text[:300]}"
            )
        except Exception as exc:
            last = exc
        if attempt + 1 < attempts:
            time.sleep(HTTP_BACKOFF_S)
    raise last or RuntimeError("DOWNLOAD_FAILED")


def fetch_inventory_bounded(
    session: requests.Session, cid: str, start: str, end: str
):
    net, sta, loc, cha = v1.parse_channel_id(cid)
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
    response = bounded_get(session, FDSN + "/station/1/query", params)
    inventory = read_inventory(io.BytesIO(response.content))
    return inventory, {
        "url": response.url,
        "bytes": len(response.content),
        "sha256": v1.sha256_bytes(response.content),
        "content_type": response.headers.get("content-type"),
        "attempt_limit": HTTP_ATTEMPTS,
        "per_attempt_timeout_s": HTTP_TIMEOUT_S,
    }


def fetch_trace_bounded(
    session: requests.Session, cid: str, start: str, end: str
):
    net, sta, loc, cha = v1.parse_channel_id(cid)
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
    response = bounded_get(session, FDSN + "/dataselect/1/query", params)
    stream = read(io.BytesIO(response.content)).select(
        network=net, station=sta, channel=cha
    )
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
        "sha256": v1.sha256_bytes(response.content),
        "npts": int(trace.stats.npts),
        "sampling_rate_hz": sampling_rate,
        "starttime": str(trace.stats.starttime),
        "endtime": str(trace.stats.endtime),
        "attempt_limit": HTTP_ATTEMPTS,
        "per_attempt_timeout_s": HTTP_TIMEOUT_S,
    }


def _to_epoch(value: str) -> float:
    return float(UTCDateTime(value))


def intervals_overlap(a: dict[str, str], b: dict[str, str]) -> bool:
    return _to_epoch(a["start_utc"]) < _to_epoch(b["end_utc"]) and _to_epoch(
        b["start_utc"]
    ) < _to_epoch(a["end_utc"])


def events_conflict(a: dict[str, Any], b: dict[str, Any]) -> bool:
    station_ids = sorted(set(a["stations"]) & set(b["stations"]))
    if not station_ids:
        raise ValueError("NO_COMMON_STATIONS_FOR_INDEPENDENCE_CHECK")
    return any(
        intervals_overlap(
            a["stations"][sid]["event_window"],
            b["stations"][sid]["event_window"],
        )
        for sid in station_ids
    )


def maximum_nonoverlapping_passing_subset(
    events: list[dict[str, Any]],
) -> list[str]:
    passing = [
        event
        for event in events
        if event.get("positive_control_replicated_both_stations") is True
    ]
    passing.sort(key=lambda event: str(event["source_time_code"]))
    for size in range(len(passing), -1, -1):
        valid: list[list[str]] = []
        for subset in itertools.combinations(passing, size):
            if any(events_conflict(a, b) for a, b in itertools.combinations(subset, 2)):
                continue
            valid.append([str(event["source_time_code"]) for event in subset])
        if valid:
            valid.sort()
            return valid[0]
    return []


def decide_control(
    events: list[dict[str, Any]], minimum_events: int
) -> tuple[str, dict[str, Any]]:
    complete_analyzed = sum(
        1
        for event in events
        if len(event.get("stations", {})) == 2
        and all(
            row.get("data_status") == "ANALYZED"
            for row in event["stations"].values()
        )
    )
    raw_candidate_ids = sorted(
        str(event["source_time_code"])
        for event in events
        if event.get("positive_control_replicated_both_stations") is True
    )
    independent_ids = maximum_nonoverlapping_passing_subset(events)
    diagnostic = {
        "complete_events_analyzed_on_both_stations": complete_analyzed,
        "raw_candidate_events_passing_both_stations": len(raw_candidate_ids),
        "raw_candidate_source_time_codes": raw_candidate_ids,
        "independent_events_passing_both_stations": len(independent_ids),
        "independent_source_time_codes": independent_ids,
        "overlap_deduplication_applied": True,
    }
    if complete_analyzed < minimum_events:
        return "BLOCKED_POSITIVE_CONTROL_DATA_ACCESS_OR_RESPONSE", diagnostic
    if len(independent_ids) >= minimum_events:
        return "PASS_HA10_INBAND_TPHASE_PIPELINE_CONTROL", diagnostic
    return "FAIL_HA10_INBAND_TPHASE_PIPELINE_CONTROL", diagnostic


def verify_frozen_contracts() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    expected = {
        PROTOCOL: EXPECTED_PROTOCOL_BLOB,
        ADDENDUM: EXPECTED_ADDENDUM_BLOB,
        WINDOWS: EXPECTED_WINDOWS_BLOB,
        FROZEN_119: EXPECTED_FROZEN_119_BLOB,
        V1_HELPER: EXPECTED_V1_HELPER_BLOB,
    }
    for path, blob in expected.items():
        if git_blob_sha1_file(path) != blob:
            raise RuntimeError(f"FROZEN_BLOB_DRIFT:{path.name}")

    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    addendum = json.loads(ADDENDUM.read_text(encoding="utf-8"))
    windows = json.loads(WINDOWS.read_text(encoding="utf-8"))
    frozen_119 = json.loads(FROZEN_119.read_text(encoding="utf-8"))
    if protocol.get("status") != "PREREGISTERED_POST_NEGATIVE_PIPELINE_CONTROL_BEFORE_CONTROL_RUN":
        raise RuntimeError("POSITIVE_CONTROL_PROTOCOL_STATUS_MISMATCH")
    if addendum.get("status") != "POSTRUN_IMPLEMENTATION_CLARIFICATION__MONOTONIC_NON_PROMOTIONAL":
        raise RuntimeError("INDEPENDENCE_ADDENDUM_STATUS_MISMATCH")
    if windows.get("status") != "WINDOWS_FROZEN_READY_FOR_FFT":
        raise RuntimeError("WINDOW_FREEZE_NOT_READY")
    if windows.get("window_freeze_sha256") != EXPECTED_WINDOW_FREEZE_SHA256:
        raise RuntimeError("WINDOW_FREEZE_SEMANTIC_HASH_DRIFT")
    if frozen_119.get("summary", {}).get("verdict") != EXPECTED_FROZEN_119_VERDICT:
        raise RuntimeError("FROZEN_119_NEGATIVE_VERDICT_DRIFT")
    if protocol["epistemic_position"].get("authority_delta_for_119hz") != 0:
        raise RuntimeError("POSITIVE_CONTROL_AUTHORITY_DELTA_NOT_ZERO")
    return protocol, windows, frozen_119


def ensure_noncanonical_output(path: Path) -> None:
    resolved = path.resolve()
    try:
        resolved.relative_to(DATA.resolve())
    except ValueError:
        return
    raise RuntimeError("CANONICAL_DATA_WRITE_FORBIDDEN_BY_POSITIVE_CONTROL_V2")


def run(output: Path) -> dict[str, Any]:
    ensure_noncanonical_output(output)
    protocol, windows, frozen_119 = verify_frozen_contracts()
    estimator = protocol["spectral_estimator"]
    thresholds = protocol["admission_thresholds"]
    nperseg = int(estimator["welch_nperseg_samples"])
    noverlap = int(estimator["welch_noverlap_samples"])

    session = requests.Session()
    session.headers["User-Agent"] = (
        "Janus-Echo-Cousteau/2.0 bounded overlap-safe in-band positive control"
    )

    inventories: dict[str, Any] = {}
    response_metadata: dict[str, Any] = {}
    inventory_errors: dict[str, str] = {}
    for channel in protocol["channels"]:
        cid = channel["id"]
        try:
            inventory, metadata = fetch_inventory_bounded(
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
                event_trace, event_metadata = fetch_trace_bounded(
                    session,
                    cid,
                    station_window["event_window"]["start_utc"],
                    station_window["event_window"]["end_utc"],
                )
                f, raw, corrected = v1.psd_with_response(
                    event_trace,
                    inventories[cid],
                    nperseg=nperseg,
                    noverlap=noverlap,
                )
                event_spectrum = v1.summarize_control_spectrum(
                    f, raw, corrected, protocol
                )

                noise_spectra: list[dict[str, Any]] = []
                noise_metadata: list[dict[str, Any]] = []
                for noise_window in station_window["noise_windows"]:
                    noise_trace, metadata = fetch_trace_bounded(
                        session,
                        cid,
                        noise_window["start_utc"],
                        noise_window["end_utc"],
                    )
                    nf, nraw, ncorrected = v1.psd_with_response(
                        noise_trace,
                        inventories[cid],
                        nperseg=nperseg,
                        noverlap=noverlap,
                    )
                    if np.max(np.abs(nf - f)) > 1e-9:
                        raise RuntimeError("NOISE_FREQUENCY_GRID_MISMATCH")
                    noise_spectra.append(
                        v1.summarize_control_spectrum(
                            nf, nraw, ncorrected, protocol
                        )
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
                raw_snr = v1.db_ratio(
                    event_spectrum["raw_integrated_power"], raw_noise
                )
                corrected_snr = v1.db_ratio(
                    event_spectrum["corrected_integrated_power_pa2"],
                    corrected_noise,
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
                        "raw_event_vs_noise_snr_db": v1.db_ratio(
                            event_sub["raw_integrated_power"], raw_sub_noise
                        ),
                        "corrected_event_vs_noise_snr_db": v1.db_ratio(
                            event_sub["corrected_integrated_power_pa2"],
                            corrected_sub_noise,
                        ),
                    }

                per_station = thresholds["per_station_event"]
                checks = {
                    "raw_broadband_snr_ge_3db": raw_snr
                    >= float(
                        per_station[
                            "minimum_raw_broadband_event_vs_noise_snr_db"
                        ]
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
        thresholds["cross_station"][
            "minimum_independent_events_passing_both_stations"
        ]
    )
    verdict, independence = decide_control(results, minimum_events)

    receipt = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-HA10-TPHASE-INBAND-ARRIVAL-POSITIVE-CONTROL-REPLAY-V2",
        "created_utc": utc_now(),
        "gate_id": protocol["gate_id"],
        "implementation_semantics": "V2_BOUNDED_DOWNLOADS_AND_NONOVERLAPPING_INDEPENDENT_EVENT_COUNT",
        "protocol_git_blob_sha1": EXPECTED_PROTOCOL_BLOB,
        "independence_addendum_git_blob_sha1": EXPECTED_ADDENDUM_BLOB,
        "window_freeze_git_blob_sha1": EXPECTED_WINDOWS_BLOB,
        "frozen_119hz_result_git_blob_sha1": EXPECTED_FROZEN_119_BLOB,
        "spectral_v1_helper_git_blob_sha1": EXPECTED_V1_HELPER_BLOB,
        "frozen_119hz_verdict": frozen_119["summary"]["verdict"],
        "frozen_119hz_negative_result_immutable": True,
        "authority_delta_for_119hz": 0,
        "target_identity": "UNCONFIRMED",
        "target_evidence_delta": "NONE_FROM_POSITIVE_CONTROL_REPLAY",
        "network_budget": {
            "attempts_per_request": HTTP_ATTEMPTS,
            "timeout_per_attempt_s": HTTP_TIMEOUT_S,
            "backoff_s": HTTP_BACKOFF_S,
        },
        "response_metadata": response_metadata,
        "inventory_errors": inventory_errors,
        "events": results,
        "independence": independence,
        "summary": {
            "selected_complete_events": len(selected_events),
            "complete_events_analyzed_on_both_stations": independence[
                "complete_events_analyzed_on_both_stations"
            ],
            "raw_candidate_events_passing_both_stations": independence[
                "raw_candidate_events_passing_both_stations"
            ],
            "independent_events_passing_both_stations": independence[
                "independent_events_passing_both_stations"
            ],
            "minimum_events_required": minimum_events,
            "blocked_station_event_pairs": len(data_errors),
            "verdict": verdict,
        },
        "data_errors": data_errors,
        "claim_ceiling": protocol["claim_ceiling"],
        "hard_rules": list(protocol["hard_rules"])
        + [
            "OVERLAPPING_RECEPTION_WINDOWS_DO_NOT_COUNT_AS_INDEPENDENT",
            "BOUNDED_NETWORK_RETRIES_MUST_ALLOW_BLOCKED_RECEIPT",
        ],
        "status": "CONTROL_REPLAY_COMPLETE",
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
