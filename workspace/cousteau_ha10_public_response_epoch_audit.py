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
from obspy import UTCDateTime, read_inventory

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cousteau"
PROTOCOL = DATA / "JANUS-ECHO-COUSTEAU-HA10-PUBLIC-RESPONSE-EPOCH-AUDIT-PROTOCOL-2026-08-22-v1.1.json"
REVERSE_AUDIT = DATA / "JANUS-ECHO-COUSTEAU-HA10-REVERSE-SPIRAL-CALIBRATION-CROSSTALK-AUDIT-2026-08-22-v1.0.json"
FROZEN_119 = DATA / "JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0.json"
EXPECTED_PROTOCOL_BLOB = "182f057811e75f9ee2024316c103fec258945641"
EXPECTED_REVERSE_AUDIT_BLOB = "775b53d17e4fd2c185ab1b47cf6e5bfe8d5d3ddf"
EXPECTED_FROZEN_119_BLOB = "eb8b48fb7f043160c057f9df6264a781412ed854"
EXPECTED_FROZEN_119_VERDICT = "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git_blob_sha1_bytes(data: bytes) -> str:
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def git_blob_sha1_file(path: Path) -> str:
    return git_blob_sha1_bytes(path.read_bytes())


def canonical_sha256(value: Any) -> str:
    raw = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def jsonable(value: Any) -> Any:
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            return str(value)
        return value
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, np.generic):
        return jsonable(value.item())
    if isinstance(value, (list, tuple)):
        return [jsonable(x) for x in value]
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in sorted(value.items(), key=lambda kv: str(kv[0]))}
    if hasattr(value, "real") and hasattr(value, "imag"):
        try:
            return {"real": float(value.real), "imag": float(value.imag)}
        except Exception:
            pass
    return str(value)


def response_fingerprint_payload(response: Any) -> dict[str, Any]:
    sensitivity = response.instrument_sensitivity
    sens = None
    if sensitivity is not None:
        sens = {
            "value": jsonable(getattr(sensitivity, "value", None)),
            "frequency": jsonable(getattr(sensitivity, "frequency", None)),
            "input_units": jsonable(getattr(sensitivity, "input_units", None)),
            "output_units": jsonable(getattr(sensitivity, "output_units", None)),
        }
    stages: list[dict[str, Any]] = []
    for stage in response.response_stages:
        public: dict[str, Any] = {"class": stage.__class__.__name__}
        for key, value in sorted(getattr(stage, "__dict__", {}).items()):
            if key.startswith("_"):
                continue
            public[key] = jsonable(value)
        stages.append(public)
    return {
        "instrument_sensitivity": sens,
        "stage_count": len(stages),
        "stages": stages,
    }


def sensitivity_summary(response: Any) -> dict[str, Any] | None:
    sensitivity = response.instrument_sensitivity
    if sensitivity is None:
        return None
    return {
        "value": float(sensitivity.value),
        "frequency_hz": float(sensitivity.frequency),
        "input_units": str(sensitivity.input_units),
        "output_units": str(sensitivity.output_units),
    }


def response_magnitude_db(response: Any, frequencies_hz: list[float], sample_rate_hz: float) -> dict[str, float]:
    freqs = np.asarray(frequencies_hz, dtype=float)
    if hasattr(response, "get_evalresp_response_for_frequencies"):
        values = response.get_evalresp_response_for_frequencies(freqs, output="DEF")
        values = np.asarray(values)
    else:
        nfft = 262144
        values_grid, freq_grid = response.get_evalresp_response(
            t_samp=1.0 / sample_rate_hz,
            nfft=nfft,
            output="DEF",
        )
        amp_grid = np.abs(np.asarray(values_grid, dtype=complex))
        db_grid = np.full_like(amp_grid, np.nan, dtype=float)
        valid = np.isfinite(amp_grid) & (amp_grid > 0)
        db_grid[valid] = 20.0 * np.log10(amp_grid[valid])
        out: dict[str, float] = {}
        for f in freqs:
            out[f"{f:g}"] = float(np.interp(f, freq_grid[valid], db_grid[valid]))
        return out
    amps = np.abs(values.astype(complex))
    out = {}
    for f, amp in zip(freqs, amps):
        if not (math.isfinite(float(amp)) and float(amp) > 0):
            raise RuntimeError(f"INVALID_RESPONSE_MAGNITUDE:{f}:{amp}")
        out[f"{float(f):g}"] = float(20.0 * math.log10(float(amp)))
    return out


def get_bytes(session: requests.Session, url: str, params: dict[str, Any], *, tries: int = 3, timeout: int = 60) -> requests.Response:
    last: Exception | None = None
    for attempt in range(tries):
        try:
            response = session.get(url, params=params, timeout=timeout)
            if response.status_code == 200 and response.content:
                return response
            last = RuntimeError(f"HTTP_{response.status_code}:{response.text[:300]}")
        except Exception as exc:
            last = exc
        if attempt + 1 < tries:
            time.sleep(1.0 * (attempt + 1))
    raise last or RuntimeError("DOWNLOAD_FAILED")


def fetch_station_inventory(session: requests.Session, protocol: dict[str, Any], station: str) -> tuple[Any, dict[str, Any]]:
    contract = protocol["public_metadata_contract"]
    params = {
        "net": contract["network"],
        "sta": station,
        "loc": contract["location"],
        "cha": contract["channel"],
        "starttime": contract["query_start_utc"],
        "endtime": contract["query_end_utc"],
        "level": contract["level"],
        "format": contract["format"],
        "nodata": contract["nodata"],
    }
    response = get_bytes(session, contract["service"], params)
    inv = read_inventory(io.BytesIO(response.content))
    return inv, {
        "url": response.url,
        "bytes": len(response.content),
        "sha256": sha256_bytes(response.content),
        "content_type": response.headers.get("content-type"),
    }


def channel_epochs(inv: Any, station: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for network in inv:
        if network.code != "IM":
            continue
        for sta in network:
            if sta.code != station:
                continue
            for channel in sta:
                if channel.code != "EDH" or (channel.location_code or "") != "":
                    continue
                rows.append({
                    "start_utc": str(channel.start_date) if channel.start_date else None,
                    "end_utc": str(channel.end_date) if channel.end_date else None,
                    "sample_rate_hz": float(channel.sample_rate) if channel.sample_rate is not None else None,
                    "latitude": float(channel.latitude) if channel.latitude is not None else None,
                    "longitude": float(channel.longitude) if channel.longitude is not None else None,
                    "depth_m": float(channel.depth) if channel.depth is not None else None,
                    "azimuth_deg": float(channel.azimuth) if channel.azimuth is not None else None,
                    "dip_deg": float(channel.dip) if channel.dip is not None else None,
                })
    rows.sort(key=lambda r: r["start_utc"] or "")
    return rows


def select_channel(inv: Any, station: str, when: str) -> Any:
    t = UTCDateTime(when)
    candidates: list[Any] = []
    for network in inv:
        if network.code != "IM":
            continue
        for sta in network:
            if sta.code != station:
                continue
            for channel in sta:
                if channel.code != "EDH" or (channel.location_code or "") != "":
                    continue
                if channel.start_date is not None and t < channel.start_date:
                    continue
                if channel.end_date is not None and t > channel.end_date:
                    continue
                candidates.append(channel)
    if len(candidates) != 1:
        raise RuntimeError(f"CHANNEL_EPOCH_SELECTION_COUNT:{station}:{when}:{len(candidates)}")
    return candidates[0]


def anchor_row(inv: Any, station: str, when: str, frequencies_hz: list[float]) -> dict[str, Any]:
    channel = select_channel(inv, station, when)
    response = channel.response
    if response is None:
        raise RuntimeError(f"NO_RESPONSE:{station}:{when}")
    payload = response_fingerprint_payload(response)
    return {
        "status": "ANALYZED",
        "time_utc": when,
        "channel_epoch_start_utc": str(channel.start_date) if channel.start_date else None,
        "channel_epoch_end_utc": str(channel.end_date) if channel.end_date else None,
        "sample_rate_hz": float(channel.sample_rate),
        "instrument_sensitivity": sensitivity_summary(response),
        "response_stage_count": len(response.response_stages),
        "response_fingerprint_sha256": canonical_sha256(payload),
        "fixed_frequency_response_magnitude_db": response_magnitude_db(response, frequencies_hz, float(channel.sample_rate)),
    }


def curve_delta_db(later: dict[str, Any], earlier: dict[str, Any]) -> dict[str, Any]:
    a = earlier["fixed_frequency_response_magnitude_db"]
    b = later["fixed_frequency_response_magnitude_db"]
    deltas = {key: float(b[key] - a[key]) for key in a if key in b}
    finite = [abs(v) for v in deltas.values() if math.isfinite(v)]
    return {
        "per_frequency_delta_db": deltas,
        "max_abs_delta_db": max(finite) if finite else None,
        "fingerprint_changed": later["response_fingerprint_sha256"] != earlier["response_fingerprint_sha256"],
        "same_channel_epoch": (
            later.get("channel_epoch_start_utc") == earlier.get("channel_epoch_start_utc")
            and later.get("channel_epoch_end_utc") == earlier.get("channel_epoch_end_utc")
        ),
    }


def near_fault_boundaries(epochs: list[dict[str, Any]], fault_date: str, window_days: int) -> list[dict[str, Any]]:
    fault = UTCDateTime(fault_date + "T00:00:00Z")
    limit = float(window_days) * 86400.0
    out: list[dict[str, Any]] = []
    for epoch in epochs:
        for field in ("start_utc", "end_utc"):
            value = epoch.get(field)
            if not value:
                continue
            delta_s = float(UTCDateTime(value) - fault)
            if abs(delta_s) <= limit:
                out.append({"field": field, "time_utc": value, "delta_days": delta_s / 86400.0})
    return out


def decide_verdict(stations: dict[str, Any], protocol: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    contract = protocol["public_metadata_contract"]
    threshold = float(protocol["classification_contract"]["numerical_curve_change_floor_db"])
    required_n = int(protocol["completeness_gate"]["required_north_channels_with_pre_and_post_response"])
    required_s = int(protocol["completeness_gate"]["required_south_channels_with_pre_and_post_response"])

    def complete(station: str) -> bool:
        row = stations.get(station, {})
        anchors = row.get("anchors", {})
        return (
            anchors.get("PRE_FAULT_DAY", {}).get("status") == "ANALYZED"
            and anchors.get("POST_SEP_2013", {}).get("status") == "ANALYZED"
            and row.get("pre_to_post_sep_2013", {}).get("max_abs_delta_db") is not None
        )

    complete_n = [s for s in contract["north_stations"] if complete(s)]
    complete_s = [s for s in contract["south_stations"] if complete(s)]
    diagnostics = {"complete_north": complete_n, "complete_south": complete_s}
    if len(complete_n) < required_n or len(complete_s) < required_s:
        return "BLOCKED_PUBLIC_RESPONSE_EPOCH_AUDIT", diagnostics

    south_same = True
    south_encodes = False
    for station in contract["south_stations"]:
        comparison = stations[station]["pre_to_post_sep_2013"]
        max_delta = float(comparison["max_abs_delta_db"])
        if comparison["fingerprint_changed"] or max_delta >= threshold:
            south_same = False
        if comparison["fingerprint_changed"] and max_delta >= threshold and stations[station]["near_fault_epoch_boundaries"]:
            south_encodes = True
    if south_same:
        return "PUBLIC_STATIONXML_DOES_NOT_ENCODE_KNOWN_FAULT_ERA_CHANGE", diagnostics
    if south_encodes:
        return "PUBLIC_STATIONXML_ENCODES_FAULT_ERA_RESPONSE_CHANGE", diagnostics
    return "MIXED_OR_PARTIAL_PUBLIC_RESPONSE_EPOCH_ENCODING", diagnostics


def run(output: Path) -> dict[str, Any]:
    if DATA in output.resolve().parents or output.resolve() == DATA.resolve():
        raise RuntimeError("CANONICAL_DATA_OUTPUT_FORBIDDEN_USE_EPHEMERAL_ARTIFACT_PATH")
    if git_blob_sha1_file(PROTOCOL) != EXPECTED_PROTOCOL_BLOB:
        raise RuntimeError("PROTOCOL_BLOB_MISMATCH")
    if git_blob_sha1_file(REVERSE_AUDIT) != EXPECTED_REVERSE_AUDIT_BLOB:
        raise RuntimeError("REVERSE_AUDIT_BLOB_MISMATCH")
    if git_blob_sha1_file(FROZEN_119) != EXPECTED_FROZEN_119_BLOB:
        raise RuntimeError("FROZEN_119_BLOB_MISMATCH")
    frozen = json.loads(FROZEN_119.read_text(encoding="utf-8"))
    if frozen.get("verdict") != EXPECTED_FROZEN_119_VERDICT:
        raise RuntimeError("FROZEN_119_VERDICT_MISMATCH")
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    frequencies_hz = [float(x) for x in protocol["frequency_probe_contract"]["frequencies_hz"]]
    stations_to_run = protocol["public_metadata_contract"]["north_stations"] + protocol["public_metadata_contract"]["south_stations"]
    anchors = protocol["frozen_anchor_times_utc"]
    fault_date = protocol["historical_anchor"]["fault_onset_utc_date"]
    window_days = int(protocol["classification_contract"]["fault_epoch_boundary_window_days"])

    session = requests.Session()
    session.headers.update({"User-Agent": "JANUS-Cousteau-Public-Response-Epoch-Audit/1.0"})
    station_results: dict[str, Any] = {}
    for station in stations_to_run:
        station_row: dict[str, Any] = {"station": station, "anchors": {}}
        try:
            inv, transport = fetch_station_inventory(session, protocol, station)
            epochs = channel_epochs(inv, station)
            station_row["stationxml"] = transport
            station_row["public_channel_epochs"] = epochs
            station_row["near_fault_epoch_boundaries"] = near_fault_boundaries(epochs, fault_date, window_days)
            for anchor_name, when in anchors.items():
                try:
                    station_row["anchors"][anchor_name] = anchor_row(inv, station, when, frequencies_hz)
                except Exception as exc:
                    station_row["anchors"][anchor_name] = {
                        "status": "UNRESOLVED",
                        "time_utc": when,
                        "error": f"{type(exc).__name__}:{exc}",
                    }
            pre = station_row["anchors"].get("PRE_FAULT_DAY", {})
            post = station_row["anchors"].get("POST_SEP_2013", {})
            later = station_row["anchors"].get("JANUS_2014_REFERENCE", {})
            if pre.get("status") == "ANALYZED" and post.get("status") == "ANALYZED":
                station_row["pre_to_post_sep_2013"] = curve_delta_db(post, pre)
            if pre.get("status") == "ANALYZED" and later.get("status") == "ANALYZED":
                station_row["pre_to_janus_2014"] = curve_delta_db(later, pre)
        except Exception as exc:
            station_row["stationxml"] = {"status": "FETCH_OR_PARSE_FAILED", "error": f"{type(exc).__name__}:{exc}"}
            station_row["public_channel_epochs"] = []
            station_row["near_fault_epoch_boundaries"] = []
        station_results[station] = station_row

    verdict, completeness = decide_verdict(station_results, protocol)
    receipt = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-HA10-PUBLIC-RESPONSE-EPOCH-AUDIT-RUN-001",
        "status": "EXECUTED_EPHEMERAL_RECEIPT",
        "started_completed_utc": utc_now(),
        "protocol_path": str(PROTOCOL.relative_to(ROOT)).replace("\\", "/"),
        "protocol_git_blob_sha1": EXPECTED_PROTOCOL_BLOB,
        "reverse_spiral_audit_git_blob_sha1": EXPECTED_REVERSE_AUDIT_BLOB,
        "frozen_119hz_result_git_blob_sha1": EXPECTED_FROZEN_119_BLOB,
        "frozen_119hz_verdict": EXPECTED_FROZEN_119_VERDICT,
        "historical_fault_anchor_utc_date": fault_date,
        "stations": station_results,
        "completeness": completeness,
        "verdict": verdict,
        "interpretation": {
            "public_stationxml_equals_ctbto_internal_calibration_database": False,
            "absence_of_public_epoch_change_proves_no_hardware_change": False,
            "presence_of_epoch_change_proves_fault_cause": False,
            "119hz_authority_delta": 0,
            "target_identity": "UNCONFIRMED",
            "next_step": "SUBMIT_EXACT_RECEIPT_TO_JANUS_DEMIURGE_REVERSE_COUNCIL",
        },
        "source_writeback": False,
        "claim_ceiling": protocol["claim_ceiling"],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise RuntimeError("OUTPUT_ALREADY_EXISTS_APPEND_ONLY_RUN_REQUIRED")
    output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    receipt = run(Path(args.output))
    print(json.dumps({"verdict": receipt["verdict"], "completeness": receipt["completeness"]}, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
