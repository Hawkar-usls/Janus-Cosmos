#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "workspace" / "cousteau_ha10_public_response_epoch_audit.py"
spec = importlib.util.spec_from_file_location("cousteau_response_epoch_v1", HELPER)
if spec is None or spec.loader is None:
    raise RuntimeError("CANNOT_LOAD_V1_HELPER")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

EXPECTED_HELPER_BLOB = "75c98757b7232b9c0c218d8d24d99ae5ad60d16e"
EXPECTED_FROZEN_119_VERDICT = "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE"


def contains_exact_scalar(value: Any, needle: str) -> bool:
    if isinstance(value, dict):
        return any(contains_exact_scalar(v, needle) for v in value.values())
    if isinstance(value, list):
        return any(contains_exact_scalar(v, needle) for v in value)
    return value == needle


def verify_frozen_inputs() -> tuple[dict[str, Any], dict[str, Any]]:
    if m.git_blob_sha1_file(HELPER) != EXPECTED_HELPER_BLOB:
        raise RuntimeError("V1_HELPER_BLOB_MISMATCH")
    if m.git_blob_sha1_file(m.PROTOCOL) != m.EXPECTED_PROTOCOL_BLOB:
        raise RuntimeError("PROTOCOL_BLOB_MISMATCH")
    if m.git_blob_sha1_file(m.REVERSE_AUDIT) != m.EXPECTED_REVERSE_AUDIT_BLOB:
        raise RuntimeError("REVERSE_AUDIT_BLOB_MISMATCH")
    if m.git_blob_sha1_file(m.FROZEN_119) != m.EXPECTED_FROZEN_119_BLOB:
        raise RuntimeError("FROZEN_119_BLOB_MISMATCH")
    frozen = json.loads(m.FROZEN_119.read_text(encoding="utf-8"))
    if not contains_exact_scalar(frozen, EXPECTED_FROZEN_119_VERDICT):
        raise RuntimeError("FROZEN_119_VERDICT_NOT_BOUND_IN_EXACT_BYTES")
    protocol = json.loads(m.PROTOCOL.read_text(encoding="utf-8"))
    return frozen, protocol


def run(output: Path) -> dict[str, Any]:
    if m.DATA in output.resolve().parents or output.resolve() == m.DATA.resolve():
        raise RuntimeError("CANONICAL_DATA_OUTPUT_FORBIDDEN_USE_EPHEMERAL_ARTIFACT_PATH")
    _frozen, protocol = verify_frozen_inputs()

    frequencies_hz = [float(x) for x in protocol["frequency_probe_contract"]["frequencies_hz"]]
    stations_to_run = (
        protocol["public_metadata_contract"]["north_stations"]
        + protocol["public_metadata_contract"]["south_stations"]
    )
    anchors = protocol["frozen_anchor_times_utc"]
    fault_date = protocol["historical_anchor"]["fault_onset_utc_date"]
    window_days = int(protocol["classification_contract"]["fault_epoch_boundary_window_days"])

    session = m.requests.Session()
    session.headers.update({"User-Agent": "JANUS-Cousteau-Public-Response-Epoch-Audit/2.0"})
    station_results: dict[str, Any] = {}
    for station in stations_to_run:
        station_row: dict[str, Any] = {"station": station, "anchors": {}}
        try:
            inv, transport = m.fetch_station_inventory(session, protocol, station)
            epochs = m.channel_epochs(inv, station)
            station_row["stationxml"] = transport
            station_row["public_channel_epochs"] = epochs
            station_row["near_fault_epoch_boundaries"] = m.near_fault_boundaries(
                epochs, fault_date, window_days
            )
            for anchor_name, when in anchors.items():
                try:
                    station_row["anchors"][anchor_name] = m.anchor_row(
                        inv, station, when, frequencies_hz
                    )
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
                station_row["pre_to_post_sep_2013"] = m.curve_delta_db(post, pre)
            if pre.get("status") == "ANALYZED" and later.get("status") == "ANALYZED":
                station_row["pre_to_janus_2014"] = m.curve_delta_db(later, pre)
        except Exception as exc:
            station_row["stationxml"] = {
                "status": "FETCH_OR_PARSE_FAILED",
                "error": f"{type(exc).__name__}:{exc}",
            }
            station_row["public_channel_epochs"] = []
            station_row["near_fault_epoch_boundaries"] = []
        station_results[station] = station_row

    verdict, completeness = m.decide_verdict(station_results, protocol)
    receipt = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-HA10-PUBLIC-RESPONSE-EPOCH-AUDIT-RUN-001",
        "status": "EXECUTED_EPHEMERAL_RECEIPT",
        "runner_version": "v2",
        "started_completed_utc": m.utc_now(),
        "protocol_path": str(m.PROTOCOL.relative_to(ROOT)).replace("\\", "/"),
        "protocol_git_blob_sha1": m.EXPECTED_PROTOCOL_BLOB,
        "runner_v1_helper_git_blob_sha1": EXPECTED_HELPER_BLOB,
        "reverse_spiral_audit_git_blob_sha1": m.EXPECTED_REVERSE_AUDIT_BLOB,
        "frozen_119hz_result_git_blob_sha1": m.EXPECTED_FROZEN_119_BLOB,
        "frozen_119hz_verdict": EXPECTED_FROZEN_119_VERDICT,
        "frozen_119hz_verdict_binding": "EXACT_SCALAR_FOUND_RECURSIVELY_INSIDE_EXACT_BOUND_FROZEN_JSON_BYTES",
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
    output.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    receipt = run(Path(args.output))
    print(
        json.dumps(
            {"verdict": receipt["verdict"], "completeness": receipt["completeness"]},
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
