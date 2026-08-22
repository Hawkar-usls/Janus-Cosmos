#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import timedelta
from pathlib import Path
from typing import Any

import requests

from workspace import cousteau_ha10_n1_s2_independent_scalar_collapse_holdout as prior
from workspace import cousteau_ha10_tphase_inband_positive_control_v2 as transport_v2

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cousteau"
PROTOCOL = DATA / "JANUS-ECHO-COUSTEAU-HA10-N1-S2-JANUARY-0600UTC-SAME-SCALAR-HOLDOUT-PROTOCOL-2026-08-22-v1.0.json"
JAN1800_SUMMARY = DATA / "JANUS-ECHO-COUSTEAU-HA10-N1-S2-OUT-OF-MONTH-SCALAR-HOLDOUT-RUN-001-SUMMARY-2026-08-22-v1.0.json"
DEC_SCALAR_SUMMARY = DATA / "JANUS-ECHO-COUSTEAU-HA10-N1-S2-INDEPENDENT-SCALAR-COLLAPSE-HOLDOUT-RUN-001-SUMMARY-2026-08-22-v1.0.json"
FROZEN_119 = DATA / "JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0.json"
PRIOR_RUNNER = ROOT / "workspace" / "cousteau_ha10_n1_s2_independent_scalar_collapse_holdout.py"
TRANSPORT_V2 = ROOT / "workspace" / "cousteau_ha10_tphase_inband_positive_control_v2.py"

EXPECTED_PROTOCOL_BLOB = "3bfd9adfd7fc305081c5d239ca100eee94c78686"
EXPECTED_JAN1800_SUMMARY_BLOB = "09ea4c580d93e88c152b81989a6316153e28d8f7"
EXPECTED_DEC_SCALAR_SUMMARY_BLOB = "6b2c075396b3fc60df3fafe41024a72c86a3afff"
EXPECTED_FROZEN_119_BLOB = "eb8b48fb7f043160c057f9df6264a781412ed854"
EXPECTED_PRIOR_RUNNER_BLOB = "5e1e6ea9886ce4c45bdf4f1a3fbb9155395fdc7f"
EXPECTED_TRANSPORT_V2_BLOB = "c3d37fe4f5a75514390862be4c1d2870c4c01fe9"
EXPECTED_JAN1800_VERDICT = "OUT_OF_MONTH_SCALAR_HOLDOUT_FAIL"
EXPECTED_DEC_SCALAR_VERDICT = "SCALAR_COLLAPSE_HOLDOUT_PASS"
EXPECTED_FROZEN_119_VERDICT = "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE"
EXPECTED_SCALAR_DB = 12.616188132659623


def git_blob_sha1_file(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def ensure_noncanonical_output(path: Path) -> None:
    try:
        path.resolve().relative_to(DATA.resolve())
    except ValueError:
        return
    raise RuntimeError("CANONICAL_DATA_WRITE_FORBIDDEN_BY_JANUARY_0600_SAME_SCALAR")


def verify_frozen_contracts() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    expected = {
        PROTOCOL: EXPECTED_PROTOCOL_BLOB,
        JAN1800_SUMMARY: EXPECTED_JAN1800_SUMMARY_BLOB,
        DEC_SCALAR_SUMMARY: EXPECTED_DEC_SCALAR_SUMMARY_BLOB,
        FROZEN_119: EXPECTED_FROZEN_119_BLOB,
        PRIOR_RUNNER: EXPECTED_PRIOR_RUNNER_BLOB,
        TRANSPORT_V2: EXPECTED_TRANSPORT_V2_BLOB,
    }
    for path, blob in expected.items():
        if git_blob_sha1_file(path) != blob:
            raise RuntimeError(f"FROZEN_BLOB_DRIFT:{path.name}")

    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    jan1800 = json.loads(JAN1800_SUMMARY.read_text(encoding="utf-8"))
    decscalar = json.loads(DEC_SCALAR_SUMMARY.read_text(encoding="utf-8"))
    frozen = json.loads(FROZEN_119.read_text(encoding="utf-8"))

    if protocol.get("status") != "PREREGISTERED_AFTER_JANUARY_1800_FAIL_BEFORE_JANUARY_0600_WINDOWS_ARE_DOWNLOADED":
        raise RuntimeError("JANUARY_0600_PROTOCOL_NOT_PREREGISTERED")
    if jan1800.get("result", {}).get("verdict") != EXPECTED_JAN1800_VERDICT:
        raise RuntimeError("JANUARY_1800_FAIL_BINDING_DRIFT")
    if decscalar.get("result", {}).get("verdict") != EXPECTED_DEC_SCALAR_VERDICT:
        raise RuntimeError("DECEMBER_SCALAR_PASS_BINDING_DRIFT")
    if frozen.get("summary", {}).get("verdict") != EXPECTED_FROZEN_119_VERDICT:
        raise RuntimeError("FROZEN_119_VERDICT_DRIFT")
    scalar = float(protocol["frozen_scalar"]["source_value_db"])
    upstream_scalar = float(decscalar["frozen_contract"]["frozen_scalar_db"])
    if scalar != EXPECTED_SCALAR_DB or scalar != upstream_scalar:
        raise RuntimeError("FROZEN_SCALAR_BINDING_DRIFT")
    if protocol["frozen_scalar"]["holdout_refit_allowed"] is not False:
        raise RuntimeError("HOLDOUT_REFIT_MUST_BE_FORBIDDEN")
    if protocol["frozen_scalar"]["alternate_scalar_search_allowed"] is not False:
        raise RuntimeError("ALTERNATE_SCALAR_SEARCH_MUST_BE_FORBIDDEN")
    if protocol["epistemic_position"]["may_erase_prior_january_1800_fail"] is not False:
        raise RuntimeError("PRIOR_FAIL_MUST_REMAIN_IMMUTABLE")
    if protocol["epistemic_position"]["authority_delta_for_119hz"] != 0:
        raise RuntimeError("AUTHORITY_DELTA_FOR_119HZ_NOT_ZERO")
    wc = protocol["holdout_window_contract"]
    if len(wc["dates"]) != 20 or wc["dates"][0] != "2015-01-01" or wc["dates"][-1] != "2015-01-20":
        raise RuntimeError("JANUARY_DATE_BLOCK_DRIFT")
    if wc["window_start_time_utc_each_day"] != "06:00:00Z":
        raise RuntimeError("JANUARY_0600_CLOCK_DRIFT")
    if int(wc["minimum_complete_paired_windows"]) != 15:
        raise RuntimeError("MINIMUM_COMPLETE_PAIRS_DRIFT")
    if protocol["frequency_contract"]["119hz_or_117_121hz_bins_may_be_used"] is not False:
        raise RuntimeError("TARGET_BAND_EXCLUSION_NOT_FROZEN")
    expected_thresholds = {
        "maximum_absolute_broadband_residual_median_db": 2.0,
        "maximum_broadband_residual_iqr_db": 3.0,
        "window_absolute_residual_tolerance_db": 3.0,
        "minimum_fraction_windows_within_tolerance": 0.8,
        "maximum_absolute_fixed_subband_residual_median_db": 2.5,
    }
    for key, value in expected_thresholds.items():
        if float(protocol["classification_thresholds"][key]) != value:
            raise RuntimeError(f"THRESHOLD_DRIFT:{key}")
    return protocol, jan1800, decscalar, frozen


def map_verdict(prior_verdict: str) -> str:
    mapping = {
        "BLOCKED_SCALAR_COLLAPSE_HOLDOUT_DATA_ACCESS": "BLOCKED_JANUARY_0600_SAME_SCALAR_DATA_ACCESS",
        "SCALAR_COLLAPSE_HOLDOUT_PASS": "JANUARY_0600_SAME_SCALAR_PASS",
        "SCALAR_COLLAPSE_HOLDOUT_FAIL": "JANUARY_0600_SAME_SCALAR_FAIL",
    }
    if prior_verdict not in mapping:
        raise RuntimeError(f"UNEXPECTED_PRIOR_DECISION:{prior_verdict}")
    return mapping[prior_verdict]


def run(output: Path) -> dict[str, Any]:
    ensure_noncanonical_output(output)
    protocol, jan1800, decscalar, frozen = verify_frozen_contracts()
    wc = protocol["holdout_window_contract"]
    estimator = protocol["spectral_estimator"]
    thresholds = protocol["classification_thresholds"]
    scalar_db = float(protocol["frozen_scalar"]["source_value_db"])
    nperseg = int(estimator["welch_nperseg_samples"])
    noverlap = int(estimator["welch_noverlap_samples"])
    duration_s = int(wc["window_duration_s"])
    minimum_pairs = int(wc["minimum_complete_paired_windows"])
    channels = list(protocol["channels"])
    if channels != ["IM.H10N1..EDH", "IM.H10S2..EDH"]:
        raise RuntimeError("CHANNEL_SET_DRIFT")

    bands: dict[str, list[float]] = {
        "broadband": [float(x) for x in protocol["frequency_contract"]["broadband_hz"]],
        **{name: [float(x) for x in band] for name, band in protocol["frequency_contract"]["fixed_subbands_hz"].items()},
    }

    session = requests.Session()
    session.headers["User-Agent"] = "Janus-Echo-Cousteau/1.0 January 06UTC same-scalar discriminator"
    first = prior.parse_start(wc["dates"][0], wc["window_start_time_utc_each_day"])
    last = prior.parse_start(wc["dates"][-1], wc["window_start_time_utc_each_day"])
    inv_start = (first - timedelta(days=1)).isoformat().replace("+00:00", "Z")
    inv_end = (last + timedelta(days=2)).isoformat().replace("+00:00", "Z")

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

    windows: list[dict[str, Any]] = []
    complete: list[dict[str, Any]] = []
    data_errors: list[dict[str, Any]] = []
    for date_text in wc["dates"]:
        start_dt = prior.parse_start(date_text, wc["window_start_time_utc_each_day"])
        end_dt = start_dt + timedelta(seconds=duration_s)
        start = start_dt.isoformat().replace("+00:00", "Z")
        end = end_dt.isoformat().replace("+00:00", "Z")
        row: dict[str, Any] = {"date": date_text, "start_utc": start, "end_utc": end, "stations": {}, "pair_complete": False}
        for cid in channels:
            try:
                if cid not in inventories:
                    raise RuntimeError("RESPONSE_UNAVAILABLE:" + inventory_errors.get(cid, "UNKNOWN"))
                trace, waveform = transport_v2.fetch_trace_bounded(session, cid, start, end)
                powers = prior.powers_for_bands(trace, inventories[cid], bands=bands, nperseg=nperseg, noverlap=noverlap)
                row["stations"][cid] = {"data_status": "ANALYZED", "waveform": waveform, "band_powers": powers}
            except Exception as exc:
                row["stations"][cid] = {"data_status": "BLOCKED", "error_type": type(exc).__name__, "error": str(exc)}
                data_errors.append({"date": date_text, "station": cid, "error_type": type(exc).__name__, "error": str(exc)})

        n1 = row["stations"].get(channels[0], {})
        s2 = row["stations"].get(channels[1], {})
        if n1.get("data_status") == "ANALYZED" and s2.get("data_status") == "ANALYZED":
            row["pair_complete"] = True
            row["broadband"] = {}
            row["subbands"] = {}
            for name in bands:
                n1p = n1["band_powers"][name]
                s2p = s2["band_powers"][name]
                raw_ratio = prior.ratio_db(float(s2p["raw_integrated_power"]), float(n1p["raw_integrated_power"]))
                corr_ratio = prior.ratio_db(float(s2p["corrected_integrated_power_pa2"]), float(n1p["corrected_integrated_power_pa2"]))
                result = {
                    "band_hz": n1p["band_hz"],
                    "raw_s2_minus_n1_db": raw_ratio,
                    "corrected_s2_minus_n1_db": corr_ratio,
                    "raw_residual_db": raw_ratio - scalar_db,
                    "corrected_residual_db": corr_ratio - scalar_db,
                }
                if name == "broadband":
                    row["broadband"] = result
                else:
                    row["subbands"][name] = result
            complete.append(row)
        windows.append(row)

    prior_verdict, broadband, subbands, gates = prior.decide(
        complete,
        minimum_pairs=minimum_pairs,
        max_abs_median_db=float(thresholds["maximum_absolute_broadband_residual_median_db"]),
        max_iqr_db=float(thresholds["maximum_broadband_residual_iqr_db"]),
        tolerance_db=float(thresholds["window_absolute_residual_tolerance_db"]),
        minimum_fraction=float(thresholds["minimum_fraction_windows_within_tolerance"]),
        max_abs_subband_median_db=float(thresholds["maximum_absolute_fixed_subband_residual_median_db"]),
    )
    verdict = map_verdict(prior_verdict)
    cross_gate = (
        "SUPPORTS_DAYPART_OR_REGIME_DEPENDENCE_OF_WINDOW_LEVEL_RESIDUAL_VARIABILITY"
        if verdict == "JANUARY_0600_SAME_SCALAR_PASS"
        else "SUPPORTS_BROADER_OUT_OF_MONTH_NONSTATIONARITY_BEYOND_ONE_DAYPART"
        if verdict == "JANUARY_0600_SAME_SCALAR_FAIL"
        else "BLOCKED_NO_DISCRIMINATION"
    )

    receipt = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-HA10-N1-S2-JANUARY-0600UTC-SAME-SCALAR-HOLDOUT-RUN",
        "created_utc": prior.utc_now(),
        "gate_id": protocol["gate_id"],
        "protocol_git_blob_sha1": EXPECTED_PROTOCOL_BLOB,
        "january_1800_fail_summary_git_blob_sha1": EXPECTED_JAN1800_SUMMARY_BLOB,
        "january_1800_verdict": jan1800["result"]["verdict"],
        "december_scalar_summary_git_blob_sha1": EXPECTED_DEC_SCALAR_SUMMARY_BLOB,
        "december_scalar_verdict": decscalar["result"]["verdict"],
        "frozen_119hz_result_git_blob_sha1": EXPECTED_FROZEN_119_BLOB,
        "frozen_119hz_verdict": frozen["summary"]["verdict"],
        "authority_delta_for_119hz": 0,
        "target_identity": "UNCONFIRMED",
        "target_evidence_delta": "NONE_FROM_JANUARY_0600_DISCRIMINATOR",
        "source_writeback": False,
        "frozen_scalar": protocol["frozen_scalar"],
        "frequency_contract": protocol["frequency_contract"],
        "classification_thresholds": thresholds,
        "inventory_metadata": inventory_meta,
        "inventory_errors": inventory_errors,
        "windows": windows,
        "summary": {
            "frozen_window_count": len(wc["dates"]),
            "complete_paired_windows": len(complete),
            "blocked_station_windows": len(data_errors),
            "verdict": verdict,
            "cross_gate_discriminator": cross_gate,
            "frozen_scalar_db": scalar_db,
            "broadband_residual_statistics": broadband,
            "fixed_subband_residual_statistics": subbands,
            "gate_checks": gates,
        },
        "data_errors": data_errors,
        "claim_ceiling": protocol["claim_ceiling"],
        "hard_rules": protocol["hard_rules"],
        "status": "DIAGNOSTIC_RUN_COMPLETE",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"NO_OVERWRITE:{output}")
    output.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
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
