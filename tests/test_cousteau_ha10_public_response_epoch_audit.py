from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "workspace" / "cousteau_ha10_public_response_epoch_audit.py"
spec = importlib.util.spec_from_file_location("cousteau_response_epoch", MODULE_PATH)
assert spec and spec.loader
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def protocol():
    return json.loads(m.PROTOCOL.read_text(encoding="utf-8"))


def analyzed_anchor(fp: str, value: float = 10.0):
    return {
        "status": "ANALYZED",
        "response_fingerprint_sha256": fp,
        "channel_epoch_start_utc": "2010-01-01T00:00:00.000000Z",
        "channel_epoch_end_utc": None,
        "fixed_frequency_response_magnitude_db": {
            "1": value,
            "10": value,
            "30": value,
            "55": value,
            "80": value,
            "100": value,
        },
    }


def station_row(pre_fp="same", post_fp="same", post_delta=0.0, near=False):
    pre = analyzed_anchor(pre_fp, 10.0)
    post = analyzed_anchor(post_fp, 10.0 + post_delta)
    return {
        "anchors": {"PRE_FAULT_DAY": pre, "POST_SEP_2013": post},
        "pre_to_post_sep_2013": m.curve_delta_db(post, pre),
        "near_fault_epoch_boundaries": ([{"field": "start_utc", "delta_days": 1.0}] if near else []),
    }


def full_station_map(south_row=None, north_row=None):
    p = protocol()
    south_row = south_row or station_row()
    north_row = north_row or station_row()
    out = {}
    for station in p["public_metadata_contract"]["north_stations"]:
        out[station] = json.loads(json.dumps(north_row))
    for station in p["public_metadata_contract"]["south_stations"]:
        out[station] = json.loads(json.dumps(south_row))
    return out


def test_exact_frozen_git_blob_bindings():
    assert m.git_blob_sha1_file(m.PROTOCOL) == m.EXPECTED_PROTOCOL_BLOB
    assert m.git_blob_sha1_file(m.REVERSE_AUDIT) == m.EXPECTED_REVERSE_AUDIT_BLOB
    assert m.git_blob_sha1_file(m.FROZEN_119) == m.EXPECTED_FROZEN_119_BLOB


def test_curve_delta_definition():
    pre = analyzed_anchor("a", 10.0)
    post = analyzed_anchor("b", 10.25)
    result = m.curve_delta_db(post, pre)
    assert result["fingerprint_changed"] is True
    assert result["max_abs_delta_db"] == pytest.approx(0.25)
    assert all(v == pytest.approx(0.25) for v in result["per_frequency_delta_db"].values())


def test_no_public_fault_epoch_encoding_when_all_south_identical():
    verdict, diag = m.decide_verdict(full_station_map(), protocol())
    assert verdict == "PUBLIC_STATIONXML_DOES_NOT_ENCODE_KNOWN_FAULT_ERA_CHANGE"
    assert len(diag["complete_south"]) == 3
    assert len(diag["complete_north"]) == 3


def test_public_fault_epoch_encoding_requires_curve_change_and_near_boundary():
    stations = full_station_map()
    stations["H10S1"] = station_row(pre_fp="pre", post_fp="post", post_delta=0.8, near=True)
    verdict, _ = m.decide_verdict(stations, protocol())
    assert verdict == "PUBLIC_STATIONXML_ENCODES_FAULT_ERA_RESPONSE_CHANGE"


def test_changed_response_without_near_fault_boundary_is_mixed():
    stations = full_station_map()
    stations["H10S1"] = station_row(pre_fp="pre", post_fp="post", post_delta=0.8, near=False)
    verdict, _ = m.decide_verdict(stations, protocol())
    assert verdict == "MIXED_OR_PARTIAL_PUBLIC_RESPONSE_EPOCH_ENCODING"


def test_missing_required_channel_blocks():
    stations = full_station_map()
    stations["H10S3"]["anchors"]["POST_SEP_2013"] = {"status": "UNRESOLVED"}
    stations["H10S3"].pop("pre_to_post_sep_2013", None)
    verdict, _ = m.decide_verdict(stations, protocol())
    assert verdict == "BLOCKED_PUBLIC_RESPONSE_EPOCH_AUDIT"


def test_canonical_data_output_is_refused_before_network_access(tmp_path):
    with pytest.raises(RuntimeError, match="CANONICAL_DATA_OUTPUT_FORBIDDEN"):
        m.run(m.DATA / "SHOULD_NOT_EXIST.json")
