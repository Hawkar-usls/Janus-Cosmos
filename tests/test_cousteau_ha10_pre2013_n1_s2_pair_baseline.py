from __future__ import annotations

import importlib.util
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "workspace" / "cousteau_ha10_pre2013_n1_s2_pair_baseline.py"
spec = importlib.util.spec_from_file_location("cousteau_pre2013_pair", MODULE_PATH)
assert spec and spec.loader
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_frozen_protocol_and_authority():
    p, ref, frozen = m.verify_inputs()
    assert p["gate_id"] == "COUSTEAU_HA10_PRE2013_N1_S2_PAIR_BASELINE_V1"
    assert p["status"] == "PREREGISTERED_BEFORE_PRE2013_PAIR_WAVEFORM_VALUE_INSPECTION"
    assert p["channels"] == ["IM.H10N1..EDH", "IM.H10S2..EDH"]
    assert p["spectral_contract"]["119hz_excluded"] is True
    assert p["authority"]["authority_delta_for_119hz"] == 0
    assert ref and frozen


def test_month_chunks_cover_interval_without_overlap():
    start = datetime(2012, 12, 15, tzinfo=timezone.utc)
    end = datetime(2013, 3, 2, tzinfo=timezone.utc)
    chunks = list(m.month_chunks(start, end))
    assert chunks[0][0] == start
    assert chunks[-1][1] == end
    for left, right in zip(chunks, chunks[1:]):
        assert left[1] == right[0]


def test_scan_horizon_is_strictly_prefault():
    p, _, _ = m.verify_inputs()
    assert p["discovery_contract"]["scan_end_utc_exclusive"] == "2013-07-19T00:00:00Z"
    assert p["discovery_contract"]["selection_uses_waveform_sample_values"] is False


def test_canonical_output_is_forbidden():
    bad = ROOT / "data" / "cousteau" / "SHOULD_NOT_WRITE_PAIR.json"
    try:
        m.v1.ensure_noncanonical_output(bad)
    except RuntimeError as exc:
        assert "CANONICAL_DATA_WRITE_FORBIDDEN" in str(exc)
    else:
        raise AssertionError("canonical output must be rejected")
