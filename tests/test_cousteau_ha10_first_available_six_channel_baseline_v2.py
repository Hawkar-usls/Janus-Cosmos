from __future__ import annotations

import importlib.util
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "workspace" / "cousteau_ha10_first_available_six_channel_baseline_v2.py"
spec = importlib.util.spec_from_file_location("cousteau_first_six_v2", MODULE_PATH)
assert spec and spec.loader
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_transport_addendum_is_bound_and_non_promotional():
    a = m.verify_transport_addendum()
    assert a["failed_execution"]["github_actions_run_id"] == 32578915494
    assert a["failed_execution"]["scientific_result_produced"] is False
    assert a["transport_only_repair"]["waveform_value_blinding_unchanged"] is True
    assert a["transport_only_repair"]["scientific_thresholds_unchanged"] is True
    assert a["epistemic_firewall"]["authority_delta_for_119hz"] == 0


def test_month_chunks_are_chronological_complete_and_nonoverlapping():
    start = datetime(2004, 9, 14, tzinfo=timezone.utc)
    end = datetime(2004, 12, 15, tzinfo=timezone.utc)
    chunks = list(m.iter_month_chunks(start, end))
    assert chunks == [
        (datetime(2004, 9, 14, tzinfo=timezone.utc), datetime(2004, 10, 1, tzinfo=timezone.utc)),
        (datetime(2004, 10, 1, tzinfo=timezone.utc), datetime(2004, 11, 1, tzinfo=timezone.utc)),
        (datetime(2004, 11, 1, tzinfo=timezone.utc), datetime(2004, 12, 1, tzinfo=timezone.utc)),
        (datetime(2004, 12, 1, tzinfo=timezone.utc), datetime(2004, 12, 15, tzinfo=timezone.utc)),
    ]
    assert chunks[0][0] == start
    assert chunks[-1][1] == end
    for left, right in zip(chunks, chunks[1:]):
        assert left[1] == right[0]


def test_naive_month_chunk_input_is_rejected():
    try:
        list(m.iter_month_chunks(datetime(2004, 1, 1), datetime(2004, 2, 1)))
    except ValueError as exc:
        assert "AWARE_UTC" in str(exc)
    else:
        raise AssertionError("naive timestamps must be rejected")
