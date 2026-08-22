from __future__ import annotations

import importlib.util
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "workspace" / "cousteau_ha10_first_available_six_channel_baseline.py"
spec = importlib.util.spec_from_file_location("cousteau_first_six", MODULE_PATH)
assert spec and spec.loader
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_frozen_inputs_and_authority_firewall():
    protocol, reference, frozen = m.verify_frozen_inputs()
    assert protocol["gate_id"] == "COUSTEAU_HA10_FIRST_AVAILABLE_SIX_CHANNEL_BASELINE_V1"
    assert protocol["status"] == "PREREGISTERED_BEFORE_AVAILABILITY_PROBES_OR_NEW_WAVEFORM_VALUE_INSPECTION"
    assert protocol["spectral_contract"]["119hz_excluded"] is True
    assert protocol["authority"]["authority_delta_for_119hz"] == 0
    assert frozen["summary"]["verdict"] == "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE"
    assert reference


def test_selection_body_is_fixed_and_value_blind():
    start = datetime(2010, 1, 2, 12, 0, 0, tzinfo=timezone.utc)
    end = datetime(2010, 1, 2, 12, 0, 1, tzinfo=timezone.utc)
    body = m.selection_body(["H10N1", "H10S2"], [(start, end)])
    assert body.startswith("quality=M\nformat=miniseed\n")
    assert "IM H10N1 -- EDH 2010-01-02T12:00:00 2010-01-02T12:00:01" in body
    assert "IM H10S2 -- EDH 2010-01-02T12:00:00 2010-01-02T12:00:01" in body
    assert "119" not in body


def test_header_coverage_requires_complete_contiguous_250hz_window():
    start = datetime(2010, 1, 2, 12, 0, 0, tzinfo=timezone.utc)
    end = datetime(2010, 1, 2, 12, 10, 0, tzinfo=timezone.utc)
    a = start.timestamp()
    b = end.timestamp()
    assert m.header_covers_window([(a, b, 250.0, 150000)], start, end)
    assert not m.header_covers_window([(a, b - 5.0, 250.0, 148750)], start, end)
    assert not m.header_covers_window([(a, b, 200.0, 120000)], start, end)


def test_predeclared_baseline_classifier_boundaries():
    ref = 12.479728560324869
    assert m.classify_baseline(12.0, ref, 6.0, 3.0) == "EARLY_LARGE_S2_N1_BASELINE_SIMILAR_TO_2014"
    assert m.classify_baseline(5.9, ref, 6.0, 3.0) == "EARLY_S2_N1_BASELINE_MATERIALLY_LOWER_THAN_2014"
    assert m.classify_baseline(8.0, ref, 6.0, 3.0) == "EARLY_S2_N1_BASELINE_DIFFERENT_OR_MIXED"


def test_output_under_canonical_data_tree_is_forbidden(tmp_path):
    canonical = ROOT / "data" / "cousteau" / "SHOULD_NOT_WRITE.json"
    try:
        m.ensure_noncanonical_output(canonical)
    except RuntimeError as exc:
        assert "CANONICAL_DATA_WRITE_FORBIDDEN" in str(exc)
    else:
        raise AssertionError("canonical output path must be rejected")
    m.ensure_noncanonical_output(tmp_path / "receipt.json")
