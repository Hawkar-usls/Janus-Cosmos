from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "workspace" / "cousteau_ha10_public_response_epoch_audit_v2.py"
spec = importlib.util.spec_from_file_location("cousteau_response_epoch_v2", MODULE_PATH)
assert spec and spec.loader
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_v1_helper_is_exactly_bound():
    assert m.m.git_blob_sha1_file(m.HELPER) == m.EXPECTED_HELPER_BLOB


def test_nested_frozen_verdict_binding_accepts_exact_scalar():
    assert m.contains_exact_scalar({"a": [{"b": m.EXPECTED_FROZEN_119_VERDICT}]}, m.EXPECTED_FROZEN_119_VERDICT)
    assert not m.contains_exact_scalar({"a": [{"b": "OTHER"}]}, m.EXPECTED_FROZEN_119_VERDICT)


def test_actual_exact_frozen_bytes_contain_expected_verdict():
    frozen = json.loads(m.m.FROZEN_119.read_text(encoding="utf-8"))
    assert m.m.git_blob_sha1_file(m.m.FROZEN_119) == m.m.EXPECTED_FROZEN_119_BLOB
    assert m.contains_exact_scalar(frozen, m.EXPECTED_FROZEN_119_VERDICT)


def test_verify_frozen_inputs_passes_without_network():
    frozen, protocol = m.verify_frozen_inputs()
    assert m.contains_exact_scalar(frozen, m.EXPECTED_FROZEN_119_VERDICT)
    assert protocol["gate_id"] == "COUSTEAU_HA10_PUBLIC_RESPONSE_EPOCH_AUDIT_V1"


def test_canonical_data_output_is_refused_before_network_access():
    with pytest.raises(RuntimeError, match="CANONICAL_DATA_OUTPUT_FORBIDDEN"):
        m.run(m.m.DATA / "SHOULD_NOT_EXIST_V2.json")
