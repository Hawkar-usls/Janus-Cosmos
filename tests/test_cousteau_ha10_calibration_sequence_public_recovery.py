from pathlib import Path
import importlib.util

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "workspace" / "cousteau_ha10_calibration_sequence_public_recovery.py"
SPEC = importlib.util.spec_from_file_location("cousteau_calibration_public_recovery", MODULE_PATH)
assert SPEC and SPEC.loader
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)


def test_frozen_contract_and_authority():
    protocol, frozen = gate.verify_contracts()
    assert protocol["gate_id"] == "COUSTEAU_HA10_CALIBRATION_SEQUENCE_PUBLIC_RECOVERY_V1"
    assert protocol["authority"]["authority_delta_for_119hz"] == 0
    assert frozen["summary"]["verdict"] == "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE"


def test_raw_sequence_binary_is_not_confused_with_presentation():
    row = {
        "status": "FETCHED",
        "requested_url": "https://example.ctbto.org/H10S_calibration_sequence.mseed",
        "final_url": "https://example.ctbto.org/H10S_calibration_sequence.mseed",
        "content_type": "application/vnd.fdsn.mseed",
        "data": b"binary",
    }
    cls, evidence = gate.classify(row)
    assert cls == "RAW_CALIBRATION_SEQUENCE_BYTES"
    assert "RAW_SAMPLE_PAYLOAD_SHAPE" in evidence


def test_numeric_response_requires_machine_readable_numeric_payload():
    text = "H10S calibration transfer function frequency_db " + " ".join(f"{i} {0.1*i}" for i in range(100))
    row = {
        "status": "FETCHED",
        "requested_url": "https://www.ctbto.org/H10S_calibration.csv",
        "final_url": "https://www.ctbto.org/H10S_calibration.csv",
        "content_type": "text/csv",
        "data": text.encode(),
    }
    cls, _ = gate.classify(row)
    assert cls == "MACHINE_READABLE_NUMERIC_RESPONSE"


def test_vdec_page_is_access_pointer_not_data():
    text = "vDEC request access to IMS hydroacoustic data under contract and confidentiality conditions"
    row = {
        "status": "FETCHED",
        "requested_url": "https://www.ctbto.org/resources/for-researchers-experts/vdec",
        "final_url": "https://www.ctbto.org/resources/for-researchers-experts/vdec",
        "content_type": "text/html",
        "data": f"<html><body>{text}</body></html>".encode(),
    }
    cls, _ = gate.classify(row)
    assert cls == "ACCESS_POINTER_ONLY"


def test_published_pdf_is_derived_not_raw(monkeypatch):
    monkeypatch.setattr(gate, "payload_text", lambda row: "H10S calibration 0.8396 dB cross-talk transfer function")
    row = {
        "status": "FETCHED",
        "requested_url": "https://presentations.copernicus.org/EGU2020/EGU2020-5481_presentation.pdf",
        "final_url": "https://presentations.copernicus.org/EGU2020/EGU2020-5481_presentation.pdf",
        "content_type": "application/pdf",
        "data": b"pdf",
    }
    cls, _ = gate.classify(row)
    assert cls == "PUBLISHED_DERIVED_CALIBRATION_FIGURE_OR_PRESENTATION"


def test_public_not_found_does_not_become_custodian_absence():
    protocol, _ = gate.verify_contracts()
    assert protocol["claim_ceiling"]["public_not_found_is_not_custodian_absence"] is True
    assert protocol["claim_ceiling"]["derived_curve_is_not_raw_sequence"] is True


def test_canonical_output_is_refused():
    with pytest.raises(RuntimeError, match="CANONICAL_DATA_WRITE_FORBIDDEN"):
        gate.run(gate.DATA / "forbidden.json")
