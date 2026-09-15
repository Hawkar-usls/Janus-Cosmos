from pathlib import Path

import pytest

from workspace import cousteau_ha10_pre_post_2013_crosstalk_triplet_gate as gate


def _protocol():
    protocol, _ = gate.verify_frozen_contracts()
    return protocol


def _row(*, south_coh, north_coh, south_corr, north_corr, s1=0.0, s2=0.0, s3=0.0):
    return {
        "south_driver_coherence": south_coh,
        "north_driver_coherence": north_coh,
        "south_driver_abs_correlation": south_corr,
        "north_driver_abs_correlation": north_corr,
        "raw_s1_minus_n1_db": s1,
        "raw_s2_minus_n1_db": s2,
        "raw_s3_minus_n1_db": s3,
        "corrected_s1_minus_n1_db": s1,
        "corrected_s2_minus_n1_db": s2,
        "corrected_s3_minus_n1_db": s3,
    }


def test_frozen_contracts_bind_protocol_and_negative():
    protocol, frozen = gate.verify_frozen_contracts()
    assert protocol["gate_id"] == "COUSTEAU_HA10_PRE_POST_2013_CROSSTALK_TRIPLET_GATE_V1"
    assert protocol["authority"]["authority_delta_for_119hz"] == 0
    assert frozen["summary"]["verdict"] == "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE"


def test_supported_requires_all_four_coupling_checks():
    protocol = _protocol()
    pre = [_row(south_coh=0.10, north_coh=0.10, south_corr=0.10, north_corr=0.10) for _ in range(12)]
    post = [_row(south_coh=0.35, north_coh=0.12, south_corr=0.35, north_corr=0.12) for _ in range(12)]
    verdict, diag = gate.decide(pre, post, protocol)
    assert verdict == "SUPPORTED_POST_2013_SOUTH_TRIPLET_COUPLING_INCREASE"
    assert all(diag["gate_checks"].values())


def test_north_control_blocks_false_common_increase():
    protocol = _protocol()
    pre = [_row(south_coh=0.10, north_coh=0.10, south_corr=0.10, north_corr=0.10) for _ in range(12)]
    post = [_row(south_coh=0.30, north_coh=0.28, south_corr=0.30, north_corr=0.28) for _ in range(12)]
    verdict, diag = gate.decide(pre, post, protocol)
    assert verdict == "NO_STRONG_POST_2013_SOUTH_TRIPLET_COUPLING_INCREASE"
    assert diag["gate_checks"]["coherence_difference_in_difference"] is False
    assert diag["gate_checks"]["abs_correlation_difference_in_difference"] is False


def test_blocked_when_epoch_has_too_few_complete_windows():
    protocol = _protocol()
    pre = [_row(south_coh=0.1, north_coh=0.1, south_corr=0.1, north_corr=0.1) for _ in range(11)]
    post = [_row(south_coh=0.4, north_coh=0.1, south_corr=0.4, north_corr=0.1) for _ in range(20)]
    verdict, diag = gate.decide(pre, post, protocol)
    assert verdict == "BLOCKED_PRE_POST_2013_TRIPLET_DATA_ACCESS"
    assert diag["pre_complete"] == 11


def test_canonical_output_is_refused():
    with pytest.raises(RuntimeError, match="CANONICAL_DATA_WRITE_FORBIDDEN"):
        gate.ensure_noncanonical_output(gate.DATA / "forbidden.json")


def test_ephemeral_output_is_allowed(tmp_path: Path):
    gate.ensure_noncanonical_output(tmp_path / "receipt.json")
