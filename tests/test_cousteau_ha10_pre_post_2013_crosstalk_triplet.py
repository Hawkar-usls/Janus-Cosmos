from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "workspace" / "cousteau_ha10_pre_post_2013_crosstalk_triplet.py"
spec = importlib.util.spec_from_file_location("triplet_gate", MODULE_PATH)
assert spec and spec.loader
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def protocol():
    return m.verify_frozen_inputs()


def row(*, contrast_shift=0.0, s12_shift=0.0, s13_shift=0.0, s2_n1_db=10.0):
    n12 = 0.10
    n13 = 0.10
    s12 = 0.10 + contrast_shift + s12_shift
    s13 = 0.10 + contrast_shift + s13_shift
    south_mean = (s12 + s13) / 2.0
    north_mean = (n12 + n13) / 2.0
    ratio = 10 ** (s2_n1_db / 10.0)
    raw_power = {
        "H10N1": 1.0,
        "H10N2": 1.1,
        "H10N3": 0.9,
        "H10S1": 1.2,
        "H10S2": ratio,
        "H10S3": 1.3,
    }
    corrected = {k: v / 4.0 for k, v in raw_power.items()}
    return {
        "analysis": {
            "pair_coherence_10_80": {
                "H10S1-H10S2": s12,
                "H10S1-H10S3": s13,
                "H10S2-H10S3": 0.10,
                "H10N1-H10N2": n12,
                "H10N1-H10N3": n13,
                "H10N2-H10N3": 0.10,
            },
            "south_minus_north_coherence_contrast": south_mean - north_mean,
            "raw_integrated_power_10_80": raw_power,
            "public_response_corrected_integrated_power_10_80": corrected,
            "raw_s2_minus_n1_power_ratio_db": s2_n1_db,
            "corrected_s2_minus_n1_power_ratio_db": s2_n1_db,
        }
    }


def test_exact_frozen_blob_bindings_and_janus_selection():
    p = protocol()
    assert m.git_blob_sha1_file(m.PROTOCOL) == m.EXPECTED_PROTOCOL_BLOB
    assert m.git_blob_sha1_file(m.RESPONSE_SUMMARY) == m.EXPECTED_RESPONSE_SUMMARY_BLOB
    assert m.git_blob_sha1_file(m.REVERSE_AUDIT) == m.EXPECTED_REVERSE_AUDIT_BLOB
    assert m.git_blob_sha1_file(m.FROZEN_119) == m.EXPECTED_FROZEN_119_BLOB
    assert p["selected_by_janus"]["top1_votes"] == "17/17"
    assert p["spectral_contract"]["primary_domain"] == "RAW_COUNTS"
    assert p["spectral_contract"]["119hz_excluded"] is True


def test_synthetic_linear_coherent_post_change_detects_signature():
    p = protocol()
    pre = [row(s2_n1_db=10.0) for _ in range(15)]
    post = [row(contrast_shift=0.08, s2_n1_db=10.5) for _ in range(15)]
    result = m.summarize_complete_rows(pre, post, p, permutations_override=5000)
    assert result["verdict"] == "LINEAR_COHERENT_SOUTH_TRIPLET_CROSSTALK_SIGNATURE_DETECTED"
    assert result["primary"]["post_minus_pre_contrast_delta_msc"] == pytest.approx(0.08)
    assert result["primary"]["one_sided_permutation_p"] <= 0.01
    assert all(result["primary"]["checks"].values())


def test_no_coherence_change_does_not_detect_signature():
    p = protocol()
    pre = [row(s2_n1_db=10.0) for _ in range(15)]
    post = [row(s2_n1_db=10.4) for _ in range(15)]
    result = m.summarize_complete_rows(pre, post, p, permutations_override=1000)
    assert result["verdict"] == "NO_PREREGISTERED_LINEAR_COHERENT_CROSSTALK_SIGNATURE"
    assert result["primary"]["post_minus_pre_contrast_delta_msc"] == pytest.approx(0.0)


def test_large_s2_n1_offset_preexistence_classification():
    p = protocol()
    assert m.classify_preexistence(12.0, 13.5, p) == "LARGE_S2_N1_OFFSET_PREEXISTS_AND_IS_STABLE_ACROSS_FAULT"
    assert m.classify_preexistence(2.0, 9.0, p) == "LARGE_S2_N1_OFFSET_EMERGES_AFTER_FAULT"
    assert m.classify_preexistence(12.0, 17.0, p) == "S2_N1_OFFSET_PRE_POST_MIXED"


def test_summary_reports_preexisting_offset_independently_of_primary_signature():
    p = protocol()
    pre = [row(s2_n1_db=11.8) for _ in range(15)]
    post = [row(contrast_shift=0.08, s2_n1_db=12.7) for _ in range(15)]
    result = m.summarize_complete_rows(pre, post, p, permutations_override=5000)
    assert result["verdict"] == "LINEAR_COHERENT_SOUTH_TRIPLET_CROSSTALK_SIGNATURE_DETECTED"
    assert result["preexistence"]["classification"] == "LARGE_S2_N1_OFFSET_PREEXISTS_AND_IS_STABLE_ACROSS_FAULT"
    assert result["preexistence"]["raw_pre_median_s2_minus_n1_db"] == pytest.approx(11.8)


def test_completeness_gate_blocks_without_fifteen_each_epoch():
    p = protocol()
    result = m.summarize_complete_rows([row() for _ in range(14)], [row() for _ in range(20)], p, permutations_override=100)
    assert result["verdict"] == "BLOCKED_PRE_POST_2013_TRIPLET_DATA"
    assert result["preexistence_classification"] == "BLOCKED_PRE_POST_2013_TRIPLET_DATA"


def test_permutation_is_deterministic_for_fixed_seed():
    pre = [0.0] * 15
    post = [0.1] * 15
    a = m.one_sided_permutation_p(pre, post, 2000, 20260822)
    b = m.one_sided_permutation_p(pre, post, 2000, 20260822)
    assert a == b


def test_band_power_uses_fixed_band_and_is_positive():
    import numpy as np
    f = np.linspace(0, 125, 4097)
    psd = np.ones_like(f)
    power = m.band_integrated_power(f, psd, 10.0, 80.0)
    assert math.isfinite(power)
    assert 69.9 < power < 70.1


def test_canonical_source_data_output_refused_before_network():
    with pytest.raises(RuntimeError, match="CANONICAL_DATA_OUTPUT_FORBIDDEN"):
        m.run(m.DATA / "SHOULD_NOT_EXIST_TRIPLET.json")
