import json
from pathlib import Path

import numpy as np
import pytest

from janus_cosmos.core import (
    GateConfig, analyze_image, block_shuffle_surrogate, bonferroni_alpha,
    directional_by_scale, minimum_test_nulls, normalize_image,
    phase_iaaft_surrogate, stable_seed,
)
from janus_cosmos.discovery import DiscoveryConfig, canonical_name, discover_target, select_products
from janus_cosmos.pipeline import load_manifest, parse_seeds


def synthetic(size=96):
    y, x = np.indices((size, size))
    g = np.exp(-((x-size/2)**2 + (y-size/2)**2)/(2*(size/7)**2))
    arm = np.exp(-((y - (0.5*x + 15))**2)/(2*2.0**2))
    return (g + 0.6*arm).astype(np.float32)


def test_directional_score_does_not_mutate():
    x = normalize_image(synthetic(), 128)
    before = x.copy()
    assert len(directional_by_scale(x)) == 3
    assert np.allclose(x, before)


def test_block_shuffle_preserves_values():
    x = normalize_image(synthetic(), 128)
    y = block_shuffle_surrogate(x, np.random.default_rng(42), 16)
    assert np.array_equal(np.sort(x.ravel()), np.sort(y.ravel()))


def test_phase_surrogate_preserves_marginal_distribution():
    x = normalize_image(synthetic(), 64)
    y = phase_iaaft_surrogate(x, np.random.default_rng(42), iterations=2)
    assert np.allclose(np.sort(x.ravel()), np.sort(y.ravel()), atol=1e-6)


def test_bonferroni_resolution_guard():
    alpha = bonferroni_alpha(0.05, filter_count=14, null_model_count=2)
    assert alpha == 0.05 / 28
    assert minimum_test_nulls(alpha) == 560
    assert 1 / 561 < alpha


def test_seed_is_stable_and_keyed():
    assert stable_seed('a', 1) == stable_seed('a', 1)
    assert stable_seed('a', 1) != stable_seed('a', 2)


def test_small_smoke_analysis_runs():
    cfg = GateConfig(image_size=64, calibration_nulls=8, iaaft_iterations=1, block_size=8)
    out = analyze_image(
        synthetic(64), target='SYNTH', filter_name='F555W', test_nulls=16,
        seeds=[1, 2], alpha=0.05, include_legacy=False, config=cfg,
    )
    assert 0.0 < out['phase_iaaft']['p_empirical'] <= 1.0
    assert 0.0 < out['block_shuffle']['p_empirical'] <= 1.0
    assert isinstance(out['robust_candidate'], bool)


def test_canonical_ngc_spacing_and_product_selection():
    assert canonical_name('NGC1425') == 'NGC 1425'
    rows = [
        {"filters":"F555W", "productFilename":"x_raw.fits", "dataURI":"mast:HST/x_raw.fits", "productType":"SCIENCE", "productSubGroupDescription":"RAW", "calib_level":1, "size":100},
        {"filters":"F555W", "productFilename":"x_drz.fits", "dataURI":"mast:HST/x_drz.fits", "productType":"SCIENCE", "productSubGroupDescription":"DRZ", "calib_level":3, "size":90},
        {"filters":"F814W", "productFilename":"i_drz.fits", "dataURI":"mast:HST/i_drz.fits", "productType":"SCIENCE", "productSubGroupDescription":"DRZ", "calib_level":3, "size":120},
    ]
    selected = select_products(rows, ["F555W", "F814W"])
    assert [x['productFilename'] for x in selected] == ['x_drz.fits', 'i_drz.fits']


class FakeObs:
    @staticmethod
    def query_criteria(**kwargs):
        assert kwargs['objectname'] == 'NGC 1425'
        return [1, 2]

    @staticmethod
    def get_unique_product_list(obs):
        return [
            {"filters":"F555W", "productFilename":"v_drz.fits", "dataURI":"mast:HST/v_drz.fits", "productType":"SCIENCE", "productSubGroupDescription":"DRZ", "calib_level":3, "size":100},
            {"filters":"F814W", "productFilename":"i_drz.fits", "dataURI":"mast:HST/i_drz.fits", "productType":"SCIENCE", "productSubGroupDescription":"DRZ", "calib_level":3, "size":100},
        ]


def test_discover_target_with_fake_backend():
    selected, status = discover_target('NGC1425', cfg=DiscoveryConfig(retries=1, filters=('F555W','F814W')), observations=FakeObs)
    assert status['status'] == 'OK'
    assert len(selected) == 2


def test_manifest_and_seed_validation(tmp_path: Path):
    p = tmp_path / 'm.json'
    p.write_text(json.dumps({"targets":[{"target":"T","filters":[{"filter":"F","url":"https://example.invalid/a.fits"}]}]}))
    assert load_manifest(p)['targets'][0]['target'] == 'T'
    assert parse_seeds('1,2,3') == [1,2,3]
    with pytest.raises(Exception):
        parse_seeds('')
