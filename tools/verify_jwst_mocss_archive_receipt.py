#!/usr/bin/env python3
"""Offline scientific authority and internal-consistency gate for MAST receipts.

No network, FITS downloads, external packages, or model updates in CI. This
checks recorded independent-replay calculations, not their astrophysical truth.
"""

import datetime
import json
import math
import re
from pathlib import Path


BASE = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent


def find(name_in_checkout, name_locally):
    path = BASE / name_in_checkout
    if not path.exists():
        path = HERE / name_locally
    return json.loads(path.read_text())


def verify():
    result = find('domains/space/JWST_MOCSS_HUDF_MAST_INDEPENDENT_RECEIPT_2026-10-16.json',
                  'MAST_HUDF_INDEPENDENT_REPLAY_RECEIPT_2026-10-16.json')
    spatial = find('domains/space/JWST_MOCSS_HUDF_ORIGIN_SPATIAL_TEST_2026-10-16.json',
                   'SPATIAL_ORIGIN_AND_DISTANCE_TEST_2026-10-16.json')
    assert result['schema'] == 'janus.independent_mast_subset_receipt.v1'
    assert result['target'] == {'ra_deg': 53.13127, 'dec_deg': -27.78908,
                                'observation_line_angstrom': 16004}
    archive = result['archive_discovery']['query_program_3290']
    assert archive['total_observation_rows'] >= 80
    assert archive['long_science_detector_rows'] == 32
    rates = result['detector_rate_fits']
    assert set(rates) == {'2023-08-24', '2023-11-19'}
    for name, rate in rates.items():
        h = rate['primary']
        assert h['DETECTOR'] == 'NRS2' and h['GRATING'] == 'G140M'
        assert h['NINTS'] == 6 and h['NGROUPS'] == 19
        assert h['EFFEXPTM'] > 8000 and rate['science_shape'] == [2048, 2048]
        assert 0.5 < rate['sample_finite_fraction'] <= 1
        assert 0 < rate['dq_sample_nonzero_fraction'] < 1
        d = rate['archive_file']
        assert re.fullmatch('[0-9a-f]{64}', d['sha256']) and d['size_bytes'] > 8e7
    chrono = sorted(result['mosaic_forced_aperture'])
    assert chrono == ['2023-10-03', '2025-10-06', '2026-10-10']
    assert datetime.date(2023, 8, 24) < datetime.date(2023, 10, 3) < datetime.date(2023, 11, 19)
    for date in chrono:
        x = result['mosaic_forced_aperture'][date]
        assert x['target_on_valid_image'] is True
        assert x['pipeline_image_unit'] == 'MJy/sr'
        assert x['blank_aperture_control_count'] >= 200
        assert abs(x['target_aperture_snr_empirical']) < 3
        assert 0 < x['blank_aperture_control_robust_sigma_njy'] < 20
        assert 0 < x['pixel_scale_arcsec'] < 0.1
        lo, hi = x['retrieved_sci_range_bytes_inclusive']
        assert hi - lo + 1 == x['science_stripe_digest']['size_bytes']
        assert re.fullmatch('[0-9a-f]{64}', x['science_stripe_digest']['sha256'])
        assert 'no PSF aperture correction' in x['interpretation']
        c = result['catalog_crossmatches'][date]
        assert c['matches_lt_2arcsec'] == 0 and c['source_count'] > 1000
    pc = result['atnf_pulsar_crossmatch']
    assert pc['cone_1p0_degree']['decoded_binary_stream_byte_counts'] == [0]
    assert pc['positive_control_crab_0p1_degree']['decoded_binary_stream_byte_counts'][0] > 0
    expected = 16004 * 4960.295 / 5008.239
    actual = result['line_companion_predictions']['if_primary_is_OIII_5008'][
        'predicted_OIII_4960_angstrom']
    assert math.isclose(actual, expected, rel_tol=0, abs_tol=0.02)
    assert result['line_companion_predictions']['science_gate'].startswith('PREDICTIONS_ONLY')
    assert 'MOCSS sky-to-detector solution and wavelength extraction' in result['not_done']
    screen = spatial['published_MOCSS_same_wavelength_catalog_screen']
    assert screen['other_reported_emitter_sources_screened'] == 23
    assert screen['representative_line_identifications_screened'] == 41
    assert screen['representative_lines_within_plusminus_100_angstrom'] == 0
    assert spatial['epistemic_gate']['PHYSICAL_DISTANCE'] == 'CONDITIONAL_ON_LINE_ID'
    assert spatial['epistemic_gate']['SPATIAL_PROPAGATION_LENGTH'] == 'UNIDENTIFIABLE'
    assert spatial['epistemic_gate']['SOURCE_IDENTITY'] == 'UNRESOLVED'
    assert spatial['world_truth'] is False
    for case in spatial['conditional_distance_scenarios']['cases']:
        assert math.isclose(case['z'], 16004 / case['rest_wavelength_angstrom' if 'rest_wavelength_angstrom' in case else 'rest_wavelength_angstrom_approx'] - 1, abs_tol=0.0001)
        assert case['luminosity_distance_gpc'] > 0
    print('PASS: archive metadata, FITS rate QA, 3 epoch apertures, [OIII] corrected math, ATNF positive control, 23-source sky catalog screen, unresolved spatial/distance gates')


if __name__ == '__main__':
    verify()
