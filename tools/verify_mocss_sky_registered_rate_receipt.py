#!/usr/bin/env python3
"""Dependency-free CI contract for the 32-FITS sky-registered RATE replay.

Validates the committed receipt's internal provenance and scientific ceilings;
does NOT pretend to rerun archival FITS or assert line-flux measurements in CI.
"""

import json
import math
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent


def main():
    p = ROOT / 'MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json'
    if not p.exists():
        p = ROOT.parent / 'domains' / 'space' / p.name
    x = json.loads(p.read_text())
    assert x['schema'] == 'janus-mocss-independent-wcs-sky-rate-probe-v1'
    assert x['finding_classification'] == 'GEOMETRY_VALIDATED_RATE_LINE_CONTRAST_UNCALIBRATED'
    assert x['crds'] == 'jwst_1413.pmap'
    assert x['jwst_pipeline'] == '1.19.1'
    assert len(x['exposures']) == 16
    assert len(x['groups']) == 4
    assert len(x['two_full_msa_reference_files']) == 2
    assert abs(x['comparison_oiii_angstrom'] -
               16004 * 4960.295 / 5008.239) < 0.005
    assert abs(x['comparison_pagamma_angstrom'] -
               16004 * 10941.1 / 12821.6) < 0.005
    sha_pat = re.compile(r'[0-9a-f]{64}')
    for source in x['two_full_msa_reference_files'].values():
        assert sha_pat.fullmatch(source['sha256'])
        assert source['bytes'] > 10_000_000
    for i, rec in enumerate(x['exposures'], start=1):
        assert rec['exposure_index'] == i
        assert rec['group'] == 'ABCD'[(i - 1) // 4]
        assert rec['source_msametid_signed16'] == (
            (rec['source_msametid_unwrapped'] + 32768) % 65536) - 32768
        assert len(rec['detectors']) == 2
        detector = 'nrs2' if i <= 8 else 'nrs1'
        assert rec['expected_science_detector_from_WCS'] == detector
        for det in ('nrs1', 'nrs2'):
            file = rec['detectors'][det]
            assert file['bytes'] == 83_966_400
            assert sha_pat.fullmatch(file['sha256'])
            assert file['rate_file'].endswith(det + '_rate.fits')
        trace = rec['detectors'][detector]['sky_registered_trace']['16004A']
        assert trace['on_detector']
        assert abs(trace['slit_x']) < 0.75
        assert abs(trace['slit_y']) < 0.75
        assert 0 <= trace['x_px'] < 2048 and 0 <= trace['y_px'] < 2048
        probes = rec['raw_rate_line_contrasts']
        assert probes['status'].startswith('RATE_CONTRAST_')
        assert 'heuristic_robust_z_unvalidated' in probes['line_probes']['16004A']
        assert math.isfinite(probes['line_probes']['16004A']['heuristic_robust_z_unvalidated'])
    for group in x['groups']:
        label = group['group']
        assert group['num_dither_frames'] == 4
        assert group['detector'] == ('nrs2' if label in 'AB' else 'nrs1')
        values = [group[f'{name}']['16004A']['provisional_center_z_not_gaussian_sigma']
                  for name in ('robustness_31pixel_median_highpass',
                               'experimental_41pixel_median_highpass_probe',
                               'robustness_61pixel_median_highpass')]
        assert all(math.isfinite(t) for t in values)
        assert all(v > 0 for v in values), (label, values)
        if label in 'ABD':
            assert min(values) > 2
        companion = group['experimental_41pixel_median_highpass_probe']['15850p8A']
        assert abs(companion['provisional_center_z_not_gaussian_sigma']) < 1
    sky = x['sky_consistency']
    assert sky['status'] == 'DESCRIPTIVE_UNVALIDATED_LOCAL_SKY_COMPARISON'
    assert sky['source_location_preselected']
    assert sky['source_field_arcsec'] == 3.0
    assert len(sky['neighbors']) == 3
    assert sky['wavelength_16004A']['combined_target_relative_robust_units'] > 0
    assert sky['companion_control_15850p8A']['combined_target_relative_robust_units'] < 0
    assert 'physical' in ' '.join(x['unresolved_work']).lower()
    print('PASS: 32 original FITS fingerprints, 16 sky-registered traces and raw-control science ceilings')
    print('INDEPENDENT_SKY_WCS_32_RATE_INPUTS=VERIFIED_RECEIPT')
    print('SPECTROPHOTOMETRIC_DETECTION_SIGNIFICANCE=NOT_VALIDATED')
    print('OIII_COMPANION_FLUX_LIMIT=NOT_CALIBRATED')
    print('PHYSICAL_SOURCE_AND_TRAIL=UNKNOWN')


if __name__ == '__main__':
    main()
