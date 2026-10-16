#!/usr/bin/env python3
"""Compile a compact, offline-auditable 32-FITS MOCSS replay evidence receipt."""

import datetime as dt
import hashlib
import json
import math
from pathlib import Path


ROOT = Path('janus_jwst')


def read(name):
    path = ROOT / name
    if not path.exists():
        path = Path('domains/space') / name
    return json.loads(path.read_text())


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def json_safe(obj):
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [json_safe(v) for v in obj]
    return obj


def main():
    inputs = read('MOCSS_16_EXPOSURE_MANIFEST_2026-10-16.json')['exposures']
    wcs = {det: read('mocss_16_exposure_' + det + '_wcs_trace_2026-10-16.json')['records']
           for det in ('nrs1', 'nrs2')}
    contrast = read('mocss_rate_16_exposure_line_contrasts_2026-10-16.json')['records']
    highpass = read('mocss_wcs_registered_rate_highpass_2026-10-16.json')['groups']
    alternative_31 = read('mocss_wcs_registered_rate_highpass_31px_2026-10-16.json')['groups']
    alternative_61 = read('mocss_wcs_registered_rate_highpass_61px_2026-10-16.json')['groups']
    rawconstant = read('mocss_wcs_registered_rate_stack_2026-10-16.json')['groups']
    spatial = read('mocss_local_sky_coherence_receipt_2026-10-16.json')
    record = {
        'schema': 'janus-mocss-independent-wcs-sky-rate-probe-v1',
        'created_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'finding_classification': 'GEOMETRY_VALIDATED_RATE_LINE_CONTRAST_UNCALIBRATED',
        'archive': 'MAST JWST GO-3290',
        'crds': 'jwst_1413.pmap', 'jwst_pipeline': '1.19.1',
        'target_ra_deg': 53.13127, 'target_dec_deg': -27.78908,
        'observed_lambda_angstrom': 16004.0,
        'comparison_oiii_angstrom': 15850.7933,
        'comparison_pagamma_angstrom': 13656.7483,
        'two_full_msa_reference_files': {},
        'exposures': [], 'groups': [], 'sky_consistency': spatial,
        'unresolved_work': [
            'Independent calibration of 32 original RATE arrays to full spectrophotometric products',
            'Physical shutter throughput, flat/1-f background, variance/DQ and slitless deblending',
            'Upper limits for the 15850.8-A and Pa-gamma comparison line with injection recovery',
            'Flux-calibrated extended emission map with independent control positions',
            'No source identity, beam direction or physical line length established']}
    for filename in ('jw03290001001_05_msa.fits', 'jw03290002001_01_msa.fits'):
        path = Path('/tmp') / filename
        record['two_full_msa_reference_files'][filename] = {
            'bytes': path.stat().st_size, 'sha256': sha256(path)}
    for exposure in inputs:
        ix = exposure['exposure_index']
        a = {'exposure_index': ix, 'group': exposure['group'],
             'orientation': exposure['orientation'],
             'source_msametid_unwrapped': exposure['msa_metadata_id'],
             'source_msametid_signed16': ((exposure['msa_metadata_id'] + 32768) % 65536) - 32768,
             'detectors': {}}
        for detector in ('nrs1', 'nrs2'):
            original = Path('/tmp') / (exposure['exposure_prefix'] + '_' + detector + '_rate.fits')
            frame = wcs[detector][ix - 1]
            a['detectors'][detector] = {
                'rate_file': original.name, 'bytes': original.stat().st_size,
                'sha256': sha256(original),
                'wcs_error': frame.get('wcs_error'),
                'selected_virtual_shutter': frame.get('best_shutter'),
                'sky_registered_trace': frame.get('trace')}
        selection = 'nrs2' if ix <= 8 else 'nrs1'
        a['expected_science_detector_from_WCS'] = selection
        a['raw_rate_line_contrasts'] = contrast[ix - 1]['detectors'][selection]
        record['exposures'].append(a)
        print('HASHED', ix, exposure['group'], flush=True)
    for idx, group in enumerate(highpass):
        orig = rawconstant[idx]
        assert group['group'] == orig['group']
        record['groups'].append({
            'group': group['group'], 'detector': group['detector'],
            'num_dither_frames': len(group['frames']),
            'constant_row_background_probe': orig['lines'],
            'experimental_41pixel_median_highpass_probe': group['lines'],
            'robustness_31pixel_median_highpass': alternative_31[idx]['lines'],
            'robustness_61pixel_median_highpass': alternative_61[idx]['lines']})
    assert len(record['exposures']) == 16
    assert len(record['groups']) == 4
    result = ROOT / 'MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json'
    result.write_text(json.dumps(json_safe(record), indent=2, allow_nan=False) + '\n')
    print('WROTE', result, result.stat().st_size,
          '32 original RATE SHA256, 2 MSA SHA256; UNCALIBRATED')


if __name__ == '__main__':
    main()
