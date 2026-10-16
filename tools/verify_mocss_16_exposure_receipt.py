#!/usr/bin/env python3
"""Offline strict science gate for the committed MOCSS input receipt.

Validates provenance consistency and prevents an archive-header inspection
from being misrepresented as an independently extracted spectral line map.
"""

import argparse
import json
from pathlib import Path


def run(manifest_path, receipt_path):
    manifest = json.loads(Path(manifest_path).read_text())
    receipt = json.loads(Path(receipt_path).read_text())
    assert manifest['num_long_exposures'] == 16
    assert manifest['num_detector_frames'] == 32
    assert len(manifest['exposures']) == 16
    assert len(receipt['per_exposure_msa_id_consistency']) == 16
    assert [e['group'] for e in manifest['exposures']] == list('AAAABBBBCCCCDDDD')
    assert all(set(e['detectors']) == {'nrs1', 'nrs2'} for e in manifest['exposures'])
    assert all(e['metadata_id_verified'] for e in receipt['per_exposure_msa_id_consistency'])
    for raw, matching in zip(manifest['exposures'], receipt['per_exposure_msa_id_consistency']):
        assert raw['exposure_index'] == matching['exposure']
        assert raw['msa_metadata_id'] == matching['rate_primary_header_msa_id']
        unsigned16 = raw['msa_metadata_id'] % 65536
        expected_signed = unsigned16 if unsigned16 < 32768 else unsigned16 - 65536
        assert matching['msa_shutter_info_signed16_id'] == expected_signed
        assert matching['matching_open_shutter_rows'] > 200_000
    coverage = receipt['replay_coverage']
    assert coverage['long_exposures_with_FITS_headers_verified'] == 16
    assert coverage['detector_FITS_headers_verified'] == 32
    assert coverage['MSA_shutter_metadata_ids_verified'] == 16
    assert coverage['independent_line_reextractions'] == 0
    assert coverage['off_source_at_16004A_sky_map_reconstructions'] == 0
    assert receipt['research_stage'].startswith('PARTIAL_INPUT_VERIFICATION')
    companion = receipt['hypothesis_diagnostics']['OIII_5008_primary_hypothesis']
    expected_companion = 16004.0 * 4960.295 / 5008.239
    assert abs(companion['derived_companion_A'] - expected_companion) < 0.015
    assert companion['paper_companion_claim_A'] == 15864
    assert companion['paper_minus_derived_A'] > 13.0
    assert companion['predicted_companion_flux_erg_s_cm2'] < companion['paper_reported_local_3sigma_limit_around_15864_erg_s_cm2']
    const = receipt['hypothesis_diagnostics']['stack_flux_constancy']
    assert const['degrees_of_freedom'] == 3
    assert 0.48 < const['null_survival_probability'] < 0.52
    assert all(not d['pipeline_stage2_s2d_or_cal_present']
               for d in receipt['mast_product_availability'].values())
    print('PASS: archive input 16/16, MSA IDs 16/16, no fake spectral-map promotion; OIII companion arithmetic audited')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', default='domains/space/MOCSS_16_EXPOSURE_MANIFEST_2026-10-16.json')
    parser.add_argument('--receipt', default='domains/space/MOCSS_16_EXPOSURE_GEOMETRY_RECEIPT_2026-10-16.json')
    args = parser.parse_args()
    run(args.manifest, args.receipt)
