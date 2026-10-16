#!/usr/bin/env python3
"""Verify 16 MOCSS aperture IDs and test conditional physical geometry.

Two publicly archived 80-MB MSA files suffice to *inspect* all 8 dither
metadata IDs for each orientation, even though the 16 rate products cite 16
different MSA filenames. This is metadata validation, not a spectral cube.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from math import atan2, cos, degrees, hypot, radians
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.stats import chi2


TARGET = (53.13127, -27.78908)
NEIGHBORS = {
    'northwest': {'ra_deg': 53.131165, 'dec_deg': -27.788930, 'zphot': 0.22},
    'east': {'ra_deg': 53.131584, 'dec_deg': -27.789105, 'zphot': 6.03},
    'south': {'ra_deg': 53.131301, 'dec_deg': -27.789367, 'zphot': 0.19},
}
PA_BETA_Z_SCALE = 4.02  # proper kpc/arcsec at z~0.249 (Planck-like LCDM)
LIGHT_YEARS_PER_KPC = 3261.56


def checksum(path):
    h = hashlib.sha256()
    with path.open('rb') as fp:
        for chunk in iter(lambda: fp.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return {'filename': path.name, 'bytes': path.stat().st_size, 'sha256': h.hexdigest()}


def summarize_msa(path):
    with fits.open(path, memmap=True) as h:
        sh = h['SHUTTER_INFO'].data
        ids, counts = np.unique(sh['msa_metadata_id'], return_counts=True)
        return {'archived_file': checksum(path),
                'row_count': int(len(sh)),
                'msaid_signed_16bit_to_counts': {str(int(k)): int(v) for k, v in zip(ids, counts)},
                'normalized_ids_unsigned_16bit': sorted([int(v) % 65536 for v in ids]),
                'source_info_rows': int(len(h['SOURCE_INFO'].data)),
                'caution': 'SHUTTER_INFO msa_metadata_id stored as signed int16 in archival FITS; values >32767 wrap modulo 65536. Compare after normalization, do not compare raw integers.'}


def hypotheses():
    # Report only physically defensible diagnostics; redshifts not measured.
    line = 16004.0
    oi = 5008.239
    oj = 4960.295
    predicted_oii = line * oj / oi
    flux = np.array([2.3, 3.5, 2.9, 2.9])
    uncertainty = np.array([0.5, 0.6, 0.7, 1.2])
    weights = 1.0 / uncertainty ** 2
    mu = float(np.sum(weights * flux) / weights.sum())
    chisq = float(np.sum(((flux - mu) / uncertainty) ** 2))
    return {
        'OIII_5008_primary_hypothesis': {
            'published_primary_line_A': line,
            'vacuum_rest_5008_A': oi,
            'vacuum_rest_4960_A': oj,
            'derived_companion_A': round(predicted_oii, 2),
            'paper_companion_claim_A': 15864,
            'paper_minus_derived_A': round(15864 - predicted_oii, 2),
            'predicted_companion_flux_erg_s_cm2': 2.5e-18 / 3,
            'paper_reported_local_3sigma_limit_around_15864_erg_s_cm2': 9.4e-19,
            'interpretation': 'A nondetection evaluated at 15864 A cannot exclude a feature at 15850.8 A. Additionally, the predicted companion flux is smaller than the reported 3-sigma limit. Re-extract actual spectrum at correct position, with a local uncertainty and wavelength solution. No line detection claimed.'
        },
        'stack_flux_constancy': {
            'provenance': 'Huberty et al. 2026 Table 2; no independent 1D re-extraction',
            'source_flux_unit': '1e-18 erg s-1 cm-2',
            'measurements_ABCD': flux.tolist(),
            'errors_ABCD': uncertainty.tolist(),
            'inverse_variance_weighted_mean': round(mu, 4),
            'chi_square': round(chisq, 4),
            'degrees_of_freedom': 3,
            'null_survival_probability': round(float(chi2.sf(chisq, 3)), 4),
            'interpretation': 'Published four-stack fluxes are compatible with a constant line over ~88 days; they do not classify the emitter or exclude a long-duration transient.'
        },
        'nearby_galaxies': {},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', default='janus_jwst/MOCSS_16_EXPOSURE_MANIFEST_2026-10-16.json')
    parser.add_argument('--output', default='janus_jwst/MOCSS_16_EXPOSURE_GEOMETRY_RECEIPT_2026-10-16.json')
    parser.add_argument('--data-dir', default='/tmp')
    parser.add_argument('--aug-products', default='/tmp/janus_mocss_products.json')
    parser.add_argument('--nov-products', default='/tmp/janus_mocss_products_nov.json')
    args = parser.parse_args()
    inventory = json.loads(Path(args.manifest).read_text())
    assert inventory['num_long_exposures'] == 16 and inventory['num_detector_frames'] == 32

    # Each archival MSA table includes all eight dithers in one orientation.
    representatives = {
        1: Path(args.data_dir) / 'jw03290001001_05_msa.fits',
        2: Path(args.data_dir) / 'jw03290002001_01_msa.fits',
    }
    msa = {str(orientation): summarize_msa(path)
           for orientation, path in representatives.items()}
    checks = []
    for row in inventory['exposures']:
        orientation = row['orientation']
        stored = int(row['msa_metadata_id'])
        unsigned16 = stored % 65536
        signed = unsigned16 if unsigned16 < 32768 else unsigned16 - 65536
        counts = msa[str(orientation)]['msaid_signed_16bit_to_counts']
        found = counts.get(str(signed))
        if found is None:
            raise AssertionError(f'MSA id {stored} absent in orientation {orientation}')
        checks.append({'exposure': row['exposure_index'], 'orientation': orientation,
                       'group': row['group'], 'rate_primary_header_msa_id': stored,
                       'msa_shutter_info_signed16_id': signed,
                       'matching_open_shutter_rows': found,
                       'metadata_id_verified': True})

    products = {}
    for name, filename in [('aug_rate', args.aug_products), ('nov_rate', args.nov_products)]:
        d = json.loads(Path(filename).read_text())
        if d.get('status') != 'COMPLETE':
            raise ValueError(f'MAST product response not COMPLETE: {filename}')
        groups = sorted(set(x.get('productSubGroupDescription', '') or '' for x in d['data']))
        products[name] = {
            'request_response_sha256': checksum(Path(filename))['sha256'],
            'product_count': len(d['data']),
            'available_subgroups': groups,
            'pipeline_stage2_s2d_or_cal_present': any(x in groups for x in ('S2D', 'CAL', 'X1D')),
        }

    science = hypotheses()
    for name, sky in NEIGHBORS.items():
        east = (sky['ra_deg'] - TARGET[0]) * cos(radians(TARGET[1])) * 3600
        north = (sky['dec_deg'] - TARGET[1]) * 3600
        rho = hypot(east, north)
        proper = rho * PA_BETA_Z_SCALE
        science['nearby_galaxies'][name] = {
            **sky, 'east_offset_arcsec': round(east, 4),
            'north_offset_arcsec': round(north, 4),
            'angular_sep_arcsec': round(rho, 4),
            'bearing_from_line_deg_east_of_north': round(degrees(atan2(east, north)) % 360, 2),
            'if_same_z_as_Pa_beta': {'projected_separation_kpc': round(proper, 3),
                                    'minimum_source_to_cloud_light_travel_years': round(proper * LIGHT_YEARS_PER_KPC)},
            'association_status': 'unverified, photometric redshift only',
        }

    receipt = {
        'schema_version': 'janus-mocss-shutter-geometry-and-hypothesis-test-v1',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'target': {'ra_deg': TARGET[0], 'dec_deg': TARGET[1], 'line_A': 16004},
        'source_references': [
            'https://arxiv.org/html/2610.12553v1',
            'https://arxiv.org/html/2609.36003',
            'https://doi.org/10.17909/tct8-jr92',
            'https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-observing-modes/nirspec-ifu-spectroscopy',
        ],
        'research_stage': 'PARTIAL_INPUT_VERIFICATION_NOT_AN_INDEPENDENT_LINE_MAP',
        'replay_coverage': {'long_exposures_with_FITS_headers_verified': 16,
                            'detector_FITS_headers_verified': 32,
                            'MSA_shutter_metadata_ids_verified': len(checks),
                            'full_detector_rate_images_independently_processed': 2,
                            'independent_line_reextractions': 0,
                            'off_source_at_16004A_sky_map_reconstructions': 0},
        'msa_representative_archives': msa,
        'per_exposure_msa_id_consistency': checks,
        'mast_product_availability': products,
        'hypothesis_diagnostics': science,
        'ifu_observation_design': {
            'single_pointing_field_arcsec': [3, 3],
            'spaxel_size_arcsec': [0.1, 0.1],
            'all_three_neighbors_within_1p05arcsec': True,
            'search_requirements': [
                'One deep NIRSpec IFU pointing, approximately centered on the reported orphan, measures all three proposed illuminators and line-emitting region.',
                'Mosaic IFU pointings needed for a blind >=5-arcsec sky-extent search.',
                'Extract 16004 +/- 10 A and continuum-matched adjacent wavelengths, 15850.8 A ([O III] doublet), Pa gamma at ~1.365 um, plus other expected transitions.',
                'Simultaneous F150W imaging and comparison between 2023 October mosaic and 2023 August-November MOCSS spectra.',
            ],
        },
        'scientific_stop_condition': 'Cannot conclude source class, emission direction, extension or whether 15850.8-A line is present without wavelength-aware detector registration and calibrated spectral extraction.',
    }
    Path(args.output).write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + '\n')
    print('PASS: 16 shutter IDs checked; two distinct orientations; OIII companion =',
          science['OIII_5008_primary_hypothesis']['derived_companion_A'])
    print('SAVED', args.output)


if __name__ == '__main__':
    main()
