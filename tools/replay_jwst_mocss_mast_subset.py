#!/usr/bin/env python3
"""Reproducible, deliberately bounded JWST GO-3290 / HUDF checks.

Inputs are untouched public MAST rate FITS, NIRCam I2D FITS and ECSV catalogs.
This is NOT a NIRSpec MOCSS 1-D spectral re-extraction: that requires the
custom MSA geometry, background treatment and wavelength/sky registration.

Run with: PYTHONPATH=/tmp/janus_astrolib python3 janus_jwst/replay_mast.py
The defaults point at the explicitly named files downloaded from MAST.
"""

import csv
import base64
import datetime as dt
import hashlib
import json
import math
import os
import random
import re
import statistics
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area, proj_plane_pixel_scales


RA = 53.13127
DEC = -27.78908
ROOT = Path('/tmp')
CATALOGS = {
    '2023-10-03': ROOT / 'janus_jwst_2023_f150w_cat.ecsv',
    '2025-10-06': ROOT / 'janus_jwst_2025_f150w_cat.ecsv',
    '2026-10-10': ROOT / 'janus_jwst_2026_f150w_cat.ecsv',
}
MOSAICS = {
    '2023-10-03': ('jw01180-o013_t011_nircam_clear-f150w_i2d.fits', ROOT / 'janus_2023_i2d_header.bin', ROOT / 'janus_2023_f150w_stripe.bin'),
    '2025-10-06': ('jw08060-o002_t002_nircam_clear-f150w_i2d.fits', ROOT / 'janus_2025_i2d_header.bin', ROOT / 'janus_2025_f150w_stripe.bin'),
    '2026-10-10': ('jw12577-o002_t002_nircam_clear-f150w_i2d.fits', ROOT / 'janus_2026_i2d_header.bin', ROOT / 'janus_2026_f150w_stripe.bin'),
}
RATES = {
    '2023-08-24': ROOT / 'jw03290001001_03109_00002_nrs2_rate.fits',
    '2023-11-19': ROOT / 'jw03290002001_11101_00002_nrs2_rate.fits',
}
OUTPUT = Path('janus_jwst/MAST_HUDF_INDEPENDENT_REPLAY_RECEIPT_2026-10-16.json')


def archive_discovery():
    rows = json.loads((ROOT / 'janus_jwst_3290_observations.json').read_text())['data']
    nearby = json.loads((ROOT / 'janus_jwst_cone.json').read_text())['data']
    nircam = [x for x in nearby if x.get('obs_collection') == 'JWST'
              and x.get('instrument_name') == 'NIRCAM/IMAGE'
              and x.get('filters') == 'F150W' and x.get('t_min')]
    obs = []
    for row in sorted(nircam, key=lambda x: x['t_min']):
        epoch = (dt.datetime(1858, 11, 17) +
                 dt.timedelta(days=row['t_min'])).strftime('%Y-%m-%d')
        if epoch in ['2023-10-03', '2025-10-06', '2026-10-10']:
            obs.append({'date': epoch, 'obsid': row['obsid'],
                        'proposal_id': row['proposal_id'], 'data_rights': row['dataRights'],
                        'archive_observation_id': row['obs_id']})
    return {
        'query_program_3290': {'service': 'Mast.Caom.Filtered',
                               'filter': 'proposal_id=3290',
                               'total_observation_rows': len(rows),
                               'long_science_detector_rows': sum(
                                   r.get('t_exptime', 0) > 8000 for r in rows)},
        'query_sky_cone': {'service': 'Mast.Caom.Cone',
                           'ra_deg': RA, 'dec_deg': DEC, 'radius_deg': 0.002,
                           'all_missions_observation_rows': len(nearby),
                           'jwst_rows': sum(r.get('obs_collection') == 'JWST' for r in nearby)},
        'retrieved_f150w_mosaic_parent_observations': obs,
        'caveat': 'MAST CAOM footprint intersection and MSA survey target IDs alone do not establish that a spectral slit exposes this precise source; FITS science WCS determines imaging coverage.',
    }


def pulsar_cone_controls():
    def inspect(path):
        root = ET.parse(path).getroot()
        ns = {'v': 'http://www.ivoa.net/xml/VOTable/v1.3'}
        status = [x.get('value') for x in root.findall('.//v:INFO', ns)
                  if x.get('name') == 'QUERY_STATUS']
        streams = root.findall('.//v:STREAM', ns)
        payload_bytes = [len(base64.b64decode(''.join((s.text or '').split())))
                         for s in streams]
        return {'query_status': status, 'decoded_binary_stream_byte_counts': payload_bytes,
                'archive_response_sha256': digest(path)['sha256']}
    return {
        'catalog': 'ATNF pulsar catalogue mirror, data.xao.ac.cn/pulsarcatalog/q/cone/scs.xml',
        'cone_0p1_degree': inspect(ROOT / 'janus_atnf_cone.xml'),
        'cone_1p0_degree': inspect(ROOT / 'janus_atnf_cone_1deg.xml'),
        'positive_control_crab_0p1_degree': inspect(ROOT / 'janus_atnf_control_crab.xml'),
        'interpretation': 'Valid zero-row ATNF cones after a non-empty Crab positive control. Only catalogued pulsars are excluded as a positional match; unknown pulsars are not ruled out.',
        'time_resolution_caveat': 'Each chosen JWST long rate frame coadds six integrations over 8315.667 s; an ordinary subsecond or second pulsar period cannot be tested in a rate image.',
    }


def predicted_companions():
    obs = 16004.0
    o3_5008 = 5008.239
    o2_3730 = 3729.875
    return {
        'if_primary_is_OIII_5008': {
            'predicted_OIII_4960_angstrom': round(obs * 4960.295 / o3_5008, 2),
            'predicted_Hbeta_angstrom': round(obs * 4862.683 / o3_5008, 2),
            'OIII_4960_to_5008_expected_flux_ratio_approximately': 1 / 3,
            'note': 'The prior COSMOS draft gave ~15864 A for 4960; that approximation is wrong. The consistent prediction is ~15850.79 A.'},
        'if_primary_is_OII_3730': {
            'predicted_OII_3727_angstrom': round(obs * 3727.092 / o2_3730, 2),
            'doublet_separation_angstrom': round(obs * (o2_3730 - 3727.092) / o2_3730, 2)},
        'science_gate': 'PREDICTIONS_ONLY__NO_COMPANION_LINE_LIMIT_WAS_REEXTRACTED',
    }


def digest(path):
    sha = hashlib.sha256()
    with open(path, 'rb') as fp:
        for block in iter(lambda: fp.read(4 * 1024 * 1024), b''):
            sha.update(block)
    return {'filename': path.name, 'size_bytes': path.stat().st_size,
            'sha256': sha.hexdigest()}


def coordinates_offset_arcsec(ra, dec):
    return math.hypot((ra - RA) * 3600 * math.cos(math.radians(DEC)),
                      (dec - DEC) * 3600)


def catalog_check(path):
    with path.open() as fp:
        lines = (line for line in fp if line.strip() and not line.startswith('#'))
        fields = next(lines).split()
        rows = []
        for line in lines:
            values = line.split()
            if len(values) != len(fields):
                continue
            row = dict(zip(fields, values))
            try:
                row['sep_arcsec'] = coordinates_offset_arcsec(
                    float(row['sky_centroid.ra']), float(row['sky_centroid.dec']))
            except (ValueError, KeyError):
                continue
            rows.append(row)
    rows.sort(key=lambda row: row['sep_arcsec'])
    return {
        'archive_file': digest(path),
        'source_count': len(rows),
        'matches_lt_1arcsec': sum(r['sep_arcsec'] < 1 for r in rows),
        'matches_lt_2arcsec': sum(r['sep_arcsec'] < 2 for r in rows),
        'nearest': [
            {'label': r['label'], 'separation_arcsec': round(r['sep_arcsec'], 4),
             'catalog_aper_total_abmag': float(r['aper_total_abmag'])}
            for r in rows[:3]
        ],
        'catalog_nonmatch_is_not_forced_photometry': True,
    }, rows


def rate_check(path):
    with fits.open(path, memmap=True, do_not_scale_image_data=True) as hdul:
        prim = hdul[0].header
        sci = hdul['SCI']
        sample = np.asarray(sci.data[::16, ::16], dtype=float)
        good = np.isfinite(sample)
        finite = sample[good]
        dq = hdul['DQ'].data if 'DQ' in hdul else None
        result = {
            'archive_file': digest(path),
            'primary': {k: prim.get(k) for k in [
                'TELESCOP', 'INSTRUME', 'DETECTOR', 'FILTER', 'GRATING',
                'EXP_TYPE', 'DATE-OBS', 'TIME-OBS', 'EXPTIME', 'EFFEXPTM',
                'NINTS', 'NGROUPS', 'READPATT', 'SUBARRAY', 'CRDS_CTX']},
            'science_shape': list(sci.data.shape),
            'science_unit': sci.header.get('BUNIT'),
            'sample_stride_pixels': 16,
            'sample_count': int(sample.size),
            'sample_finite_fraction': round(float(good.mean()), 6),
            'sample_finite_median': float(np.median(finite)),
            'sample_finite_p01_p99': [float(v) for v in np.percentile(finite, [1, 99])],
        }
        if dq is not None:
            # FITS unsigned integer convention stores detector DQ as signed
            # INT32 with BZERO=2**31. Inspect physical unsigned values.
            dsample = np.asarray(dq[::16, ::16], dtype=np.int64) + int(
                hdul['DQ'].header.get('BZERO', 0))
            result['dq_sample_nonzero_fraction'] = round(
                float(np.count_nonzero(dsample) / dsample.size), 6)
        return result


def aperture_photometry(image, weight, wcs, sky_ra, sky_dec, y_origin=0):
    cx, cy = wcs.all_world2pix([[sky_ra, sky_dec]], 0)[0]
    cy -= y_origin
    if not (math.isfinite(cx) and math.isfinite(cy)):
        return None
    pix_scale = float(np.mean(proj_plane_pixel_scales(wcs.celestial)) * 3600)
    r_ap = 0.12 / pix_scale
    r_in = 0.24 / pix_scale
    r_out = 0.38 / pix_scale
    if cx < r_out + 2 or cy < r_out + 2 or cx > image.shape[1] - r_out - 2 or cy > image.shape[0] - r_out - 2:
        return None
    x0, x1 = math.floor(cx - r_out - 1), math.ceil(cx + r_out + 2)
    y0, y1 = math.floor(cy - r_out - 1), math.ceil(cy + r_out + 2)
    ys, xs = np.mgrid[y0:y1, x0:x1]
    rsq = (xs - cx) ** 2 + (ys - cy) ** 2
    good = np.isfinite(image[y0:y1, x0:x1])
    if weight is not None:
        good &= (weight[y0:y1, x0:x1] > 0)
    ap = good & (rsq <= r_ap ** 2)
    sky = good & (rsq >= r_in ** 2) & (rsq <= r_out ** 2)
    if ap.sum() < 0.8 * math.pi * r_ap ** 2 or sky.sum() < 35:
        return None
    background = float(np.median(image[y0:y1, x0:x1][sky]))
    pixels_mjysr = float((image[y0:y1, x0:x1][ap] - background).sum())
    pixel_area_sr = float(proj_plane_pixel_area(wcs.celestial)) * (math.pi / 180) ** 2
    return {'aperture_njy': pixels_mjysr * pixel_area_sr * 1e15,
            'x': float(cx), 'y': float(cy), 'aperture_pixels': int(ap.sum()),
            'local_background_mjy_sr': background, 'pixel_scale_arcsec': pix_scale}


def mosaic_check(product, header_file, stripe_file, catalog_rows):
    with header_file.open('rb') as fp:
        fits.Header.fromfile(fp, endcard=True, padding=True)
        header = fits.Header.fromfile(fp, endcard=True, padding=True)
        data_start = fp.tell()
    assert header.get('EXTNAME') == 'SCI'
    assert header.get('BUNIT', '').replace(' ', '').lower() in ('mjy/sr', 'mjy/steradian')
    wcs = WCS(header)
    _, target_y = wcs.all_world2pix([[RA, DEC]], 0)[0]
    y0 = int(target_y) - 450
    nrows = 901
    width = header['NAXIS1']
    start, end = data_start + y0 * width * 4, data_start + (y0 + nrows) * width * 4 - 1
    assert stripe_file.stat().st_size == end - start + 1, 'MAST byte-range response was not complete'
    image = np.memmap(stripe_file, dtype='>f4', mode='r', shape=(nrows, width))
    weight = None
    target = aperture_photometry(image, weight, wcs, RA, DEC, y_origin=y0)
    if target is None:
        return {'archive_file': digest(stripe_file), 'target_on_valid_image': False}
    random_gen = random.Random(3290)
    empty_control_fluxes = []
    for _ in range(1100):
        dx, dy = random_gen.uniform(-12, 12), random_gen.uniform(-12, 12)
        if math.hypot(dx, dy) < 1.5:
            continue
        ra = RA + dx / (3600 * math.cos(math.radians(DEC)))
        dec = DEC + dy / 3600
        too_close_to_catalog_source = any(
            math.hypot((ra - float(r['sky_centroid.ra'])) * 3600 * math.cos(math.radians(DEC)),
                       (dec - float(r['sky_centroid.dec'])) * 3600) < 0.8
            for r in catalog_rows if r['sep_arcsec'] < 19)
        if too_close_to_catalog_source:
            continue
        c = aperture_photometry(image, weight, wcs, ra, dec, y_origin=y0)
        if c is not None:
            empty_control_fluxes.append(c['aperture_njy'])
        if len(empty_control_fluxes) >= 250:
            break
    med = statistics.median(empty_control_fluxes)
    sigma = 1.4826 * statistics.median(abs(x - med) for x in empty_control_fluxes)
    return {
        'archive_product': product,
        'mast_data_uri': 'mast:JWST/product/' + product,
        'retrieved_sci_range_bytes_inclusive': [start, end],
        'science_stripe_digest': digest(stripe_file),
        'fits_primary_plus_sci_header_digest': digest(header_file),
        'science_image_naxis_xy': [width, header['NAXIS2']],
        'science_stripe_row_start_zero_indexed': y0,
        'science_stripe_nrows': nrows,
        'target_on_valid_image': True,
        'pipeline_image_unit': header.get('BUNIT'),
        'aperture_radius_arcsec': 0.12,
        'annulus_inner_outer_arcsec': [0.24, 0.38],
        'target_wcs_pixel_xy_zero_indexed': [round(target['x'], 3), round(target['y'] + y0, 3)],
        'target_aperture_pixels': target['aperture_pixels'],
        'pixel_scale_arcsec': target['pixel_scale_arcsec'],
        'target_background_subtracted_aperture_njy': round(target['aperture_njy'], 5),
        'blank_aperture_control_count': len(empty_control_fluxes),
        'blank_aperture_control_median_njy': round(med, 5),
        'blank_aperture_control_robust_sigma_njy': round(sigma, 5),
        'target_aperture_snr_empirical': round((target['aperture_njy'] - med) / sigma, 4),
        'local_aperture_3sigma_conservative_bound_njy': round(max(0., target['aperture_njy'] - med) + 3 * sigma, 5),
        'interpretation': 'Diagnostic 0.12-arcsec aperture only, 901-row SCI cutout from exact public FITS byte range; no PSF aperture correction, injected-source completeness, or correlated-noise calibration. Do not treat this as a certified total-flux limit.'
    }


def main():
    cats, rows = {}, {}
    for epoch, path in CATALOGS.items():
        cats[epoch], rows[epoch] = catalog_check(path)
    result = {
        'schema': 'janus.independent_mast_subset_receipt.v1',
        'target': {'ra_deg': RA, 'dec_deg': DEC, 'observation_line_angstrom': 16004},
        'archive': 'https://mast.stsci.edu/api/v0/invoke',
        'program': 'JWST GO-3290',
        'analysis_scope': 'Archive enumeration; two independent detector-rate FITS quality summaries; three NIRCam catalog crossmatches; image-space forced aperture checks for available public mosaics.',
        'not_done': ['MOCSS sky-to-detector solution and wavelength extraction',
                     'Independently measured 1.6004 micron line flux and companion-line limits',
                     'Exposure time-domain periodogram capable of ruling out pulsar',
                     'Contemporaneous (same-day) F150W imaging and G140M spectroscopy'],
        'archive_discovery': archive_discovery(),
        'atnf_pulsar_crossmatch': pulsar_cone_controls(),
        'line_companion_predictions': predicted_companions(),
        'catalog_crossmatches': cats,
        'detector_rate_fits': {epoch: rate_check(path) for epoch, path in RATES.items()},
        'mosaic_forced_aperture': {epoch: mosaic_check(*files, rows[epoch]) for epoch, files in MOSAICS.items() if files[2].is_file()},
    }
    OUTPUT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
