#!/usr/bin/env python3
"""Sky-coordinate quick-look grid of WCS-registered MOCSS RATE line cutouts.

Uses independent NIRSpec WCS at four group anchors and high-pass detector
coadds. This is not a flux-calibrated NIRSpec IFU map, not physical extension
or a statistically validated detection. Nearby galaxies are *test positions*.
"""

import argparse
import json
import logging
import math
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates
from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.nirspec import nrs_wcs_set_input


RA, DEC = 53.13127, -27.78908
NEIGHBORS = {
    'NW_low_z_photo': (53.131165, -27.788930),
    'S_low_z_photo': (53.131301, -27.789367),
    'E_high_z_photo': (53.131584, -27.789105),
}
HALF_X, HALF_Y = 76, 26
SPACING, HALF_FIELD = 0.1, 1.5
GRID = np.arange(-HALF_FIELD, HALF_FIELD + 0.5 * SPACING, SPACING)
EN_EAST, EN_NORTH = np.meshgrid(GRID, GRID, indexing='xy')


def group_wcs_info(group, records):
    index = {'A': 1, 'B': 8, 'C': 10, 'D': 13}[group]
    detector = 'nrs2' if group in 'AB' else 'nrs1'
    entry = records[(index, detector)]
    source = ('/tmp/janus_mocss_wcs_all16_nrs2/all16_input_receipt.json'
              if detector == 'nrs2' else '/tmp/janus_mocss_wcs_all16/all16_input_receipt.json')
    inputs = json.loads(Path(source).read_text())['records']
    probe = next(x for x in inputs if x['exposure_index'] == index)
    m = AssignWcsStep.call(probe['probe_rate'], save_results=False)
    w = nrs_wcs_set_input(m, entry['best_shutter']['slit_id'])
    lam = 1.6004
    def get_xy(east_arcsec, north_arcsec):
        x, y = w.invert(RA + east_arcsec / (3600 * math.cos(math.radians(DEC))),
                        DEC + north_arcsec / 3600, lam)
        return np.array([float(x), float(y)])
    h = 0.20
    reference_xy = get_xy(0, 0)
    col_east = (get_xy(h, 0) - get_xy(-h, 0)) / (2 * h)
    col_north = (get_xy(0, h) - get_xy(0, -h)) / (2 * h)
    j = np.column_stack((col_east, col_north))
    m.close()
    assert np.all(np.isfinite(j)), (group, j)
    return {'anchor_exposure_index': index, 'detector': detector,
            'jacobian_detector_pixels_per_arcsec': j.tolist(),
            'anchor_target_x_y': reference_xy.tolist()}


def map_one(coadd, jacobian):
    j = np.asarray(jacobian)
    dx = j[0, 0] * EN_EAST + j[0, 1] * EN_NORTH
    dy = j[1, 0] * EN_EAST + j[1, 1] * EN_NORTH
    valid = np.isfinite(coadd).astype(float)
    smoothed_data = gaussian_filter(np.nan_to_num(coadd, nan=0.),
                                    sigma=(1.3, 1.0))
    smoothed_valid = gaussian_filter(valid, sigma=(1.3, 1.0))
    mapped = map_coordinates(smoothed_data, [dy + HALF_Y, dx + HALF_X],
                             order=1, mode='constant', cval=0)
    weights = map_coordinates(smoothed_valid, [dy + HALF_Y, dx + HALF_X],
                              order=1, mode='constant', cval=0)
    return np.divide(mapped, weights, out=np.full_like(mapped, np.nan),
                     where=weights > 0.8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--coadds', default='/tmp/mocss_rate_registered_highpass_four_groups.npz')
    ap.add_argument('--output', default='janus_jwst/mocss_16004A_local_sky_projection_2026-10-16.json')
    ap.add_argument('--map-output', default='/tmp/mocss_16004_local_sky_maps.npz')
    args = ap.parse_args()
    logging.disable(logging.INFO)
    input_records = {}
    for detector in ['nrs1', 'nrs2']:
        source = Path('janus_jwst/mocss_16_exposure_' + detector + '_wcs_trace_2026-10-16.json')
        for r in json.loads(source.read_text())['records']:
            input_records[(r['exposure_index'], detector)] = r
    coadds = np.load(args.coadds)
    output = {'status': 'WCS_LOCAL_SKY_RATE_CONTRAST_NOT_SPECTRAL_IFU_MAP',
              'target_ra_dec': [RA, DEC],
              'grid_half_field_arcsec': HALF_FIELD,
              'grid_step_arcsec': SPACING,
              'wavelength_primary_angstrom': 16004,
              'wavelength_oiii_companion_angstrom': 15850.7933,
              'neighbors': {}, 'groups': {},
              'important_limitations': [
                  'Two-PA sky coordinates reconstructed from open-shutter MOS; no spectrophotometric calibration',
                  'Each map is a model of 16004-A photon arrival locations, not evidence for an emission beam',
                  '2D RATE background was high-pass filtered; artifacts/systematics require independent science-grade calibration',
                  'No corrected source intensity, extension limit or line detection significance']}
    for name, (ra, dec) in NEIGHBORS.items():
        output['neighbors'][name] = {
            'east_arcsec': (ra - RA) * 3600 * math.cos(math.radians(DEC)),
            'north_arcsec': (dec - DEC) * 3600}
    products = {'east_arcsec': EN_EAST, 'north_arcsec': EN_NORTH}
    for group in 'ABCD':
        info = group_wcs_info(group, input_records)
        j = info['jacobian_detector_pixels_per_arcsec']
        group_products = {}
        for specname in ('16004A', '15850p8A'):
            source = coadds[group + '_' + specname]
            sky = map_one(source, j)
            products[group + '_' + specname] = sky
            points = {'target': float(sky[len(GRID)//2, len(GRID)//2])}
            for key, coords in output['neighbors'].items():
                e,n=coords['east_arcsec'], coords['north_arcsec']
                k_e = int(round((e + HALF_FIELD) / SPACING))
                k_n = int(round((n + HALF_FIELD) / SPACING))
                points[key] = float(sky[k_n, k_e])
            group_products[specname] = {'sampled_rate_contrast_dn_per_s': points}
        output['groups'][group] = dict(info, line_maps=group_products)
        print('SKY-MAPPED GROUP', group, 'WCS pixels/arcsec',
              np.round(np.asarray(j), 2).tolist(), flush=True)
    np.savez_compressed(args.map_output, **products)
    Path(args.output).write_text(json.dumps(output, indent=2) + '\n')
    print('WROTE', args.output, args.map_output,
          'GEOMETRY+RAW_RATE CONTROL ONLY')


if __name__ == '__main__':
    main()
