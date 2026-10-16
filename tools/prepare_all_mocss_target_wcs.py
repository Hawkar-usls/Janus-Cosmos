#!/usr/bin/env python3
"""Prepare targeted, read-only-reproducible WCS probes for all 16 GO-3290 visits.

Local affine sky maps from two full JWST WCS anchors predict the shutter
location in each dither using its SCI RA_REF/DEC_REF. The subsequent official
assign_wcs run must independently test and refine the prediction.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits
from prepare_mocss_wcs_probe import create_probe


DEC = -27.78908


def affine_anchor(probe_record, quadrant):
    p = [x for x in probe_record['points']
         if x['projects_on_detector'] and x['quadrant'] == quadrant]
    xy = np.array([[1, x['row'], x['column']] for x in p])
    en = np.array([[x['offset_east_arcsec'], x['offset_north_arcsec']] for x in p])
    model, _, _, _ = np.linalg.lstsq(xy, en, rcond=None)
    return model


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', default='domains/space/MOCSS_16_EXPOSURE_MANIFEST_2026-10-16.json')
    parser.add_argument('--wcs-grids', default='janus_jwst/wcs_probe_16exposure_2026-10-16.json')
    parser.add_argument('--half-width', type=int, default=4)
    parser.add_argument('--output-dir', default='/tmp/janus_mocss_wcs_all16')
    parser.add_argument('--detector', choices=['nrs1', 'nrs2'], default='nrs1')
    args = parser.parse_args()
    records = json.loads(Path(args.manifest).read_text())['exposures']
    anchors = json.loads(Path(args.wcs_grids).read_text())['records']
    anchor_2023aug = anchors[0]
    anchor_2023nov = anchors[1]
    models = {1: affine_anchor(anchor_2023aug, 2),
              2: affine_anchor(anchor_2023nov, 4)}
    anchor_rate_paths = {1: Path('/tmp/jw03290001001_03109_00002_' + args.detector + '_rate.fits'),
                         2: Path('/tmp/jw03290002001_11101_00002_' + args.detector + '_rate.fits')}
    reference = {}
    for orientation, path in anchor_rate_paths.items():
        h = fits.getheader(path, 'SCI')
        reference[orientation] = (h['RA_REF'], h['DEC_REF'])
    output = []
    for exposure in records:
        orientation = exposure['orientation']
        q = 2 if orientation == 1 else 4
        rate = Path('/tmp/' + exposure['exposure_prefix'] + '_' + args.detector + '_rate.fits')
        msa = Path('/tmp/jw03290001001_05_msa.fits' if orientation == 1
                   else '/tmp/jw03290002001_01_msa.fits')
        hdr = fits.getheader(rate, 'SCI')
        old_ra, old_dec = reference[orientation]
        de = (hdr['RA_REF'] - old_ra) * 3600 * math.cos(math.radians(DEC))
        dn = (hdr['DEC_REF'] - old_dec) * 3600
        model = models[orientation]
        center = np.linalg.solve(model[1:3].T, -(model[0] + np.array([de, dn])))
        guess = {'quadrant': q, 'row': float(center[0]), 'column': float(center[1]),
                 'prediction_method': 'reference-RADEC-shift-affine-WCS',
                 'anchor_exposure': anchor_rate_paths[orientation].name,
                 'delta_east_arcsec': de, 'delta_north_arcsec': dn}
        assert 8 <= center[0] <= 357 and 8 <= center[1] <= 163, guess
        result = create_probe(rate, msa, Path(args.output_dir),
                              grid=3, target=guess, half_width=args.half_width)
        result['exposure_index'] = exposure['exposure_index']
        result['group'] = exposure['group']
        result['orientation'] = orientation
        output.append(result)
        print('PREPARED', exposure['exposure_index'], exposure['group'],
              'Q', q, 'row/col', *(round(float(x), 2) for x in center), flush=True)
    destination = Path(args.output_dir) / 'all16_input_receipt.json'
    destination.write_text(json.dumps({'status': '16_WCS_TEST_INPUTS_ONLY',
                                       'records': output}, indent=2) + '\n')
    print('PREPARED 16/16 virtual-shutter WCS test models, NOT SCIENCE EXTRACTION')


if __name__ == '__main__':
    main()
