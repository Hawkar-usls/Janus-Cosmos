#!/usr/bin/env python3
"""Fit local MSA row/column to sky-plane offsets in two sampled orientations.

The affine model is an INTERPOLATION proposal, not a final astrometric
solution. Its output should seed a finer independent assign_wcs probe.
"""

import json
import numpy as np
from pathlib import Path


def solve(record):
    candidates = []
    for quadrant in range(1, 5):
        points = [x for x in record['points']
                  if x['projects_on_detector'] and x['quadrant'] == quadrant]
        if len(points) < 5:
            continue
        xy = np.array([[1, x['row'], x['column']] for x in points], dtype=float)
        en = np.array([[x['offset_east_arcsec'], x['offset_north_arcsec']]
                       for x in points], dtype=float)
        model, _, _, _ = np.linalg.lstsq(xy, en, rcond=None)
        center = np.linalg.solve(model[1:3].T, -model[0])
        rms = np.sqrt(np.mean(np.sum((xy @ model - en)**2, axis=1)))
        candidates.append({'quadrant': quadrant, 'row': float(center[0]),
                           'column': float(center[1]),
                           'rms_fit_arcsec': float(rms),
                           'num_reference_slits': len(points),
                           'in_shutter_array': bool(1 <= center[0] <= 365
                                                    and 1 <= center[1] <= 171)})
    return {'source_rate': record['source_rate'],
            'candidates': candidates}


if __name__ == '__main__':
    filename = Path('janus_jwst/wcs_probe_16exposure_2026-10-16.json')
    original = json.loads(filename.read_text())
    result = {'status': 'AFFINE_APPROXIMATION_FOR_TARGET_SEARCH',
              'records': [solve(x) for x in original['records']]}
    outfile = Path('janus_jwst/shutter_proposals_2026-10-16.json')
    outfile.write_text(json.dumps(result, indent=2) + '\n')
    for record in result['records']:
        print(record['source_rate'])
        for row in record['candidates']:
            print('Q', row['quadrant'], 'row', round(row['row'], 2),
                  'col', round(row['column'], 2), 'rms',
                  round(row['rms_fit_arcsec'], 3),
                  'in array', row['in_shutter_array'])
