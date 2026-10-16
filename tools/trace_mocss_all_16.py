#!/usr/bin/env python3
"""Official 1.19.1 WCS inversion for all sixteen GO-3290 NRS1 RATE frames.

Output is detector/sky geometry only: it cannot establish an emission line.
Each frame receives its own MSA metadata and SCI pointing; virtual shutter
IDs are never treated as real source detections.
"""

import argparse
import json
import logging
import math
from pathlib import Path

import numpy as np
from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.nirspec import nrs_wcs_set_input


RA, DEC = 53.13127, -27.78908
WAVELENGTHS = {'16004A': 1.6004, '15850p8A': 1.58507933,
               '13656p75A_Pagamma': 1.36567483}


def coords(w, wl):
    r = {'wavelength_um': wl}
    try:
        x, y = w.invert(RA, DEC, wl)
        sx, sy, lam = w.get_transform('world', 'slit_frame')(RA, DEC, wl)
        r.update(x_px=float(x), y_px=float(y),
                 slit_x=float(sx), slit_y=float(sy),
                 slit_wavelength_m=float(lam),
                 on_detector=bool(np.isfinite(x) and np.isfinite(y)
                                  and 0 <= x < 2048 and 0 <= y < 2048),
                 inside_shutter=bool(abs(sx) <= 0.5 and abs(sy) <= 0.5))
    except (ValueError, ArithmeticError, RuntimeError) as exc:
        r['error'] = str(exc)
    return r


def analyze_frame(rec):
    model = AssignWcsStep.call(rec['probe_rate'], save_results=False)
    valid_ids = set(map(int, model.meta.wcs.get_transform('gwa', 'slit_frame').slit_ids))
    candidates = []
    for virtual in rec['slit_positions']:
        if virtual['slit_id'] not in valid_ids:
            continue
        w = nrs_wcs_set_input(model, virtual['slit_id'])
        ra, dec, _ = w.get_transform('slit_frame', 'world')(0., 0., 1.6004e-6)
        dist = math.hypot((ra - RA) * 3600 * math.cos(math.radians(DEC)),
                          (dec - DEC) * 3600)
        if dist > 0.75:
            continue
        row = dict(virtual)
        row['shutter_center_sep_arcsec'] = dist
        row['wcs_16004A'] = coords(w, WAVELENGTHS['16004A'])
        candidates.append(row)
    candidates.sort(key=lambda x: x['shutter_center_sep_arcsec'])
    covered = [c for c in candidates if c['wcs_16004A'].get('on_detector')
               and c['wcs_16004A'].get('inside_shutter')]
    best = covered[0] if covered else (candidates[0] if candidates else None)
    result = {'exposure_index': rec['exposure_index'], 'group': rec['group'],
              'orientation': rec['orientation'], 'rate': rec['source_rate'],
              'num_virtual_slits': len(rec['slit_positions']),
              'num_projected_slits': len(valid_ids),
              'nearest_shutters': candidates[:8],
              'physically_covered_by_selected_shutter': bool(covered)}
    if best:
        w = nrs_wcs_set_input(model, best['slit_id'])
        result['best_shutter'] = {key: best[key] for key in (
            'slit_id', 'quadrant', 'row', 'column', 'shutter_center_sep_arcsec')}
        result['trace'] = {name: coords(w, wavelength)
                           for name, wavelength in WAVELENGTHS.items()}
    model.close()
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--inputs', default='/tmp/janus_mocss_wcs_all16/all16_input_receipt.json')
    p.add_argument('--output', default='janus_jwst/mocss_16_exposure_nrs1_wcs_trace_2026-10-16.json')
    args = p.parse_args()
    logging.disable(logging.INFO)
    inputs = json.loads(Path(args.inputs).read_text())['records']
    receipt = {'status': 'WCS_DETECTOR_COORDINATES_NOT_LINE_DETECTION',
               'crds_context': 'jwst_1413.pmap',
               'pipeline_version': '1.19.1',
               'target_icrs_deg': [RA, DEC],
               'wavelengths_um': WAVELENGTHS,
               'records': []}
    out = Path(args.output)
    for i, rec in enumerate(inputs, 1):
        try:
            r = analyze_frame(rec)
        except Exception as exc:
            r = {'exposure_index': rec['exposure_index'],
                 'group': rec['group'], 'rate': rec['source_rate'],
                 'wcs_error': f'{type(exc).__name__}: {exc}'}
        receipt['records'].append(r)
        out.write_text(json.dumps(receipt, indent=2) + '\n')
        trace = r.get('trace', {}).get('16004A', {})
        print('WCS', i, rec['group'], 'on detector', trace.get('on_detector'),
              'in shutter', trace.get('inside_shutter'),
              'x/y', round(trace.get('x_px', float('nan')), 2),
              round(trace.get('y_px', float('nan')), 2),
              'error', r.get('wcs_error', ''), flush=True)
    assert len(receipt['records']) == 16
    print('WROTE', out, 'WCS GEOMETRY ONLY')


if __name__ == '__main__':
    main()
