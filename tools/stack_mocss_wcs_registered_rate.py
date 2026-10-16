#!/usr/bin/env python3
"""Exploratory WCS-registered 2D RATE cutouts for four MOCSS dither stacks.

NRS2: A/B, NRS1: C/D, as independently solved from archival MSA+CRDS WCS.
The output is background-subtracted in raw DN/s, not a `calwebb_spec2`
flux-calibrated spectrum. Detection statistics are deliberately provisional.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.ndimage import map_coordinates, median_filter


HALF_X, HALF_Y = 76, 26
XX = np.arange(-HALF_X, HALF_X + 1)
YY = np.arange(-HALF_Y, HALF_Y + 1)
GRID_Y, GRID_X = np.meshgrid(YY, XX, indexing='ij')


def sample_stamp(filename, x0, y0, background_method, median_width):
    with fits.open('/tmp/' + filename, memmap=False) as h:
        science = np.asarray(h['SCI'].data, dtype=float)
        dq = h['DQ'].data
        valid = np.isfinite(science) & ((dq & 3) == 0)
        sample_coords = [GRID_Y + y0, GRID_X + x0]
        numer = map_coordinates(np.where(valid, science, 0), sample_coords,
                                order=1, mode='constant', cval=0)
        denom = map_coordinates(valid.astype(float), sample_coords,
                                order=1, mode='constant', cval=0)
        interp = np.divide(numer, denom, out=np.full_like(numer, np.nan),
                           where=denom > 0.9)
    # Each detector row has its own pedestal and smooth wavelength-dependent
    # residual. Median filtering across a 41-pixel x-window is a CONTROL,
    # not a replacement for the published MOCSS row-wise sky model + flats.
    if background_method == 'local_median_41':
        row_median = np.nanmedian(interp, axis=1)
        filled = np.where(np.isfinite(interp), interp, row_median[:, None])
        pedestal = median_filter(filled, size=(1, median_width), mode='nearest')
    else:
        mask = np.abs(XX) > 12
        pedestal = np.nanmedian(interp[:, mask], axis=1)[:, None]
    residual = interp - pedestal
    return residual


def empirical_score(image, px, py):
    x = np.arange(image.shape[1]) - HALF_X
    y = np.arange(image.shape[0]) - HALF_Y
    weights_x = np.exp(-0.5 * ((x - px) / 1.0) ** 2)
    weights_y = np.exp(-0.5 * ((y - py) / 1.3) ** 2)
    weighted = weights_y[:, None] * weights_x[None, :]
    okay = np.isfinite(image)
    total = np.sum(np.where(okay, weighted, 0))
    return float(np.nansum(image * weighted) / total) if total > 0 else float('nan')


def summarize_line(image, predicted_dx, predicted_dy=0):
    on = empirical_score(image, predicted_dx, predicted_dy)
    control = [empirical_score(image, dx, predicted_dy)
               for dx in np.arange(-56, 57, 1)
               if abs(dx - predicted_dx) >= 12]
    control = np.asarray(control)
    baseline = float(np.nanmedian(control))
    mad = float(1.4826 * np.nanmedian(np.abs(control - baseline)))
    # The local +/-4x+/-3y peak search has a significant look-elsewhere effect.
    search = [(empirical_score(image, predicted_dx + dx, predicted_dy + dy), dx, dy)
              for dx in range(-4, 5) for dy in range(-3, 4)]
    peak, ddx, ddy = max(search, key=lambda p: p[0] if np.isfinite(p[0]) else -1e99)
    return {'expected_center_dn_per_s': on,
            'local_control_median_dn_per_s': baseline,
            'empirical_mad_dn_per_s': mad,
            'provisional_center_z_not_gaussian_sigma': (on - baseline) / mad if mad > 0 else None,
            'searched_peak_dn_per_s': peak,
            'searched_peak_delta_x_px': ddx,
            'searched_peak_delta_y_px': ddy,
            'searched_peak_z_uncorrected': (peak - baseline) / mad if mad > 0 else None,
            'warning': 'Potential peak search is 63 trials; heuristic z is not formal significance'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output', default='janus_jwst/mocss_wcs_registered_rate_stack_2026-10-16.json')
    p.add_argument('--image-output', default='/tmp/mocss_rate_registered_four_groups.npz')
    p.add_argument('--background-method', choices=['row_constant', 'local_median_41'],
                   default='row_constant')
    p.add_argument('--median-width', type=int, default=41)
    args = p.parse_args()
    assert 21 <= args.median_width <= 101 and args.median_width % 2 == 1
    geometry = {}
    for detector in ['nrs1', 'nrs2']:
        path = Path('janus_jwst/mocss_16_exposure_' + detector + '_wcs_trace_2026-10-16.json')
        for entry in json.loads(path.read_text())['records']:
            geometry[(entry['exposure_index'], detector)] = entry
    outputs = {'status': 'EXPLORATORY_WCS_REGISTERED_RATE_STACK_NOT_CALIBRATED_FLUX',
               'wcs_context': 'JWST 1.19.1 / jwst_1413.pmap',
               'input_nrs1_nrs2_frames': 32,
               'background_method': args.background_method,
               'median_filter_width_detector_pixels': args.median_width,
               'limitations': [
                   'No full standard stage 2 calibration, 2D MOCSS background or variance propagation',
                   'No shutter bar transmission correction; half of dithers near shutter edges',
                   'No validated source significance, neighboring-sky mosaics or trial correction',
                   'Subpixel resampling can correlate noise; DQ only rejects DO_NOT_USE/SATURATED'],
               'groups': []}
    images = {}
    for group, indices in [('A', range(1, 5)), ('B', range(5, 9)),
                           ('C', range(9, 13)), ('D', range(13, 17))]:
        detector = 'nrs2' if group in ('A', 'B') else 'nrs1'
        tiles = {}
        entry_records = []
        for idx in indices:
            rec = geometry[(idx, detector)]
            coords = rec.get('trace', {})
            label = {'exposure_index': idx, 'rate': rec['rate'],
                     'strictly_inside_shutter': rec.get('physically_covered_by_selected_shutter')}
            for wavelength in ('16004A', '15850p8A', '13656p75A_Pagamma'):
                c = coords.get(wavelength, {})
                if c.get('on_detector'):
                    tiles.setdefault(wavelength, []).append(
                        sample_stamp(rec['rate'], c['x_px'], c['y_px'],
                                     args.background_method, args.median_width))
            entry_records.append(label)
        group_report = {'group': group, 'detector': detector,
                        'frames': entry_records, 'lines': {}}
        for label, stacks in tiles.items():
            if not stacks:
                continue
            with np.errstate(invalid='ignore'):
                coadd = np.nanmean(stacks, axis=0)
            images[group + '_' + label] = coadd.astype('float32')
            report = summarize_line(coadd, 0.)
            report['stacked_num_input_frames'] = len(stacks)
            group_report['lines'][label] = report
            print('STACK', group, label, len(stacks),
                  'center raw heuristic', round(report['provisional_center_z_not_gaussian_sigma'], 2),
                  'scan peak dx/dy', report['searched_peak_delta_x_px'],
                  report['searched_peak_delta_y_px'], flush=True)
        outputs['groups'].append(group_report)
    np.savez_compressed(args.image_output, **images)
    Path(args.output).write_text(json.dumps(outputs, indent=2) + '\n')
    print('WROTE', args.output, args.image_output,
          'UNVALIDATED RATE-SPACE 2D COADDS ONLY')


if __name__ == '__main__':
    main()
