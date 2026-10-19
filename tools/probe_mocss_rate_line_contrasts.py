#!/usr/bin/env python3
"""Conservative detector-DN/s checks at sky-registered MOCSS wavelengths.

This is an exploratory RATE *contrast*, not a pipeline-calibrated spectrum,
not an absolute line flux, and not a detection/significance assessment.
The local controls sample nearby detector columns, which are not independent
sky sources; detector systematics and spatially extended contamination remain.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits


def mad(a):
    x = np.asarray(a, dtype=float)
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if len(x) else float('nan')


def sample_line(sci, dq, xc, yc, sigma_x=1.0, sigma_y=1.3):
    """Spatial/spectral matched-filter contrast against local per-row bias."""
    xi, yi = int(round(xc)), int(round(yc))
    if xi - 72 < 0 or xi + 72 >= sci.shape[1] or yi - 8 < 0 or yi + 8 >= sci.shape[0]:
        return {'reason': 'stamp_outside_detector'}
    xx = np.arange(xi - 70, xi + 71)
    yy = np.arange(yi - 5, yi + 6)
    patch = np.asarray(sci[np.ix_(yy, xx)], dtype=float)
    flags = np.asarray(dq[np.ix_(yy, xx)])
    patch[((flags & 3) != 0) | (~np.isfinite(patch))] = np.nan
    baseline = np.nanmedian(patch[:, np.abs(xx - xc) > 9], axis=1)
    residual = patch - baseline[:, None]
    weights_y = np.exp(-0.5 * ((yy - yc) / sigma_y) ** 2)
    weights_x = np.exp(-0.5 * ((xx - xc) / sigma_x) ** 2)

    def score_at(x):
        w = weights_y[:, None] * np.exp(-0.5 * ((xx - x) / sigma_x) ** 2)[None, :]
        ok = np.isfinite(residual)
        return float(np.nansum(np.where(ok, residual * w, np.nan)) /
                     np.sum(np.where(ok, w, 0.))) if np.any(ok) else float('nan')

    on = score_at(xc)
    controls = [score_at(float(x)) for x in range(xi - 58, xi + 59)
                if abs(x - xc) > 12]
    scatter = mad(controls)
    return {'line_contrast_dn_per_s': on,
            'local_control_median_dn_per_s': float(np.nanmedian(controls)),
            'empirical_local_mad_dn_per_s': scatter,
            'heuristic_robust_z_unvalidated':
                (on - float(np.nanmedian(controls))) / scatter if scatter > 0 else None,
            'n_local_control_columns': len(controls),
            'num_valid_stamp_pixels': int(np.isfinite(patch).sum()),
            'source_detector_x': xc, 'source_detector_y': yc}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--wcs-nrs1', default='janus_jwst/mocss_16_exposure_nrs1_wcs_trace_2026-10-16.json')
    p.add_argument('--wcs-nrs2', default='janus_jwst/mocss_16_exposure_nrs2_wcs_trace_2026-10-16.json')
    p.add_argument('--output', default='janus_jwst/mocss_rate_16_exposure_line_contrasts_2026-10-16.json')
    args = p.parse_args()
    geometry = {}
    for detector, source in [('nrs1', args.wcs_nrs1), ('nrs2', args.wcs_nrs2)]:
        d = json.loads(Path(source).read_text())
        assert len(d['records']) == 16, (source, len(d['records']))
        for r in d['records']:
            geometry[(r['exposure_index'], detector)] = r
    results = []
    for i in range(1, 17):
        frame = {'exposure_index': i, 'detectors': {}}
        for det in ('nrs1', 'nrs2'):
            meta = geometry[(i, det)]
            geometry_16004 = meta.get('trace', {}).get('16004A', {})
            edge = (geometry_16004.get('on_detector')
                    and abs(geometry_16004.get('slit_x', float('inf'))) < 0.75
                    and abs(geometry_16004.get('slit_y', float('inf'))) < 0.75)
            if not edge:
                frame['detectors'][det] = {
                    'status': 'NOT_RELIABLY_COVERED_BY_SELECTED_WCS_SHUTTER',
                    'coordinate': geometry_16004,
                    'wcs_error': meta.get('wcs_error')}
                continue
            target_rate = Path('/tmp') / meta['rate']
            with fits.open(target_rate, memmap=False) as hdul:
                sci, dq = hdul['SCI'].data, hdul['DQ'].data
                lines = {}
                for name, coord in meta['trace'].items():
                    if coord.get('on_detector') and (abs(coord.get('slit_x', float('inf'))) < 0.75
                                                      and abs(coord.get('slit_y', float('inf'))) < 0.75):
                        lines[name] = sample_line(sci, dq, coord['x_px'], coord['y_px'])
                    else:
                        lines[name] = {'reason': 'target_not_on_detector_or_inside_shutter'}
            frame['detectors'][det] = {'status': ('RATE_CONTRAST_MEASURED_STRICT_SHUTTER'
                                                 if geometry_16004['inside_shutter'] else
                                                 'RATE_CONTRAST_EDGE_SHUTTER_UNCERTAIN'),
                                       'rate': meta['rate'],
                                       'shutter': meta['best_shutter'],
                                       'line_probes': lines}
        results.append(frame)
        statuses = [det for det in ('nrs1', 'nrs2')
                    if frame['detectors'][det]['status'].startswith('RATE_CONTRAST_')]
        vals = [(det,
                 frame['detectors'][det]['line_probes'].get('16004A', {}).get(
                     'heuristic_robust_z_unvalidated')) for det in statuses]
        print('FRAME', i, 'detectors/raw-line-contrasts', vals, flush=True)
    report = {'status': 'UNVALIDATED_UNCALIBRATED_DETECTOR_CONTRAST',
              'units': 'original NRS RATE science array DN/s',
              'not_significance': 'heuristic z values are not Gaussian sigmas; no pipeline 1/f/flat/pathloss/error correction',
              'limitations': ['No standard JWST calwebb_spec2 science calibration',
                              'No row-based MOCSS background map or forward PSF fit',
                              'No differential MSA shutter throughput correction',
                              'No trials, injection recovery, independent sky controls'],
              'records': results}
    Path(args.output).write_text(json.dumps(report, indent=2) + '\n')
    print('WROTE', args.output, 'HEURISTIC CONTRAST ONLY')


if __name__ == '__main__':
    main()
