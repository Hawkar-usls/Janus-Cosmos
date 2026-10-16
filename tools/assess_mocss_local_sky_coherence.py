#!/usr/bin/env python3
"""Descriptive sky-location consistency check, with strong non-detection caveats."""

import json
from pathlib import Path

import numpy as np


def rank_diagnostic(maps, name):
    xx = np.asarray(maps['east_arcsec'])
    yy = np.asarray(maps['north_arcsec'])
    mask = np.hypot(xx, yy) > 0.45
    groups = {}
    normalized = []
    for group in 'ABCD':
        z = np.asarray(maps[group + '_' + name], dtype=float)
        base = np.nanmedian(z[mask])
        mad = 1.4826 * np.nanmedian(np.abs(z[mask] - base))
        norm = (z - base) / mad if mad > 0 else z * float('nan')
        normalized.append(norm)
        groups[group] = {'background_median_dn_per_s': float(base),
                         'background_mad_dn_per_s': float(mad),
                         'target_relative_robust_units': float(norm[len(z)//2,len(z)//2])}
    full = np.nanmean(normalized, axis=0)
    valid = np.isfinite(full)
    center = float(full[len(full)//2, len(full)//2])
    center_rank_fraction_above = float(np.mean(full[valid] >= center))
    highest = np.unravel_index(np.nanargmax(full), full.shape)
    return {'groups': groups,
            'combined_target_relative_robust_units': center,
            'fraction_of_map_cells_at_least_as_high_as_target_not_pvalue': center_rank_fraction_above,
            'max_map_east_north_arcsec': [float(xx[highest]), float(yy[highest])],
            'max_map_relative_robust_units': float(full[highest]),
            'combined_map': full}


def main():
    f = np.load('/tmp/mocss_16004_local_sky_maps.npz')
    line = rank_diagnostic(f, '16004A')
    control = rank_diagnostic(f, '15850p8A')
    receipt = json.loads(Path('janus_jwst/mocss_16004A_local_sky_projection_2026-10-16.json').read_text())
    xx = np.asarray(f['east_arcsec'])
    yy = np.asarray(f['north_arcsec'])
    neighbors = {}
    for name, coord in receipt['neighbors'].items():
        position = np.argmin((xx - coord['east_arcsec'])**2 + (yy - coord['north_arcsec'])**2)
        row, col = np.unravel_index(position, xx.shape)
        neighbors[name] = {'16004A_relative_robust_units': float(line['combined_map'][row,col]),
                           '15850p8A_relative_robust_units': float(control['combined_map'][row,col])}
    for obj in (line, control):
        del obj['combined_map']
    output = {'status': 'DESCRIPTIVE_UNVALIDATED_LOCAL_SKY_COMPARISON',
              'source_field_arcsec': 3.0, 'source_location_preselected': True,
              'wavelength_16004A': line, 'companion_control_15850p8A': control,
              'neighbors': neighbors,
              'interpretation_guardrails': [
                  'A rank fraction is NOT a false-alarm probability or significance',
                  'WCS registration was tested but no standard spectrophotometric calibration',
                  'Spatial cells and detector residuals are highly correlated',
                  'No physical emission length or propagation direction can be measured',
                  'Validation requires full flat, source extraction, simulated injection and independent epoch or IFU']}
    output_path = Path('janus_jwst/mocss_local_sky_coherence_receipt_2026-10-16.json')
    output_path.write_text(json.dumps(output, indent=2) + '\n')
    print('LINE rank', line['fraction_of_map_cells_at_least_as_high_as_target_not_pvalue'],
          'target', line['combined_target_relative_robust_units'],
          'max offset', line['max_map_east_north_arcsec'])
    print('CONTROL rank', control['fraction_of_map_cells_at_least_as_high_as_target_not_pvalue'],
          'target', control['combined_target_relative_robust_units'])
    print('NEIGHBORS', neighbors)
    print('WROTE', output_path, 'NOT SCIENCE SIGNIFICANCE')


if __name__ == '__main__':
    main()
