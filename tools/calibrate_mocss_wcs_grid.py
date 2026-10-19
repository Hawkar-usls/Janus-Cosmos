#!/usr/bin/env python3
"""CRDS-pinned sky positions of MOCSS virtual shutters; NO line detection.

Requires `prepare_mocss_wcs_probe.py` generated temporary FITS and official
JWST 1.19.1 pipeline with a locally synced jwst_1413.pmap CRDS context.
"""

import argparse
import json
import logging
import math
from pathlib import Path

from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.nirspec import nrs_wcs_set_input
import jwst


TARGET_RA = 53.13127
TARGET_DEC = -27.78908
TARGET_WAVELENGTH_M = 1.6004e-6


def offset_arcsec(ra, dec):
    return ((ra - TARGET_RA) * 3600 * math.cos(math.radians(TARGET_DEC)),
            (dec - TARGET_DEC) * 3600)


def calibrate_one(record):
    model = AssignWcsStep.call(record['probe_rate'], save_results=False)
    assert model.meta.cal_step.assign_wcs == 'COMPLETE'
    found = {int(n) for n in model.meta.wcs.get_transform(
        'gwa', 'slit_frame').slit_ids}
    points = []
    for item in record['slit_positions']:
        entry = dict(item)
        entry['projects_on_detector'] = item['slit_id'] in found
        if entry['projects_on_detector']:
            w = nrs_wcs_set_input(model, item['slit_id'])
            transform = w.get_transform('slit_frame', 'world')
            ra, dec, lam = transform(0.0, 0.0, TARGET_WAVELENGTH_M)
            dx, dy = offset_arcsec(float(ra), float(dec))
            entry.update(sky_ra_deg=float(ra), sky_dec_deg=float(dec),
                         wavelength_m=float(lam), offset_east_arcsec=dx,
                         offset_north_arcsec=dy,
                         angular_distance_arcsec=math.hypot(dx, dy),
                         detector_bbox={k: [float(v.lower), float(v.upper)]
                                        for k, v in w.bounding_box.intervals.items()})
        points.append(entry)
    model.close()
    nearest = sorted((p for p in points if p['projects_on_detector']),
                     key=lambda p: p['angular_distance_arcsec'])[:10]
    return {'source_rate': record['source_rate'],
            'msametid_original': record['original_msa_metadata_id'],
            'msametid_int16': record['normalized_signed16_id'],
            'total_virtual_slits': len(points),
            'projected_slits': len(found),
            'nearest_probes': nearest,
            'points': points}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default='/tmp/janus_mocss_wcs_probe/probe_input_nrs2_receipt.json')
    parser.add_argument('--output', default='janus_jwst/wcs_probe_16exposure_2026-10-16.json')
    args = parser.parse_args()
    logging.getLogger('jwst').setLevel(logging.WARNING)
    receipt = json.loads(Path(args.input).read_text())
    output = {'status': 'GEOMETRIC_WCS_PROBE_NOT_SCIENCE_EXTRACTION',
              'pipeline_version': jwst.__version__,
              'crds_context': 'jwst_1413.pmap',
              'target_icrs_ra_deg': TARGET_RA,
              'target_icrs_dec_deg': TARGET_DEC,
              'line_wavelength_angstrom': 16004.0,
              'probe_centers_wavelength_input_m': TARGET_WAVELENGTH_M,
              'scientific_limitations': [
                  'Only 2 of 16 exposures, NRS2 only, virtual reference shutters',
                  'Test metadata contains selected virtual shutters, not actual extracted spectra',
                  'No background subtraction, line detection, spectral or imaging association',
                  'Original FITS files left unchanged'],
              'records': [calibrate_one(rec) for rec in receipt['records']]}
    Path(args.output).write_text(json.dumps(output, indent=2) + '\n')
    for row in output['records']:
        print(row['source_rate'], 'projected', row['projected_slits'],
              'nearest arcsec', [round(p['angular_distance_arcsec'], 3)
                                for p in row['nearest_probes'][:3]], flush=True)
    print('WROTE', args.output, 'NO LINE MAP', flush=True)


if __name__ == '__main__':
    main()
