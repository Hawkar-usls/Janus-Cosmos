#!/usr/bin/env python3
"""Create limited, reversible WCS probe inputs from original JWST MOCSS FITS.

The products have TEST_ prefixes and do not represent science extractions.
Each grid cell is a different virtual background slit and is not a detected
source. Original archival rate and MSA data are read, never modified.

For high MSAMETID (e.g. 37338), FITS SHUTTER_INFO stores the ID as signed
int16; the current JWST 1.19.1 compare uses the full primary-header integer.
Normalize only the test-rate metadata ID to signed16. Preserve the source ID
and the modulo mapping in the machine-readable provenance receipt.
"""

import argparse
from pathlib import Path
import json
import os
import shutil

import numpy as np
from astropy.io import fits


def create_probe(rate_path, msa_path, output_directory, grid, target=None, half_width=3):
    output_directory.mkdir(parents=True, exist_ok=True)
    with fits.open(rate_path, memmap=True) as src:
        original_id = int(src[0].header['MSAMETID'])
        dither = int(src[0].header['PATT_NUM'])
        signed16 = (original_id + 32768) % 65536 - 32768
        with fits.open(msa_path, memmap=True) as msa:
            data = msa['SHUTTER_INFO'].data
            subset = data[(data['msa_metadata_id'] == signed16) &
                          (data['dither_point_index'] == dither)]
            assert len(subset) > 200000
            if target:
                desired_rows = np.arange(round(target['row']) - half_width,
                                         round(target['row']) + half_width + 1)
                desired_cols = np.arange(round(target['column']) - half_width,
                                         round(target['column']) + half_width + 1)
                assert desired_rows.min() >= 1 and desired_rows.max() <= 365
                assert desired_cols.min() >= 1 and desired_cols.max() <= 171
            else:
                desired_rows = np.linspace(26, 335, grid, dtype=int)
                desired_cols = np.linspace(16, 156, grid, dtype=int)
            selected = []
            for quadrant in ([target['quadrant']] if target else range(1, 5)):
                avail = subset[subset['shutter_quadrant'] == quadrant]
                for row in desired_rows:
                    for col in desired_cols:
                        near = (np.abs(avail['shutter_row'].astype(int) - row) +
                                np.abs(avail['shutter_column'].astype(int) - col))
                        selected.append(tuple(avail[int(np.argmin(near))]))
            entries = np.array(selected, dtype=subset.dtype)
            assert len({(int(e['shutter_quadrant']), int(e['shutter_row']),
                        int(e['shutter_column'])) for e in entries}) == len(entries)
            for i, row in enumerate(entries):
                row['slitlet_id'] = i + 1
                row['background'] = 'Y'
                row['primary_source'] = 'N'
                row['source_id'] = 0
                row['estimated_source_in_shutter_x'] = np.nan
                row['estimated_source_in_shutter_y'] = np.nan
            suffix = f'target_h{half_width}' if target else f'grid{grid}'
            mini_msa = output_directory / f'TEST_{rate_path.stem}_MSA_{suffix}.fits'
            original_msa = msa['SHUTTER_INFO']
            hdu = fits.BinTableHDU(data=entries, header=original_msa.header,
                                   name='SHUTTER_INFO')
            out_hdus = fits.HDUList([fits.PrimaryHDU(header=msa[0].header),
                                    fits.ImageHDU(data=msa['SHUTTER_IMAGE'].data,
                                                  header=msa['SHUTTER_IMAGE'].header,
                                                  name='SHUTTER_IMAGE'),
                                    hdu,
                                    fits.BinTableHDU(data=msa['SOURCE_INFO'].data[:0],
                                                     header=msa['SOURCE_INFO'].header,
                                                     name='SOURCE_INFO')])
            out_hdus.writeto(mini_msa, overwrite=True, output_verify='silentfix')

        mini_rate = output_directory / f'TEST_{rate_path.stem}_{suffix}.fits'
        shutil.copyfile(rate_path, mini_rate)
        with fits.open(mini_rate, mode='update', memmap=False) as probe:
            probe[0].header['MSAMETFL'] = mini_msa.name
            probe[0].header['MSAMETID'] = signed16
        return {'source_rate': rate_path.name, 'source_msa': msa_path.name,
                'virtual_probe': True, 'science_line_extracted': False,
                'original_msa_metadata_id': original_id,
                'normalized_signed16_id': signed16, 'dither': dither,
                'grid_dimension_per_quadrant': grid if not target else None,
                'target_shutter_proposal': target,
                'num_virtual_slits': len(entries),
                'probe_msa': str(mini_msa), 'probe_rate': str(mini_rate),
                'slit_positions': [
                    {'slit_id': i + 1, 'quadrant': int(row['shutter_quadrant']),
                     'row': int(row['shutter_row']), 'column': int(row['shutter_column'])}
                    for i, row in enumerate(entries)
                ]}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--data-dir', default='/tmp')
    p.add_argument('--output-dir', default='/tmp/janus_mocss_wcs_probe')
    p.add_argument('--grid', type=int, default=3)
    p.add_argument('--target-proposals', default=None,
                   help='JSON output of locate_target_shutters.py')
    p.add_argument('--half-width', type=int, default=3)
    p.add_argument('--detector', choices=['nrs1', 'nrs2'], default='nrs2')
    args = p.parse_args()
    root = Path(args.data_dir)
    out = Path(args.output_dir)
    examples = [
        (f'jw03290001001_03109_00002_{args.detector}_rate.fits', 'jw03290001001_05_msa.fits'),
        (f'jw03290002001_11101_00002_{args.detector}_rate.fits', 'jw03290002001_01_msa.fits'),
    ]
    proposals = None
    if args.target_proposals:
        proposals = json.loads(Path(args.target_proposals).read_text())
    receipts = []
    for rate, msa in examples:
        target = None
        if proposals:
            matches = [rec for rec in proposals['records']
                       if rec['source_rate'].replace('nrs2', args.detector) == rate]
            assert len(matches) == 1
            valid = [x for x in matches[0]['candidates'] if x['in_shutter_array']]
            assert len(valid) == 1
            target = valid[0]
        receipts.append(create_probe(root / rate, root / msa, out,
                                     args.grid, target=target,
                                     half_width=args.half_width))
    receipt_name = ('probe_target_' + args.detector + '_receipt.json' if proposals
                    else 'probe_input_' + args.detector + '_receipt.json')
    (out / receipt_name).write_text(
        json.dumps({'status': 'WCS_PROBE_INPUT_PREPARED_ONLY',
                    'records': receipts}, indent=2) + '\n')
    print('PASS: Prepared', len(receipts), 'probe inputs with',
          ((2 * args.half_width + 1)**2 if proposals else args.grid**2 * 4),
          'virtual slits each. NO LINE MAP.')


if __name__ == '__main__':
    main()
