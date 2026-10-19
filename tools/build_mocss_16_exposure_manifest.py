#!/usr/bin/env python3
"""Enumerate every archival GO-3290 long MOCSS exposure and its metadata.

This is an input/provenance manifest, *not* an extracted 16004-Angstrom map.
Queries the official MAST CAOM snapshot downloaded by the JANUS subset fetcher
and reads the primary FITS header of both detectors for every exposure.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from io import BytesIO
import hashlib
import json
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

from astropy.io import fits


FILE_BASE = 'https://mast.stsci.edu/api/v0.1/Download/file?uri='
PRIMARY_KEYS = (
    'DETECTOR', 'EXP_TYPE', 'GRATING', 'FILTER', 'APERNAME', 'MSAMETFL',
    'MSAMETID', 'PROGRAM', 'OBSERVTN', 'VISIT', 'ACT_ID', 'EXPOSURE',
    'PATT_NUM', 'DATE-OBS', 'TIME-OBS', 'TARG_RA', 'TARG_DEC', 'PA_APER',
    'V3I_YANG', 'ROLL_REF', 'RA_REF', 'DEC_REF', 'EXPTIME', 'EFFEXPTM',
    'NINTS', 'NGROUPS', 'SUBARRAY',
)


def read_primary_header(obs_id, target_directory):
    name = obs_id + '_rate.fits'
    saved_file = target_directory / name
    if saved_file.is_file() and saved_file.stat().st_size > 100_000:
        with fits.open(saved_file, memmap=True) as hdus:
            header = hdus[0].header
        return {'source': 'local_full_rate_file', 'bytes': saved_file.stat().st_size,
                'primary_header': {key: header.get(key) for key in PRIMARY_KEYS}}
    uri = quote('mast:JWST/product/' + name, safe='')
    req = Request(FILE_BASE + uri, headers={'Range': 'bytes=0-32767',
                                          'User-Agent': 'JANUS-MOCSS-manifest/1.0'})
    with urlopen(req, timeout=90) as response:
        packet = response.read()
        code = response.status
        content_range = response.headers.get('Content-Range', '')
    if code != 206 or not content_range.startswith('bytes 0-32767/') or len(packet) != 32768:
        raise ValueError(f'FITS Range protocol violated for {name}: HTTP={code}; {content_range}')
    header = fits.Header.fromstring(packet.decode('ascii', errors='replace'), sep='')
    return {'source': 'MAST_range_0_32767', 'bytes': 32768,
            'header_range_sha256': hashlib.sha256(packet).hexdigest(),
            'primary_header': {key: header.get(key) for key in PRIMARY_KEYS}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--archive', default='/tmp/janus_jwst_3290_observations.json')
    parser.add_argument('--data-dir', default='/tmp')
    parser.add_argument('--output', default='janus_jwst/MOCSS_16_EXPOSURE_MANIFEST_2026-10-16.json')
    parser.add_argument('--threads', type=int, default=6)
    args = parser.parse_args()

    snapshot_bytes = Path(args.archive).read_bytes()
    archive = json.loads(snapshot_bytes)
    assert archive['status'] == 'COMPLETE', 'MAST query was not complete'
    long_rows = [row for row in archive['data'] if
                 row.get('proposal_id') == '3290' and
                 row.get('t_exptime', 0) > 8000 and
                 row.get('obs_id', '').startswith('jw03290') and
                 row.get('obs_id', '').split('_')[-1] in {'nrs1', 'nrs2'}]
    exposure_pairs = {}
    for row in long_rows:
        key = '_'.join(row['obs_id'].split('_')[:-1])
        exposure_pairs.setdefault(key, {})[row['obs_id'].split('_')[-1]] = row
    if len(exposure_pairs) != 16:
        raise ValueError(f'Expected 16 long exposures, got {len(exposure_pairs)}')
    for key, pair in exposure_pairs.items():
        if set(pair) != {'nrs1', 'nrs2'}:
            raise ValueError(f'One detector missing for {key}: {list(pair)}')

    obs_ids = [row['obs_id'] for row in long_rows]
    headers = {}
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {pool.submit(read_primary_header, oid, Path(args.data_dir)): oid for oid in obs_ids}
        for future in as_completed(futures):
            oid = futures[future]
            headers[oid] = future.result()
            print('HEADER', oid, headers[oid]['primary_header'].get('MSAMETFL'), flush=True)
    manifest = []
    for index, (key, pair) in enumerate(sorted(exposure_pairs.items(),
                                               key=lambda kv: min(x['t_min'] for x in kv[1].values())), 1):
        detectors = {which: {'obs_id': row['obs_id'],
                             'mast_obsid': row['obsid'],
                             'mast_t_min': row['t_min'],
                             'mast_exptime_seconds': row['t_exptime'],
                             **headers[row['obs_id']]}
                     for which, row in sorted(pair.items())}
        assert detectors['nrs1']['primary_header']['MSAMETID'] == detectors['nrs2']['primary_header']['MSAMETID']
        assert detectors['nrs1']['primary_header']['MSAMETFL'] == detectors['nrs2']['primary_header']['MSAMETFL']
        assert all(detectors[n]['primary_header']['GRATING'] == 'G140M' for n in detectors)
        manifest.append({'exposure_index': index,
                         'group': 'ABCD'[(index - 1) // 4],
                         'orientation': 1 if index <= 8 else 2,
                         'exposure_prefix': key,
                         'msa_file': detectors['nrs1']['primary_header']['MSAMETFL'],
                         'msa_metadata_id': detectors['nrs1']['primary_header']['MSAMETID'],
                         'detectors': detectors})

    out = {
        'schema_version': 'janus-mocss-16-exposure-input-inventory-v1',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'project': 'JWST GO-3290 MOCSS HUDF',
        'archive_endpoint': 'https://mast.stsci.edu/api/v0/invoke',
        'mast_snapshot_sha256': hashlib.sha256(snapshot_bytes).hexdigest(),
        'num_long_exposures': len(manifest),
        'num_detector_frames': len(long_rows),
        'num_distinct_msa_files': len(set(item['msa_file'] for item in manifest)),
        'target': {'ra_deg': 53.13127, 'dec_deg': -27.78908,
                   'observed_line_A': 16004, 'search_half_width_A': 10},
        'scope': 'input inventory and primary-FITS-header validation; NOT 2D/1D WCS-extracted line map',
        'scientific_blockers': [
            'Each MOCSS all-open MSA geometry must be filtered by msa_metadata_id and shutter sky position.',
            'JWST assign_wcs with CRDS jwst_1413.pmap and custom slitlet IDs must locate 16004-A traces.',
            'Row-based sky subtraction, within-group alignment, variance and dq combination still required.',
            'Independent local ±5-arcsec spatial and off-line control tests, injections and trials correction pending.',
        ],
        'exposures': manifest,
    }
    Path(args.output).write_text(json.dumps(out, indent=2, ensure_ascii=False) + '\n')
    print('SAVED', args.output, '16 exposures; MSA files:', out['num_distinct_msa_files'])


if __name__ == '__main__':
    main()
