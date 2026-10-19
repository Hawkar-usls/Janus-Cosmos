#!/usr/bin/env python3
"""Download a small, exact, provenance-checkable MAST subset for JWST HUDF.

The three mosaics total several GB, so this retrieves 901-row SCI byte ranges
covering the target and nearby blank-sky controls (~107 MB), NOT full mosaics.
Two original NIRSpec rate frames are included (~168 MB). Public archival data
only; no credentials or third-party mirrors are required.
"""

import argparse
import json
import os
from pathlib import Path
import tempfile
import time
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen


BASE = 'https://mast.stsci.edu/api/v0.1/Download/file?uri='
ROOT_URI = 'mast:JWST/product/'
PRODUCTS = [
    ('jw03290001001_03109_00002_nrs2_rate.fits', 'jw03290001001_03109_00002_nrs2_rate.fits', None),
    ('jw03290002001_11101_00002_nrs2_rate.fits', 'jw03290002001_11101_00002_nrs2_rate.fits', None),
    ('jw01180-o013_t011_nircam_clear-f150w_cat.ecsv', 'janus_jwst_2023_f150w_cat.ecsv', None),
    ('jw08060-o002_t002_nircam_clear-f150w_cat.ecsv', 'janus_jwst_2025_f150w_cat.ecsv', None),
    ('jw12577-o002_t002_nircam_clear-f150w_cat.ecsv', 'janus_jwst_2026_f150w_cat.ecsv', None),
    ('jw01180-o013_t011_nircam_clear-f150w_i2d.fits', 'janus_2023_i2d_header.bin', (0, 131071)),
    ('jw08060-o002_t002_nircam_clear-f150w_i2d.fits', 'janus_2025_i2d_header.bin', (0, 131071)),
    ('jw12577-o002_t002_nircam_clear-f150w_i2d.fits', 'janus_2026_i2d_header.bin', (0, 65535)),
    ('jw01180-o013_t011_nircam_clear-f150w_i2d.fits', 'janus_2023_f150w_stripe.bin', (93813672, 129875295)),
    ('jw08060-o002_t002_nircam_clear-f150w_i2d.fits', 'janus_2025_f150w_stripe.bin', (86662400, 122269919)),
    ('jw12577-o002_t002_nircam_clear-f150w_i2d.fits', 'janus_2026_f150w_stripe.bin', (83896000, 119503519)),
]
ATNF_CONES = [
    ('janus_atnf_cone.xml', '53.13127', '-27.78908', '0.1'),
    ('janus_atnf_cone_1deg.xml', '53.13127', '-27.78908', '1.0'),
    ('janus_atnf_control_crab.xml', '83.633', '22.0145', '0.1'),
]


def read_and_store(url, path, expected_range=None, force=False):
    if path.exists() and not force and (
        expected_range is None or path.stat().st_size == expected_range[1] - expected_range[0] + 1
    ):
        print('REUSE', path.name, path.stat().st_size)
        return
    headers = {'User-Agent': 'JANUS-SPACE-source-frozen-MAST-replay/1.0'}
    if expected_range:
        headers['Range'] = f'bytes={expected_range[0]}-{expected_range[1]}'
    for attempt in range(3):
        tmp = None
        try:
            req = Request(url, headers=headers)
            with urlopen(req, timeout=210) as response:
                if expected_range:
                    actual = response.headers.get('Content-Range', '')
                    expected = f'bytes {expected_range[0]}-{expected_range[1]}/'
                    if response.status != 206 or not actual.startswith(expected):
                        raise RuntimeError(f'Remote service ignored exact FITS byte range: {actual}')
                with tempfile.NamedTemporaryFile(dir=path.parent, prefix='mast_partial_', delete=False) as out:
                    tmp = Path(out.name)
                    while True:
                        chunk = response.read(4 * 1024 * 1024)
                        if not chunk:
                            break
                        out.write(chunk)
            if expected_range and tmp.stat().st_size != expected_range[1] - expected_range[0] + 1:
                raise RuntimeError(f'Truncated range {path.name}: {tmp.stat().st_size} bytes')
            os.replace(tmp, path)
            print('FETCHED', path.name, path.stat().st_size)
            return
        except Exception:
            if tmp and tmp.exists():
                tmp.unlink()
            if attempt == 2:
                raise
            time.sleep(1 + attempt)


def archive_query(service, params, pagesize, out, force):
    query = json.dumps({'service': service, 'params': params,
                        'format': 'json', 'pagesize': pagesize}, separators=(',', ':'))
    url = 'https://mast.stsci.edu/api/v0/invoke?' + urlencode({'request': query})
    read_and_store(url, out, force=force)
    snapshot = json.loads(out.read_text())
    if snapshot.get('status') != 'COMPLETE':
        raise RuntimeError(f'MAST query failed: {snapshot.get("msg")}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', default='/tmp')
    parser.add_argument('--force', action='store_true', help='Re-download existing inputs')
    opts = parser.parse_args()
    root = Path(opts.output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    archive_query('Mast.Caom.Filtered', {
        'columns': 'obsid,obs_id,proposal_id,t_min,t_max,obs_collection,dataproduct_type,instrument_name,t_exptime,s_ra,s_dec',
        'filters': [{'paramName': 'proposal_id', 'values': ['3290']}]},
        300, root / 'janus_jwst_3290_observations.json', opts.force)
    archive_query('Mast.Caom.Cone', {'ra': 53.13127, 'dec': -27.78908, 'radius': 0.002},
                  2000, root / 'janus_jwst_cone.json', opts.force)
    for name, ra, dec, sr in ATNF_CONES:
        read_and_store('https://data.xao.ac.cn/pulsarcatalog/q/cone/scs.xml?' +
                       urlencode({'RA': ra, 'DEC': dec, 'SR': sr}), root / name, force=opts.force)
    for name, out, interval in PRODUCTS:
        read_and_store(BASE + quote(ROOT_URI + name, safe=''), root / out, interval, opts.force)
    print('READY: run JANUS_DATA_ROOT=' + str(root) +
          ' python3 tools/replay_jwst_mocss_mast_subset.py')


if __name__ == '__main__':
    main()
