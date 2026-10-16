#!/usr/bin/env python3
"""Fetch all 16 official NRS1 or NRS2 RATE arrays from the GO-3290 manifest.

No analysis or user credentials are involved. Output is transient; science
receipts, not multi-GB raw FITS files, belong in the Git research branch.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
from urllib.parse import quote

from fetch_jwst_mocss_mast_subset import read_and_store, BASE, ROOT_URI


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', default='domains/space/MOCSS_16_EXPOSURE_MANIFEST_2026-10-16.json')
    parser.add_argument('--dest', default='/tmp')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--detector', choices=['nrs1', 'nrs2'], default='nrs1')
    args = parser.parse_args()
    records = json.loads(Path(args.manifest).read_text())['exposures']
    assert len(records) == 16
    downloads = []
    for record in records:
        filename = record['exposure_prefix'] + f'_{args.detector}_rate.fits'
        downloads.append((BASE + quote(ROOT_URI + filename, safe=''),
                          Path(args.dest) / filename))
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(read_and_store, url, path): path
                   for url, path in downloads}
        for future in as_completed(futures):
            future.result()
    print('PASS: 16/16 MAST', args.detector.upper(),
          'RATE inputs present; NO SPECTRAL ANALYSIS')


if __name__ == '__main__':
    main()
