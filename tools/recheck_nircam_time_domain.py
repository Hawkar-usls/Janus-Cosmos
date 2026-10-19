#!/usr/bin/env python3
"""Independent aperture time series from public NIRCam F150W i2d FITS byte ranges.

This intentionally does not infer physical source classes from a non-detection.
The fixed sky coordinates, aperture, background annulus, and random controls are
the same for all epochs. The FITS SCI bytes are retrieved directly from MAST.
"""

import argparse
import datetime as dt
import hashlib
import io
import json
import math
import random
import statistics
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from replay_jwst_mocss_mast_subset import RA, DEC, aperture_photometry

MAST = 'https://mast.stsci.edu/api/v0.1/Download/file?uri='
INPUTS = Path('/tmp/janus_mast_exact_cone_20261016.json')
DATA = Path('/tmp/janus_f150w_time_domain')
OUTPUT = Path('/tmp/janus_nircam_f150w_time_domain_20261016.json')


def fetch_byte_range(product, offset, end, destination):
    expected = end - offset + 1
    if destination.is_file() and destination.stat().st_size == expected:
        return
    url = MAST + quote('mast:JWST/product/' + product, safe='')
    headers = {'Range': f'bytes={offset}-{end}', 'User-Agent': 'JANUS-COSMOS-NIRCam-time-domain-v1'}
    request = Request(url, headers=headers)
    with urlopen(request, timeout=210) as response:
        content_range = response.headers.get('Content-Range', '')
        if response.status != 206 or not content_range.startswith(f'bytes {offset}-{end}/'):
            raise RuntimeError(f'Range not honored for {product}: {response.status}, {content_range}')
        payload = response.read()
    if len(payload) != expected:
        raise RuntimeError(f'Truncated range for {product}: {len(payload)} != {expected}')
    destination.write_bytes(payload)


def parse_headers(blob):
    fp = io.BytesIO(blob)
    first = fits.Header.fromfile(fp, endcard=True, padding=True)
    second = fits.Header.fromfile(fp, endcard=True, padding=True)
    if second.get('EXTNAME') != 'SCI' or second.get('BITPIX') != -32:
        raise RuntimeError('Unexpected FITS file header / SCI data layout')
    return first, second, fp.tell()


def recheck(row, nrows):
    name = row['obs_id'] + '_i2d.fits'
    stem = row['obs_id']
    header_file = DATA / f'{stem}_header.bin'
    fetch_byte_range(name, 0, 262143, header_file)
    first, sci, data_start = parse_headers(header_file.read_bytes())
    wcs = WCS(sci)
    cx, cy = wcs.all_world2pix([[RA, DEC]], 0)[0]
    width, height = sci['NAXIS1'], sci['NAXIS2']
    y0 = int(cy) - nrows // 2
    if not (30 < cx < width-30 and 0 < y0 and y0+nrows < height):
        return {'product': name, 'coverage': False, 'x': float(cx), 'y': float(cy)}
    byte_start = data_start + y0*width*4
    byte_end = byte_start + nrows*width*4-1
    data_file = DATA / f'{stem}_sci_{nrows}rows.bin'
    fetch_byte_range(name, byte_start, byte_end, data_file)
    image = np.memmap(data_file, dtype='>f4', mode='r', shape=(nrows, width))
    target = aperture_photometry(image, None, wcs, RA, DEC, y_origin=y0)
    if target is None:
        return {'product': name, 'coverage': False, 'x': float(cx), 'y': float(cy),
                'reason': 'target pixels masked/invalid or local background incomplete'}
    rand = random.Random(16004)
    controls = []
    nearest_coords = [(53.131165, -27.788930), (53.131584, -27.789105), (53.131301, -27.789367)]
    for _ in range(3500):
        dx, dy = rand.uniform(-6, 6), rand.uniform(-6, 6)
        if math.hypot(dx, dy) < 1.6:
            continue
        control_ra = RA + dx/(3600*math.cos(math.radians(DEC)))
        control_dec = DEC + dy/3600
        if any(math.hypot((control_ra-r)*3600*math.cos(math.radians(DEC)),(control_dec-d)*3600) < 0.7 for r,d in nearest_coords):
            continue
        phot = aperture_photometry(image,None,wcs,control_ra,control_dec,y_origin=y0)
        if phot is not None:
            controls.append(phot['aperture_njy'])
        if len(controls) >= 350:
            break
    if len(controls) < 100:
        return {'product': name, 'coverage': False, 'reason': 'insufficient nearby sky controls', 'n': len(controls)}
    med = statistics.median(controls)
    sigma = 1.4826 * statistics.median(abs(x-med) for x in controls)
    fixed_flux = target['aperture_njy']-med
    with data_file.open('rb') as fp:
        image_digest=hashlib.file_digest(fp,'sha256').hexdigest()
    return {
        'product': name, 'date_time_utc': (dt.datetime(1858,11,17)+dt.timedelta(days=row['t_min'])).isoformat(),
        'proposal_id': row['proposal_id'], 't_exptime_s': row['t_exptime'],
        'coverage': True, 'fixed_sky_ra_deg': RA, 'fixed_sky_dec_deg': DEC,
        'x_pixel': round(float(cx),3),'y_pixel': round(float(cy),3),
        'aperture_arcsec':0.12,'background_annulus_arcsec':[0.24,0.38],
        'aperture_pixels':target['aperture_pixels'],
        'flux_aperture_njy':round(float(fixed_flux),4),
        'control_noise_1sigma_njy':round(float(sigma),4),
        'fixed_aperture_snr':round(float(fixed_flux/sigma),3),
        'local_sky_control_n':len(controls),
        'primary_date_obs':first.get('DATE-OBS'),
        'header_sha256':hashlib.sha256(header_file.read_bytes()).hexdigest(),
        'sci_subset_sha256':image_digest,
        'sci_byte_range_inclusive':[byte_start,byte_end],
        'sci_subset_rows':nrows,
        'method_note':'Not an aperture-corrected flux; no PSF injection or exact correlated noise model. Nearby blank controls may include astrophysical sources.',
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--rows',type=int,default=501)
    parser.add_argument('--start',default='2022-01-01')
    parser.add_argument('--end',default='2024-12-31')
    opts=parser.parse_args()
    DATA.mkdir(parents=True, exist_ok=True)
    archive=json.loads(INPUTS.read_text())
    rows=[x for x in archive['data'] if x.get('obs_collection')=='JWST'
          and x.get('instrument_name')=='NIRCAM/IMAGE'
          and (x.get('filters') or '')=='F150W' and x.get('t_min')
          and x.get('obs_id','').startswith('jw') and '-o' in x.get('obs_id','')]
    rows.sort(key=lambda x:x['t_min'])
    rows=[x for x in rows if opts.start <= (dt.datetime(1858,11,17)+dt.timedelta(days=x['t_min'])).strftime('%Y-%m-%d') <= opts.end]
    results=[]
    for row in rows:
        print('CHECK',row['obs_id'],flush=True)
        try:
            out=recheck(row,opts.rows)
        except Exception as exc:
            out={'product':row['obs_id'],'error':str(exc)}
        print('RESULT',json.dumps(out,ensure_ascii=False),flush=True)
        results.append(out)
        OUTPUT.write_text(json.dumps({'source':'MAST CAOM exact 1 arcsec cone plus native JWST i2d SCI FITS byte ranges',
            'ra_deg':RA,'dec_deg':DEC,'results':results},ensure_ascii=False,indent=2)+'\n')
    print('RECEIPT',OUTPUT)


if __name__=='__main__':
    main()
