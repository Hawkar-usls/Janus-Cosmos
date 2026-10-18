#!/usr/bin/env python3
"""Control-source and astrometric-offset stress tests for NIRCam time series.

Read only original archived FITS byte ranges downloaded by
recheck_nircam_time_domain.py. Test nearby JADES source positions as positive
controls and forced-flux variations within an allowed positional search.
"""

import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from recheck_nircam_time_domain import DATA, OUTPUT, parse_headers
from replay_jwst_mocss_mast_subset import RA, DEC, aperture_photometry

REPORT = Path('/tmp/janus_nircam_f150w_astrometry_controls_20261016.json')
NEIGHBORS = {
    'NW': {'ra':53.131165,'dec':-27.788930,'published_total_njy':4.75},
    'East': {'ra':53.131584,'dec':-27.789105,'published_total_njy':4.69},
    'South': {'ra':53.131301,'dec':-27.789367,'published_total_njy':3.57},
}


def main():
    receipt=json.loads(OUTPUT.read_text())
    entries=[x for x in receipt['results'] if x.get('coverage') and x.get('date_time_utc','').startswith('2023-10')]
    per_epoch=[]
    local_grid={}
    for x in entries:
        stem=x['product'].removesuffix('_i2d.fits')
        _,sci,offset=parse_headers((DATA/f'{stem}_header.bin').read_bytes())
        wcs=WCS(sci)
        # Use the byte-range offset difference exactly (SCI header is padded).
        y0=(x['sci_byte_range_inclusive'][0]-offset)//(4*sci['NAXIS1'])
        shape=(x['sci_subset_rows'],sci['NAXIS1'])
        image=np.memmap(DATA/f'{stem}_sci_{x["sci_subset_rows"]}rows.bin',dtype='>f4',mode='r',shape=shape)
        neighbors={}
        for key,point in NEIGHBORS.items():
            p=aperture_photometry(image,None,wcs,point['ra'],point['dec'],y_origin=y0)
            neighbors[key]=None if p is None else round(p['aperture_njy'],4)
        per_epoch.append({'date':x['date_time_utc'],'product':x['product'],
                          'target_flux_aperture_njy':x['flux_aperture_njy'],
                          'one_sigma_controls_njy':x['control_noise_1sigma_njy'],
                          'neighbor_apertures_njy':neighbors})
        for ox in range(-6,7):
            for oy in range(-6,7):
                dx,dy=ox*0.025,oy*0.025
                if math.hypot(dx,dy)>0.151:
                    continue
                key=f'{ox},{oy}'
                ra=RA+dx/(3600*math.cos(math.radians(DEC)))
                dec=DEC+dy/3600
                p=aperture_photometry(image,None,wcs,ra,dec,y_origin=y0)
                if p is None:
                    continue
                local_grid.setdefault(key,{'offset_east_arcsec':dx,'offset_north_arcsec':dy,'measurements':[]})['measurements'].append({'flux':p['aperture_njy'],'sigma':x['control_noise_1sigma_njy']})
    weighted_neighbors={}
    for key in NEIGHBORS:
        pts=[(e['neighbor_apertures_njy'][key],e['one_sigma_controls_njy']) for e in per_epoch if e['neighbor_apertures_njy'][key] is not None]
        weights=[1/s**2 for _,s in pts]
        weighted_neighbors[key]={'weighted_mean_aperture_njy':round(sum(v*w for (v,_),w in zip(pts,weights))/sum(weights),4),
                                  'nominal_sigma_njy':round(1/math.sqrt(sum(weights)),4),
                                  'published_total_flux_2022_njy':NEIGHBORS[key]['published_total_njy']}
    grid=[]
    for x in local_grid.values():
        pts=x['measurements']
        if len(pts)!=len(per_epoch): continue
        ww=[1/(p['sigma']**2) for p in pts]
        f=sum(p['flux']*w for p,w in zip(pts,ww))/sum(ww)
        ss=1/math.sqrt(sum(ww))
        grid.append({'offset_east_arcsec':x['offset_east_arcsec'],'offset_north_arcsec':x['offset_north_arcsec'],
                     'mean_njy':round(f,4),'nominal_sigma_njy':round(ss,4),'snr':round(f/ss,3)})
    grid.sort(key=lambda x:x['snr'],reverse=True)
    out={'science_status':'diagnostic astrometric and positive-control validation, NOT definitive deblended PSF photometry',
         'source_ra':RA,'source_dec':DEC,'epochs':per_epoch,
         'neighbor_positive_controls':weighted_neighbors,
         'offset_grid_spacing_arcsec':0.025,'offset_grid_radius_arcsec':0.15,
         'number_of_offset_trials':len(grid),
         'offset_peak_UNCORRECTED_for_trials':grid[:5],
         'at_zero_offset':next(x for x in grid if x['offset_east_arcsec']==0 and x['offset_north_arcsec']==0),
         'caveats':['Nearby galaxies may contaminate some sky annuli and are extended, not PSF standards.',
                    'Grid peak is a post-selection statistic and cannot be cited as an independent significant detection.',
                    'Per-image mosaic noise and mosaicking covariance have not been cross-validated by injection/recovery.']}
    REPORT.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(out,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
