#!/usr/bin/env python3
"""Crossmatch HUDF orphan position against public JADES DR4 NIRSpec and MAST.

Official archive inputs are not included in git: fetch the DR4 spectroscopic
catalogue from https://jades.herts.ac.uk/DR4/ and query Mast.Caom.Cone.
This is a target-selection catalogue comparison, not a spectral re-extraction.
"""

import datetime as dt
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits

RA, DEC = 53.13127, -27.78908
CATALOG = Path('/tmp/janus_jades_dr4_spectral_catalog_20261016.fits')
MAST = Path('/tmp/janus_mast_exact_cone_20261016.json')
OUTPUT = Path('/tmp/janus_jades_dr4_mast_coverage_20261016.json')


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    with fits.open(CATALOG,memmap=True) as hdul:
        x=hdul['Obs_info'].data
        sep=np.hypot((x['RA_TARG']-RA)*math.cos(math.radians(DEC)),x['Dec_TARG']-DEC)*3600
        nearest=[]
        for i in np.argsort(sep)[:12]:
            nearest.append({'separation_arcsec':round(float(sep[i]),4),
                            'ra_deg':float(x['RA_TARG'][i]),'dec_deg':float(x['Dec_TARG'][i]),
                            'unique_id':str(x['Unique_ID'][i]),
                            'z_spec':float(x['z_Spec'][i]),'z_phot':float(x['z_phot'][i]),
                            'z_spec_flag':str(x['z_Spec_flag'][i]),'pid':str(x['PID'][i])})
        targets=len(x)
        n_within_5=int(np.count_nonzero(sep<=5))
        n_within_2=int(np.count_nonzero(sep<=2))
    caom=json.loads(MAST.read_text())
    rows=caom['data']
    jwst_spectra=[r for r in rows if r.get('obs_collection')=='JWST' and r.get('dataproduct_type')=='spectrum'
                  and r.get('t_min') and (r.get('em_min') or 9999)<1600.4<(r.get('em_max') or -1)]
    groups={}
    for r in jwst_spectra:
        key=(str(r['proposal_id']),r['instrument_name'],r['filters'])
        groups.setdefault(key,[]).append(r)
    relevant=[{'proposal':k[0],'instrument':k[1],'mode':k[2],
               'number_caom_rows':len(v),
               'earliest_utc':(dt.datetime(1858,11,17)+dt.timedelta(days=min(z['t_min'] for z in v))).isoformat(timespec='minutes'),
               'footprint_only_not_proof_of_slit_overlap':True} for k,v in sorted(groups.items())]
    out={
      'source_ra_deg':RA,'source_dec_deg':DEC,
      'jades_dr4':{'url':'https://jades.herts.ac.uk/DR4/Combined_DR4_external_v1.2.1.fits',
                   'sha256':sha(CATALOG),'n_target_rows':targets,
                   'n_target_centers_within_2arcsec':n_within_2,
                   'n_target_centers_within_5arcsec':n_within_5,
                   'closest_12':nearest},
      'mast_cone':{'service':'Mast.Caom.Cone','radius_arcsec':1.008,
                    'all_mission_rows':len(rows),'jwst_rows':sum(r.get('obs_collection')=='JWST' for r in rows),
                    'raw_query_sha256':sha(MAST),
                    'jwst_archived_spectrum_caom_rows_covering_1600nm':len(jwst_spectra),
                    'individual_programs':relevant,
                    'source_point_is_not_necessarily_inside_any_MSA_slit':True,
                    'MOCSS_go3290_rate_images_indexed_as_images_not_spectra':True},
      'science_conclusion':'No known JADES DR4 NIRSpec targeted object within 5 arcsec; do NOT attribute any catalog redshift to the orphan. Other footprint-covering NIRSpec PRISM rows do not establish an open slit at this coordinate.',
      'gate':'CATALOG_CROSSMATCH_ONLY__NO_NEW_ORPHAN_SPECTRUM_OR_SPECTROSCOPIC_DISTANCE'
    }
    OUTPUT.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(out,ensure_ascii=False,indent=2))


if __name__=='__main__': main()
