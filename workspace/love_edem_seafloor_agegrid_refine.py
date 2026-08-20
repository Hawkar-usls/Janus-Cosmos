#!/usr/bin/env python3
import json
from pathlib import Path
from love_edem_seafloor_tectonic_provenance import sample_nc, LAT, LON

BASE='https://earthbyte.org/webdav/ftp/earthbyte/agegrid/2020/Grids/'
OUT=Path('data/love/JANUS-LOVE-EDEM-SEAFLOOR-AGEGRID-REFINEMENT-RUN-001-RECEIPT.json')

def s(name,stem):
    try: return sample_nc(BASE+name,stem)
    except Exception as e: return {'url':BASE+name,'resolved':False,'error':repr(e)}

def main():
    grids={
      'age_GeeK2007_2m':s('age.2020.1.GeeK2007.2m.nc','janus_age_geek_2m.nc'),
      'age_GTS2012_2m':s('age.2020.1.GTS2012.2m.nc','janus_age_gts_2m.nc'),
      'age_GTS2012_6m':s('age.2020.1.GTS2012.6m.nc','janus_age_gts_6m.nc'),
      'confidence_6m':s('conf.2020.1.GeeK2007.6m.nc','janus_conf_6m.nc'),
      'age_misfit_6m':s('age_misfit.2020.1.GeeK2007.6m.nc','janus_misfit_6m.nc'),
      'asymmetry_6m':s('asym.2020.1.GeeK2007.6m.nc','janus_asym_6m.nc'),
      'obliquity_6m':s('obliq.2020.1.GeeK2007.6m.nc','janus_obliq_6m.nc'),
      'full_rate_band_6m':s('full_rate_bands.2020.1.GeeK2007.6m.nc','janus_rateband_6m.nc')
    }
    a2=grids['age_GeeK2007_2m'].get('value'); g2=grids['age_GTS2012_2m'].get('value')
    receipt={
      'artifact_id':'JANUS-LOVE-EDEM-SEAFLOOR-AGEGRID-REFINEMENT-RUN-001-2026-08-20-v1.0',
      'schema':'janus.cosmos.love_edem.seafloor.agegrid_refinement.v1',
      'frozen_point':{'lat_deg':LAT,'lon_deg_east':LON,'no_recenter':True},
      'grids':grids,
      'derived':{
        'GeeK2007_2m_minus_prior_6m_ma': None if a2 is None else a2-96.09878540039062,
        'GTS2012_minus_GeeK2007_2m_ma': None if a2 is None or g2 is None else g2-a2,
        'inside_Cretaceous_Normal_Superchron_by_age': None if a2 is None else bool(83.0 <= a2 <= 120.6)
      },
      'confidence_boundary':'Seton et al. 2020 assigns lower confidence to crust within the Cretaceous Normal Superchron and to other poorly constrained settings. Confidence-grid numeric value is exported raw; interpretation must follow the published code/legend rather than guessing its scale.',
      'claim_ceiling':'AGE_GRID_REPLICATION_AND_PARAMETER_AUDIT_ONLY__NO_ARTIFACT_OR_RELAY_CLAIM'
    }
    receipt['status']='REFINEMENT_COMPLETE' if a2 is not None else 'REFINEMENT_PARTIAL'
    OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'age_GeeK2007_2m':a2,'age_GTS2012_2m':g2,'confidence':grids['confidence_6m'].get('value'),'misfit':grids['age_misfit_6m'].get('value'),'asym':grids['asymmetry_6m'].get('value'),'obliq':grids['obliquity_6m'].get('value')},indent=2))
if __name__=='__main__': main()
