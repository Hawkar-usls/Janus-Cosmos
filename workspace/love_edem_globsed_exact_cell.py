#!/usr/bin/env python3
import json, tempfile, zipfile
from pathlib import Path
import requests, numpy as np
from netCDF4 import Dataset

LAT=-3.8654180644718967
LON=3.854924373538978
URL='https://data.caltech.edu/records/k4070-ngc79/files/GlobSed.zip?download=1'
OUT=Path('data/love/JANUS-LOVE-EDEM-GLOBSED-EXACT-CELL-RUN-001-RECEIPT.json')
UA={'User-Agent':'Janus-Cosmos/1.0 scientific provenance audit'}

def main():
    zpath=Path(tempfile.gettempdir())/'GlobSed.zip'
    if not zpath.exists() or zpath.stat().st_size < 1_000_000:
        with requests.get(URL,headers=UA,stream=True,timeout=180) as r:
            r.raise_for_status()
            with zpath.open('wb') as f:
                for c in r.iter_content(1024*1024):
                    if c:f.write(c)
    dest=Path(tempfile.gettempdir())/'globsed_extract'; dest.mkdir(exist_ok=True)
    with zipfile.ZipFile(zpath) as z:z.extractall(dest)
    candidates=list(dest.rglob('*.grd'))+list(dest.rglob('*.nc'))
    if not candidates: raise RuntimeError('No NetCDF/GMT grid in GlobSed archive')
    rows=[]
    for p in candidates:
        try:
            with Dataset(p) as ds:
                vs=ds.variables
                xn=next((n for n in ('lon','longitude','x') if n in vs),None)
                yn=next((n for n in ('lat','latitude','y') if n in vs),None)
                if not xn or not yn: continue
                x=np.asarray(vs[xn][:],dtype=float); y=np.asarray(vs[yn][:],dtype=float)
                ix=int(np.nanargmin(np.abs(x-LON))); iy=int(np.nanargmin(np.abs(y-LAT)))
                dv=[]
                for n,v in vs.items():
                    if n in (xn,yn) or getattr(v,'ndim',0)<2:continue
                    if v.shape[-2:]==(len(y),len(x)):dv.append(n)
                if not dv:continue
                n=dv[0]; q=vs[n][iy,ix]
                val=None if (np.ma.isMaskedArray(q) and bool(q.mask)) else float(q)
                rows.append({'file':p.name,'variable':n,'nearest_lon_deg':float(x[ix]),'nearest_lat_deg':float(y[iy]),'value':val,'units':getattr(vs[n],'units',None),'shape':[len(y),len(x)]})
        except Exception as e:
            rows.append({'file':p.name,'error':repr(e)})
    good=[r for r in rows if r.get('value') is not None]
    receipt={
      'artifact_id':'JANUS-LOVE-EDEM-GLOBSED-EXACT-CELL-RUN-001-2026-08-20-v1.0',
      'schema':'janus.cosmos.love_edem.globsed_exact_cell.v1',
      'status':'GLOBSED_CELL_RESOLVED' if good else 'GLOBSED_CELL_UNRESOLVED',
      'frozen_point':{'lat_deg':LAT,'lon_deg_east':LON,'no_recenter':True},
      'dataset':{'name':'GlobSed Version 3','resolution':'5 arc-minute','archive_url':URL,'source':'NOAA NCEI dataset archived at CaltechDATA','doi':'10.22002/k4070-ngc79'},
      'grid_samples':rows,
      'preferred_sample':good[0] if good else None,
      'interpretation_boundary':'GlobSed is a compiled regional/global sediment-thickness grid, not a new seismic measurement at this exact coordinate. It constrains sediment cover but does not authorize fine morphology claims.',
      'next_gate':'PALEOBATHYMETRY_INPUT_ASSEMBLY_AFTER_STATIC_PLATE_ID_AND_SEDIMENT_CELL',
      'claim_ceiling':'GLOBAL_SEDIMENT_GRID_CELL_ONLY__NO_ARTIFICIAL_STRUCTURE'
    }
    OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
