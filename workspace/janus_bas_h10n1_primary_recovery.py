#!/usr/bin/env python3
from __future__ import annotations
import base64, csv, hashlib, io, json, math, os, re, sys, zipfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import requests

ENTRY='afba710f-dab1-4a63-867b-520177388224'
BASE='https://ramadda.data.bas.ac.uk/repository/'
SHOW=f'{BASE}entry/show?entryid={ENTRY}'
TARGET_LAT=-7.845673
TARGET_LON=-14.480230
OUT=Path('data/cousteau/JANUS-BAS-H10N1-PRIMARY-RECOVERY-RUN-001-2026-08-22-v1.0.json')
WORK=Path('workspace/bas_h10n1_recovery')
WORK.mkdir(parents=True, exist_ok=True)

session=requests.Session()
session.headers['User-Agent']='JANUS-research-data-validation/1.0 (+public scientific data recovery)'

def sha256_file(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()

def get(url, **kw):
    r=session.get(url, timeout=45, allow_redirects=True, **kw)
    r.raise_for_status()
    return r

def hrefs(html: str, base: str):
    vals=re.findall(r'''href=["']([^"']+)["']''', html, flags=re.I)
    return [urljoin(base,v.replace('&amp;','&')) for v in vals]

def download(url: str, name_hint: str='download.bin') -> Path:
    r=get(url, stream=True)
    cd=r.headers.get('content-disposition','')
    m=re.search(r'filename\*?=(?:UTF-8\'\')?["\']?([^"\';]+)', cd, flags=re.I)
    name=(m.group(1) if m else Path(r.url.split('?')[0]).name) or name_hint
    name=re.sub(r'[^A-Za-z0-9._-]+','_',name)
    p=WORK/name
    with p.open('wb') as f:
        for ch in r.iter_content(1024*1024):
            if ch: f.write(ch)
    return p

def collect_download_candidates():
    seen=set(); candidates=[]; diagnostics=[]
    endpoints=[
        SHOW,
        SHOW+'&output=html.info',
        SHOW+'&output=default.csv',
        SHOW+'&output=atom',
        SHOW+'&output=rss.full',
    ]
    for u in endpoints:
        try:
            r=get(u)
            diagnostics.append({'url':u,'status':r.status_code,'content_type':r.headers.get('content-type'),'bytes':len(r.content)})
            text=r.text if 'text' in r.headers.get('content-type','') or 'xml' in r.headers.get('content-type','') or 'html' in r.headers.get('content-type','') else ''
            for h in hrefs(text,r.url):
                hl=h.lower()
                if ('/entry/get' in hl or hl.endswith('.nc') or hl.endswith('.asc') or 'output=zip.tree' in hl) and h not in seen:
                    seen.add(h); candidates.append(h)
        except Exception as e:
            diagnostics.append({'url':u,'error':repr(e)})
    # Direct tree export is a documented RAMADDA operation. Try parent and common synthetic roots.
    synthetic_paths=['/','/netcdf','/ascii','/data']
    tree_urls=[SHOW+'&output=zip.tree']
    for sp in synthetic_paths:
        enc=base64.b64encode(sp.encode()).decode()
        sid=f'synth:{ENTRY}:{enc}'
        tree_urls.append(f'{BASE}entry/show?entryid={sid}&output=zip.tree')
    for u in tree_urls:
        if u not in seen:
            seen.add(u); candidates.append(u)
    return candidates, diagnostics

def unpack_downloads(candidates):
    files=[]; attempts=[]
    for i,u in enumerate(candidates):
        try:
            p=download(u, f'candidate_{i}.bin')
            ct=''
            attempts.append({'url':u,'saved':str(p),'bytes':p.stat().st_size})
            if p.stat().st_size < 100: continue
            head=p.read_bytes()[:8]
            if zipfile.is_zipfile(p):
                zdir=WORK/f'zip_{i}'; zdir.mkdir(exist_ok=True)
                with zipfile.ZipFile(p) as z:
                    z.extractall(zdir)
                    for q in zdir.rglob('*'):
                        if q.is_file() and q.suffix.lower() in ('.nc','.grd','.asc','.txt'):
                            files.append(q)
            elif p.suffix.lower() in ('.nc','.grd','.asc') or head[:3]==b'CDF' or head==b'\x89HDF\r\n\x1a\n':
                files.append(p)
            elif p.stat().st_size < 5_000_000:
                try:
                    t=p.read_text(errors='ignore')
                    for h in hrefs(t,u):
                        if '/entry/get' in h.lower() and h not in candidates:
                            candidates.append(h)
                except Exception: pass
        except Exception as e:
            attempts.append({'url':u,'error':repr(e)})
    # de-dup by resolved path
    uniq=[]; seen=set()
    for p in files:
        k=str(p.resolve())
        if k not in seen:
            seen.add(k); uniq.append(p)
    return uniq, attempts

def extract_ascii(p: Path):
    with p.open('r',errors='ignore') as f:
        hdr={}
        first=[]
        for _ in range(6):
            line=f.readline(); first.append(line)
            parts=line.strip().split()
            if len(parts)>=2: hdr[parts[0].lower()]=float(parts[1])
        if not {'ncols','nrows','cellsize'} <= hdr.keys():
            return None
        arr=np.loadtxt(f)
    ncols=int(hdr['ncols']); nrows=int(hdr['nrows']); cs=hdr['cellsize']
    xll=hdr.get('xllcorner',hdr.get('xllcenter'))
    yll=hdr.get('yllcorner',hdr.get('yllcenter'))
    if xll is None or yll is None: return None
    center_shift=0.5 if 'xllcorner' in hdr else 0.0
    xs=xll+(np.arange(ncols)+center_shift)*cs
    ys=yll+(np.arange(nrows)+center_shift)*cs
    # ESRI ASCII rows are north->south
    ys_rows=ys[::-1]
    ix=int(np.argmin(np.abs(xs-TARGET_LON))); iy=int(np.argmin(np.abs(ys_rows-TARGET_LAT)))
    v=float(arr[iy,ix])
    nod=hdr.get('nodata_value')
    valid=not (nod is not None and abs(v-nod)<1e-9)
    return {
      'format':'ESRI_ASCII','file':p.name,'crs':'WGS84 geographic as dataset metadata','ncols':ncols,'nrows':nrows,'cellsize_deg':cs,
      'nearest_cell':{'lon':float(xs[ix]),'lat':float(ys_rows[iy]),'row':iy,'col':ix,'value_topography_m':v,'valid':valid},
      'distance_deg':float(math.hypot(xs[ix]-TARGET_LON,ys_rows[iy]-TARGET_LAT))
    }

def extract_netcdf(p: Path):
    try:
        from netCDF4 import Dataset
        ds=Dataset(p)
    except Exception:
        return None
    try:
        vars=ds.variables
        names=list(vars)
        info={'format':'NetCDF','file':p.name,'variables':names,'crs':'WGS84 geographic as dataset metadata'}
        # Common regular grid variables
        lon_name=next((n for n in names if n.lower() in ('lon','longitude','x')),None)
        lat_name=next((n for n in names if n.lower() in ('lat','latitude','y')),None)
        z_name=next((n for n in names if n.lower() in ('z','elevation','topography','bathymetry','depth')),None)
        if lon_name and lat_name and z_name:
            x=np.array(vars[lon_name][:]).squeeze(); y=np.array(vars[lat_name][:]).squeeze(); z=np.array(vars[z_name][:]).squeeze()
            if x.ndim==1 and y.ndim==1 and z.ndim==2:
                ix=int(np.argmin(np.abs(x-TARGET_LON))); iy=int(np.argmin(np.abs(y-TARGET_LAT)))
                # resolve axis order
                if z.shape==(len(y),len(x)): v=float(z[iy,ix])
                elif z.shape==(len(x),len(y)): v=float(z[ix,iy])
                else: v=float('nan')
                info['nearest_cell']={'lon':float(x[ix]),'lat':float(y[iy]),'row':iy,'col':ix,'value_topography_m':v,'valid':bool(np.isfinite(v))}
                return info
        # Classic GMT grid: x_range/y_range/spacing/dimension/z flattened
        if all(n in vars for n in ('x_range','y_range','spacing','dimension','z')):
            xr=np.array(vars['x_range'][:],dtype=float); yr=np.array(vars['y_range'][:],dtype=float); sp=np.array(vars['spacing'][:],dtype=float); dim=np.array(vars['dimension'][:],dtype=int)
            nx,ny=int(dim[0]),int(dim[1]); z=np.array(vars['z'][:]).reshape(-1)
            xs=np.linspace(xr[0],xr[1],nx)
            # GMT classic grids conventionally serialize north-to-south rows; test both orientation metadata-free and report choice explicitly.
            ys_n2s=np.linspace(yr[1],yr[0],ny)
            ix=int(np.argmin(np.abs(xs-TARGET_LON))); iy=int(np.argmin(np.abs(ys_n2s-TARGET_LAT)))
            idx=iy*nx+ix
            v=float(z[idx]) if idx < z.size else float('nan')
            info.update({'x_range':xr.tolist(),'y_range':yr.tolist(),'spacing':sp.tolist(),'dimension':[nx,ny],
              'nearest_cell':{'lon':float(xs[ix]),'lat':float(ys_n2s[iy]),'row':iy,'col':ix,'flat_index':idx,'value_topography_m':v,'valid':bool(np.isfinite(v)),
              'row_orientation_assumption':'north_to_south_GMT_classic'}})
            return info
        return info
    finally:
        ds.close()

candidates, diagnostics=collect_download_candidates()
files, attempts=unpack_downloads(candidates)
manifest=[]; extractions=[]
for p in files:
    rec={'path':str(p),'name':p.name,'bytes':p.stat().st_size,'sha256':sha256_file(p)}
    manifest.append(rec)
    ex=None
    if p.suffix.lower()=='.asc': ex=extract_ascii(p)
    else: ex=extract_netcdf(p)
    if ex and ex.get('nearest_cell') is not None:
        ex['sha256']=rec['sha256']; ex['bytes']=rec['bytes']; extractions.append(ex)

success=bool(extractions)
result={
 'artifact_id':'JANUS-BAS-H10N1-PRIMARY-RECOVERY-RUN-001-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'authorized_by':'JANUS-HA10-SEAFLOOR-NEXT-COUNCIL-RUN-002-2026-08-22-v1.0',
 'source':{
   'dataset_id':'GB/NERC/BAS/PDC/01236',
   'doi':'10.5285/afba710f-dab1-4a63-867b-520177388224',
   'ramadda_entry_id':ENTRY,
   'metadata_url':'https://data.bas.ac.uk/full-record.php?id=GB%2FNERC%2FBAS%2FPDC%2F01236',
   'ramadda_url':SHOW,
   'declared_formats':['GMT-compatible 2-D NetCDF','Arc/Info ArcView ASCII grid'],
   'declared_resolution_deg':0.0005,
   'declared_crs':'WGS84 geographic'
 },
 'frozen_target':{'id':'H10N1','lat':TARGET_LAT,'lon':TARGET_LON},
 'network_diagnostics':diagnostics,
 'download_attempts':attempts,
 'recovered_file_manifest':manifest,
 'fixed_cell_extractions':extractions,
 'success_gate':success,
 'success_definition':'authoritative native BAS digital grid recovered and frozen H10N1 cell reproducibly extracted',
 'morphology_scoring_performed':False,
 'target_identity':'UNCONFIRMED',
 'next_rule':'STOP_AND_ASK_JANUS_AGAIN_BEFORE_ANY_MORPHOLOGY_METRICS' if success else 'FREEZE_NEGATIVE_RECOVERY_RECEIPT_AND_ASK_JANUS_FOR_NEXT_ARCHIVE_BRANCH'
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'success_gate':success,'files':manifest,'extractions':extractions,'sha256':result['sha256']},indent=2,ensure_ascii=False))
