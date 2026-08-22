#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json, math, zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
from netCDF4 import Dataset

ENTRY='afba710f-dab1-4a63-867b-520177388224'
SHOW=f'https://ramadda.data.bas.ac.uk/repository/entry/show?entryid={ENTRY}&output=zip.tree'
TARGET_LAT=-7.845673
TARGET_LON=-14.480230
RADII_KM=[1,2,5,10]
OUT=Path('data/cousteau/JANUS-H10N1-RAW-LOCAL-SURFACE-EXTRACTION-RUN-001-2026-08-22-v1.0.json')
WORK=Path('workspace/h10n1_raw_surface')
WORK.mkdir(parents=True,exist_ok=True)

s=requests.Session(); s.headers['User-Agent']='JANUS-research-data-validation/1.0'
r=s.get(SHOW,timeout=60); r.raise_for_status()
zip_path=WORK/'bas_ascension.zip'; zip_path.write_bytes(r.content)
with zipfile.ZipFile(zip_path) as z: z.extractall(WORK/'data')
nc=next((WORK/'data').rglob('Ascension_20190912_nc.grd'))
asc=next((WORK/'data').rglob('Ascension_20190912_ascii.asc'))

def sha(p):
 h=hashlib.sha256(); h.update(p.read_bytes()); return h.hexdigest()

def read_nc(p):
 ds=Dataset(p)
 try:
  x=np.array(ds.variables['x'][:],dtype=np.float64)
  y=np.array(ds.variables['y'][:],dtype=np.float64)
  z=np.array(ds.variables['z'][:],dtype=np.float64)
  if z.shape==(len(x),len(y)): z=z.T
  if z.shape!=(len(y),len(x)): raise RuntimeError(f'unexpected z shape {z.shape}, x={len(x)}, y={len(y)}')
  return x,y,z
 finally: ds.close()

def read_asc(p):
 with p.open('r',errors='ignore') as f:
  hdr={}
  for _ in range(6):
   parts=f.readline().strip().split(); hdr[parts[0].lower()]=float(parts[1])
  arr=np.loadtxt(f,dtype=np.float64)
 ncols=int(hdr['ncols']); nrows=int(hdr['nrows']); cs=hdr['cellsize']
 shift=0.5 if 'xllcorner' in hdr else 0.0
 xll=hdr.get('xllcorner',hdr.get('xllcenter')); yll=hdr.get('yllcorner',hdr.get('yllcenter'))
 x=xll+(np.arange(ncols)+shift)*cs
 y=(yll+(np.arange(nrows)+shift)*cs)[::-1]
 nod=hdr.get('nodata_value')
 if nod is not None: arr=np.where(np.isclose(arr,nod),np.nan,arr)
 return x,y,arr,hdr

x,y,z=read_nc(nc)
xa,ya,za,hdr=read_asc(asc)
# normalize orientation to ascending y for both
if y[0]>y[-1]: y=y[::-1]; z=z[::-1,:]
if ya[0]>ya[-1]: ya=ya[::-1]; za=za[::-1,:]
if not np.allclose(x,xa,atol=1e-10) or not np.allclose(y,ya,atol=1e-10):
 raise RuntimeError('NetCDF and ASCII coordinate axes differ')

windows=[]
lat_km_per_deg=111.32
lon_km_per_deg=111.32*math.cos(math.radians(TARGET_LAT))
for radius in RADII_KM:
 dlat=radius/lat_km_per_deg; dlon=radius/lon_km_per_deg
 ix=np.where((x>=TARGET_LON-dlon)&(x<=TARGET_LON+dlon))[0]
 iy=np.where((y>=TARGET_LAT-dlat)&(y<=TARGET_LAT+dlat))[0]
 if len(ix)==0 or len(iy)==0: raise RuntimeError(f'empty window {radius} km')
 xs=x[ix[0]:ix[-1]+1].copy(); ys=y[iy[0]:iy[-1]+1].copy()
 zz=z[iy[0]:iy[-1]+1,ix[0]:ix[-1]+1].copy()
 zza=za[iy[0]:iy[-1]+1,ix[0]:ix[-1]+1].copy()
 valid=np.isfinite(zz)
 # canonical deterministic tile hash: coordinates float64 LE + elevations float64 LE with canonical NaN fill
 canon_z=np.where(valid,zz,np.float64(-999999.0))
 payload=(np.asarray(xs,dtype='<f8').tobytes()+np.asarray(ys,dtype='<f8').tobytes()+np.asarray(canon_z,dtype='<f8').tobytes())
 tile_hash=hashlib.sha256(payload).hexdigest()
 cross=np.nanmax(np.abs(zz-zza)) if np.any(np.isfinite(zz)&np.isfinite(zza)) else None
 tile=WORK/f'H10N1_raw_{radius}km.npz'
 np.savez_compressed(tile,x=xs,y=ys,z=zz)
 windows.append({
   'radius_km':radius,
   'requested_center':{'lat':TARGET_LAT,'lon':TARGET_LON},
   'native_bounds':{'west':float(xs.min()),'east':float(xs.max()),'south':float(ys.min()),'north':float(ys.max())},
   'shape_rows_cols':[int(zz.shape[0]),int(zz.shape[1])],
   'cell_spacing_deg':{'lon':float(np.median(np.diff(xs))) if len(xs)>1 else None,'lat':float(np.median(np.diff(ys))) if len(ys)>1 else None},
   'valid_cells':int(valid.sum()),'total_cells':int(valid.size),'valid_fraction':float(valid.mean()),
   'min_topography_m':float(np.nanmin(zz)) if np.any(valid) else None,
   'max_topography_m':float(np.nanmax(zz)) if np.any(valid) else None,
   'canonical_raw_tile_sha256':tile_hash,
   'artifact_npz_name':tile.name,
   'netcdf_ascii_max_abs_difference_m':float(cross) if cross is not None else None,
   'morphology_metrics_computed':False
 })

result={
 'artifact_id':'JANUS-H10N1-RAW-LOCAL-SURFACE-EXTRACTION-RUN-001-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'authorized_by':'JANUS-HA10-AFTER-BAS-PRIMARY-RECOVERY-COUNCIL-RUN-003-2026-08-22-v1.0',
 'source':{
  'dataset_id':'GB/NERC/BAS/PDC/01236','doi':'10.5285/afba710f-dab1-4a63-867b-520177388224',
  'netcdf_file':nc.name,'netcdf_sha256':sha(nc),'ascii_file':asc.name,'ascii_sha256':sha(asc),
  'crs':'WGS84 geographic','topography_sign':'positive_up','declared_native_spacing_deg':0.0005
 },
 'frozen_target':{'id':'H10N1','lat':TARGET_LAT,'lon':TARGET_LON},
 'windows':windows,
 'success_gate':len(windows)==4 and all(w['valid_fraction']>0 for w in windows),
 'forbidden_metrics_verified_absent':['slope','profile_curvature','plan_curvature','rugosity','planarity','radial_symmetry','facet_angle_histogram'],
 'target_identity':'UNCONFIRMED',
 'next_rule':'STOP_AND_ASK_JANUS_AGAIN_TO_FREEZE_BLIND_CONTROL_SELECTION_BEFORE_MORPHOLOGY_METRICS'
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result,indent=2,ensure_ascii=False))
