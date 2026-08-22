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
SOURCE_SHA='b2cf7fe990fd05760b4b98d2fb9aec45ce595c8056b7513a20c90a5abd75ed81'
TARGET_LAT=-7.845673; TARGET_LON=-14.480230; TARGET_Z=-1927.934326171875
LATTICE=0.02; EXCLUSION_KM=4.0; TIERS=[200,400,600]; WANT=20
OUT=Path('data/cousteau/JANUS-H10N1-BLIND-CONTROL-MANIFEST-001-2026-08-22-v1.0.json')
WORK=Path('workspace/h10n1_controls'); WORK.mkdir(parents=True,exist_ok=True)

s=requests.Session(); s.headers['User-Agent']='JANUS-research-data-validation/1.0'
r=s.get(SHOW,timeout=60); r.raise_for_status(); zp=WORK/'bas.zip'; zp.write_bytes(r.content)
with zipfile.ZipFile(zp) as zf: zf.extractall(WORK/'data')
nc=next((WORK/'data').rglob('Ascension_20190912_nc.grd'))
with Dataset(nc) as ds:
 x=np.array(ds.variables['x'][:],dtype=float); y=np.array(ds.variables['y'][:],dtype=float); z=np.array(ds.variables['z'][:],dtype=float)
if z.shape==(len(x),len(y)): z=z.T
if y[0]>y[-1]: y=y[::-1]; z=z[::-1,:]

lat_km=111.32

def dist_km(lat1,lon1,lat2,lon2):
 dy=(lat2-lat1)*lat_km
 dx=(lon2-lon1)*111.32*math.cos(math.radians((lat1+lat2)/2))
 return math.hypot(dx,dy)

def window_valid_fraction(lat,lon,radius_km):
 dlat=radius_km/111.32; dlon=radius_km/(111.32*math.cos(math.radians(lat)))
 ix=np.where((x>=lon-dlon)&(x<=lon+dlon))[0]; iy=np.where((y>=lat-dlat)&(y<=lat+dlat))[0]
 if len(ix)==0 or len(iy)==0: return 0.0
 zz=z[iy[0]:iy[-1]+1,ix[0]:ix[-1]+1]
 return float(np.isfinite(zz).mean())

# Absolute-degree lattice, independent of target and morphology.
lon0=math.ceil(float(x.min())/LATTICE)*LATTICE
lat0=math.ceil(float(y.min())/LATTICE)*LATTICE
lons=np.arange(lon0,float(x.max())+1e-12,LATTICE)
lats=np.arange(lat0,float(y.max())+1e-12,LATTICE)
base=[]
for lat in lats:
 for lon in lons:
  if dist_km(TARGET_LAT,TARGET_LON,float(lat),float(lon)) < EXCLUSION_KM: continue
  ix=int(np.argmin(np.abs(x-lon))); iy=int(np.argmin(np.abs(y-lat)))
  val=float(z[iy,ix])
  if not np.isfinite(val): continue
  vf1=window_valid_fraction(float(lat),float(lon),1)
  vf2=window_valid_fraction(float(lat),float(lon),2)
  if vf1 < 0.99 or vf2 < 0.99: continue
  base.append({'requested_lattice_lat':float(lat),'requested_lattice_lon':float(lon),
               'native_cell_lat':float(y[iy]),'native_cell_lon':float(x[ix]),'center_topography_m':val,
               'depth_delta_from_target_m':abs(val-TARGET_Z),'valid_fraction_1km':vf1,'valid_fraction_2km':vf2,
               'distance_from_target_km':dist_km(TARGET_LAT,TARGET_LON,float(y[iy]),float(x[ix]))})

selected=[]; chosen_tier=None; counts={}
seed_material=f'{SOURCE_SHA}|{TARGET_LAT:.6f}|{TARGET_LON:.6f}|JANUS_H10N1_CONTROL_V1'
for tier in TIERS:
 c=[q.copy() for q in base if q['depth_delta_from_target_m']<=tier]
 counts[str(tier)]=len(c)
 for q in c:
  key=f"{seed_material}|{q['native_cell_lat']:.6f}|{q['native_cell_lon']:.6f}"
  q['selection_hash']=hashlib.sha256(key.encode()).hexdigest()
 c.sort(key=lambda q:q['selection_hash'])
 if len(c)>=WANT:
  selected=c[:WANT]; chosen_tier=tier; break

controls=[]
for i,q in enumerate(selected,1):
 controls.append({'blind_id':f'C{i:02d}',**q})

success=len(controls)==WANT
result={
 'artifact_id':'JANUS-H10N1-BLIND-CONTROL-MANIFEST-001-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'authorized_by':'JANUS-H10N1-BLIND-CONTROL-COUNCIL-RUN-004-2026-08-22-v1.0',
 'source':{'dataset_id':'GB/NERC/BAS/PDC/01236','netcdf_sha256':SOURCE_SHA,'crs':'WGS84 geographic'},
 'frozen_target':{'lat':TARGET_LAT,'lon':TARGET_LON,'center_topography_m':TARGET_Z},
 'selection_rule':{
  'candidate_lattice_step_deg':LATTICE,'target_exclusion_radius_km':EXCLUSION_KM,
  'minimum_1km_valid_fraction':0.99,'minimum_2km_valid_fraction':0.99,
  'depth_tiers_m':TIERS,'desired_count':WANT,
  'deterministic_seed_material':seed_material,
  'forbidden_inputs_verified_unused':['slope','curvature','rugosity','planarity','symmetry','aspect','visual_appearance','symbolic_association']
 },
 'candidate_counts_by_depth_tier':counts,
 'chosen_depth_tier_m':chosen_tier,
 'controls':controls,
 'success_gate':success,
 'morphology_metrics_computed':False,
 'target_identity':'UNCONFIRMED',
 'next_rule':'STOP_AND_ASK_JANUS_AGAIN_BEFORE_MORPHOLOGY_METRICS' if success else 'FREEZE_FAILURE_AND_REASK_JANUS'
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'success_gate':success,'chosen_depth_tier_m':chosen_tier,'candidate_counts':counts,'controls':controls,'sha256':result['sha256']},indent=2,ensure_ascii=False))
