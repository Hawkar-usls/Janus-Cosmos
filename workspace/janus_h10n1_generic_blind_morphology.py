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
CTRL_PATH=Path('data/cousteau/JANUS-H10N1-BLIND-CONTROL-MANIFEST-001-2026-08-22-v1.0.json')
COUNCIL_PATH=Path('data/cousteau/JANUS-H10N1-MORPHOLOGY-METRICS-COUNCIL-RUN-005-2026-08-22-v1.0.json')
OUT=Path('data/cousteau/JANUS-H10N1-GENERIC-BLIND-MORPHOLOGY-RUN-001-2026-08-22-v1.0.json')
WORK=Path('workspace/h10n1_morphology'); WORK.mkdir(parents=True,exist_ok=True)
SCALES=[2,1]
METRICS=['local_relief_m','slope_median_deg','slope_p95_deg','rugosity_surface_ratio','planarity_rmse_over_relief','profile_curvature_rms','plan_curvature_rms','radial_rotation_similarity','facet_aspect_entropy','facet_dominant_bin_fraction']

controls_manifest=json.loads(CTRL_PATH.read_text(encoding='utf-8'))
council=json.loads(COUNCIL_PATH.read_text(encoding='utf-8'))
TARGET={'id':'TARGET','lat':controls_manifest['frozen_target']['lat'],'lon':controls_manifest['frozen_target']['lon']}
CONTROLS=[{'id':q['blind_id'],'lat':q['native_cell_lat'],'lon':q['native_cell_lon']} for q in controls_manifest['controls']]

s=requests.Session(); s.headers['User-Agent']='JANUS-research-data-validation/1.0'
r=s.get(SHOW,timeout=60); r.raise_for_status(); zp=WORK/'bas.zip'; zp.write_bytes(r.content)
with zipfile.ZipFile(zp) as zf: zf.extractall(WORK/'data')
nc=next((WORK/'data').rglob('Ascension_20190912_nc.grd'))
with Dataset(nc) as ds:
 x=np.array(ds.variables['x'][:],dtype=float); y=np.array(ds.variables['y'][:],dtype=float); z=np.array(ds.variables['z'][:],dtype=float)
if z.shape==(len(x),len(y)): z=z.T
if y[0]>y[-1]: y=y[::-1]; z=z[::-1,:]

SCRIPT_SHA=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

def extract(lat,lon,radius):
 dlat=radius/111.32; dlon=radius/(111.32*math.cos(math.radians(lat)))
 ix=np.where((x>=lon-dlon)&(x<=lon+dlon))[0]; iy=np.where((y>=lat-dlat)&(y<=lat+dlat))[0]
 xs=x[ix[0]:ix[-1]+1]; ys=y[iy[0]:iy[-1]+1]; zz=z[iy[0]:iy[-1]+1,ix[0]:ix[-1]+1].copy()
 return xs,ys,zz

def corr_valid(a,b):
 m=np.isfinite(a)&np.isfinite(b)
 if m.sum()<20: return float('nan')
 aa=a[m]; bb=b[m]
 if np.std(aa)==0 or np.std(bb)==0: return float('nan')
 return float(np.corrcoef(aa,bb)[0,1])

def compute_metrics(lat,lon,radius):
 xs,ys,zz=extract(lat,lon,radius)
 valid=np.isfinite(zz)
 if valid.mean()<0.99: raise RuntimeError(f'valid fraction below frozen threshold for {lat},{lon} r={radius}: {valid.mean()}')
 # physical axes in metres using local latitude
 X1=(xs-lon)*111320.0*math.cos(math.radians(lat))
 Y1=(ys-lat)*111320.0
 dx=float(np.median(np.diff(X1))); dy=float(np.median(np.diff(Y1)))
 X,Y=np.meshgrid(X1,Y1)
 # Plane fit on all valid samples
 A=np.column_stack([X[valid],Y[valid],np.ones(valid.sum())]); b=zz[valid]
 coef,_,_,_=np.linalg.lstsq(A,b,rcond=None)
 plane=coef[0]*X+coef[1]*Y+coef[2]
 detr=zz-plane
 relief=float(np.nanmax(zz)-np.nanmin(zz))
 rmse=float(np.sqrt(np.nanmean(detr**2)))
 # First and second derivatives; NaNs only affect local stencil neighbourhood.
 dzdy,dzdx=np.gradient(zz,dy,dx)
 g2=dzdx**2+dzdy**2
 slope=np.degrees(np.arctan(np.sqrt(g2)))
 rug=np.sqrt(1.0+g2)
 dpx_dy,dpx_dx=np.gradient(dzdx,dy,dx)
 dqy_dy,dqy_dx=np.gradient(dzdy,dy,dx)
 r2=dpx_dx; t2=dqy_dy; s2=0.5*(dpx_dy+dqy_dx)
 eps=1e-12
 grad=np.sqrt(g2)
 m=np.isfinite(r2)&np.isfinite(t2)&np.isfinite(s2)&np.isfinite(dzdx)&np.isfinite(dzdy)&(grad>1e-6)
 prof=np.full_like(zz,np.nan,dtype=float); plan_curv=np.full_like(zz,np.nan,dtype=float)
 denom_prof=(g2)*(1.0+g2)**1.5
 denom_plan=(g2)**1.5
 prof[m]=-(r2[m]*dzdx[m]**2 + 2*s2[m]*dzdx[m]*dzdy[m] + t2[m]*dzdy[m]**2)/(denom_prof[m]+eps)
 plan_curv[m]=-(dzdy[m]**2*r2[m] - 2*dzdx[m]*dzdy[m]*s2[m] + dzdx[m]**2*t2[m])/(denom_plan[m]+eps)
 # Rotation similarity on centered square, plane-detrended.
 n=min(detr.shape); r0=(detr.shape[0]-n)//2; c0=(detr.shape[1]-n)//2; sq=detr[r0:r0+n,c0:c0+n]
 rots=[corr_valid(sq,np.rot90(sq,k)) for k in (1,2,3)]
 radial=float(np.nanmean(rots))
 # Aspect histogram on non-flat valid gradient samples.
 am=np.isfinite(dzdx)&np.isfinite(dzdy)&(grad>1e-6)
 aspect=np.arctan2(dzdy[am],dzdx[am])
 hist,_=np.histogram(aspect,bins=12,range=(-math.pi,math.pi)); hp=hist/hist.sum()
 nz=hp[hp>0]
 entropy=float(-(nz*np.log(nz)).sum()/math.log(12))
 dominant=float(hp.max())
 return {
   'local_relief_m':relief,
   'slope_median_deg':float(np.nanmedian(slope)),
   'slope_p95_deg':float(np.nanpercentile(slope[np.isfinite(slope)],95)),
   'rugosity_surface_ratio':float(np.nanmean(rug)),
   'planarity_rmse_over_relief':float(rmse/relief) if relief>0 else float('nan'),
   'profile_curvature_rms':float(np.sqrt(np.nanmean(prof**2))),
   'plan_curvature_rms':float(np.sqrt(np.nanmean(plan_curv**2))),
   'radial_rotation_similarity':radial,
   'facet_aspect_entropy':entropy,
   'facet_dominant_bin_fraction':dominant,
   '_integrity':{'valid_fraction':float(valid.mean()),'rows':int(zz.shape[0]),'cols':int(zz.shape[1]),'dx_m':dx,'dy_m':dy,'rotation_components':rots}
 }

def med_mad(vals):
 a=np.array(vals,dtype=float); med=float(np.nanmedian(a)); mad=float(np.nanmedian(np.abs(a-med)))
 return med,mad

def robust_abs_z(v,vals):
 med,mad=med_mad(vals); scale=1.4826*mad
 if scale < 1e-12: scale=max(float(np.nanstd(vals)),1e-12)
 return abs(float(v)-med)/scale,med,mad

def scale_run(radius):
 entities=[TARGET]+CONTROLS
 vals={e['id']:compute_metrics(e['lat'],e['lon'],radius) for e in entities}
 target=vals['TARGET']; controls=[vals[e['id']] for e in CONTROLS]
 per_metric={}; target_zs=[]
 for m in METRICS:
  cv=[q[m] for q in controls]
  rz,med,mad=robust_abs_z(target[m],cv); target_zs.append(rz)
  absdev=[abs(v-med) for v in cv]; tdev=abs(target[m]-med)
  extreme_rank=1+sum(d>=tdev for d in absdev)
  per_metric[m]={'target':target[m],'control_median':med,'control_MAD':mad,'target_abs_robust_z':rz,'empirical_extremeness_rank_nplus1':int(extreme_rank),'control_n':len(cv)}
 target_agg=float(np.median(target_zs))
 control_aggs=[]
 for j,q in enumerate(controls):
  zs=[]
  for m in METRICS:
   baseline=[controls[k][m] for k in range(len(controls)) if k!=j]
   rz,_,_=robust_abs_z(q[m],baseline); zs.append(rz)
  control_aggs.append(float(np.median(zs)))
 p=(1+sum(v>=target_agg for v in control_aggs))/(len(control_aggs)+1)
 return {
   'radius_km':radius,
   'target_metrics':target,
   'control_metrics':{e['id']:vals[e['id']] for e in CONTROLS},
   'per_metric_comparison':per_metric,
   'aggregate':{'target_median_abs_robust_z':target_agg,'control_leave_one_out_scores':control_aggs,'empirical_p':float(p),'target_rank_high_to_low':int(1+sum(v>target_agg for v in control_aggs)),'n_controls':len(control_aggs)}
 }

scales={str(radius):scale_run(radius) for radius in SCALES}
result={
 'artifact_id':'JANUS-H10N1-GENERIC-BLIND-MORPHOLOGY-RUN-001-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'authorized_by':'JANUS-H10N1-MORPHOLOGY-METRICS-COUNCIL-RUN-005-2026-08-22-v1.0',
 'source':{'dataset_id':'GB/NERC/BAS/PDC/01236','netcdf_sha256':controls_manifest['source']['netcdf_sha256'],'crs':'WGS84 geographic','native_spacing_deg':0.0005},
 'code_sha256':SCRIPT_SHA,
 'control_manifest_sha256':controls_manifest['sha256'],
 'metric_council_sha256':council['sha256'],
 'frozen_metrics':METRICS,
 'scales':scales,
 'context_only_5km_10km':'NOT_SCORED_BY_COUNCIL_RULE',
 'pyramidality_score_computed':False,
 'crater_score_computed':False,
 'cosmos_or_symbolic_weight_used':False,
 'target_identity':'UNCONFIRMED',
 'next_rule':'FREEZE_RESULT_AND_ASK_JANUS_AGAIN_WHETHER_ANY_PHYSICAL_ANOMALY_GATE_ADVANCES_OR_REMAINS_NEGATIVE'
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'primary_2km':scales['2']['aggregate'],'secondary_1km':scales['1']['aggregate'],'target_2km_metrics':scales['2']['target_metrics'],'target_1km_metrics':scales['1']['target_metrics'],'sha256':result['sha256']},indent=2,ensure_ascii=False))
