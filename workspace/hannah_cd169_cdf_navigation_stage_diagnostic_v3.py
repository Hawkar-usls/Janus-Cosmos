#!/usr/bin/env python3
"""Post-run diagnostic after V3B formal-rule false positive.

This is NOT a promotion gate. It tests whether archived raw2prism CDF headers are
pre-navigation, and whether ship_latlon reflects a deterministic signed-degrees
conversion defect. No sonar image/intensity arrays are read and no parameters are fit.
"""
from __future__ import annotations
import argparse, ftplib, hashlib, json, math, os, statistics, tempfile
from datetime import datetime, timezone
from pathlib import Path

HOST='livftp.noc.ac.uk';BASE='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/TOBI'
CDF_FILES=[f'cd169p{i}.cdf' for i in range(1,12)];NAV=BASE+'/cd169.nav';CABLE=BASE+'/cd169.cable'
TARGET=datetime(2005,2,28,1,7,25,tzinfo=timezone.utc)

def ftp():
 f=ftplib.FTP(timeout=90);f.connect(HOST,21);f.login('anonymous','janus-probe@example.invalid');f.voidcmd('TYPE I');return f

def get(path,dest=None):
 f=ftp();h=hashlib.sha256();parts=[];n=0;out=open(dest,'wb') if dest else None
 try:
  def cb(b):
   nonlocal n;n+=len(b);h.update(b);out.write(b) if out else parts.append(b)
  f.retrbinary('RETR '+path,cb,1048576)
 finally:
  if out:out.close()
  try:f.quit()
  except Exception:
   try:f.close()
   except Exception:pass
 return {'bytes':None if dest else b''.join(parts),'sha256':h.hexdigest(),'size_bytes':n}

def parse_series(raw,key):
 out=[]
 for line in raw.decode('utf-8','replace').splitlines():
  p=line.split()
  if len(p)<4:continue
  try:
   ds,hm=p[1],p[2];t=datetime(2000+int(ds[:2]),int(ds[2:4]),int(ds[4:6]),int(hm[:2]),int(hm[2:4]),tzinfo=timezone.utc)
   nums=[float(x) for x in p[3:]]
  except Exception:continue
  if key=='nav' and len(nums)>=2:out.append((t,nums[0],nums[1]))
  elif key=='cable' and nums:out.append((t,nums[-1]))
 return sorted(out)

def interp(rows,t,cols):
 if not rows:return None
 q=t.timestamp()
 if q<rows[0][0].timestamp() or q>rows[-1][0].timestamp():return None
 lo=0;hi=len(rows)-1
 while hi-lo>1:
  m=(lo+hi)//2
  if rows[m][0].timestamp()<=q:lo=m
  else:hi=m
 a,b=rows[lo],rows[hi];den=b[0].timestamp()-a[0].timestamp();f=0 if den==0 else (q-a[0].timestamp())/den
 return tuple(a[i]+f*(b[i]-a[i]) for i in range(1,cols+1))

def hav(a,b):
 r=6371008.8;p1,l1=map(math.radians,a);p2,l2=map(math.radians,b);dp=p2-p1;dl=l2-l1
 h=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
 return 2*r*math.asin(min(1,math.sqrt(h)))

def dtrow(d,t,s):
 try:return datetime(int(d[0]),int(d[1]),int(d[2]),int(t[0]),int(t[1]),int(math.floor(float(s))),tzinfo=timezone.utc)
 except Exception:return None

def undo_signed_degree_plus_minutes_bug(x):
 # If raw degree D<0 and old converter used D + minutes/60, observed x lies in [D,D+1).
 # floor(x)=D, so correct D-minutes/60 = 2*D-x. No fitted constants.
 return 2.0*math.floor(float(x))-float(x) if float(x)<0 else float(x)

def pct(xs,p):
 if not xs:return None
 s=sorted(xs);q=(len(s)-1)*p;k=int(q);f=q-k
 return s[k] if k==len(s)-1 else s[k]*(1-f)+s[k+1]*f

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
 out={'schema':'janus.cosmos.cousteau.hannah_cd169.cdf_navigation_stage_diagnostic.v3','status':'STARTED','role':'POST_RUN_DIAGNOSTIC_NOT_PROMOTION_GATE','trigger':'V3B formal checker passed relative dominance despite physically invalid absolute coordinates','image_or_intensity_arrays_read':False,'parameters_fitted':False,'files':[]}
 try:
  import numpy as np
  from netCDF4 import Dataset
  navr=get(NAV);cabr=get(CABLE);nav=parse_series(navr['bytes'],'nav');cable=parse_series(cabr['bytes'],'cable')
  ship_orig=[];ship_fix=[];allrows=[];latlon_total=0;latlon_zero=0;latlon_valid_nonzero=0
  with tempfile.TemporaryDirectory() as td:
   for fn in CDF_FILES:
    local=os.path.join(td,fn);dr=get(BASE+'/'+fn,local);ds=Dataset(local,'r')
    try:
     fr={'file':fn,'sha256':dr['sha256'],'creation_program':str(getattr(ds,'creation_program','')),'latlon_semantics':str(getattr(ds,'latlon','')),'ship_latlon_semantics':str(getattr(ds,'ship_latlon',''))}
     date=np.asarray(ds.variables['date'][:]);tim=np.asarray(ds.variables['time'][:]);sec=np.asarray(ds.variables['seconds'][:]);ll=np.asarray(ds.variables['latlon'][:]);sl=np.asarray(ds.variables['ship_latlon'][:]);sat=np.asarray(ds.variables['ss_attributes_tobi'][:]) if 'ss_attributes_tobi' in ds.variables else None
     n=min(len(date),len(tim),len(sec),len(ll),len(sl));fz=0;fnz=0;valid=0;target_candidate=None
     for i in range(n):
      dt=dtrow(date[i],tim[i],sec[i]);
      try:la,lo=float(ll[i][0]),float(ll[i][1]);sla,slo=float(sl[i][0]),float(sl[i][1])
      except Exception:continue
      latlon_total+=1
      if abs(la)<1e-12 and abs(lo)<1e-12:latlon_zero+=1;fz+=1
      elif math.isfinite(la) and math.isfinite(lo) and -90<=la<=90 and -180<=lo<=180:latlon_valid_nonzero+=1;fnz+=1
      if dt is None:continue
      ext=interp(nav,dt,2)
      fla,flo=undo_signed_degree_plus_minutes_bug(sla),undo_signed_degree_plus_minutes_bug(slo)
      row={'file':fn,'row_index':i,'utc':dt,'sonar_lat':la,'sonar_lon':lo,'ship_lat_raw2prism':sla,'ship_lon_raw2prism':slo,'ship_lat_signbug_inverse':fla,'ship_lon_signbug_inverse':flo}
      if sat is not None and i<len(sat):row['tobi_header_water_pressure_cable']=[float(x) for x in np.ravel(sat[i])[:2]]
      if ext:
       ro=hav((sla,slo),ext);rf=hav((fla,flo),ext);ship_orig.append(ro);ship_fix.append(rf);row['raw_ship_to_external_nav_m']=ro;row['signbug_inverse_to_external_nav_m']=rf
      allrows.append(row);valid+=1
      if abs((dt-TARGET).total_seconds())<=5 and (target_candidate is None or abs((dt-TARGET).total_seconds())<abs((target_candidate['utc']-TARGET).total_seconds())):target_candidate=row
     fr.update(rows=n,latlon_zero_rows=fz,latlon_valid_nonzero_rows=fnz,valid_timestamp_rows=valid)
     if target_candidate:
      tc=dict(target_candidate);tc['utc']=tc['utc'].isoformat().replace('+00:00','Z');fr['target_candidate']=tc
     out['files'].append(fr)
    finally:ds.close()
  out['latlon_population']={'rows_total':latlon_total,'exact_zero_zero_rows':latlon_zero,'valid_nonzero_rows':latlon_valid_nonzero,'zero_fraction':latlon_zero/latlon_total if latlon_total else None}
  out['ship_latlon_sign_bug_diagnostic']={'formula':'for negative observed x: corrected = 2*floor(x)-x','basis':'tests old converter hypothesis D + minutes/60 instead of D - minutes/60 for signed negative degree D','matched_rows':len(ship_fix),'raw_ship_latlon_to_external_nav_median_m':statistics.median(ship_orig) if ship_orig else None,'corrected_to_external_nav_median_m':statistics.median(ship_fix) if ship_fix else None,'corrected_p90_m':pct(ship_fix,.9),'improvement_factor_median':(statistics.median(ship_orig)/statistics.median(ship_fix)) if ship_fix and statistics.median(ship_fix)>0 else None,'no_fit':True}
  tr=min(allrows,key=lambda r:abs((r['utc']-TARGET).total_seconds())) if allrows else None
  if tr:
   ext=interp(nav,tr['utc'],2);cab=interp(cable,tr['utc'],1);tc=dict(tr);tc['utc']=tr['utc'].isoformat().replace('+00:00','Z');tc['delta_to_target_s']=(tr['utc']-TARGET).total_seconds();tc['external_ship_nav']=None if ext is None else {'lat':ext[0],'lon':ext[1]};tc['external_cable_m']=None if cab is None else cab[0]
   if ext:tc['corrected_ship_to_external_nav_m']=hav((tr['ship_lat_signbug_inverse'],tr['ship_lon_signbug_inverse']),ext)
   out['target']=tc
  allzero=(latlon_total>0 and latlon_valid_nonzero==0)
  strongfix=bool(ship_fix and statistics.median(ship_fix)<statistics.median(ship_orig)/100)
  out['diagnosis']={
   'cdf_latlon_vehicle_channel_population':'UNPOPULATED_ZERO_ZERO' if allzero else 'HAS_NONZERO_VALUES',
   'ship_latlon_signed_degree_conversion_bug':'STRONGLY_SUPPORTED' if strongfix else 'NOT_ESTABLISHED',
   'cdf_stage':'PRE_MRGNAV_OR_NAVIGATION_NOT_MERGED' if allzero else 'UNRESOLVED',
   'V3B_prior_formal_pass_scientific_disposition':'REJECT_FALSE_POSITIVE_UNDER_SPECIFIED_ABSOLUTE_SANITY',
   'G3B_CALIBRATED_NAV_PRODUCT':'FAIL_CDF_NOT_CALIBRATED_VEHICLE_NAV' if allzero else 'HOLD',
   'native_veh_nav_recovered':False
  }
  out['status']='CDF_PRE_MRGNAV_STAGE_DIAGNOSED' if allzero and strongfix else 'CDF_NAV_STAGE_DIAGNOSTIC_HOLD'
 except Exception as e:
  out['status']='CDF_NAV_STAGE_DIAGNOSTIC_FAILED';out['error_type']=type(e).__name__;out['error']=str(e)
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
 print(json.dumps({'status':out['status'],'latlon_population':out.get('latlon_population'),'signbug':out.get('ship_latlon_sign_bug_diagnostic'),'target':out.get('target'),'diagnosis':out.get('diagnosis')},indent=2));return 0 if out['status'] in {'CDF_PRE_MRGNAV_STAGE_DIAGNOSED','CDF_NAV_STAGE_DIAGNOSTIC_HOLD'} else 2
if __name__=='__main__':raise SystemExit(main())
