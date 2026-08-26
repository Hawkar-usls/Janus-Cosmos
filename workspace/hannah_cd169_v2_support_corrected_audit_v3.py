#!/usr/bin/env python3
"""Support-corrected audit of the unchanged V2 model.

Fixes only the evaluation data-selection bug caused by hardcoding sd11285. Uses
raw2prism CDF header values (time, altitude, gyro heading, pressure) across the
sd11285->sd11286 boundary. No imagery/intensity and no V2 parameter changes.
"""
from __future__ import annotations
import argparse, ftplib, hashlib, json, math, os, statistics, tempfile
from datetime import datetime, timezone
from pathlib import Path
import hannah_cd169_ground_fixed_nav_reconstruction as v1
import hannah_cd169_veh_nav_directional_repair_v2 as v2

HOST='livftp.noc.ac.uk';BASE='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/TOBI'
FILES=['cd169p5.cdf','cd169p6.cdf']
ANCHORS=[
 ('00:00:00',-3.907,-12.153,6481),('00:15:00',-3.895,-12.167,6277),('00:29:00',-3.886,-12.151,6262),
 ('01:00:00',-3.869,-12.145,6032),('01:29:00',-3.855,-12.135,6155),('02:00:00',-3.835,-12.128,6050),
 ('02:30:00',-3.814,-12.128,6016),('03:00:00',-3.795,-12.131,6112),('03:30:00',-3.774,-12.134,6314),
 ('04:00:00',-3.752,-12.136,6308),('04:30:00',-3.734,-12.140,6634),('05:00:00',-3.715,-12.146,6586),
 ('05:30:00',-3.683,-12.181,6304),('06:00:00',-3.670,-12.151,6110),('06:31:00',-3.649,-12.154,6200)
]
DISCOVERY={'01:00:00','01:29:00'}

def ftp():
 f=ftplib.FTP(timeout=90);f.connect(HOST,21);f.login('anonymous','janus-probe@example.invalid');f.voidcmd('TYPE I');return f

def get(path,dest):
 f=ftp();h=hashlib.sha256();n=0
 try:
  with open(dest,'wb') as out:
   def cb(b):
    nonlocal n;n+=len(b);h.update(b);out.write(b)
   f.retrbinary('RETR '+path,cb,1048576)
 finally:
  try:f.quit()
  except Exception:
   try:f.close()
   except Exception:pass
 return {'sha256':h.hexdigest(),'size_bytes':n}

def dtrow(d,t,s):
 try:return datetime(int(d[0]),int(d[1]),int(d[2]),int(t[0]),int(t[1]),int(math.floor(float(s))),tzinfo=timezone.utc)
 except:return None

def pct(xs,p):
 if not xs:return None
 s=sorted(xs);q=(len(s)-1)*p;k=int(math.floor(q));f=q-k
 return s[k] if k==len(s)-1 else s[k]*(1-f)+s[k+1]*f

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);a=ap.parse_args();out={'schema':'janus.cosmos.cousteau.hannah_cd169.v2_support_corrected_audit.v3','status':'STARTED','model':'UNCHANGED_V2_RECIPROCAL_RAW_TOBI_GYRO_LAYBACK','model_parameters_changed':False,'image_or_intensity_read':False,'reason':'Correct V2 evaluation support across sd11285->sd11286; old V2 code hardcoded sd11285.'}
 try:
  import numpy as np
  from netCDF4 import Dataset
  nr=v1.get(v1.NAV);cr=v1.get(v1.CABLE);nav=v1.parse_nav(nr);cable=v1.parse_cable(cr);v1.add_arc(nav);headers=[];out['cdf_files']=[]
  with tempfile.TemporaryDirectory() as td:
   for fn in FILES:
    local=os.path.join(td,fn);dr=get(BASE+'/'+fn,local);out['cdf_files'].append({'file':fn,**dr});ds=Dataset(local,'r')
    try:
     date=np.asarray(ds.variables['date'][:]);tim=np.asarray(ds.variables['time'][:]);sec=np.asarray(ds.variables['seconds'][:]);ssa=np.asarray(ds.variables['ss_attributes'][:]);sat=np.asarray(ds.variables['ss_attributes_tobi'][:]);n=min(len(date),len(tim),len(sec),len(ssa),len(sat))
     for i in range(n):
      dt=dtrow(date[i],tim[i],sec[i])
      if dt is None:continue
      headers.append({'utc':dt,'index':i,'cdf_file':fn,'altitude_m':float(ssa[i][0]),'heading_deg':float(ssa[i][1])%360,'pressure_m_equiv':float(sat[i][0])})
    finally:ds.close()
  headers.sort(key=lambda x:x['utc']);rows=[]
  for tim,lat,lon,logged_lay in ANCHORS:
   hh,mm,ss=map(int,tim.split(':'));t=datetime(2005,2,28,hh,mm,ss,tzinfo=timezone.utc);h=min(headers,key=lambda x:abs((x['utc']-t).total_seconds()));rec=v2.reconstruct_v2(nav,cable,h);res=v1.hav_m((rec['vehicle']['lat'],rec['vehicle']['lon']),(lat,lon)) if rec.get('valid') else None
   rows.append({'anchor_utc':t.isoformat().replace('+00:00','Z'),'split':'DISCOVERY' if tim in DISCOVERY else 'HELDOUT','cdf_file':h['cdf_file'],'cdf_row_index':h['index'],'header_utc':h['utc'].isoformat().replace('+00:00','Z'),'header_delta_s':(h['utc']-t).total_seconds(),'heading_deg':h['heading_deg'],'pressure_m_equiv':h['pressure_m_equiv'],'v2_vehicle':rec.get('vehicle'),'v2_residual_m':res,'v2_horizontal_layback_m':rec.get('horizontal_layback_m'),'layback_delta_vs_logged_m':None if not rec.get('valid') else rec['horizontal_layback_m']-logged_lay})
  held=[r['v2_residual_m'] for r in rows if r['split']=='HELDOUT' and r['v2_residual_m'] is not None];out['anchors']=rows;out['heldout_support_corrected_stats']={'n':len(held),'median_m':statistics.median(held),'p90_m':pct(held,.9),'max_m':max(held)}
  old={'n':13,'median_m':116.23842045021642,'p90_m':1695.6104609129748,'max_m':4147.310164134658};out['old_v2_stats_from_sd11285_hardcode']=old
  r0631=next(r for r in rows if r['anchor_utc'].endswith('06:31:00Z'));out['support_correction']={'06:31_old_raw_time':'2005-02-28T06:22:16Z','06:31_correct_header_time':r0631['header_utc'],'06:31_correct_cdf_file':r0631['cdf_file'],'06:31_correct_residual_m':r0631['v2_residual_m'],'classification':'OLD_06_31_RESIDUAL_SUPERSEDED_BY_CORRECT_SUPPORT'}
  out['interpretation']='Audit changes only per-anchor header support, not V2 physics. If aggregate metrics change, this is an evaluation correction, not V4 fitting.';out['status']='V2_SUPPORT_CORRECTED_AUDIT_READY'
 except Exception as e:out['status']='V2_SUPPORT_CORRECTED_AUDIT_FAILED';out['error_type']=type(e).__name__;out['error']=str(e)
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8');print(json.dumps({'status':out['status'],'stats':out.get('heldout_support_corrected_stats'),'support_correction':out.get('support_correction')},indent=2));return 0 if out['status']=='V2_SUPPORT_CORRECTED_AUDIT_READY' else 2
if __name__=='__main__':raise SystemExit(main())
