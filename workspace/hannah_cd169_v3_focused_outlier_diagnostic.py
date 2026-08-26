#!/usr/bin/env python3
from __future__ import annotations
import argparse, ftplib, hashlib, json, math, os, tempfile
from datetime import datetime, timezone
from pathlib import Path

HOST='livftp.noc.ac.uk';BASE='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/TOBI'
RAW6='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11286/TOBI.DAT';BLOCK=40960
CASES={
 '00:15':{'utc':'2005-02-28T00:15:00Z','file':'cd169p5.cdf','ship':(-3.849,-12.133),'tobi':(-3.895,-12.167)},
 '05:30':{'utc':'2005-02-28T05:30:00Z','file':'cd169p5.cdf','ship':(-3.636,-12.149),'tobi':(-3.683,-12.181)},
 '06:31':{'utc':'2005-02-28T06:31:00Z','file':'cd169p6.cdf','ship':(-3.592,-12.155),'tobi':(-3.649,-12.154)}
}

def ftp():
 f=ftplib.FTP(timeout=90);f.connect(HOST,21);f.login('anonymous','janus-probe@example.invalid');f.voidcmd('TYPE I');return f

def get(path,dest=None):
 f=ftp();parts=[];h=hashlib.sha256();n=0;out=open(dest,'wb') if dest else None
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

def parse_nav(raw):
 out=[]
 for line in raw.decode('utf-8','replace').splitlines():
  p=line.split()
  if len(p)<5:continue
  try:
   ds,hm=p[1],p[2];t=datetime(2000+int(ds[:2]),int(ds[2:4]),int(ds[4:6]),int(hm[:2]),int(hm[2:4]),tzinfo=timezone.utc);lat=float(p[3]);lon=float(p[4])
  except Exception:continue
  out.append((t,lat,lon))
 return sorted(out)

def interp(nav,t):
 q=t.timestamp();lo=0;hi=len(nav)-1
 if q<=nav[0][0].timestamp():return nav[0][1:]
 if q>=nav[-1][0].timestamp():return nav[-1][1:]
 while hi-lo>1:
  m=(lo+hi)//2
  if nav[m][0].timestamp()<=q:lo=m
  else:hi=m
 a,b=nav[lo],nav[hi];f=(q-a[0].timestamp())/(b[0].timestamp()-a[0].timestamp());return (a[1]+f*(b[1]-a[1]),a[2]+f*(b[2]-a[2]))

def hav(a,b):
 r=6371008.8;p1,l1=map(math.radians,a);p2,l2=map(math.radians,b);dp=p2-p1;dl=l2-l1;h=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2;return 2*r*math.asin(min(1,math.sqrt(h)))

def bearing(a,b):
 p1=math.radians(a[0]);p2=math.radians(b[0]);dl=math.radians(b[1]-a[1]);y=math.sin(dl)*math.cos(p2);x=math.cos(p1)*math.sin(p2)-math.sin(p1)*math.cos(p2)*math.cos(dl);return (math.degrees(math.atan2(y,x))+360)%360

def adiff(a,b):return abs((a-b+180)%360-180)

def fix(x):return 2*math.floor(float(x))-float(x) if float(x)<0 else float(x)

def dtrow(d,t,s):
 try:return datetime(int(d[0]),int(d[1]),int(d[2]),int(t[0]),int(t[1]),int(math.floor(float(s))),tzinfo=timezone.utc)
 except:return None

def raw_sd11286_0631():
 # sd11286 first valid is 06:22:24; 06:31 is nominally index 129. Read a small sequential prefix only.
 f=ftp();buf=bytearray();want=140*BLOCK
 try:
  sock=f.transfercmd('RETR '+RAW6)
  try:
   while len(buf)<want:
    b=sock.recv(min(1048576,want-len(buf)))
    if not b:break
    buf.extend(b)
  finally:
   try:sock.close()
   except:pass
 finally:
  try:f.close()
  except:pass
 import struct,statistics
 best=None
 for i in range(len(buf)//BLOCK):
  b=bytes(buf[i*BLOCK:(i+1)*BLOCK]);tm,dt,alt=struct.unpack_from('<HHH',b,0x32)
  try:stamp=datetime(1980+((dt>>9)&127),(dt>>5)&15,dt&31,(tm>>11)&31,(tm>>5)&63,(tm&31)*2,tzinfo=timezone.utc)
  except:continue
  target=datetime(2005,2,28,6,31,0,tzinfo=timezone.utc);delta=abs((stamp-target).total_seconds())
  if best is None or delta<best['abs_delta_s']:
   gyro=struct.unpack_from('<8h',b,0xD4);heading_vals=[x/10.0-10.1 for x in gyro];s=sum(math.sin(math.radians(x)) for x in heading_vals);c=sum(math.cos(math.radians(x)) for x in heading_vals);heading=(math.degrees(math.atan2(s,c))+360)%360
   best={'record_index':i,'utc':stamp.isoformat().replace('+00:00','Z'),'abs_delta_s':delta,'altitude_raw':alt,'heading_deg':heading,'block_sha256':hashlib.sha256(b).hexdigest(),'bytes_read':len(buf),'source':'sd11286/TOBI.DAT'}
 return best

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);a=ap.parse_args();out={'schema':'janus.cosmos.cousteau.hannah_cd169.v3_focused_outlier_diagnostic.v1','status':'STARTED','role':'DIAGNOSTIC_ONLY','image_or_intensity_read':False,'parameters_fitted':False,'outliers_dropped':False}
 try:
  import numpy as np
  from netCDF4 import Dataset
  navr=get(BASE+'/cd169.nav');nav=parse_nav(navr['bytes']);out['cd169_nav_sha256']=navr['sha256'];rows={}
  with tempfile.TemporaryDirectory() as td:
   cache={}
   for key,cfg in CASES.items():
    fn=cfg['file']
    if fn not in cache:
     local=os.path.join(td,fn);dr=get(BASE+'/'+fn,local);cache[fn]=(local,dr)
    local,dr=cache[fn];ds=Dataset(local,'r')
    try:
     date=np.asarray(ds.variables['date'][:]);tim=np.asarray(ds.variables['time'][:]);sec=np.asarray(ds.variables['seconds'][:]);sl=np.asarray(ds.variables['ship_latlon'][:]);ssa=np.asarray(ds.variables['ss_attributes'][:]);target=datetime.fromisoformat(cfg['utc'].replace('Z','+00:00'));best=None
     for i in range(min(len(date),len(tim),len(sec),len(sl),len(ssa))):
      dt=dtrow(date[i],tim[i],sec[i])
      if dt is None:continue
      dd=abs((dt-target).total_seconds())
      if best is None or dd<best[0]:best=(dd,i,dt)
     dd,i,dt=best;raw_ship=(float(sl[i][0]),float(sl[i][1]));corrected=(fix(raw_ship[0]),fix(raw_ship[1]));ext=interp(nav,dt);heading=float(ssa[i][1]);logged_b=bearing(cfg['ship'],cfg['tobi']);recip=(heading+180)%360
     tags=[];source_res=hav(ext,cfg['ship']);corr_log_res=hav(corrected,cfg['ship']);bearing_mis=adiff(recip,logged_b)
     if source_res>=1000:tags.append('SHIP_NAV_SOURCE_DISAGREEMENT_CANDIDATE')
     if source_res<=250 and bearing_mis>=20:tags.append('HEADING_OR_TOW_DYNAMICS_CANDIDATE')
     rows[key]={'anchor_utc':cfg['utc'],'cdf_file':fn,'cdf_sha256':dr['sha256'],'cdf_row_index':i,'cdf_utc':dt.isoformat().replace('+00:00','Z'),'cdf_time_delta_s':(dt-target).total_seconds(),'science_log_ship':list(cfg['ship']),'science_log_tobi':list(cfg['tobi']),'external_cd169_nav_ship':list(ext),'external_nav_vs_science_log_ship_m':source_res,'raw2prism_ship_latlon':list(raw_ship),'signbug_inverse_ship_latlon':list(corrected),'signbug_inverse_vs_science_log_ship_m':corr_log_res,'cdf_heading_deg':heading,'reciprocal_heading_deg':recip,'logged_ship_to_tobi_distance_m':hav(cfg['ship'],cfg['tobi']),'logged_ship_to_tobi_bearing_deg':logged_b,'bearing_mismatch_deg':bearing_mis,'diagnostic_tags':tags}
    finally:ds.close()
  raw6=raw_sd11286_0631();rows['06:31']['raw_support']=raw6
  if raw6 and raw6['abs_delta_s']<=10:rows['06:31']['diagnostic_tags'].append('RAW_TIME_SUPPORT_GOOD_SD11286')
  out['cases']=rows
  out['diagnosis']={
   '00:15':'SHIP_NAV_SOURCE_DISAGREEMENT_DOMINANT_CANDIDATE' if 'SHIP_NAV_SOURCE_DISAGREEMENT_CANDIDATE' in rows['00:15']['diagnostic_tags'] and rows['00:15']['bearing_mismatch_deg']<5 else 'UNRESOLVED',
   '05:30':'HEADING_OR_TOW_DYNAMICS_DOMINANT_CANDIDATE' if 'HEADING_OR_TOW_DYNAMICS_CANDIDATE' in rows['05:30']['diagnostic_tags'] else 'UNRESOLVED',
   '06:31':'RAW_SUPPORT_CONFIRMED_IN_SD11286' if raw6 and raw6['abs_delta_s']<=10 else 'RAW_SUPPORT_HOLD',
   'authority':'POSTHOC_DIAGNOSTIC_NOT_BLINDED_MODEL_VALIDATION'
  }
  out['status']='V3_FOCUSED_OUTLIER_DIAGNOSTIC_READY'
 except Exception as e:out['status']='V3_FOCUSED_OUTLIER_DIAGNOSTIC_FAILED';out['error_type']=type(e).__name__;out['error']=str(e)
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8');print(json.dumps({'status':out['status'],'diagnosis':out.get('diagnosis'),'cases':out.get('cases')},indent=2));return 0 if out['status']=='V3_FOCUSED_OUTLIER_DIAGNOSTIC_READY' else 2
if __name__=='__main__':raise SystemExit(main())
