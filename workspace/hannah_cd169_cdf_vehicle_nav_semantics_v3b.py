#!/usr/bin/env python3
from __future__ import annotations
import argparse, ftplib, hashlib, json, math, os, statistics, tempfile
from datetime import datetime, timezone
from pathlib import Path

HOST='livftp.noc.ac.uk'
BASE='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/TOBI'
CDF_FILES=[f'cd169p{i}.cdf' for i in range(1,12)]
NAV_PATH=BASE+'/cd169.nav'
TARGET=datetime(2005,2,28,1,7,25,tzinfo=timezone.utc)
FROZEN=(-3.8654180644718967,-12.142441475)
RAW_SAMPLE_M=0.75
ANCHORS=[
 ('00:00:00',-3.907,-12.153),('00:15:00',-3.895,-12.167),('00:29:00',-3.886,-12.151),
 ('01:00:00',-3.869,-12.145),('01:29:00',-3.855,-12.135),('02:00:00',-3.835,-12.128),
 ('02:30:00',-3.814,-12.128),('03:00:00',-3.795,-12.131),('03:30:00',-3.774,-12.134),
 ('04:00:00',-3.752,-12.136),('04:30:00',-3.734,-12.140),('05:00:00',-3.715,-12.146),
 ('05:30:00',-3.683,-12.181),('06:00:00',-3.670,-12.151),('06:31:00',-3.649,-12.154)
]

def ftp():
    f=ftplib.FTP(timeout=90);f.connect(HOST,21);f.login('anonymous','janus-probe@example.invalid');f.voidcmd('TYPE I');return f

def get(path,dest=None):
    f=ftp();h=hashlib.sha256();parts=[];n=0
    try:
        out=open(dest,'wb') if dest else None
        def cb(b):
            nonlocal n;n+=len(b);h.update(b)
            if out:out.write(b)
            else:parts.append(b)
        f.retrbinary('RETR '+path,cb,1048576)
        if out:out.close()
    finally:
        try:f.quit()
        except Exception:
            try:f.close()
            except Exception:pass
    return {'bytes':b''.join(parts) if dest is None else None,'size_bytes':n,'sha256':h.hexdigest()}

def parse_nav(raw):
    rows=[]
    for line in raw.decode('utf-8','replace').splitlines():
        p=line.split()
        if len(p)<5:continue
        try:
            ds,hm=p[1],p[2]
            d=datetime(2000+int(ds[:2]),int(ds[2:4]),int(ds[4:6]),int(hm[:2]),int(hm[2:4]),tzinfo=timezone.utc)
            lat=float(p[3]);lon=float(p[4])
        except Exception:continue
        if -90<=lat<=90 and -180<=lon<=180:rows.append((d,lat,lon))
    return sorted(rows)

def interp_nav(nav,t):
    if not nav:return None
    q=t.timestamp()
    if q<=nav[0][0].timestamp():return {'lat':nav[0][1],'lon':nav[0][2],'extrapolated':True}
    if q>=nav[-1][0].timestamp():return {'lat':nav[-1][1],'lon':nav[-1][2],'extrapolated':True}
    lo=0;hi=len(nav)-1
    while hi-lo>1:
        m=(lo+hi)//2
        if nav[m][0].timestamp()<=q:lo=m
        else:hi=m
    a,b=nav[lo],nav[hi];f=(q-a[0].timestamp())/(b[0].timestamp()-a[0].timestamp())
    return {'lat':a[1]+f*(b[1]-a[1]),'lon':a[2]+f*(b[2]-a[2]),'extrapolated':False}

def hav_m(a,b):
    r=6371008.8;p1,l1=map(math.radians,a);p2,l2=map(math.radians,b)
    dp=p2-p1;dl=l2-l1;h=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*r*math.asin(min(1,math.sqrt(h)))

def en_offset(origin,target):
    lat0=math.radians((origin[0]+target[0])/2);north=(target[0]-origin[0])*111195.08;east=(target[1]-origin[1])*111195.08*math.cos(lat0)
    return east,north

def finite2(x):
    try:
        a=float(x[0]);b=float(x[1]);return math.isfinite(a) and math.isfinite(b) and -90<=a<=90 and -180<=b<=180
    except Exception:return False

def make_dt(date_row,time_row,sec):
    try:
        y,mo,da=[int(x) for x in date_row[:3]];hh,mm=[int(x) for x in time_row[:2]];s=float(sec)
        whole=int(math.floor(s));micro=int(round((s-whole)*1e6))
        return datetime(y,mo,da,hh,mm,whole,micro,tzinfo=timezone.utc)
    except Exception:return None

def row_serial(r):
    o=dict(r);o['utc']=r['utc'].isoformat().replace('+00:00','Z');return o

def percentile(xs,p):
    if not xs:return None
    s=sorted(xs);q=(len(s)-1)*p;k=int(math.floor(q));f=q-k
    return s[k] if k==len(s)-1 else s[k]*(1-f)+s[k+1]*f

def nearest(rows,t):
    if not rows:return None
    return min(rows,key=lambda r:abs((r['utc']-t).total_seconds()))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
    out={'schema':'janus.cosmos.cousteau.hannah_cd169.v3b_cdf_vehicle_nav_semantics.run.v1','contract':'JANUS-HANNAH-CD169-V3B-CDF-VEHICLE-NAV-SEMANTICS-CONTRACT-2026-08-26-v1.0','status':'STARTED','image_or_intensity_arrays_read':False,'navigation_parameters_fitted':False,'recentered':False,'files':[]}
    try:
        import numpy as np
        from netCDF4 import Dataset
        navres=get(NAV_PATH);nav=parse_nav(navres['bytes'])
        out['external_ship_nav']={'path':'sd11281/TOBI/cd169.nav','sha256':navres['sha256'],'size_bytes':navres['size_bytes'],'rows':len(nav)}
        allrows=[]
        with tempfile.TemporaryDirectory() as td:
            for fn in CDF_FILES:
                path=BASE+'/'+fn;local=os.path.join(td,fn);dr=get(path,local)
                fr={'file':fn,'sha256':dr['sha256'],'size_bytes':dr['size_bytes'],'rows_read':0,'global_semantics':{}}
                ds=Dataset(local,'r')
                try:
                    fr['creation_program']=str(getattr(ds,'creation_program',''))
                    fr['global_semantics']['latlon']=str(getattr(ds,'latlon',''))
                    fr['global_semantics']['ship_latlon']=str(getattr(ds,'ship_latlon',''))
                    required=['date','time','seconds','latlon','ship_latlon']
                    if not all(k in ds.variables for k in required):
                        fr['status']='MISSING_REQUIRED_HEADER_VARIABLE';out['files'].append(fr);continue
                    date=np.asarray(ds.variables['date'][:]);tim=np.asarray(ds.variables['time'][:]);sec=np.asarray(ds.variables['seconds'][:]);ll=np.asarray(ds.variables['latlon'][:]);sl=np.asarray(ds.variables['ship_latlon'][:])
                    ssa=np.asarray(ds.variables['ss_attributes'][:]) if 'ss_attributes' in ds.variables else None
                    sat=np.asarray(ds.variables['ss_attributes_tobi'][:]) if 'ss_attributes_tobi' in ds.variables else None
                    dep=np.asarray(ds.variables['depths'][:]) if 'depths' in ds.variables else None
                    n=min(len(date),len(tim),len(sec),len(ll),len(sl));fr['rows_read']=int(n);valid=0
                    first=None;last=None
                    for i in range(n):
                        dt=make_dt(date[i],tim[i],sec[i])
                        if dt is None or not finite2(ll[i]) or not finite2(sl[i]):continue
                        sonar=(float(ll[i][0]),float(ll[i][1]));ship=(float(sl[i][0]),float(sl[i][1]))
                        row={'file':fn,'row_index':i,'utc':dt,'sonar_lat':sonar[0],'sonar_lon':sonar[1],'ship_lat':ship[0],'ship_lon':ship[1],'sonar_ship_separation_m':hav_m(sonar,ship)}
                        if ssa is not None and i<len(ssa):
                            vals=[float(x) for x in np.ravel(ssa[i])[:5]];row['ss_attributes']=vals
                            if len(vals)>=2:row['fish_altitude_m']=vals[0];row['heading_deg']=vals[1]
                        if sat is not None and i<len(sat):row['ss_attributes_tobi']=[float(x) for x in np.ravel(sat[i])[:4]]
                        if dep is not None and i<len(dep):row['depths']=[float(x) for x in np.ravel(dep[i])[:4]]
                        allrows.append(row);valid+=1;first=dt if first is None else min(first,dt);last=dt if last is None else max(last,dt)
                    fr['valid_navigation_rows']=valid;fr['first_utc']=first.isoformat().replace('+00:00','Z') if first else None;fr['last_utc']=last.isoformat().replace('+00:00','Z') if last else None;fr['status']='NAV_ROWS_READY';out['files'].append(fr)
                finally:ds.close()
        out['cdf_row_count']=len(allrows)
        # Direct metadata semantics, frozen before value inspection.
        direct=all(f['global_semantics'].get('latlon','').lower().startswith('sonar position') and f['global_semantics'].get('ship_latlon','').lower().startswith('ship position') for f in out['files'] if f.get('status')=='NAV_ROWS_READY')
        out['S1_direct_metadata_semantics_pass']=direct
        # Match every CDF row to independent ship nav and compute two competing residuals.
        shipctl=[]
        for r in allrows:
            nv=interp_nav(nav,r['utc'])
            if nv and not nv['extrapolated']:
                shipctl.append({'file':r['file'],'row_index':r['row_index'],'utc':r['utc'],'ship_to_external_m':hav_m((r['ship_lat'],r['ship_lon']),(nv['lat'],nv['lon'])),'sonar_to_external_ship_m':hav_m((r['sonar_lat'],r['sonar_lon']),(nv['lat'],nv['lon'])),'sonar_ship_separation_m':r['sonar_ship_separation_m']})
        ares=[x['ship_to_external_m'] for x in shipctl];bres=[x['sonar_to_external_ship_m'] for x in shipctl];sep=[x['sonar_ship_separation_m'] for x in shipctl]
        s2=bool(ares and bres and statistics.median(ares)<statistics.median(bres))
        out['S2_ship_control']={'matched_rows':len(shipctl),'ship_latlon_to_cd169_nav_median_m':statistics.median(ares) if ares else None,'ship_latlon_to_cd169_nav_p90_m':percentile(ares,.9),'latlon_to_cd169_nav_median_m':statistics.median(bres) if bres else None,'latlon_to_cd169_nav_p10_m':percentile(bres,.1),'dominance_pass':s2}
        out['S3_sonar_ship_separation']={'n':len(sep),'median_m':statistics.median(sep) if sep else None,'p10_m':percentile(sep,.1),'p90_m':percentile(sep,.9),'numerically_identical':bool(sep and max(sep)<0.01)}
        # Historical TOBI science-log checks; nearest CDF row selected by time only.
        checks=[]
        for tim,lat,lon in ANCHORS:
            hh,mm,ss=map(int,tim.split(':'));t=datetime(2005,2,28,hh,mm,ss,tzinfo=timezone.utc);r=nearest(allrows,t)
            if r is None:continue
            nv=interp_nav(nav,r['utc'])
            checks.append({'anchor_utc':t.isoformat().replace('+00:00','Z'),'cdf':row_serial(r),'cdf_time_delta_s':(r['utc']-t).total_seconds(),'latlon_to_logged_tobi_m':hav_m((r['sonar_lat'],r['sonar_lon']),(lat,lon)),'ship_latlon_to_logged_tobi_m':hav_m((r['ship_lat'],r['ship_lon']),(lat,lon)),'ship_latlon_to_external_nav_m':None if not nv else hav_m((r['ship_lat'],r['ship_lon']),(nv['lat'],nv['lon'])),'external_nav_extrapolated':None if not nv else nv['extrapolated']})
        out['S4_historical_tobi_checks']={'role':'HISTORICAL_CONSISTENCY_NOT_INDEPENDENT_ABSOLUTE_GROUND_TRUTH','anchors':checks}
        # Frozen target by timestamp only.
        tr=nearest(allrows,TARGET)
        if tr:
            sonar=(tr['sonar_lat'],tr['sonar_lon']);ship=(tr['ship_lat'],tr['ship_lon']);e,n=en_offset(sonar,FROZEN)
            heading=tr.get('heading_deg');alt=tr.get('fish_altitude_m');mapping=None
            if heading is not None and alt is not None and math.isfinite(heading) and math.isfinite(alt) and alt>0:
                hr=math.radians(heading);along=e*math.sin(hr)+n*math.cos(hr);cross=e*math.cos(hr)-n*math.sin(hr);slant=math.sqrt(alt*alt+cross*cross);mapping={'east_m':e,'north_m':n,'heading_deg':heading,'fish_altitude_m':alt,'along_track_miss_m':along,'cross_track_m':cross,'ground_range_m':abs(cross),'slant_range_m':slant,'raw_sample_index_round':int(round(slant/RAW_SAMPLE_M)),'side':'stbd' if cross>=0 else 'port','sample_spacing_m':RAW_SAMPLE_M}
            out['target']={'frozen_utc':TARGET.isoformat().replace('+00:00','Z'),'cdf_row':row_serial(tr),'cdf_time_delta_s':(tr['utc']-TARGET).total_seconds(),'archived_sonar_position_to_frozen_m':hav_m(sonar,FROZEN),'archived_ship_position_to_frozen_m':hav_m(ship,FROZEN),'geometry_mapping_without_intensity':mapping}
        g3b=bool(direct and s2 and tr is not None and abs((tr['utc']-TARGET).total_seconds())<=10)
        out['gate']={'G3B_CALIBRATED_NAV_PRODUCT':'PASS_ARCHIVED_PRISM_CDF_SONAR_NAVIGATION' if g3b else 'HOLD','native_intermediate_veh_nav_recovered':False,'authority':'ARCHIVED_PRISM_CDF_PER_PING_SONAR_POSITION_NOT_NATIVE_VEH_NAV_BYTES' if g3b else 'METADATA_DIAGNOSTIC_ONLY','no_fit':True,'no_recenter':True,'image_or_intensity_arrays_read':False}
        out['status']='V3B_CDF_CALIBRATED_SONAR_NAV_PASS' if g3b else 'V3B_CDF_NAV_HOLD'
    except Exception as e:
        out['status']='V3B_CDF_NAV_FAILED';out['error_type']=type(e).__name__;out['error']=str(e)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'S1':out.get('S1_direct_metadata_semantics_pass'),'S2':out.get('S2_ship_control'),'S3':out.get('S3_sonar_ship_separation'),'target':out.get('target'),'gate':out.get('gate')},indent=2))
    return 0 if out['status'] in {'V3B_CDF_CALIBRATED_SONAR_NAV_PASS','V3B_CDF_NAV_HOLD'} else 2
if __name__=='__main__':raise SystemExit(main())
