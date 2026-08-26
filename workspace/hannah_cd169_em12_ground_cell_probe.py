#!/usr/bin/env python3
from __future__ import annotations
import argparse, ftplib, hashlib, json, math, statistics
from pathlib import Path

HOST='livftp.noc.ac.uk'
ROOT='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/EM12'
TARGET=(-3.8654180644718967,-12.142441475)  # lat, lon
R_EARTH=6371008.8


def ftp():
    f=ftplib.FTP(timeout=120);f.connect(HOST,21);f.login('anonymous','janus-probe@example.invalid');f.voidcmd('TYPE I');return f

def list_candidates():
    f=ftp()
    try:
        names=f.nlst(ROOT)
        out=[]
        for p in names:
            name=p.rstrip('/').split('/')[-1]
            if name.lower().endswith('.xyz.ascii'):
                full=p if p.startswith('/') else ROOT+'/'+name
                try:sz=f.size(full)
                except Exception:sz=None
                out.append((name,full,sz))
        return sorted(out)
    finally:
        try:f.quit()
        except Exception:
            try:f.close()
            except Exception:pass

def parse_xyz_line(raw):
    s=raw.decode('ascii','ignore').strip()
    if not s or s.startswith(('#',';','!')):return None
    p=s.replace(',',' ').split()
    if len(p)<3:return None
    try:lon=float(p[0]);lat=float(p[1]);z=float(p[2])
    except Exception:return None
    if not (math.isfinite(lon) and math.isfinite(lat) and math.isfinite(z)):return None
    if not (-180<=lon<=180 and -90<=lat<=90):return None
    return lon,lat,z

def hav_m(lat,lon):
    p1=math.radians(TARGET[0]);p2=math.radians(lat);dp=p2-p1;dl=math.radians(lon-TARGET[1])
    h=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R_EARTH*math.asin(min(1.0,math.sqrt(h)))

def en_m(lat,lon):
    north=math.radians(lat-TARGET[0])*R_EARTH
    east=math.radians(lon-TARGET[1])*R_EARTH*math.cos(math.radians((lat+TARGET[0])/2))
    return east,north

def scan_xy(name,path,size):
    f=ftp();h=hashlib.sha256();n=0;line_no=0;nearest=None;lon_min=lon_max=lat_min=lat_max=None;c50=c250=c1000=0
    try:
        sock=f.transfercmd('RETR '+path)
        stream=sock.makefile('rb')
        try:
            for raw in stream:
                line_no+=1;h.update(raw)
                q=parse_xyz_line(raw)
                if q is None:continue
                lon,lat,_=q;n+=1
                lon_min=lon if lon_min is None else min(lon_min,lon);lon_max=lon if lon_max is None else max(lon_max,lon)
                lat_min=lat if lat_min is None else min(lat_min,lat);lat_max=lat if lat_max is None else max(lat_max,lat)
                d=hav_m(lat,lon)
                if d<=50:c50+=1
                if d<=250:c250+=1
                if d<=1000:c1000+=1
                if nearest is None or d<nearest['distance_m']:
                    nearest={'distance_m':d,'line_number':line_no,'longitude':lon,'latitude':lat}
        finally:
            try:stream.close()
            except Exception:pass
            try:sock.close()
            except Exception:pass
    finally:
        try:f.close()
        except Exception:pass
    d=nearest['distance_m'] if nearest else None
    if d is None:cls='NO_PARSED_XY'
    elif d<=50:cls='EXACT_LOCAL'
    elif d<=250:cls='LOCAL_PATCH'
    elif d<=1000:cls='CORRIDOR_ONLY'
    else:cls='NOT_LOCAL'
    return {'file':name,'path':path,'reported_size_bytes':size,'sha256':h.hexdigest(),'parsed_xy_count':n,'line_count':line_no,'longitude_min_max':[lon_min,lon_max],'latitude_min_max':[lat_min,lat_max],'nearest_xy':nearest,'count_xy_within_50m':c50,'count_xy_within_250m':c250,'count_xy_within_1000m':c1000,'coverage_class':cls}

def pct(xs,p):
    if not xs:return None
    s=sorted(xs);q=(len(s)-1)*p;k=int(math.floor(q));f=q-k
    return s[k] if k==len(s)-1 else s[k]*(1-f)+s[k+1]*f

def medmad(xs):
    if not xs:return None,None
    m=statistics.median(xs);return m,statistics.median(abs(x-m) for x in xs)

def solve3(a,b):
    m=[list(map(float,row))+[float(v)] for row,v in zip(a,b)]
    for col in range(3):
        piv=max(range(col,3),key=lambda r:abs(m[r][col]))
        if abs(m[piv][col])<1e-12:return None
        m[col],m[piv]=m[piv],m[col]
        v=m[col][col]
        for j in range(col,4):m[col][j]/=v
        for r in range(3):
            if r==col:continue
            f=m[r][col]
            for j in range(col,4):m[r][j]-=f*m[col][j]
    return [m[i][3] for i in range(3)]

def fit_plane(points):
    if len(points)<3:return None,None
    n=len(points);sx=sy=sz=sxx=syy=sxy=sxz=syz=0.0
    for x,y,z,_ in points:
        sx+=x;sy+=y;sz+=z;sxx+=x*x;syy+=y*y;sxy+=x*y;sxz+=x*z;syz+=y*z
    coef=solve3([[n,sx,sy],[sx,sxx,sxy],[sy,sxy,syy]],[sz,sxz,syz])
    if coef is None:return None,None
    a,b,c=coef;res=[z-(a+b*x+c*y) for x,y,z,_ in points]
    return {'intercept':a,'east_slope_per_m':b,'north_slope_per_m':c,'slope_magnitude':math.sqrt(b*b+c*c)},res

def scan_depth(name,path,expected_sha):
    f=ftp();h=hashlib.sha256();local=[];control=[];nearest=None;line_no=0
    try:
        sock=f.transfercmd('RETR '+path);stream=sock.makefile('rb')
        try:
            for raw in stream:
                line_no+=1;h.update(raw);q=parse_xyz_line(raw)
                if q is None:continue
                lon,lat,z=q;d=hav_m(lat,lon);e,n=en_m(lat,lon)
                if nearest is None or d<nearest['distance_m']:nearest={'distance_m':d,'line_number':line_no,'longitude':lon,'latitude':lat,'field3_z':z}
                if d<=250:local.append((e,n,z,d))
                elif 500<=d<=1500:control.append((e,n,z,d))
        finally:
            try:stream.close()
            except Exception:pass
            try:sock.close()
            except Exception:pass
    finally:
        try:f.close()
        except Exception:pass
    sha=h.hexdigest()
    if sha!=expected_sha:raise RuntimeError(f'{name}: phase B SHA mismatch')
    lz=[p[2] for p in local];cz=[p[2] for p in control];lm,lmad=medmad(lz);cm,cmad=medmad(cz);plane,res=fit_plane(local);rm,rmad=medmad(res or [])
    return {'file':name,'sha256_reverified':sha,'nearest_sounding':nearest,'local_radius_m':250,'local_n':len(local),'local_depth_median':lm,'local_depth_mad':lmad,'local_depth_p05_p95_relief':None if not lz else pct(lz,0.95)-pct(lz,0.05),'local_depth_p05':pct(lz,0.05),'local_depth_p95':pct(lz,0.95),'local_min_max_depth':None if not lz else [min(lz),max(lz)],'local_linear_plane_fit_coefficients':plane,'local_plane_residual_median':rm,'local_plane_residual_mad':rmad,'local_plane_residual_p05_p95_span':None if not res else pct(res,0.95)-pct(res,0.05),'control_annulus_m':[500,1500],'control_annulus_n':len(control),'control_depth_median':cm,'control_depth_mad':cmad}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
    out={'schema':'janus.cosmos.cousteau.cd169.g7b_em12_ground_cell_probe.v1','status':'STARTED','contract':'JANUS-HANNAH-CD169-G7B-EM12-COVERAGE-CONTRACT-2026-08-26-v1.0','target':{'latitude':TARGET[0],'longitude':TARGET[1]},'phase_A_depth_values_used':False,'depth_based_file_selection':False,'recentered':False,'geometry_classification':False,'pyramid_test':False,'anomaly_threshold':None}
    try:
        candidates=list_candidates();out['candidate_count']=len(candidates);rows=[]
        for name,path,size in candidates:rows.append(scan_xy(name,path,size))
        rows.sort(key=lambda r:(float('inf') if not r['nearest_xy'] else r['nearest_xy']['distance_m'],r['file']))
        winners=[r for r in rows if r['nearest_xy'] and r['nearest_xy']['distance_m']<=250]
        out['phase_A_xy_only']={'files':rows,'coverage_pass':bool(winners),'winning_files_geometry_only':[r['file'] for r in winners],'nearest_overall':rows[0] if rows else None,'coverage_rule':'nearest_xy_distance_m <= 250'}
        if winners:
            byname={name:(path,size) for name,path,size in candidates};depth=[]
            for r in winners:
                path,_=byname[r['file']];depth.append(scan_depth(r['file'],path,r['sha256']))
            out['phase_B_depth_after_coverage_pass']={'triggered':True,'files':depth,'scientific_role':'INDEPENDENT_EARTH_FIXED_BATHYMETRY_AND_LOCAL_TERRAIN_CONSTRAINT_ONLY'}
            out['status']='G7B_EM12_LOCAL_COVERAGE_PASS__BATHYMETRY_CONSTRAINT_READY'
        else:
            out['phase_B_depth_after_coverage_pass']={'triggered':False,'reason':'NO_XY_WITHIN_250M'}
            out['status']='G7B_EM12_LOCAL_COVERAGE_FAIL__DEPTH_NOT_READ'
    except Exception as e:
        out['status']='G7B_EM12_PROBE_FAILED';out['error_type']=type(e).__name__;out['error']=str(e)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    ph=out.get('phase_A_xy_only',{});print(json.dumps({'status':out['status'],'candidate_count':out.get('candidate_count'),'coverage_pass':ph.get('coverage_pass'),'winners':ph.get('winning_files_geometry_only'),'nearest':None if not ph.get('nearest_overall') else {'file':ph['nearest_overall']['file'],'nearest_xy':ph['nearest_overall']['nearest_xy'],'class':ph['nearest_overall']['coverage_class']},'depth':out.get('phase_B_depth_after_coverage_pass')},indent=2));return 0 if out['status'].startswith('G7B_EM12_') else 2
if __name__=='__main__':raise SystemExit(main())
