#!/usr/bin/env python3
from __future__ import annotations
import argparse, ftplib, hashlib, json, math, statistics
from pathlib import Path

HOST='livftp.noc.ac.uk'
ROOT='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/EM12'
TARGETS=[
 {"id":"KN19207_CAND_002","lat":-3.9727527956056825,"lon":-12.272824298723462},
 {"id":"KN19207_CAND_003","lat":-4.015075679897318,"lon":-12.29915403590432},
]
R=6371008.8

def ftp():
    f=ftplib.FTP(timeout=120); f.connect(HOST,21); f.login('anonymous','janus-kusto@example.invalid'); f.voidcmd('TYPE I'); return f

def list_xyz():
    f=ftp()
    try:
        out=[]
        for p in f.nlst(ROOT):
            name=p.rstrip('/').split('/')[-1]
            if name.lower().endswith('.xyz.ascii'):
                full=p if p.startswith('/') else ROOT+'/'+name
                try: sz=f.size(full)
                except Exception: sz=None
                out.append((name,full,sz))
        return sorted(out)
    finally:
        try:f.quit()
        except: f.close()

def parse_xy(raw):
    s=raw.decode('ascii','ignore').strip()
    if not s or s[0] in '#;!': return None
    p=s.replace(',',' ').split()
    if len(p)<2:return None
    try: lon=float(p[0]); lat=float(p[1])
    except:return None
    if not (math.isfinite(lon) and math.isfinite(lat) and -180<=lon<=180 and -90<=lat<=90):return None
    return lon,lat

def parse_xyz(raw):
    s=raw.decode('ascii','ignore').strip()
    if not s or s[0] in '#;!': return None
    p=s.replace(',',' ').split()
    if len(p)<3:return None
    try: lon=float(p[0]); lat=float(p[1]); z=float(p[2])
    except:return None
    if not all(map(math.isfinite,(lon,lat,z))):return None
    return lon,lat,z

def hav(t,lat,lon):
    p1=math.radians(t['lat']); p2=math.radians(lat)
    dp=p2-p1; dl=math.radians(lon-t['lon'])
    a=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R*math.asin(min(1,math.sqrt(a)))

def en(t,lat,lon):
    n=math.radians(lat-t['lat'])*R
    e=math.radians(lon-t['lon'])*R*math.cos(math.radians((lat+t['lat'])/2))
    return e,n

def pct(x,p):
    s=sorted(x)
    if not s:return None
    q=(len(s)-1)*p; i=int(math.floor(q)); f=q-i
    return s[i] if i==len(s)-1 else s[i]*(1-f)+s[i+1]*f

def solve3(A,b):
    m=[list(map(float,r))+[float(v)] for r,v in zip(A,b)]
    for c in range(3):
        piv=max(range(c,3),key=lambda r:abs(m[r][c]))
        if abs(m[piv][c])<1e-12:return None
        m[c],m[piv]=m[piv],m[c]
        v=m[c][c]
        for j in range(c,4):m[c][j]/=v
        for r in range(3):
            if r==c:continue
            f=m[r][c]
            for j in range(c,4):m[r][j]-=f*m[c][j]
    return [m[i][3] for i in range(3)]

def plane(points):
    if len(points)<3:return None,None
    n=len(points); sx=sy=sz=sxx=syy=sxy=sxz=syz=0.0
    for x,y,z in points:
        sx+=x;sy+=y;sz+=z;sxx+=x*x;syy+=y*y;sxy+=x*y;sxz+=x*z;syz+=y*z
    c=solve3([[n,sx,sy],[sx,sxx,sxy],[sy,sxy,syy]],[sz,sxz,syz])
    if c is None:return None,None
    a,b,d=c; res=[z-(a+b*x+d*y) for x,y,z in points]
    return {'intercept':a,'east_slope':b,'north_slope':d,'gradient_magnitude':math.hypot(b,d),'gradient_azimuth_deg':(math.degrees(math.atan2(b,d))+360)%360},res

def phaseA(fileinfo):
    name,path,size=fileinfo
    f=ftp(); h=hashlib.sha256(); per={t['id']:{'nearest':None,'counts':{50:0,100:0,250:0,500:0,1000:0}} for t in TARGETS}; n=0
    try:
        sock=f.transfercmd('RETR '+path); stream=sock.makefile('rb')
        try:
            for raw in stream:
                h.update(raw); q=parse_xy(raw)
                if q is None:continue
                n+=1; lon,lat=q
                for t in TARGETS:
                    d=hav(t,lat,lon); r=per[t['id']]
                    for rad in r['counts']:
                        if d<=rad:r['counts'][rad]+=1
                    if r['nearest'] is None or d<r['nearest']['distance_m']:
                        r['nearest']={'distance_m':d,'lon':lon,'lat':lat}
        finally:
            stream.close(); sock.close()
    finally:f.close()
    return {'file':name,'path':path,'size':size,'sha256':h.hexdigest(),'parsed_xy':n,'targets':per}

def phaseB(t,winners):
    points=[]; files=[]
    for w in winners:
        f=ftp(); h=hashlib.sha256(); name=w['file']; path=w['path']; n=0
        try:
            sock=f.transfercmd('RETR '+path); stream=sock.makefile('rb')
            try:
                for raw in stream:
                    h.update(raw); q=parse_xyz(raw)
                    if q is None:continue
                    lon,lat,z=q; d=hav(t,lat,lon)
                    if d<=1200:
                        e,north=en(t,lat,lon); points.append((e,north,z,d)); n+=1
            finally:stream.close();sock.close()
        finally:f.close()
        if h.hexdigest()!=w['sha256']: raise RuntimeError(name+' sha mismatch')
        files.append({'file':name,'sha256':h.hexdigest(),'points_within_1200m':n})
    out={'files':files,'radii':{}}
    nearest=min(points,key=lambda p:p[3]) if points else None
    out['nearest'] = None if nearest is None else {'distance_m':nearest[3],'depth':nearest[2]}
    for rad in [250,500,1000]:
        pp=[p for p in points if p[3]<=rad]
        z=[p[2] for p in pp]
        pl,res=plane([(p[0],p[1],p[2]) for p in pp])
        out['radii'][str(rad)]={
          'n':len(pp),
          'median_depth':statistics.median(z) if z else None,
          'p05':pct(z,.05),'p95':pct(z,.95),
          'p05_p95_relief':None if not z else pct(z,.95)-pct(z,.05),
          'plane':pl,
          'plane_rmse':None if not res else math.sqrt(sum(x*x for x in res)/len(res))
        }
    return out,points

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True,type=Path); a=ap.parse_args()
    files=list_xyz()
    phase_a=[phaseA(fi) for fi in files]
    result={'artifact_id':'JANUS-KUSTO-CD169-EM12-CAND002003-COVERAGE-DEPTH-RUN-2026-09-22-v1.0','candidate_file_count':len(files),'phase_A_depth_read':False,'targets':{}}
    for t in TARGETS:
        rows=[]
        for r in phase_a:
            tr=r['targets'][t['id']]
            rows.append({'file':r['file'],'path':r['path'],'size':r['size'],'sha256':r['sha256'],'nearest':tr['nearest'],'counts':tr['counts']})
        rows.sort(key=lambda x:1e99 if x['nearest'] is None else x['nearest']['distance_m'])
        winners=[x for x in rows if x['nearest'] and x['nearest']['distance_m']<=1000]
        support=bool(rows and rows[0]['nearest'] and rows[0]['nearest']['distance_m']<=100)
        item={'target':t,'phase_A':{'nearest_overall':rows[0] if rows else None,'winning_files_within_1000m':[x['file'] for x in winners],'support_pass_le100m':support}}
        if support:
            b,_=phaseB(t,winners); item['phase_B']=b
        else:item['phase_B']={'triggered':False}
        result['targets'][t['id']]=item
    result['status']='COMPLETE'
    result['claim_ceiling']='CD169_EM12_2005_MEASURED_SUPPORT_AND_LOCAL_DEPTH_ONLY__CROSS_YEAR_COMPARISON_SEPARATE'
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
