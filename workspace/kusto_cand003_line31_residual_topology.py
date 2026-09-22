#!/usr/bin/env python3
from __future__ import annotations
import concurrent.futures, ftplib, hashlib, json, math, re, struct, time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import requests
import shapely
from shapely.geometry import LineString, Point
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

HOST="livftp.noc.ac.uk"
ROOT="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281"
NAV_PATH=ROOT+"/Nav/cd169leg1_1min.listit"
PRODUCT_PATH=ROOT+"/EM12/B1-81-1_Acceptl28-33.xyz.ascii"
PRODUCT_SHA="0e8cac049b53dd6ba0e5013e8aaea2bfb48347db887edc3c4e3e9f697e96dcf0"
C={"id":"KN19207_CAND_003","lat":-4.015075679897318,"lon":-12.29915403590432}
OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
M=111320.0;LAT0=-4.0;COS0=math.cos(math.radians(LAT0))
LOCAL_R=1000.0;SELECT_R=1600.0;CELL=100.0;MATCH_R=100.0
UA={"User-Agent":"JANUS-KUSTO-CAND003-local-residual-topology/1.0"}

SURVEYS={
 "KN192-07_2008":"https://data.ngdc.noaa.gov/platforms/ocean/ships/knorr/KN192-07/multibeam/data/version1/MB/generated/",
 "KNOX15RR_2008":"https://data.ngdc.noaa.gov/platforms/ocean/ships/roger_revelle/KNOX15RR/multibeam/data/version1/MB/generated/"
}
WINDOWS={
 28:("2005-02-25T07:30:00Z","2005-02-25T21:10:00Z"),
 29:("2005-02-25T21:10:00Z","2005-02-26T07:34:00Z"),
 30:("2005-02-26T08:22:00Z","2005-02-26T09:35:00Z"),
 31:("2005-02-26T19:11:00Z","2005-02-27T14:40:00Z"),
 32:("2005-02-27T14:40:00Z","2005-02-28T07:31:00Z"),
 33:("2005-02-28T07:31:00Z","2005-02-28T14:59:00Z")
}

def xy(lon,lat):return lon*M*COS0,lat*M
def ll(x,y):return y/M,x/(M*COS0)
def parse_iso(s):return datetime.fromisoformat(s.replace("Z","+00:00")).timestamp()
CX,CY=xy(C["lon"],C["lat"])

def ftp_fetch(path):
    f=ftplib.FTP(timeout=120);f.connect(HOST,21);f.login("anonymous","janus-kusto@example.invalid");f.voidcmd("TYPE I");chunks=[]
    try:f.retrbinary("RETR "+path,chunks.append)
    finally:
        try:f.quit()
        except:f.close()
    return b"".join(chunks)

def fetch(url,timeout=120):
    last=None
    for attempt in range(7):
        try:
            r=requests.get(url,headers=UA,timeout=timeout)
            if r.status_code==429:
                last=RuntimeError(f"429 {url}");time.sleep(min(2**attempt,32));continue
            r.raise_for_status();return r.content
        except requests.exceptions.RequestException as e:
            last=e
            if attempt==6:raise
            time.sleep(min(2**attempt,32))
    raise last

def list_files(base,ext):
    html=fetch(base,60).decode("latin1","replace")
    return sorted(set(re.findall(r'href="([^"]+\.'+re.escape(ext)+r')"',html,re.I)))

def parse_fnv(text):
    out=[]
    for line in text.splitlines():
        p=line.split()
        if len(p)<19:continue
        try:out.append((float(p[15]),float(p[16]),float(p[17]),float(p[18])))
        except:pass
    return out

def segdist(px,py,ax,ay,bx,by):
    vx,vy=bx-ax,by-ay;wx,wy=px-ax,py-ay;den=vx*vx+vy*vy
    t=0 if den==0 else max(0,min(1,(wx*vx+wy*vy)/den))
    return math.hypot(px-(ax+t*vx),py-(ay+t*vy))

def fnv_min(rows):
    best=float("inf")
    for plon,plat,slon,slat in rows:
        ax,ay=xy(plon,plat);bx,by=xy(slon,slat)
        best=min(best,segdist(CX,CY,ax,ay,bx,by))
    return best

def parse_fbt(data,name):
    off=0;chunks=[];nrec=0
    ids={b"V4":(4,90),b"V5":(5,98),b"cc":("comment",2),b"##":("oldcomment",30)}
    while off+2<=len(data):
        tag=data[off:off+2]
        if tag not in ids:raise RuntimeError(f"{name} bad tag {tag!r} @{off}")
        ver,hs=ids[tag]
        if isinstance(ver,str):off+=hs+128;continue
        h=data[off:off+hs]
        t,lon,lat,sd,alt=struct.unpack_from(">5d",h,2)
        heading,speed,roll,pitch,heave,bxw,blw=struct.unpack_from(">7f",h,42)
        if ver==4:
            nb,na,ns,head=struct.unpack_from(">4h",h,70);dscale,xscale=struct.unpack_from(">2f",h,78)
        else:
            nb,na,ns,head=struct.unpack_from(">4i",h,70);dscale,xscale=struct.unpack_from(">2f",h,86)
        p=off+hs
        flags=np.frombuffer(data,dtype=np.uint8,count=nb,offset=p).copy();p+=nb
        def arr(n):
            nonlocal p
            a=np.frombuffer(data,dtype=">i2",count=n,offset=p).astype(float);p+=2*n;return a
        bath=arr(nb);across=arr(nb);along=arr(nb);p+=2*na+6*ns
        gx,gy=xy(lon,lat);hr=math.radians(heading);sh=math.sin(hr);ch=math.cos(hr)
        xt=xscale*across;lt=xscale*along
        xx=gx+lt*sh+xt*ch;yy=gy+lt*ch-xt*sh;zz=dscale*bath+sd
        good=np.where(flags==0)[0]
        if len(good):chunks.append(np.column_stack([xx[good],yy[good],zz[good]]))
        off=p;nrec+=1
    return (np.vstack(chunks) if chunks else np.empty((0,3))),nrec

# line geometry
navb=ftp_fetch(NAV_PATH);base=datetime(2005,1,1,tzinfo=timezone.utc);nav=[]
for raw in navb.decode("ascii","ignore").splitlines():
    p=raw.split()
    if len(p)<7:continue
    try:
        jd=int(p[1]);hh,mm,ss=map(int,p[2].split(":"));lat=float(p[3]);lon=float(p[5])
    except:continue
    dt=base+timedelta(days=jd-1,hours=hh,minutes=mm,seconds=ss)
    nav.append((dt.timestamp(),)+xy(lon,lat))
line_geom={}
for line,(a,b) in WINDOWS.items():
    ta,tb=parse_iso(a),parse_iso(b)
    coords=[(x,y) for t,x,y in nav if ta<=t<=tb]
    if len(coords)<2:raise RuntimeError(f"line {line} insufficient nav")
    line_geom[line]=LineString(coords)

# 2005 line31 points
prod=ftp_fetch(PRODUCT_PATH);sha=hashlib.sha256(prod).hexdigest()
if sha!=PRODUCT_SHA:raise RuntimeError("product SHA mismatch")
xs=[];ys=[];zs=[]
for raw in prod.decode("ascii","ignore").splitlines():
    p=raw.strip().replace(","," ").split()
    if len(p)<3:continue
    try:lon=float(p[0]);lat=float(p[1]);z=float(p[2])
    except:continue
    if not all(map(math.isfinite,(lon,lat,z))):continue
    x,y=xy(lon,lat);xs.append(x);ys.append(y);zs.append(z)
x=np.asarray(xs,float);y=np.asarray(ys,float);z=np.asarray(zs,float)
pts=shapely.points(x,y);lines=sorted(WINDOWS);dist=np.empty((len(lines),len(x)),float)
for j,line in enumerate(lines):dist[j,:]=shapely.distance(pts,line_geom[line])
assigned=np.asarray([lines[i] for i in np.argmin(dist,axis=0)],dtype=int)
line31=np.column_stack([x[assigned==31],y[assigned==31],z[assigned==31]])

def thin_local(p3):
    rr=np.hypot(p3[:,0]-CX,p3[:,1]-CY);p=p3[rr<=LOCAL_R];cells={}
    for row in p:
        dx=row[0]-CX;dy=row[1]-CY
        ix=math.floor((dx+LOCAL_R)/CELL);iy=math.floor((dy+LOCAL_R)/CELL)
        xc=-LOCAL_R+(ix+.5)*CELL;yc=-LOCAL_R+(iy+.5)*CELL
        if math.hypot(xc,yc)>LOCAL_R:continue
        g=(dx-xc)**2+(dy-yc)**2;cand=(g,float(row[0]),float(row[1]),float(row[2]))
        if (ix,iy) not in cells or cand[:3]<cells[(ix,iy)][:3]:cells[(ix,iy)]=cand
    return np.asarray([(v[1],v[2],v[3]) for k,v in sorted(cells.items())],float)

ref=thin_local(line31)

# 2008 data near frozen target
loaded={};manifest={}
for lab,baseurl in SURVEYS.items():
    names=list_files(baseurl,"fnv");selected=[];fnvman=[]
    def gf(n):
        b=fetch(urljoin(baseurl,n),90);rows=parse_fnv(b.decode("utf-8","replace"));return n,hashlib.sha256(b).hexdigest(),rows
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
        for n,h,rows in (f.result() for f in concurrent.futures.as_completed([ex.submit(gf,n) for n in names])):
            d=fnv_min(rows);fnvman.append({"file":n,"sha256":h,"min_target_m":d})
            if d<=SELECT_R:selected.append(n[:-4]+".fbt")
    selected=sorted(set(selected));parts=[];fbtman=[]
    def gb(n):
        b=fetch(urljoin(baseurl,n),180);p,nrec=parse_fbt(b,n);return n,hashlib.sha256(b).hexdigest(),len(b),p,nrec
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        for n,h,size,p,nrec in (f.result() for f in concurrent.futures.as_completed([ex.submit(gb,n) for n in selected])):
            parts.append(p);fbtman.append({"file":n,"sha256":h,"bytes":size,"records":nrec,"good_beams":int(len(p))})
    allp=np.vstack(parts) if parts else np.empty((0,3))
    if len(allp):allp=allp[np.hypot(allp[:,0]-CX,allp[:,1]-CY)<=SELECT_R]
    loaded[lab]=allp;manifest[lab]={"selected_fbt":selected,"fnv_manifest":sorted(fnvman,key=lambda q:q["file"]),"fbt_manifest":sorted(fbtman,key=lambda q:q["file"])}

trees={k:cKDTree(v[:,:2]) for k,v in loaded.items()}
joint=[]
for row in ref:
    vals={}
    ok=True
    for lab,p in loaded.items():
        idx=trees[lab].query_ball_point(row[:2],MATCH_R)
        if not idx:ok=False;break
        vals[lab]=float(np.median(p[np.asarray(idx,dtype=int),2]))
    if ok:joint.append((row[0],row[1],row[2],vals["KN192-07_2008"],vals["KNOX15RR_2008"]))
J=np.asarray(joint,float)
if len(J)<30:raise RuntimeError(f"insufficient joint support {len(J)}")

xr=J[:,0]-CX;yr=J[:,1]-CY
z05=J[:,2];zkn=J[:,3];zkx=J[:,4];zcon=np.median(np.column_stack([zkn,zkx]),axis=1)

def resid(z):
    A=np.column_stack([np.ones(len(z)),xr,yr]);coef=np.linalg.lstsq(A,z,rcond=None)[0];return coef,z-A@coef
c05,r05=resid(z05);ckn,rkn=resid(zkn);ckx,rkx=resid(zkx);cco,rco=resid(zcon)
delta=r05-rco

def corr(a,b):
    if len(a)<3 or np.std(a)==0 or np.std(b)==0:return {"pearson":None,"spearman":None}
    return {"pearson":float(np.corrcoef(a,b)[0,1]),"spearman":float(spearmanr(a,b).statistic)}
def trip(a,b):
    d=corr(a,b);d["sign_concordance"]=float(np.mean(np.sign(a)==np.sign(b)));return d
def rmse(a):return float(np.sqrt(np.mean(np.asarray(a)**2)))

line31g=line_geom[31];target_proj=float(line31g.project(Point(CX,CY)))
cross=np.asarray([Point(xx,yy).distance(line31g) for xx,yy in J[:,:2]],float)
along=np.asarray([line31g.project(Point(xx,yy))-target_proj for xx,yy in J[:,:2]],float)

A=np.column_stack([np.ones(len(delta)),xr,yr]);coef=np.linalg.lstsq(A,delta,rcond=None)[0];pred=A@coef
tss=float(np.sum((delta-np.mean(delta))**2));rss=float(np.sum((delta-pred)**2));r2=None if tss==0 else 1-rss/tss

quads={}
for name,mask in {
 "NE":(xr>=0)&(yr>=0),"NW":(xr<0)&(yr>=0),"SE":(xr>=0)&(yr<0),"SW":(xr<0)&(yr<0)
}.items():
    quads[name]={"n":int(np.sum(mask)),"median_delta_m":None if not np.any(mask) else float(np.median(delta[mask])),"median_abs_delta_m":None if not np.any(mask) else float(np.median(np.abs(delta[mask])))}

order=np.argsort(-np.abs(delta))[:10];top=[]
for i in order:
    la,lo=ll(J[i,0],J[i,1])
    top.append({"lat":float(la),"lon":float(lo),"delta_m":float(delta[i]),"abs_delta_m":float(abs(delta[i])),"cross_track_m":float(cross[i]),"along_track_relative_m":float(along[i]),"z2005":float(z05[i]),"z2008_consensus":float(zcon[i])})

metrics={
 "joint_point_count":int(len(J)),
 "residual_2005_vs_consensus":trip(r05,rco),
 "residual_KN192_vs_KNOX":trip(rkn,rkx),
 "rmse_2005_minus_consensus_residual_m":rmse(delta),
 "rmse_KN192_minus_KNOX_residual_m":rmse(rkn-rkx),
 "rmse_ratio_2005_error_over_2008_internal":rmse(delta)/rmse(rkn-rkx) if rmse(rkn-rkx)>0 else None,
 "delta_vs_cross_track":corr(delta,cross),
 "abs_delta_vs_cross_track":corr(np.abs(delta),cross),
 "delta_vs_along_track":corr(delta,along),
 "abs_delta_vs_along_track":corr(np.abs(delta),along),
 "affine_delta_surface":{"coef_intercept_dx_dy":[float(v) for v in coef],"R2":r2},
 "quadrants":quads,
 "top_10_absolute_delta":top,
 "cross_track_m":{"min":float(np.min(cross)),"median":float(np.median(cross)),"max":float(np.max(cross))},
 "along_track_relative_m":{"min":float(np.min(along)),"median":float(np.median(along)),"max":float(np.max(along))}
}
out={
 "artifact_id":"JANUS-KUSTO-CAND003-LINE31-LOCAL-RESIDUAL-TOPOLOGY-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CAND003-LINE31-LOCAL-RESIDUAL-TOPOLOGY-PREREG-2026-09-23-v1.0.json",
 "target":C,
 "support":{"raw_line31_within_1000m":int(np.sum(np.hypot(line31[:,0]-CX,line31[:,1]-CY)<=LOCAL_R)),"thinned_2005":int(len(ref)),"joint_points":int(len(J))},
 "metrics":metrics,
 "manifests":manifest,
 "interpretation":"Descriptive local error topology only; no causal mechanism is promoted by this run.",
 "claim_ceiling":"LOCAL_SPATIAL_TOPOLOGY_OF_CAND003_LINE31_2005_VS_2008_DISCORDANCE"
}
p=OUT/"JANUS-KUSTO-CAND003-LINE31-LOCAL-RESIDUAL-TOPOLOGY-RUN-2026-09-23-v1.0.json";p.write_text(json.dumps(out,indent=2))
print(json.dumps({"support":out["support"],"metrics":metrics},indent=2))
