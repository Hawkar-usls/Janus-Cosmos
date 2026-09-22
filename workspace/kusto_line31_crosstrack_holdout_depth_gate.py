#!/usr/bin/env python3
from __future__ import annotations
import concurrent.futures, ftplib, hashlib, json, math, re, struct, time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import requests
import shapely
from shapely.geometry import LineString
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

HOST="livftp.noc.ac.uk"
ROOT="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281"
NAV_PATH=ROOT+"/Nav/cd169leg1_1min.listit"
PRODUCT_PATH=ROOT+"/EM12/B1-81-1_Acceptl28-33.xyz.ascii"
PRODUCT_SHA="0e8cac049b53dd6ba0e5013e8aaea2bfb48347db887edc3c4e3e9f697e96dcf0"
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
M=111320.0; LAT0=-4.0; COS0=math.cos(math.radians(LAT0))
UA={"User-Agent":"JANUS-KUSTO-line31-holdout-depth-gate/1.0"}

HOLDOUTS={
 "CENTER_1":{"band":"CENTER","lat":-4.0467114,"lon":-12.2704123},
 "CENTER_2":{"band":"CENTER","lat":-4.1151018,"lon":-12.2740956},
 "CENTER_3":{"band":"CENTER","lat":-4.0084801,"lon":-12.2137473},
 "EDGE_1":{"band":"EDGE","lat":-3.9956817,"lon":-12.2701268},
 "EDGE_2":{"band":"EDGE","lat":-4.1394583,"lon":-12.2434622},
 "EDGE_3":{"band":"EDGE","lat":-4.110448,"lon":-12.2339197}
}
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
SELECT_R=1600.0; LOCAL_R=1000.0; CELL=100.0; MATCH_R=100.0

def xy(lon,lat): return lon*M*COS0,lat*M
def parse_iso(s): return datetime.fromisoformat(s.replace("Z","+00:00")).timestamp()

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
    rows=[]
    for line in text.splitlines():
        p=line.split()
        if len(p)<19:continue
        try:rows.append((float(p[15]),float(p[16]),float(p[17]),float(p[18])))
        except:pass
    return rows

def segdist(px,py,ax,ay,bx,by):
    vx,vy=bx-ax,by-ay;wx,wy=px-ax,py-ay;den=vx*vx+vy*vy
    t=0 if den==0 else max(0,min(1,(wx*vx+wy*vy)/den))
    return math.hypot(px-(ax+t*vx),py-(ay+t*vy))

def fnv_min_any_holdout(rows):
    best=float("inf")
    targets=[xy(h["lon"],h["lat"]) for h in HOLDOUTS.values()]
    for plon,plat,slon,slat in rows:
        ax,ay=xy(plon,plat);bx,by=xy(slon,slat)
        for px,py in targets:
            d=segdist(px,py,ax,ay,bx,by)
            if d<best:best=d
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

# Navigation / line geometries.
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

# 2005 product, now depth is allowed because holdouts are frozen.
prod=ftp_fetch(PRODUCT_PATH);sha=hashlib.sha256(prod).hexdigest()
if sha!=PRODUCT_SHA:raise RuntimeError(f"product SHA mismatch {sha}")
xs=[];ys=[];zs=[]
for raw in prod.decode("ascii","ignore").splitlines():
    s=raw.strip()
    if not s or s[0] in "#;!":continue
    p=s.replace(","," ").split()
    if len(p)<3:continue
    try:lon=float(p[0]);lat=float(p[1]);z=float(p[2])
    except:continue
    if not all(map(math.isfinite,(lon,lat,z))):continue
    x,y=xy(lon,lat);xs.append(x);ys.append(y);zs.append(z)
x=np.asarray(xs,float);y=np.asarray(ys,float);z=np.asarray(zs,float)
points=shapely.points(x,y);lines=sorted(WINDOWS)
dist=np.empty((len(lines),len(x)),float)
for j,line in enumerate(lines):dist[j,:]=shapely.distance(points,line_geom[line])
assigned=np.asarray([lines[i] for i in np.argmin(dist,axis=0)],dtype=int)
line31=np.column_stack([x[assigned==31],y[assigned==31],z[assigned==31]])

# Load only 2008 files whose FNV swath geometry comes near any frozen holdout.
loaded={};manifests={}
for label,baseurl in SURVEYS.items():
    names=list_files(baseurl,"fnv");selected=[];fnvman=[]
    def getfnv(name):
        b=fetch(urljoin(baseurl,name),90);rows=parse_fnv(b.decode("utf-8","replace"))
        return name,hashlib.sha256(b).hexdigest(),rows
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
        futs=[ex.submit(getfnv,n) for n in names]
        for fut in concurrent.futures.as_completed(futs):
            name,h,rows=fut.result();d=fnv_min_any_holdout(rows);fnvman.append({"file":name,"sha256":h,"min_any_holdout_m":d})
            if d<=SELECT_R:selected.append(name[:-4]+".fbt")
    selected=sorted(set(selected));parts=[];fbtman=[]
    def getfbt(name):
        b=fetch(urljoin(baseurl,name),180);pts,nrec=parse_fbt(b,name)
        return name,hashlib.sha256(b).hexdigest(),len(b),pts,nrec
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        futs=[ex.submit(getfbt,n) for n in selected]
        for fut in concurrent.futures.as_completed(futs):
            name,h,size,pts,nrec=fut.result();parts.append(pts);fbtman.append({"file":name,"sha256":h,"bytes":size,"records":nrec,"good_beams":int(len(pts))})
    allpts=np.vstack(parts) if parts else np.empty((0,3))
    # Keep only points near at least one holdout to reduce later tree size.
    if len(allpts):
        keep=np.zeros(len(allpts),dtype=bool)
        for h in HOLDOUTS.values():
            cx,cy=xy(h["lon"],h["lat"]);keep|=(np.hypot(allpts[:,0]-cx,allpts[:,1]-cy)<=SELECT_R)
        allpts=allpts[keep]
    loaded[label]=allpts
    manifests[label]={"selected_fbt":selected,"fnv_manifest":sorted(fnvman,key=lambda q:q["file"]),"fbt_manifest":sorted(fbtman,key=lambda q:q["file"])}

def thin_local(points3,cx,cy):
    rr=np.hypot(points3[:,0]-cx,points3[:,1]-cy);p=points3[rr<=LOCAL_R]
    cells={}
    for row in p:
        dx=row[0]-cx;dy=row[1]-cy
        ix=math.floor((dx+LOCAL_R)/CELL);iy=math.floor((dy+LOCAL_R)/CELL)
        xc=-LOCAL_R+(ix+.5)*CELL;yc=-LOCAL_R+(iy+.5)*CELL
        if math.hypot(xc,yc)>LOCAL_R:continue
        g=(dx-xc)**2+(dy-yc)**2
        cand=(g,float(row[0]),float(row[1]),float(row[2]))
        if (ix,iy) not in cells or cand[:3]<cells[(ix,iy)][:3]:cells[(ix,iy)]=cand
    return np.asarray([(v[1],v[2],v[3]) for k,v in sorted(cells.items())],float) if cells else np.empty((0,3))

def fit(x,y,z):
    A=np.column_stack([np.ones(len(x)),x,y]);coef=np.linalg.lstsq(A,z,rcond=None)[0];return coef,z-A@coef

def compare(ref3,target3,cx,cy):
    out={}
    if len(ref3)<30:return {"usable":False,"gate_pass":False,"verdict":"INSUFFICIENT_2005_THINNED_SUPPORT","matched_point_count":0}
    tree=cKDTree(target3[:,:2]) if len(target3) else None
    rr=[];tz=[]
    if tree is not None:
        for row in ref3:
            idx=tree.query_ball_point(row[:2],MATCH_R)
            if idx:
                rr.append(row);tz.append(float(np.median(target3[np.asarray(idx,dtype=int),2])))
    if len(rr)<30:return {"usable":False,"gate_pass":False,"verdict":"INSUFFICIENT_MATCHED_SUPPORT","matched_point_count":len(rr)}
    ref=np.asarray(rr,float);tar=np.asarray(tz,float)
    xx=ref[:,0]-cx;yy=ref[:,1]-cy;zr=ref[:,2]
    _,a=fit(xx,yy,zr);_,b=fit(xx,yy,tar)
    pear=float(np.corrcoef(a,b)[0,1]);spear=float(spearmanr(a,b).statistic);sign=float(np.mean(np.sign(a)==np.sign(b)))
    p05r,p95r=np.percentile(zr,[5,95]);p05t,p95t=np.percentile(tar,[5,95])
    gate=pear>=.60 and spear>=.60 and sign>=.65
    return {
      "usable":True,"gate_pass":bool(gate),"matched_point_count":len(ref),
      "residual_pearson":pear,"residual_spearman":spear,"residual_sign_concordance":sign,
      "median_depth_offset_target_minus_2005_m":float(np.median(tar-zr)),
      "relief_2005_m":float(p95r-p05r),"relief_2008_m":float(p95t-p05t),
      "relief_ratio_2008_over_2005":float((p95t-p05t)/(p95r-p05r)) if p95r!=p05r else None,
      "verdict":"PASS_LOCAL_TEMPORAL_MORPHOLOGY" if gate else "FAIL_LOCAL_TEMPORAL_MORPHOLOGY"
    }

results={}
for name,h in HOLDOUTS.items():
    cx,cy=xy(h["lon"],h["lat"])
    ref=thin_local(line31,cx,cy)
    comps={}
    for lab,pts in loaded.items():
        local=pts[np.hypot(pts[:,0]-cx,pts[:,1]-cy)<=SELECT_R] if len(pts) else pts
        comps[lab]=compare(ref,local,cx,cy)
    usable=len(ref)>=30 and all(c["usable"] for c in comps.values())
    passboth=usable and all(c["gate_pass"] for c in comps.values())
    results[name]={
      "band":h["band"],"lat":h["lat"],"lon":h["lon"],
      "raw_line31_2005_within_1000m":int(np.sum(np.hypot(line31[:,0]-cx,line31[:,1]-cy)<=LOCAL_R)),
      "thinned_2005_support":int(len(ref)),
      "comparisons":comps,
      "usable":bool(usable),
      "pass_both_2008_lineages":bool(passboth)
    }

agg={}
for band in ["CENTER","EDGE"]:
    rr=[r for r in results.values() if r["band"]==band and r["usable"]]
    agg[band]={"usable":len(rr),"pass":sum(r["pass_both_2008_lineages"] for r in rr),"fail":sum(not r["pass_both_2008_lineages"] for r in rr)}
if agg["CENTER"]["usable"]<2 or agg["EDGE"]["usable"]<2:
    verdict="INSUFFICIENT_HOLDOUT_SUPPORT"
elif agg["CENTER"]["pass"]>=2 and agg["EDGE"]["fail"]>=2:
    verdict="SUPPORT_EDGE_SPECIFIC_DEGRADATION"
else:
    verdict="NO_SUPPORT_FOR_EDGE_SPECIFIC_DEGRADATION"

out={
 "artifact_id":"JANUS-KUSTO-CD169-EM12-LINE31-CROSSTRACK-HOLDOUT-DEPTH-GATE-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CD169-EM12-LINE31-CROSSTRACK-HOLDOUT-DEPTH-GATE-PREREG-2026-09-23-v1.0.json",
 "processed_product":{"sha256":sha,"line31_points":int(len(line31))},
 "results":results,
 "aggregate":agg,
 "verdict":verdict,
 "thresholds_retuned":False,
 "holdouts_recentered":False,
 "manifests":manifests,
 "claim_ceiling":"HOLDOUT_REPLICATION_OF_LINE31_CROSSTRACK_DEPENDENT_TEMPORAL_REPRODUCIBILITY"
}
p=OUT/"JANUS-KUSTO-CD169-EM12-LINE31-CROSSTRACK-HOLDOUT-DEPTH-GATE-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({"results":results,"aggregate":agg,"verdict":verdict},indent=2))
