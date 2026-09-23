#!/usr/bin/env python3
from __future__ import annotations
import concurrent.futures, json, math, re, struct, time
from pathlib import Path
from urllib.parse import urljoin, urlparse

import numpy as np
import requests
from scipy.stats import spearmanr
from shapely.geometry import Point, LineString
from shapely.ops import unary_union
from shapely.prepared import prep

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace/kusto_global_anomaly_out"; OUT.mkdir(parents=True,exist_ok=True)
PR=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-GLOBAL-STAGE2-HELDOUT-CROSSSURVEY-REPLAY-PREREG-2026-09-23-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-global-stage2-heldout-replay/1.0"}
M=111320.0
LAT0=float(np.mean([c["lat"] for c in PR["candidates"]]))
COS0=math.cos(math.radians(LAT0))
EXTREME=5500.0
SURVEYS={
 "FK140307":"https://data.ngdc.noaa.gov/platforms/ocean/ships/falkor/FK140307/multibeam/data/",
 "ZHNG09RR":"https://data.ngdc.noaa.gov/platforms/ocean/ships/roger_revelle/ZHNG09RR/multibeam/data/"
}

def get(u):
    last=None
    for a in range(7):
        try:
            r=requests.get(u,headers=UA,timeout=120)
            if r.status_code==429:
                time.sleep(min(30,2**a));continue
            r.raise_for_status();return r
        except Exception as e:
            last=e
            if a==6:raise
            time.sleep(min(30,2**a))
    raise last

def hrefs(u):
    return sorted(set(urljoin(u,h) for h in re.findall(r'href=["\']([^"\']+)["\']',get(u).text,re.I)))

def list_fnv(base,maxdepth=4):
    seen=set();out=[]
    def walk(u,d):
        if u in seen or d>maxdepth:return
        seen.add(u)
        for x in hrefs(u):
            if not x.startswith(base) or x.rstrip("/")==u.rstrip("/"):continue
            if x.endswith("/"):walk(x,d+1)
            elif x.lower().endswith(".fnv"):out.append(x)
    walk(base,0);return sorted(set(out))

def xy(lon,lat):return lon*M*COS0,lat*M
def ll(x,y):return x/(M*COS0),y/M

def parse_fnv(text):
    rows=[]
    for line in text.splitlines():
        p=line.split()
        if len(p)<19:continue
        try:rows.append((float(p[15]),float(p[16]),float(p[17]),float(p[18])))
        except:pass
    return rows

# Candidate buffer union is file-selection only, padded 25%.
bufs=[]
for c in PR["candidates"]:
    deg=1.25*c["radius_m"]/M
    bufs.append(Point(c["lon"],c["lat"]).buffer(deg))
selgeom=unary_union(bufs);psel=prep(selgeom)

def scan_fnv(u):
    rows=parse_fnv(get(u).text)
    for plon,plat,slon,slat in rows:
        if psel.intersects(LineString([(plon,plat),(slon,slat)])):
            return u,True,len(rows)
    return u,False,len(rows)

def fbt_from_fnv(u):return u[:-4]+".fbt"

def parse_fbt(data,name):
    off=0;parts=[];recs=0
    ids={b"V4":(4,90),b"V5":(5,98),b"cc":("comment",2),b"##":("oldcomment",30)}
    minx,miny,maxx,maxy=selgeom.bounds
    while off+2<=len(data):
        tag=data[off:off+2]
        if tag not in ids:raise RuntimeError(f"{name}: bad tag {tag!r}@{off}")
        ver,hs=ids[tag]
        if isinstance(ver,str):off+=hs+128;continue
        h=data[off:off+hs]
        t,lon,lat,sd,alt=struct.unpack_from(">5d",h,2)
        heading,speed,roll,pitch,heave,bx,bl=struct.unpack_from(">7f",h,42)
        if ver==4:
            nb,na,ns,shd=struct.unpack_from(">4h",h,70);dscale,xscale=struct.unpack_from(">2f",h,78)
        else:
            nb,na,ns,shd=struct.unpack_from(">4i",h,70);dscale,xscale=struct.unpack_from(">2f",h,86)
        p=off+hs
        flags=np.frombuffer(data,dtype=np.uint8,count=nb,offset=p).copy();p+=nb
        def arr(n):
            nonlocal p
            a=np.frombuffer(data,dtype=">i2",count=n,offset=p).astype(float);p+=2*n;return a
        bath=arr(nb);across=arr(nb);along=arr(nb);p+=2*na+6*ns
        navx,navy=xy(lon,lat);hr=math.radians(heading);sh=math.sin(hr);ch=math.cos(hr)
        xt=xscale*across;lt=xscale*along
        east=navx+lt*sh+xt*ch;north=navy+lt*ch-xt*sh;depth=dscale*bath+sd
        good=np.where(flags==0)[0]
        rows=[]
        for i in good:
            lonb,latb=ll(float(east[i]),float(north[i]))
            if lonb<minx or lonb>maxx or latb<miny or latb>maxy:continue
            if not psel.intersects(Point(lonb,latb)):continue
            rows.append((float(east[i]),float(north[i]),float(depth[i]),abs(float(xt[i]))))
        if rows:parts.append(np.asarray(rows,float))
        off=p;recs+=1
    return np.vstack(parts) if parts else np.empty((0,4)),recs

def load_survey(sid):
    fnvs=list_fnv(SURVEYS[sid]);hits=[];meta=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as ex:
        futs=[ex.submit(scan_fnv,u) for u in fnvs]
        for k,f in enumerate(concurrent.futures.as_completed(futs),1):
            u,hit,n=f.result()
            if hit:hits.append(u);meta.append({"fnv":u,"rows":n})
            if k%200==0:print(sid,"FNV",k,"/",len(fnvs),flush=True)
    hits=sorted(hits);parts=[];manifest=[]
    def one(u):
        fu=fbt_from_fnv(u);b=get(fu).content;pts,recs=parse_fbt(b,fu.rsplit("/",1)[-1]);return fu,len(b),pts,recs
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
        futs=[ex.submit(one,u) for u in hits]
        for k,f in enumerate(concurrent.futures.as_completed(futs),1):
            fu,n,pts,recs=f.result();parts.append(pts);manifest.append({"fbt":fu,"bytes":n,"records":recs,"selected_points":int(len(pts))})
            if k%20==0:print(sid,"FBT",k,"/",len(hits),flush=True)
    return (np.vstack(parts) if parts else np.empty((0,4))),{"fnv_total":len(fnvs),"fnv_hits":len(hits),"fbt_manifest":manifest}

def grid_cells(points,c,cell):
    cx,cy=xy(c["lon"],c["lat"]);r=c["radius_m"]
    dx=points[:,0]-cx;dy=points[:,1]-cy;rr=np.hypot(dx,dy)
    p=points[rr<=r]
    d={}
    for row in p:
        ix=math.floor((row[0]-cx+r)/cell);iy=math.floor((row[1]-cy+r)/cell)
        xc=-r+(ix+.5)*cell;yc=-r+(iy+.5)*cell
        if math.hypot(xc,yc)>r:continue
        d.setdefault((ix,iy),[]).append(float(row[2]))
    return {k:{"x":-r+(k[0]+.5)*cell,"y":-r+(k[1]+.5)*cell,"z":float(np.median(v)),"n":len(v)} for k,v in d.items()}

def fit_res(x,y,z):
    A=np.column_stack([np.ones(len(x)),x,y]);coef=np.linalg.lstsq(A,z,rcond=None)[0];return z-A@coef

print("Loading discovery FK140307",flush=True)
A,manA=load_survey("FK140307")
print("Loading HELDOUT ZHNG09RR -- first depth read",flush=True)
B,manB=load_survey("ZHNG09RR")
print("points",len(A),len(B),flush=True)

results=[]
for c in PR["candidates"]:
    cx,cy=xy(c["lon"],c["lat"]);r=c["radius_m"];cell=max(25.0,r/10.0)
    ra=np.hypot(A[:,0]-cx,A[:,1]-cy);rb=np.hypot(B[:,0]-cx,B[:,1]-cy)
    pa=A[ra<=r];pb=B[rb<=r]
    fa=float(np.mean(pa[:,3]>=EXTREME)) if len(pa) else None
    fb=float(np.mean(pb[:,3]>=EXTREME)) if len(pb) else None
    rec={**c,"cell_m":cell,"n_good_beams_discovery":int(len(pa)),"n_good_beams_heldout":int(len(pb)),
         "extreme_beam_fraction_discovery":fa,"extreme_beam_fraction_heldout":fb}
    if len(pa)<100 or len(pb)<100:
        rec.update({"verdict":"INSUFFICIENT_COMMON_SUPPORT","common_cells":0});results.append(rec);continue
    if fb is not None and fb>0.50:
        rec.update({"verdict":"NUISANCE_DOMINATED__NO_PROMOTION","common_cells":0});results.append(rec);continue
    ga=grid_cells(A,c,cell);gb=grid_cells(B,c,cell);common=sorted(set(ga)&set(gb))
    rec["common_cells"]=len(common)
    if len(common)<30:
        rec["verdict"]="INSUFFICIENT_COMMON_SUPPORT";results.append(rec);continue
    x=np.array([ga[k]["x"] for k in common]);y=np.array([ga[k]["y"] for k in common])
    za=np.array([ga[k]["z"] for k in common]);zb=np.array([gb[k]["z"] for k in common])
    resa=fit_res(x,y,za);resb=fit_res(x,y,zb)
    pear=float(np.corrcoef(resa,resb)[0,1]) if np.std(resa)>0 and np.std(resb)>0 else float("nan")
    spear=float(spearmanr(resa,resb).statistic)
    sign=float(np.mean(np.sign(resa)==np.sign(resb)))
    p05a,p95a=np.percentile(za,[5,95]);p05b,p95b=np.percentile(zb,[5,95])
    gate=bool(np.isfinite(pear) and np.isfinite(spear) and pear>=.60 and spear>=.60 and sign>=.65)
    rec.update({
      "residual_pearson":pear,"residual_spearman":spear,"residual_sign_concordance":sign,
      "median_depth_offset_heldout_minus_discovery_m":float(np.median(zb-za)),
      "relief_discovery_m":float(p95a-p05a),"relief_heldout_m":float(p95b-p05b),
      "verdict":"CROSS_SURVEY_PERSISTENT" if gate else "CROSS_SURVEY_FAIL"
    })
    results.append(rec)

out={
 "artifact_id":"JANUS-KUSTO-GLOBAL-STAGE2-HELDOUT-CROSSSURVEY-REPLAY-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-GLOBAL-STAGE2-HELDOUT-CROSSSURVEY-REPLAY-PREREG-2026-09-23-v1.0.json",
 "discovery_survey":"FK140307","heldout_survey":"ZHNG09RR",
 "heldout_depth_first_read_in_this_run":True,
 "results":results,
 "counts":{
   "candidate_count":len(results),
   "cross_survey_persistent":sum(r["verdict"]=="CROSS_SURVEY_PERSISTENT" for r in results),
   "cross_survey_fail":sum(r["verdict"]=="CROSS_SURVEY_FAIL" for r in results),
   "insufficient_support":sum(r["verdict"]=="INSUFFICIENT_COMMON_SUPPORT" for r in results),
   "nuisance_dominated":sum(r["verdict"]=="NUISANCE_DOMINATED__NO_PROMOTION" for r in results)
 },
 "source_manifest":{"FK140307":manA,"ZHNG09RR":manB},
 "claim_ceiling":"CROSS_SURVEY_PERSISTENT_SEAFLOOR_MORPHOLOGY_ONLY__NO_IDENTITY__NO_NOVELTY__NO_ARTIFICIALITY"
}
p=OUT/"JANUS-KUSTO-GLOBAL-STAGE2-HELDOUT-CROSSSURVEY-REPLAY-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({"counts":out["counts"],"results":results},indent=2))
