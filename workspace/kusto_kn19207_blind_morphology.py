#!/usr/bin/env python3
import json, math, struct
from pathlib import Path
import numpy as np
import requests

TARGET_LAT=-3.865418
TARGET_LON=-12.14244
TARGET_T=1199695295.051087
BEAM_INDEX=90
OFFSETS=[-600,-480,-360,-240,-120,120,240,360,480,600]
RADII=[250,500,1000,2000]
BASE="https://data.ngdc.noaa.gov/platforms/ocean/ships/knorr/KN192-07/multibeam/data/version1/MB/generated/"
NAME="sb20080107081434.xse.mb94.fbt"
OUT=Path("workspace/kusto_open_seafloor_out"); OUT.mkdir(parents=True,exist_ok=True)

def target_xy(lon,lat):
    return ((lon-TARGET_LON)*111320.0*math.cos(math.radians(TARGET_LAT)),
            (lat-TARGET_LAT)*111320.0)

def parse_fbt(data):
    off=0; recs=[]
    ids={b"V4":(4,90),b"V5":(5,98),b"cc":("comment",2),b"##":("oldcomment",30)}
    while off+2<=len(data):
        tag=data[off:off+2]
        if tag not in ids: raise RuntimeError(f"bad tag {tag!r} @{off}")
        ver,hs=ids[tag]
        if isinstance(ver,str):
            off+=hs+128; continue
        h=data[off:off+hs]
        t,lon,lat,sd,alt=struct.unpack_from(">5d",h,2)
        heading,speed,roll,pitch,heave,bx,bl=struct.unpack_from(">7f",h,42)
        if ver==4:
            nb,na,ns,shd=struct.unpack_from(">4h",h,70)
            dscale,xscale=struct.unpack_from(">2f",h,78)
        else:
            nb,na,ns,shd=struct.unpack_from(">4i",h,70)
            dscale,xscale=struct.unpack_from(">2f",h,86)
        p=off+hs
        flags=np.frombuffer(data,dtype=np.uint8,count=nb,offset=p).copy(); p+=nb
        def arr(n):
            nonlocal p
            a=np.frombuffer(data,dtype=">i2",count=n,offset=p).astype(float); p+=2*n; return a
        bath=arr(nb); across=arr(nb); along=arr(nb)
        p+=2*na+6*ns
        navx,navy=target_xy(lon,lat)
        hr=math.radians(heading); sh=math.sin(hr); ch=math.cos(hr)
        xtrack=xscale*across; ltrack=xscale*along
        east=navx+ltrack*sh+xtrack*ch
        north=navy+ltrack*ch-xtrack*sh
        depth=dscale*bath+sd
        recs.append({"t":t,"lon":lon,"lat":lat,"heading":heading,"flags":flags,
                     "east":east,"north":north,"depth":depth,"nb":nb})
        off=p
    return recs

def plane_residual(x,y,z):
    A=np.column_stack([np.ones(len(x)),x,y])
    coef=np.linalg.lstsq(A,z,rcond=None)[0]
    pred=A@coef
    return z-pred,float(np.sqrt(np.mean((z-pred)**2))),coef

def metrics(points,cx,cy,radius):
    x=points[:,0]-cx; y=points[:,1]-cy; z=points[:,2]
    r=np.hypot(x,y); m=r<=radius
    x=x[m]; y=y[m]; z=z[m]; r=r[m]
    out={"n_good_beams":int(len(z))}
    if len(z)<20:return out
    p05,p95=np.percentile(z,[5,95])
    out.update({
      "depth_min_m":float(z.min()),"depth_max_m":float(z.max()),
      "depth_p05_m":float(p05),"depth_p95_m":float(p95),
      "p05_p95_relief_m":float(p95-p05)
    })
    resid,prmse,_=plane_residual(x,y,z)
    out["plane_rmse_m"]=prmse
    Q=np.column_stack([np.ones(len(x)),x,y,x*x,x*y,y*y])
    qc=np.linalg.lstsq(Q,z,rcond=None)[0]
    qres=z-Q@qc
    out["quadratic_rmse_m"]=float(np.sqrt(np.mean(qres*qres)))
    H=np.array([[2*qc[3],qc[4]],[qc[4],2*qc[5]]])
    eig=np.linalg.eigvalsh(H)
    out["quadratic_hessian_eigenvalues_per_m"]=[float(v) for v in eig]
    ab=np.abs(eig); mn=float(ab.min()); mx=float(ab.max())
    out["curvature_anisotropy_abs_eigenvalue_ratio"]=None if mn<1e-12 else mx/mn
    # radial quadratic fit on plane residual
    R=np.column_stack([np.ones(len(r)),r,r*r])
    rc=np.linalg.lstsq(R,resid,rcond=None)[0]
    pred=R@rc
    ssr=float(np.sum((resid-pred)**2)); sst=float(np.sum((resid-resid.mean())**2))
    out["radial_quadratic_r2_after_plane"]=None if sst<=0 else 1-ssr/sst
    if radius>=1000:
        ann=(r>=250)&(r<=1000)
        if ann.sum()>=20:
            theta=np.arctan2(y[ann],x[ann]); rr=resid[ann]
            amp=2*abs(np.mean(rr*np.exp(-4j*theta)))
            rms=float(np.sqrt(np.mean(rr*rr)))
            out["angular_m4_normalized_amplitude_250_1000m"]=None if rms==0 else float(amp/rms)
            out["angular_m4_n"]=int(ann.sum())
    return out

data=requests.get(BASE+NAME,timeout=120,headers={"User-Agent":"JANUS-KUSTO-blind-morphology/1.0"}).content
recs=parse_fbt(data)

pts=[]
for ri,r in enumerate(recs):
    good=np.where(r["flags"]==0)[0]
    for bi in good:
        pts.append((r["east"][bi],r["north"][bi],r["depth"][bi],ri,int(bi),r["t"]))
pts=np.asarray(pts,float)

centers=[{"id":"TARGET","cx":0.0,"cy":0.0,"requested_offset_s":0,"record_time":TARGET_T,"record_index":None}]
for off in OFFSETS:
    idx=int(np.argmin([abs(r["t"]-(TARGET_T+off)) for r in recs]))
    r=recs[idx]
    if BEAM_INDEX>=r["nb"] or r["flags"][BEAM_INDEX]!=0:
        centers.append({"id":f"CONTROL_{off:+d}s","available":False,"requested_offset_s":off,
                        "reason":"beam90_not_good_or_absent","record_index":idx,"record_time":r["t"]})
    else:
        centers.append({"id":f"CONTROL_{off:+d}s","available":True,"requested_offset_s":off,
                        "actual_offset_s":r["t"]-TARGET_T,"record_index":idx,"record_time":r["t"],
                        "cx":float(r["east"][BEAM_INDEX]),"cy":float(r["north"][BEAM_INDEX])})

for c in centers:
    if c.get("available",True):
        c["metrics"]={str(rad):metrics(pts,c["cx"],c["cy"],rad) for rad in RADII}

# Rank target among available target+controls for scalar metrics at each radius.
scalar=["p05_p95_relief_m","plane_rmse_m","quadratic_rmse_m",
        "curvature_anisotropy_abs_eigenvalue_ratio","radial_quadratic_r2_after_plane",
        "angular_m4_normalized_amplitude_250_1000m"]
ranks={}
for rad in RADII:
    rr={}
    for key in scalar:
        vals=[]
        for c in centers:
            if not c.get("available",True):continue
            v=c.get("metrics",{}).get(str(rad),{}).get(key)
            if v is not None and np.isfinite(v): vals.append((c["id"],float(v)))
        vals.sort(key=lambda z:z[1])
        tv=next((v for i,v in vals if i=="TARGET"),None)
        tr=next((j+1 for j,(i,v) in enumerate(vals) if i=="TARGET"),None)
        rr[key]={"target_value":tv,"ascending_rank":tr,"n":len(vals),"ordered":vals}
    ranks[str(rad)]=rr

out={
 "artifact_id":"JANUS-KUSTO-KN19207-BLIND-LOCAL-MORPHOLOGY-RUN-2026-09-21-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-KN19207-BLIND-LOCAL-MORPHOLOGY-PREREG-2026-09-21-v1.0.json",
 "source_fbt":NAME,
 "point_count_good":int(len(pts)),
 "centers":centers,
 "target_ranks":ranks,
 "interpretation_ceiling":"BLIND_MORPHOLOGY_SCREEN_ONLY__NO_ARTIFICIALITY_OR_IDENTITY_CLAIM",
 "automatic_promotion":"FORBIDDEN__REQUIRES_POSTRUN_REVIEW_AGAINST_PREREG_WITHOUT_NEW_METRICS"
}
p=OUT/"JANUS-KUSTO-KN19207-BLIND-LOCAL-MORPHOLOGY-RUN-2026-09-21-v1.0.json"
p.write_text(json.dumps(out,indent=2),encoding="utf-8")
print(json.dumps({
 "point_count_good":out["point_count_good"],
 "controls":[{k:v for k,v in c.items() if k not in ("metrics","cx","cy")} for c in centers[1:]],
 "target_metrics":centers[0]["metrics"],
 "target_ranks":ranks
},indent=2))
