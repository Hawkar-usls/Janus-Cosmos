#!/usr/bin/env python3
import concurrent.futures, hashlib, json, math, re, struct
from pathlib import Path
from urllib.parse import urljoin
import numpy as np
import requests
from scipy.spatial import cKDTree

BASE="https://data.ngdc.noaa.gov/platforms/ocean/ships/knorr/KN192-07/multibeam/data/version1/MB/generated/"
OUT=Path("workspace/kusto_open_seafloor_out")
OUT.mkdir(parents=True,exist_ok=True)
LAT0=-4.5
COS0=math.cos(math.radians(LAT0))
M_PER_DEG=111320.0
TARGET=(-3.865418,-12.14244)
EXCLUDE_M=5000.0
RADII=[500,1000]
MIN_N={500:100,1000:300}
FAMILIES={
 "ROUGHNESS_FAMILY":["p05_p95_relief_m","plane_rmse_m","quadratic_rmse_m"],
 "CURVATURE_FAMILY":["curvature_anisotropy_abs_eigenvalue_ratio"],
 "RADIAL_FAMILY":["radial_quadratic_r2_after_plane"],
 "ANGULAR_FAMILY":["angular_m4_normalized_amplitude_250_1000m"]
}
TAIL=0.005
UA={"User-Agent":"JANUS-KUSTO-survey-wide-discovery/1.0"}

def xy_from_lonlat(lon,lat):
    return lon*M_PER_DEG*COS0, lat*M_PER_DEG

def lonlat_from_xy(x,y):
    return x/(M_PER_DEG*COS0), y/M_PER_DEG

def parse_fbt(data,name):
    off=0; points=[]; candidates=[]; rec_idx=0
    ids={b"V4":(4,90),b"V5":(5,98),b"cc":("comment",2),b"##":("oldcomment",30)}
    while off+2<=len(data):
        tag=data[off:off+2]
        if tag not in ids: raise RuntimeError(f"{name}: bad tag {tag!r} @{off}")
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
            a=np.frombuffer(data,dtype=">i2",count=n,offset=p).astype(np.float64); p+=2*n; return a
        bath=arr(nb); across=arr(nb); along=arr(nb)
        p+=2*na+6*ns
        navx,navy=xy_from_lonlat(lon,lat)
        hr=math.radians(heading); sh=math.sin(hr); ch=math.cos(hr)
        xtrack=xscale*across; ltrack=xscale*along
        east=navx+ltrack*sh+xtrack*ch
        north=navy+ltrack*ch-xtrack*sh
        depth=dscale*bath+sd
        good=np.where(flags==0)[0]
        if len(good):
            points.append(np.column_stack([east[good],north[good],depth[good]]))
        if rec_idx>=10 and (rec_idx-10)%20==0 and nb>90 and flags[90]==0:
            cx=float(east[90]); cy=float(north[90])
            lonc,latc=lonlat_from_xy(cx,cy)
            candidates.append({
              "source_file":name,"record_index":rec_idx,"epoch":float(t),
              "beam_index":90,"x":cx,"y":cy,"lon":lonc,"lat":latc
            })
        off=p; rec_idx+=1
    pts=np.vstack(points) if points else np.empty((0,3),dtype=float)
    return pts,candidates,rec_idx

def fetch_parse(url):
    r=requests.get(url,headers=UA,timeout=120); r.raise_for_status()
    b=r.content; name=url.rsplit("/",1)[-1]
    pts,cands,nrec=parse_fbt(b,name)
    return name,hashlib.sha256(b).hexdigest(),len(b),pts,cands,nrec

def plane_residual(x,y,z):
    A=np.column_stack([np.ones(len(x)),x,y])
    coef=np.linalg.lstsq(A,z,rcond=None)[0]
    res=z-A@coef
    return res,float(np.sqrt(np.mean(res*res)))

def metric_patch(pts,cx,cy,rad):
    x=pts[:,0]-cx; y=pts[:,1]-cy; z=pts[:,2]
    r=np.hypot(x,y)
    if len(z)<MIN_N[rad]: return None
    p05,p95=np.percentile(z,[5,95])
    resid,prmse=plane_residual(x,y,z)
    Q=np.column_stack([np.ones(len(x)),x,y,x*x,x*y,y*y])
    qc=np.linalg.lstsq(Q,z,rcond=None)[0]
    qres=z-Q@qc
    H=np.array([[2*qc[3],qc[4]],[qc[4],2*qc[5]]])
    eig=np.linalg.eigvalsh(H); ab=np.abs(eig)
    anis=None if float(ab.min())<1e-12 else float(ab.max()/ab.min())
    R=np.column_stack([np.ones(len(r)),r,r*r])
    rc=np.linalg.lstsq(R,resid,rcond=None)[0]
    pred=R@rc
    ssr=float(np.sum((resid-pred)**2)); sst=float(np.sum((resid-resid.mean())**2))
    radial=None if sst<=0 else float(1-ssr/sst)
    out={
      "n_good_beams":int(len(z)),
      "p05_p95_relief_m":float(p95-p05),
      "plane_rmse_m":prmse,
      "quadratic_rmse_m":float(np.sqrt(np.mean(qres*qres))),
      "curvature_anisotropy_abs_eigenvalue_ratio":anis,
      "radial_quadratic_r2_after_plane":radial
    }
    if rad==1000:
        ann=(r>=250)&(r<=1000)
        if ann.sum()>=20:
            theta=np.arctan2(y[ann],x[ann]); rr=resid[ann]
            amp=2*abs(np.mean(rr*np.exp(-4j*theta)))
            rms=float(np.sqrt(np.mean(rr*rr)))
            out["angular_m4_normalized_amplitude_250_1000m"]=None if rms==0 else float(amp/rms)
        else:
            out["angular_m4_normalized_amplitude_250_1000m"]=None
    else:
        out["angular_m4_normalized_amplitude_250_1000m"]=None
    return out

def empirical_percentiles(vals):
    order=np.argsort(vals,kind="mergesort")
    ranks=np.empty(len(vals),float)
    ranks[order]=(np.arange(len(vals))+0.5)/len(vals)
    return ranks

html=requests.get(BASE,headers=UA,timeout=60).text
names=sorted(set(re.findall(r'href="([^"]+\.fbt)"',html,re.I)))
urls=[urljoin(BASE,n) for n in names]
print("FBT discovered",len(urls))

parts=[]; candidates=[]; manifest=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=20) as ex:
    futs=[ex.submit(fetch_parse,u) for u in urls]
    for k,fut in enumerate(concurrent.futures.as_completed(futs),1):
        name,h,b,pts,cands,nrec=fut.result()
        parts.append(pts); candidates.extend(cands)
        manifest.append({"file":name,"sha256":h,"bytes":b,"records":nrec,"good_points":int(len(pts)),"candidate_centers":len(cands)})
        if k%50==0: print("downloaded",k,"/",len(urls))
points=np.vstack(parts)
manifest.sort(key=lambda x:x["file"])
print("good beam points",len(points),"candidate centers raw",len(candidates))

tx,ty=xy_from_lonlat(TARGET[1],TARGET[0])
eligible=[]
for c in candidates:
    if math.hypot(c["x"]-tx,c["y"]-ty)<EXCLUDE_M: continue
    eligible.append(c)
print("candidate centers after target exclusion",len(eligible))

tree=cKDTree(points[:,:2])
for ci,c in enumerate(eligible):
    c["metrics"]={}
    for rad in RADII:
        idx=tree.query_ball_point([c["x"],c["y"]],rad)
        patch=points[np.asarray(idx,dtype=int)] if idx else np.empty((0,3))
        c["metrics"][str(rad)]=metric_patch(patch,c["x"],c["y"],rad)
    if (ci+1)%500==0: print("metrics",ci+1,"/",len(eligible))

thresholds={}
for rad in RADII:
    rr=str(rad); thresholds[rr]={}
    keys=sorted({m for fam in FAMILIES.values() for m in fam})
    for key in keys:
        vals=[]; refs=[]
        for i,c in enumerate(eligible):
            m=c["metrics"].get(rr)
            if m is not None and m.get(key) is not None and np.isfinite(m[key]):
                vals.append(float(m[key])); refs.append(i)
        if not vals:
            thresholds[rr][key]={"n":0}; continue
        arr=np.asarray(vals,float)
        lo=float(np.quantile(arr,TAIL)); hi=float(np.quantile(arr,1-TAIL))
        pct=empirical_percentiles(arr)
        thresholds[rr][key]={"n":len(arr),"q005":lo,"q995":hi}
        for j,i in enumerate(refs):
            m=eligible[i]["metrics"][rr]
            m.setdefault("_percentiles",{})[key]=float(pct[j])
            m.setdefault("_extremes",{})[key]=bool(arr[j]<=lo or arr[j]>=hi)

flags=[]
for c in eligible:
    c["radius_screen"]={}
    for rad in RADII:
        rr=str(rad); m=c["metrics"].get(rr)
        if m is None:
            c["radius_screen"][rr]={"eligible":False}; continue
        votes=[]; details={}
        for fam,keys in FAMILIES.items():
            ex=[k for k in keys if m.get("_extremes",{}).get(k,False)]
            if ex:
                votes.append(fam); details[fam]=ex
        structural=any(f in votes for f in ["CURVATURE_FAMILY","RADIAL_FAMILY","ANGULAR_FAMILY"])
        follow=len(votes)>=2 and structural
        pcts=list(m.get("_percentiles",{}).values())
        max_ext=max([abs(p-0.5)*2 for p in pcts],default=0.0)
        c["radius_screen"][rr]={
          "eligible":True,"family_votes":votes,"family_vote_count":len(votes),
          "extreme_metrics_by_family":details,"has_structural_vote":structural,
          "followup_flag":follow,"max_empirical_tail_extremeness":max_ext
        }
        if follow:
            flags.append((c,rad))

# Flatten per-radius flagged candidates, sort exactly by prereg order, then spatially suppress.
flat=[]
for c,rad in flags:
    s=c["radius_screen"][str(rad)]
    flat.append({
      "source_file":c["source_file"],"record_index":c["record_index"],"epoch":c["epoch"],
      "lon":c["lon"],"lat":c["lat"],"x":c["x"],"y":c["y"],"radius_m":rad,
      "family_votes":s["family_votes"],"family_vote_count":s["family_vote_count"],
      "extreme_metrics_by_family":s["extreme_metrics_by_family"],
      "max_empirical_tail_extremeness":s["max_empirical_tail_extremeness"],
      "metrics":c["metrics"][str(rad)]
    })
flat.sort(key=lambda z:(-z["family_vote_count"],-z["max_empirical_tail_extremeness"],z["source_file"],z["record_index"],z["radius_m"]))
kept=[]
for z in flat:
    if any(math.hypot(z["x"]-q["x"],z["y"]-q["y"])<2000 for q in kept): continue
    kept.append(z)
    if len(kept)>=50: break

# Remove internal XY from user-facing candidates; retain lon/lat.
for z in kept:
    z.pop("x",None); z.pop("y",None)

out={
 "artifact_id":"JANUS-KUSTO-KN19207-SURVEY-WIDE-BLIND-DISCOVERY-RUN-2026-09-21-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-KN19207-SURVEY-WIDE-BLIND-DISCOVERY-PREREG-2026-09-21-v1.0.json",
 "survey":"KN192-07","fbt_files_discovered":len(names),"fbt_files_processed":len(manifest),
 "good_beam_point_count":int(len(points)),
 "raw_candidate_center_count":len(candidates),
 "candidate_center_count_after_5km_target_exclusion":len(eligible),
 "thresholds":thresholds,
 "pre_dedup_followup_flag_count":len(flat),
 "spatially_deduplicated_followup_count":len(kept),
 "followup_candidates":kept,
 "file_manifest":manifest,
 "claim_ceiling":"BLIND_SURVEY_WIDE_EXPLORATORY_DISCOVERY_ONLY__NO_IDENTITY__NO_ARTIFICIALITY",
 "independent_confirmation_required":True
}
p=OUT/"JANUS-KUSTO-KN19207-SURVEY-WIDE-BLIND-DISCOVERY-RUN-2026-09-21-v1.0.json"
p.write_text(json.dumps(out,indent=2),encoding="utf-8")
print(json.dumps({
 "fbt_files_processed":out["fbt_files_processed"],
 "good_beam_point_count":out["good_beam_point_count"],
 "candidate_center_count":out["candidate_center_count_after_5km_target_exclusion"],
 "pre_dedup_followup_flag_count":out["pre_dedup_followup_flag_count"],
 "spatially_deduplicated_followup_count":out["spatially_deduplicated_followup_count"],
 "top_candidates":kept[:20]
},indent=2))
