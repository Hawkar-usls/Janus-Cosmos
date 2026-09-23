#!/usr/bin/env python3
from __future__ import annotations
import gc, hashlib, io, json, math, zipfile
from pathlib import Path

import numpy as np
import requests
import rasterio
from rasterio.io import MemoryFile
from scipy import ndimage

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-HMAS-CANBERRA-HC0-KNOWN-WRECK-RECOVERY-PREREG-2026-09-23-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-HMAS-Canberra-HC1-blind/1.0"}
URL=PRE["dataset"]["source_url"]
CELL=float(PRE["discovery"]["native_cell_m"])
RADII=list(PRE["discovery"]["radii_m"])
STRIDE=float(PRE["discovery"]["sampling_stride_m"])
TAIL=float(PRE["discovery"]["tail_quantile"])
DEDUP=float(PRE["discovery"]["spatial_dedup_m"])
MAX_CAND=int(PRE["discovery"]["maximum_candidates"])
FAMILIES={
 "RELIEF_ROUGHNESS":["local_relief_m","local_std_m"],
 "CURVATURE":["abs_laplacian"],
 "RADIAL":["abs_center_surround_dog"],
 "LINEARITY":["structure_tensor_anisotropy"],
 "TERRACE_STEP":["gradient_magnitude"],
 "TOPOLOGIC":["isolated_extremum_persistence"]
}

def get(url):
    r=requests.get(url,headers=UA,timeout=180,allow_redirects=True);r.raise_for_status();return r.content

def robust(v):
    med=float(np.median(v));mad=float(np.median(np.abs(v-med)))
    scale=max(1e-12,1.4826*mad)
    return np.abs((v-med)/scale),med,mad

blob=get(URL)
zip_sha=hashlib.sha256(blob).hexdigest()
with zipfile.ZipFile(io.BytesIO(blob)) as zf:
    names=zf.namelist()
    tifs=[n for n in names if n.lower().endswith((".tif",".tiff")) and not n.startswith("__MACOSX/") and "/._" not in n]
    if len(tifs)!=1:raise RuntimeError(f"expected one GeoTIFF, found {tifs}")
    tif=zf.read(tifs[0])

with MemoryFile(tif) as mf:
    with mf.open() as src:
        ma=src.read(1,masked=True)
        z=np.asarray(ma.filled(np.nan),dtype=np.float32)
        valid=(~np.ma.getmaskarray(ma)) & np.isfinite(z)
        transform=src.transform
        crs=str(src.crs) if src.crs is not None else None
        nodata=src.nodata
        dtype=str(src.dtypes[0])

res=(abs(float(transform.a)),abs(float(transform.e)))
if max(abs(res[0]-CELL),abs(res[1]-CELL))/CELL>0.01:
    raise RuntimeError(f"native resolution outside 1% prereg tolerance: {res}")
CELL_EFF=float(sum(res)/2.0)

nrows,ncols=z.shape
if not np.any(valid):raise RuntimeError("no valid bathymetry")
fill=float(np.median(z[valid]));zf=z.copy();zf[~valid]=fill
maxpx=int(round(max(RADII)/CELL_EFF))
rr,cc=np.indices(z.shape,dtype=np.int32)
edge=np.minimum.reduce([rr,cc,nrows-1-rr,ncols-1-cc]);del rr,cc;gc.collect()
if np.all(valid):dist=np.full(z.shape,maxpx+1,dtype=np.float32)
else:dist=ndimage.distance_transform_edt(valid).astype(np.float32)
inner=valid & (edge>=maxpx) & (dist>=maxpx);del edge,dist;gc.collect()

step=max(1,int(round(STRIDE/CELL_EFF)))
rows=np.arange(maxpx,nrows-maxpx,step,dtype=np.int32)
cols=np.arange(maxpx,ncols-maxpx,step,dtype=np.int32)
R,C=np.meshgrid(rows,cols,indexing="ij")
sr=R.ravel();sc=C.ravel();ok=inner[sr,sc]
sr=sr[ok];sc=sc[ok];del R,C,inner,ok;gc.collect()
if len(sr)<100:raise RuntimeError(f"too few sample centers {len(sr)}")

def xy(row,col):
    x,y=rasterio.transform.xy(transform,row,col,offset="center");return float(x),float(y)

all_flags=[];scale_summaries=[]
for radius in RADII:
    print("HC1 radius",radius,flush=True)
    rpx=max(1,int(round(radius/CELL_EFF)));size=2*rpx+1;sigma=max(1.0,rpx/3.0)
    vals={}
    mean=ndimage.uniform_filter(zf,size=size,mode="nearest")
    sq=zf*zf;mean2=ndimage.uniform_filter(sq,size=size,mode="nearest")
    std=np.sqrt(np.maximum(0.0,mean2-mean*mean))
    vals["local_std_m"]=std[sr,sc].astype(np.float64)
    persist=np.abs(zf-mean)/(std+1e-6)
    vals["isolated_extremum_persistence"]=persist[sr,sc].astype(np.float64)
    del sq,mean2,std,persist;gc.collect()
    zmax=ndimage.maximum_filter(zf,size=size,mode="nearest");zmin=ndimage.minimum_filter(zf,size=size,mode="nearest")
    vals["local_relief_m"]=(zmax[sr,sc]-zmin[sr,sc]).astype(np.float64);del zmax,zmin;gc.collect()
    smooth=ndimage.gaussian_filter(zf,sigma=sigma,mode="nearest",truncate=3.0)
    lap=np.abs(ndimage.laplace(smooth,mode="nearest"))/(CELL_EFF*CELL_EFF)
    vals["abs_laplacian"]=lap[sr,sc].astype(np.float64);del lap;gc.collect()
    gy,gx=np.gradient(smooth,CELL_EFF,CELL_EFF);grad=np.hypot(gx,gy)
    vals["gradient_magnitude"]=grad[sr,sc].astype(np.float64);del grad,smooth;gc.collect()
    ts=max(1.0,rpx/2.0)
    jxx=ndimage.gaussian_filter(gx*gx,sigma=ts,mode="nearest",truncate=3.0)
    jyy=ndimage.gaussian_filter(gy*gy,sigma=ts,mode="nearest",truncate=3.0)
    jxy=ndimage.gaussian_filter(gx*gy,sigma=ts,mode="nearest",truncate=3.0)
    anis=np.sqrt((jxx-jyy)**2+4*jxy*jxy)/(jxx+jyy+1e-12)
    vals["structure_tensor_anisotropy"]=anis[sr,sc].astype(np.float64);del gx,gy,jxx,jyy,jxy,anis;gc.collect()
    g1=ndimage.gaussian_filter(zf,sigma=max(.8,rpx/4.0),mode="nearest",truncate=3.0)
    g2=ndimage.gaussian_filter(zf,sigma=max(1.2,rpx/1.5),mode="nearest",truncate=3.0)
    vals["abs_center_surround_dog"]=np.abs(g1[sr,sc]-g2[sr,sc]).astype(np.float64);del g1,g2,mean;gc.collect()

    thresholds={};rz={}
    for k,v in vals.items():
        good=np.isfinite(v)
        if not np.all(good):
            if not np.any(good):raise RuntimeError(f"{radius} {k}: all invalid")
            v=np.where(good,v,np.nanmedian(v[good]));vals[k]=v
        q=float(np.quantile(v,TAIL));zv,med,mad=robust(v)
        thresholds[k]={"q995":q,"median":med,"mad":mad};rz[k]=zv
    stack=np.column_stack([rz[k] for k in sorted(rz)])
    generic=.5*np.max(stack,axis=1)+.5*np.mean(stack,axis=1);gq=float(np.quantile(generic,TAIL))
    npre=0
    for i in range(len(sr)):
        votes=[];detail={}
        for fam,keys in FAMILIES.items():
            ex=[k for k in keys if vals[k][i]>=thresholds[k]["q995"]]
            if ex:votes.append(fam);detail[fam]=ex
        if len(votes)>=2 or generic[i]>=gq:
            row=int(sr[i]);col=int(sc[i]);x,y=xy(row,col)
            all_flags.append({"row":row,"col":col,"x_m":x,"y_m":y,"radius_m":radius,
              "family_votes":votes,"extreme_metrics_by_family":detail,"generic_robust_score":float(generic[i]),
              "metrics":{k:float(vals[k][i]) for k in vals}})
            npre+=1
    scale_summaries.append({"radius_m":radius,"sample_centers":int(len(sr)),"generic_q995":gq,
      "metric_thresholds":thresholds,"followup_count_pre_dedup":npre})
    del vals,rz,stack,generic;gc.collect()

all_flags.sort(key=lambda q:(-len(q["family_votes"]),-q["generic_robust_score"],q["radius_m"],q["row"],q["col"]))
kept=[]
for q in all_flags:
    if any(math.hypot(q["x_m"]-p["x_m"],q["y_m"]-p["y_m"])<DEDUP for p in kept):continue
    q=dict(q);q["candidate_id"]=f"HC_BATHY_C{len(kept)+1:03d}";q["blind_rank"]=len(kept)+1
    kept.append(q)
    if len(kept)>=MAX_CAND:break

out={
 "artifact_id":"JANUS-KUSTO-HMAS-CANBERRA-HC1-BLIND-BATHYMETRY-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "input":{"source_url":URL,"zip_sha256":zip_sha,"zip_bytes":len(blob),"zip_members":names,"tif_name":tifs[0],
          "shape":list(z.shape),"valid_cells":int(np.sum(valid)),"dtype":dtype,"nodata":nodata,"crs":crs,
          "transform":list(transform)[:6],"resolution_m":list(res),"effective_cell_m":CELL_EFF},
 "truth_firewall":{"official_dive_site_anchor_read":False,"wreck_identity_used_for_candidate_selection":False},
 "blind_sampling":{"stride_m":STRIDE,"sample_centers":int(len(sr)),"radii_m":RADII,"tail_quantile":TAIL},
 "scale_summaries":scale_summaries,
 "pre_dedup_followup_total":len(all_flags),
 "frozen_candidate_count":len(kept),
 "frozen_candidates":kept,
 "next_gate":"HC2_TRUTH_UNLOCK_CASE_LEVEL_RECOVERY_SCORE",
 "claim_ceiling":"BLIND_KNOWN_OBJECT_POSITIVE_CONTROL_CANDIDATE_GENERATION_ONLY"
}
p=OUT/"JANUS-KUSTO-HMAS-CANBERRA-HC1-BLIND-BATHYMETRY-RUN-2026-09-23-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(json.dumps({"artifact_id":out["artifact_id"],"shape":out["input"]["shape"],"valid_cells":out["input"]["valid_cells"],
 "sample_centers":out["blind_sampling"]["sample_centers"],"pre_dedup":len(all_flags),"candidates":len(kept),
 "candidate_json_sha256":hashlib.sha256(raw.encode()).hexdigest(),"truth_anchor_read":False},indent=2))
