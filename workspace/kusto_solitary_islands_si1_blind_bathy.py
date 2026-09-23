#!/usr/bin/env python3
from __future__ import annotations

import gc, hashlib, io, json, math, zipfile
from pathlib import Path

import numpy as np
import requests
import rasterio
from affine import Affine
from rasterio.io import MemoryFile
from scipy import ndimage

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOLITARY-ISLANDS-SI0-UNCHANGED-V2-PREREG-2026-09-24-v1.0.json").read_text())
SRC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOLITARY-ISLANDS-SI0E-EXACT-MEMBER-BINDING-RECEIPT-2026-09-24-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-SolitaryIslands-SI1-blind-bathy/1.0"}
CELL=5.0
RADII=[20,40,80,160]
STRIDE=40.0
TAIL=.995
DEDUP=250.0
MAX_CAND=100
FAMILIES={
 "RELIEF_ROUGHNESS":["local_relief_m","local_std_m"],
 "CURVATURE":["abs_laplacian"],
 "RADIAL":["abs_center_surround_dog"],
 "LINEARITY":["structure_tensor_anisotropy"],
 "TERRACE_STEP":["gradient_magnitude"],
 "TOPOLOGIC":["isolated_extremum_persistence"]
}

def get(url):
    r=requests.get(url,headers=UA,timeout=180);r.raise_for_status();return r.content

def world_affine(text):
    vals=[float(x.strip()) for x in text.splitlines() if x.strip()]
    if len(vals)!=6:raise RuntimeError(f"world file requires 6 values, got {vals}")
    A,D,B,E,C,F=vals
    return vals,Affine(A,B,C-.5*A-.5*B,D,E,F-.5*D-.5*E)

def load_bathy():
    blob=get(SRC["source_archive"]["url"])
    sha=hashlib.sha256(blob).hexdigest()
    if sha.lower()!=SRC["source_archive"]["sha256"].lower():
        raise RuntimeError("source archive SHA256 mismatch")
    if len(blob)!=int(SRC["source_archive"]["bytes"]):
        raise RuntimeError("source archive byte length mismatch")
    member=SRC["discovery_member"]["name"]
    heldout=SRC["heldout_member"]["name"]
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        names=zf.namelist()
        if member not in names:
            raise RuntimeError(f"frozen discovery member missing: {member}")
        if heldout not in names:
            raise RuntimeError(f"frozen heldout member unexpectedly missing: {heldout}")
        zi=zf.getinfo(member)
        if int(zi.file_size)!=int(SRC["discovery_member"]["uncompressed_bytes"]):
            raise RuntimeError("discovery member size mismatch")
        if int(zi.CRC)!=int(SRC["discovery_member"]["crc32"]):
            raise RuntimeError("discovery member CRC mismatch")
        # Scientific firewall: only the frozen 5 m bathymetry member is decompressed.\n        # The held-out 5 m backscatter member is never opened/read/hashed in SI1.
        tif=zf.read(member)
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
    if transform.is_identity:
        raise RuntimeError("embedded transform is identity")
    if max(abs(res[0]-CELL),abs(res[1]-CELL))/CELL>0.001:
        raise RuntimeError(f"effective resolution outside frozen 0.1% nominal tolerance: {res}")
    return blob,z,valid,transform,{
      "zip_sha256":sha,"zip_members":names,
      "tif_name":member,"tfw_name":None,
      "shape":list(z.shape),"dtype":dtype,"nodata":nodata,"crs":crs,
      "resolution_m":list(res),"georef_source":"EMBEDDED_TRANSFORM",
      "embedded_transform":list(transform)[:6],"effective_transform":list(transform)[:6],
      "world_file_values":None,
      "archive_contains_heldout_backscatter_member":True,
      "heldout_backscatter_member_name":heldout,
      "heldout_backscatter_member_decompressed":False,
      "heldout_backscatter_values_decoded":False
    }

def robust(v):
    med=float(np.median(v));mad=float(np.median(np.abs(v-med)))
    scale=max(1e-12,1.4826*mad)
    return np.abs((v-med)/scale),med,mad

def xy(transform,row,col):
    x,y=rasterio.transform.xy(transform,row,col,offset="center")
    return float(x),float(y)

blob,z,valid,transform,input_meta=load_bathy()
CELL_EFF=float(sum(input_meta["resolution_m"])/2.0)
input_meta["effective_cell_for_geometry_m"]=CELL_EFF
nrows,ncols=z.shape
if not np.any(valid):raise RuntimeError("no valid bathymetry cells")
input_meta["zip_bytes"]=len(blob);input_meta["valid_cells"]=int(np.sum(valid))

fill=float(np.median(z[valid]));zf=z.copy();zf[~valid]=fill
maxpx=int(round(max(RADII)/CELL_EFF))
rr,cc=np.indices(z.shape,dtype=np.int32)
edge=np.minimum.reduce([rr,cc,nrows-1-rr,ncols-1-cc])
del rr,cc;gc.collect()
if np.all(valid):dist=np.full(z.shape,maxpx+1,dtype=np.float32)
else:dist=ndimage.distance_transform_edt(valid).astype(np.float32)
inner=valid & (edge>=maxpx) & (dist>=maxpx)
del edge,dist;gc.collect()

step=max(1,int(round(STRIDE/CELL_EFF)))
rows=np.arange(maxpx,nrows-maxpx,step,dtype=np.int32)
cols=np.arange(maxpx,ncols-maxpx,step,dtype=np.int32)
R,C=np.meshgrid(rows,cols,indexing="ij")
sr=R.ravel();sc=C.ravel();ok=inner[sr,sc]
sr=sr[ok];sc=sc[ok]
del R,C,inner,ok;gc.collect()
if len(sr)<1000:raise RuntimeError(f"too few sample centers {len(sr)}")

all_flags=[];scale_summaries=[]

for radius in RADII:
    print("S1 radius",radius,flush=True)
    rpx=max(1,int(round(radius/CELL_EFF)));size=2*rpx+1;sigma=max(1.0,rpx/3.0)
    vals={}

    mean=ndimage.uniform_filter(zf,size=size,mode="nearest")
    sq=zf*zf
    mean2=ndimage.uniform_filter(sq,size=size,mode="nearest")
    std=np.sqrt(np.maximum(0.0,mean2-mean*mean))
    vals["local_std_m"]=std[sr,sc].astype(np.float64)
    persist=np.abs(zf-mean)/(std+1e-6)
    vals["isolated_extremum_persistence"]=persist[sr,sc].astype(np.float64)
    del sq,mean2,std,persist;gc.collect()

    zmax=ndimage.maximum_filter(zf,size=size,mode="nearest")
    zmin=ndimage.minimum_filter(zf,size=size,mode="nearest")
    vals["local_relief_m"]=(zmax[sr,sc]-zmin[sr,sc]).astype(np.float64)
    del zmax,zmin;gc.collect()

    smooth=ndimage.gaussian_filter(zf,sigma=sigma,mode="nearest",truncate=3.0)
    lap=np.abs(ndimage.laplace(smooth,mode="nearest"))/(CELL_EFF*CELL_EFF)
    vals["abs_laplacian"]=lap[sr,sc].astype(np.float64)
    del lap;gc.collect()

    gy,gx=np.gradient(smooth,CELL_EFF,CELL_EFF)
    grad=np.hypot(gx,gy)
    vals["gradient_magnitude"]=grad[sr,sc].astype(np.float64)
    del grad,smooth;gc.collect()

    ts=max(1.0,rpx/2.0)
    jxx=ndimage.gaussian_filter(gx*gx,sigma=ts,mode="nearest",truncate=3.0)
    jyy=ndimage.gaussian_filter(gy*gy,sigma=ts,mode="nearest",truncate=3.0)
    jxy=ndimage.gaussian_filter(gx*gy,sigma=ts,mode="nearest",truncate=3.0)
    anis=np.sqrt((jxx-jyy)**2+4*jxy*jxy)/(jxx+jyy+1e-12)
    vals["structure_tensor_anisotropy"]=anis[sr,sc].astype(np.float64)
    del gx,gy,jxx,jyy,jxy,anis;gc.collect()

    g1=ndimage.gaussian_filter(zf,sigma=max(.8,rpx/4.0),mode="nearest",truncate=3.0)
    g2=ndimage.gaussian_filter(zf,sigma=max(1.2,rpx/1.5),mode="nearest",truncate=3.0)
    vals["abs_center_surround_dog"]=np.abs(g1[sr,sc]-g2[sr,sc]).astype(np.float64)
    del g1,g2,mean;gc.collect()

    thresholds={};rz={}
    for k,v in vals.items():
        good=np.isfinite(v)
        if not np.all(good):
            if not np.any(good):raise RuntimeError(f"{radius} {k}: all invalid")
            v=np.where(good,v,np.nanmedian(v[good]));vals[k]=v
        q=float(np.quantile(v,TAIL));zv,med,mad=robust(v)
        thresholds[k]={"q995":q,"median":med,"mad":mad};rz[k]=zv

    stack=np.column_stack([rz[k] for k in sorted(rz)])
    generic=.5*np.max(stack,axis=1)+.5*np.mean(stack,axis=1)
    gq=float(np.quantile(generic,TAIL))
    npre=0
    for i in range(len(sr)):
        votes=[];detail={}
        for fam,keys in FAMILIES.items():
            ex=[k for k in keys if vals[k][i]>=thresholds[k]["q995"]]
            if ex:votes.append(fam);detail[fam]=ex
        if len(votes)>=2 or generic[i]>=gq:
            row=int(sr[i]);col=int(sc[i]);x,y=xy(transform,row,col)
            all_flags.append({
              "row":row,"col":col,"x_m":x,"y_m":y,"radius_m":radius,
              "family_votes":votes,"extreme_metrics_by_family":detail,
              "generic_robust_score":float(generic[i]),
              "metrics":{k:float(vals[k][i]) for k in vals}
            });npre+=1
    scale_summaries.append({"radius_m":radius,"sample_centers":int(len(sr)),"generic_q995":gq,"metric_thresholds":thresholds,"followup_count_pre_dedup":npre})
    del vals,rz,stack,generic;gc.collect()

all_flags.sort(key=lambda q:(-len(q["family_votes"]),-q["generic_robust_score"],q["radius_m"],q["row"],q["col"]))
kept=[]
for q in all_flags:
    if any(math.hypot(q["x_m"]-p["x_m"],q["y_m"]-p["y_m"])<DEDUP for p in kept):continue
    q=dict(q);q["candidate_id"]=f"SI_BATHY_C{len(kept)+1:03d}";q["blind_rank"]=len(kept)+1
    kept.append(q)
    if len(kept)>=MAX_CAND:break

out={
 "artifact_id":"JANUS-KUSTO-SOLITARY-ISLANDS-SI1-BLIND-BATHYMETRY-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],"source_receipt":SRC["artifact_id"],
 "input":input_meta,
 "truth_firewall":{"transport_archive_contains_backscatter_member":True,"backscatter_member_decompressed":False,"backscatter_values_read":False,"landform_or_sediment_labels_read":False,"towed_video_or_stills_read":False},
 "blind_sampling":{"stride_m":STRIDE,"sample_centers":int(len(sr)),"radii_m":RADII,"tail_quantile":TAIL},
 "scale_summaries":scale_summaries,
 "pre_dedup_followup_total":len(all_flags),
 "frozen_candidate_count":len(kept),
 "frozen_candidates":kept,
 "next_gate":"SI2_FROZEN_CANDIDATE_BACKSCATTER_REPLAY_V2",
 "claim_ceiling":"UNSEEN_SOLITARY_ISLANDS_MARINE_BLIND_BATHYMETRIC_ANOMALY_CANDIDATE_GENERATION_ONLY"
}
p=OUT/"JANUS-KUSTO-SOLITARY-ISLANDS-SI1-BLIND-BATHYMETRY-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"shape":input_meta["shape"],"valid_cells":input_meta["valid_cells"],
 "sample_centers":out["blind_sampling"]["sample_centers"],
 "pre_dedup":out["pre_dedup_followup_total"],"candidates":len(kept),
 "candidate_json_sha256":hashlib.sha256(raw.encode()).hexdigest(),
 "backscatter_values_read":False
},indent=2))
