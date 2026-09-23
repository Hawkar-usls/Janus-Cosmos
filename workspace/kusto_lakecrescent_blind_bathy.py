#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import io
import json
import math
import zipfile
from pathlib import Path

import numpy as np
import requests
from scipy import ndimage

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"
OUT.mkdir(parents=True,exist_ok=True)

PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-LAKE-CRESCENT-3M-BLIND-BATHYMETRY-PREREG-2026-09-23-v1.0.json").read_text())
URL="https://www.sciencebase.gov/catalog/file/get/5eea74cb82ce3bd58d8572af?name=LakeCrescent_bathy_3m_UTM10_NAD83_NAVD88.zip"
UA={"User-Agent":"JANUS-KUSTO-LakeCrescent-blind-bathy/1.0"}

RADII=[12,24,48,96]
STRIDE_M=24
TAIL=0.995
DEDUP_M=100
MAX_CAND=100

FAMILIES={
    "RELIEF_ROUGHNESS":["local_relief_m","local_std_m"],
    "CURVATURE":["abs_laplacian"],
    "RADIAL":["abs_center_surround_dog"],
    "LINEARITY":["structure_tensor_anisotropy"],
    "TERRACE_STEP":["gradient_magnitude"],
    "TOPOLOGIC":["isolated_extremum_persistence"]
}

def download():
    r=requests.get(URL,headers=UA,timeout=180)
    r.raise_for_status()
    return r.content

def read_asc_from_zip(blob):
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        names=zf.namelist()
        asc=[n for n in names if n.lower().endswith(".asc")]
        if len(asc)!=1:
            raise RuntimeError(f"Expected exactly one ASC raster, found {asc}")
        with zf.open(asc[0]) as f:
            header={}
            for _ in range(6):
                line=f.readline().decode("utf-8",errors="replace").strip()
                k,v=line.split(None,1)
                header[k.lower()]=float(v)
            arr=np.loadtxt(f,dtype=np.float32)
        return asc[0],header,arr,names

def robust_z(v):
    med=float(np.median(v))
    mad=float(np.median(np.abs(v-med)))
    scale=max(1e-12,1.4826*mad)
    return np.abs((v-med)/scale),med,mad

blob=download()
zip_sha256=hashlib.sha256(blob).hexdigest()
asc_name,header,z,names=read_asc_from_zip(blob)

ncols=int(header["ncols"]); nrows=int(header["nrows"])
if z.shape!=(nrows,ncols):
    raise RuntimeError(f"Raster shape {z.shape} != header {(nrows,ncols)}")
cell=float(header["cellsize"])
if abs(cell-3.0)>1e-9:
    raise RuntimeError(f"Frozen native cell mismatch: {cell}")
nodata=float(header.get("nodata_value",-9999))
valid=np.isfinite(z) & (z!=nodata)
if not np.any(valid):
    raise RuntimeError("No valid raster cells")

fill=float(np.median(z[valid]))
zf=z.copy()
zf[~valid]=fill

max_r_px=int(round(max(RADII)/cell))
rr,cc=np.indices(z.shape)
edge=np.minimum.reduce([rr,cc,nrows-1-rr,ncols-1-cc])
# distance to nodata in pixels; all-valid raster gets a large constant cheaply.
if np.all(valid):
    dist_valid=np.full(z.shape,max_r_px+1,dtype=np.float32)
else:
    dist_valid=ndimage.distance_transform_edt(valid).astype(np.float32)
inner=valid & (edge>=max_r_px) & (dist_valid>=max_r_px)

stride_px=max(1,int(round(STRIDE_M/cell)))
rows=np.arange(max_r_px,nrows-max_r_px,stride_px,dtype=int)
cols=np.arange(max_r_px,ncols-max_r_px,stride_px,dtype=int)
R,C=np.meshgrid(rows,cols,indexing="ij")
sample_r=R.ravel(); sample_c=C.ravel()
sample_ok=inner[sample_r,sample_c]
sample_r=sample_r[sample_ok]; sample_c=sample_c[sample_ok]
if len(sample_r)<1000:
    raise RuntimeError(f"Too few blind sample centers: {len(sample_r)}")

all_flags=[]
scale_summaries=[]

for radius_m in RADII:
    rpx=max(1,int(round(radius_m/cell)))
    size=2*rpx+1
    sigma=max(1.0,rpx/3.0)

    mean=ndimage.uniform_filter(zf,size=size,mode="nearest")
    mean2=ndimage.uniform_filter(zf*zf,size=size,mode="nearest")
    std=np.sqrt(np.maximum(0.0,mean2-mean*mean))
    zmax=ndimage.maximum_filter(zf,size=size,mode="nearest")
    zmin=ndimage.minimum_filter(zf,size=size,mode="nearest")
    relief=zmax-zmin

    smooth=ndimage.gaussian_filter(zf,sigma=sigma,mode="nearest",truncate=3.0)
    lap=np.abs(ndimage.laplace(smooth,mode="nearest"))/(cell*cell)
    g1=ndimage.gaussian_filter(zf,sigma=max(0.8,rpx/4.0),mode="nearest",truncate=3.0)
    g2=ndimage.gaussian_filter(zf,sigma=max(1.2,rpx/1.5),mode="nearest",truncate=3.0)
    dog=np.abs(g1-g2)

    gy,gx=np.gradient(smooth,cell,cell)
    grad=np.hypot(gx,gy)

    tensor_sigma=max(1.0,rpx/2.0)
    jxx=ndimage.gaussian_filter(gx*gx,sigma=tensor_sigma,mode="nearest",truncate=3.0)
    jyy=ndimage.gaussian_filter(gy*gy,sigma=tensor_sigma,mode="nearest",truncate=3.0)
    jxy=ndimage.gaussian_filter(gx*gy,sigma=tensor_sigma,mode="nearest",truncate=3.0)
    anis=np.sqrt((jxx-jyy)**2+4*jxy*jxy)/(jxx+jyy+1e-12)
    persistence=np.abs(zf-mean)/(std+1e-6)

    metric_grids={
        "local_relief_m":relief,
        "local_std_m":std,
        "abs_laplacian":lap,
        "abs_center_surround_dog":dog,
        "structure_tensor_anisotropy":anis,
        "gradient_magnitude":grad,
        "isolated_extremum_persistence":persistence
    }

    vals={k:g[sample_r,sample_c].astype(np.float64) for k,g in metric_grids.items()}
    thresholds={}
    robusts={}
    for k,v in vals.items():
        good=np.isfinite(v)
        if not np.all(good):
            v=np.where(good,v,np.nanmedian(v[good]))
            vals[k]=v
        q=float(np.quantile(v,TAIL))
        rz,med,mad=robust_z(v)
        thresholds[k]={"q995":q,"median":med,"mad":mad}
        robusts[k]=rz

    zstack=np.column_stack([robusts[k] for k in sorted(robusts)])
    generic=0.5*np.max(zstack,axis=1)+0.5*np.mean(zstack,axis=1)
    generic_q=float(np.quantile(generic,TAIL))

    fam_votes=[]
    for i in range(len(sample_r)):
        votes=[]
        detail={}
        for fam,keys in FAMILIES.items():
            ex=[k for k in keys if vals[k][i]>=thresholds[k]["q995"]]
            if ex:
                votes.append(fam); detail[fam]=ex
        follow=(len(votes)>=2) or (generic[i]>=generic_q)
        if follow:
            row=int(sample_r[i]); col=int(sample_c[i])
            # ESRI ASCII uses xllcorner/yllcorner with first array row at the north edge.
            x=float(header.get("xllcorner",header.get("xllcenter",0.0))+(col+0.5)*cell)
            y=float(header.get("yllcorner",header.get("yllcenter",0.0))+(nrows-row-0.5)*cell)
            fam_votes.append({
                "row":row,"col":col,"x_utm10_nad83_m":x,"y_utm10_nad83_m":y,
                "radius_m":radius_m,
                "family_votes":votes,
                "extreme_metrics_by_family":detail,
                "generic_robust_score":float(generic[i]),
                "metrics":{k:float(vals[k][i]) for k in vals}
            })

    all_flags.extend(fam_votes)
    scale_summaries.append({
        "radius_m":radius_m,
        "sample_centers":int(len(sample_r)),
        "generic_q995":generic_q,
        "metric_thresholds":thresholds,
        "followup_count_pre_dedup":len(fam_votes)
    })

# Blind sorting only by frozen morphology evidence.
all_flags.sort(key=lambda q:(-len(q["family_votes"]),-q["generic_robust_score"],q["radius_m"],q["row"],q["col"]))
kept=[]
for q in all_flags:
    if any(math.hypot(q["x_utm10_nad83_m"]-p["x_utm10_nad83_m"],q["y_utm10_nad83_m"]-p["y_utm10_nad83_m"])<DEDUP_M for p in kept):
        continue
    q=dict(q)
    q["candidate_id"]=f"LC_BATHY_C{len(kept)+1:03d}"
    kept.append(q)
    if len(kept)>=MAX_CAND:
        break

out={
    "artifact_id":"JANUS-KUSTO-LAKE-CRESCENT-3M-BLIND-BATHYMETRY-RUN-2026-09-23-v1.0",
    "prereg":PRE["artifact_id"],
    "input":{
        "url":URL,
        "zip_bytes":len(blob),
        "zip_sha256":zip_sha256,
        "asc_name":asc_name,
        "zip_members":names,
        "header":header,
        "shape":[nrows,ncols],
        "valid_cells":int(np.sum(valid))
    },
    "truth_firewall":{
        "backscatter_values_read":False,
        "feature_labels_read":False,
        "target_coordinates_read":False
    },
    "blind_sampling":{
        "stride_m":STRIDE_M,
        "sample_centers":int(len(sample_r)),
        "radii_m":RADII,
        "tail_quantile":TAIL
    },
    "scale_summaries":scale_summaries,
    "pre_dedup_followup_total":len(all_flags),
    "frozen_candidate_count":len(kept),
    "frozen_candidates":kept,
    "next_gate":"SEAL_CANDIDATES_IN_REPOSITORY_THEN_BACKSCATTER_REPLAY_WITHOUT_RECENTERING",
    "claim_ceiling":"BLIND_BATHYMETRIC_ANOMALY_CANDIDATE_GENERATION_ONLY"
}
p=OUT/"JANUS-KUSTO-LAKE-CRESCENT-3M-BLIND-BATHYMETRY-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
    "artifact_id":out["artifact_id"],
    "zip_sha256":zip_sha256,
    "shape":out["input"]["shape"],
    "valid_cells":out["input"]["valid_cells"],
    "sample_centers":out["blind_sampling"]["sample_centers"],
    "scale_followups":[{"radius_m":x["radius_m"],"pre_dedup":x["followup_count_pre_dedup"]} for x in scale_summaries],
    "pre_dedup_followup_total":out["pre_dedup_followup_total"],
    "frozen_candidate_count":out["frozen_candidate_count"],
    "backscatter_values_read":False,
    "feature_labels_read":False
},indent=2))
