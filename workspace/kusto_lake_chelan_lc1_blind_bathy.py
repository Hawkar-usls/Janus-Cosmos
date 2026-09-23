#!/usr/bin/env python3
from __future__ import annotations

import gc
import hashlib
import io
import json
import math
import zipfile
from pathlib import Path

import numpy as np
import requests
import rasterio
from affine import Affine
from rasterio.io import MemoryFile
from scipy import ndimage

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"
OUT.mkdir(parents=True,exist_ok=True)

PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-LAKE-CHELAN-UNSEEN-DUALCHANNEL-PREREG-2026-09-23-v1.0.json").read_text())
SRC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-LAKE-CHELAN-LC0-SOURCE-RESOLVER-RECEIPT-2026-09-23-v1.0.json").read_text())
VAL=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-LAKE-CHELAN-LC1-RASTER-VALIDITY-ADDENDUM-2026-09-23-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-LakeChelan-LC1-blind-bathy/1.0"}
RADII=[12,24,48,96]
STRIDE_M=24.0
TAIL=0.995
DEDUP_M=100.0
MAX_CAND=100
CELL_EXPECT=3.0

FAMILIES={
    "RELIEF_ROUGHNESS":["local_relief_m","local_std_m"],
    "CURVATURE":["abs_laplacian"],
    "RADIAL":["abs_center_surround_dog"],
    "LINEARITY":["structure_tensor_anisotropy"],
    "TERRACE_STEP":["gradient_magnitude"],
    "TOPOLOGIC":["isolated_extremum_persistence"]
}

def download(url):
    r=requests.get(url,headers=UA,timeout=180)
    r.raise_for_status()
    return r.content

def read_world_file(text):
    vals=[float(x.strip()) for x in text.splitlines() if x.strip()]
    if len(vals)!=6:
        raise RuntimeError(f"World file must contain six numeric lines, got {vals}")
    A,D,B,E,C,F=vals
    return vals,Affine(A,B,C-0.5*A-0.5*B,D,E,F-0.5*D-0.5*E)

def read_tile(blob):
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        names=zf.namelist()
        tifs=[n for n in names if n.lower().endswith((".tif",".tiff")) and not n.startswith("__MACOSX/") and "/._" not in n]
        tfws=[n for n in names if n.lower().endswith((".tfw",".tifw",".wld")) and not n.startswith("__MACOSX/") and "/._" not in n]
        if len(tifs)!=1:
            raise RuntimeError(f"Expected one TIFF, found {tifs}")
        if len(tfws)>1:
            raise RuntimeError(f"Expected at most one TFW, found {tfws}")
        tif_bytes=zf.read(tifs[0])
        tfw_text=zf.read(tfws[0]).decode("utf-8",errors="replace") if tfws else None

    with MemoryFile(tif_bytes) as mf:
        with mf.open() as src:
            ma=src.read(1,masked=True)
            data=np.asarray(ma.filled(np.nan),dtype=np.float32)
            native_mask=np.ma.getmaskarray(ma)
            valid=(~native_mask) & np.isfinite(data)
            embedded_transform=src.transform
            embedded_crs=str(src.crs) if src.crs is not None else None
            nodata=src.nodata
            dtype=str(src.dtypes[0])
            height,width=src.height,src.width

    embedded_res=(abs(float(embedded_transform.a)),abs(float(embedded_transform.e)))
    embedded_good=(not embedded_transform.is_identity and max(abs(embedded_res[0]-CELL_EXPECT),abs(embedded_res[1]-CELL_EXPECT))<=1e-6)
    world_vals=None
    if embedded_good:
        transform=embedded_transform
        georef_source="EMBEDDED_TRANSFORM"
    else:
        if tfw_text is None:
            raise RuntimeError(f"Embedded transform not usable ({embedded_transform}) and no world file")
        world_vals,transform=read_world_file(tfw_text)
        georef_source="WORLD_FILE"

    res=(abs(float(transform.a)),abs(float(transform.e)))
    if max(abs(res[0]-CELL_EXPECT),abs(res[1]-CELL_EXPECT))>1e-6:
        raise RuntimeError(f"Unexpected effective cell size {res}")
    return {
      "data":data,"valid":valid,"transform":transform,
      "tif_name":tifs[0],"tfw_name":tfws[0] if tfws else None,
      "zip_members":names,"world_values":world_vals,
      "embedded_transform":list(embedded_transform)[:6],
      "effective_transform":list(transform)[:6],
      "georef_source":georef_source,"embedded_crs":embedded_crs,
      "nodata":nodata,"dtype":dtype,"shape":[height,width],"resolution_m":list(res)
    }

def robust_z(v):
    med=float(np.median(v))
    mad=float(np.median(np.abs(v-med)))
    scale=max(1e-12,1.4826*mad)
    return np.abs((v-med)/scale),med,mad

def candidate_xy(transform,row,col):
    x,y=rasterio.transform.xy(transform,row,col,offset="center")
    return float(x),float(y)

tile_results={}
all_candidates=[]

for tile in sorted(SRC["files"]):
    spec=SRC["files"][tile]["bathy"]
    print("TILE",tile,"download",flush=True)
    blob=download(spec["url"])
    md5=hashlib.md5(blob).hexdigest()
    sha256=hashlib.sha256(blob).hexdigest()
    if md5.lower()!=spec["md5"].lower():
        raise RuntimeError(f"{tile}: bathy MD5 mismatch {md5} != {spec['md5']}")
    if len(blob)!=int(spec["bytes"]):
        raise RuntimeError(f"{tile}: bathy byte length mismatch {len(blob)} != {spec['bytes']}")

    rr=read_tile(blob)
    z=rr["data"];valid=rr["valid"];transform=rr["transform"]
    nrows,ncols=z.shape
    if not np.any(valid):
        raise RuntimeError(f"{tile}: no valid raster cells")

    fill=float(np.median(z[valid]))
    zf=z.copy()
    zf[~valid]=fill
    max_r_px=int(round(max(RADII)/CELL_EXPECT))
    iy,ix=np.indices(z.shape)
    edge=np.minimum.reduce([iy,ix,nrows-1-iy,ncols-1-ix])
    if np.all(valid):
        dist_valid=np.full(z.shape,max_r_px+1,dtype=np.float32)
    else:
        dist_valid=ndimage.distance_transform_edt(valid).astype(np.float32)
    inner=valid & (edge>=max_r_px) & (dist_valid>=max_r_px)

    stride_px=max(1,int(round(STRIDE_M/CELL_EXPECT)))
    rows=np.arange(max_r_px,nrows-max_r_px,stride_px,dtype=int)
    cols=np.arange(max_r_px,ncols-max_r_px,stride_px,dtype=int)
    R,C=np.meshgrid(rows,cols,indexing="ij")
    sample_r=R.ravel();sample_c=C.ravel()
    ok=inner[sample_r,sample_c]
    sample_r=sample_r[ok];sample_c=sample_c[ok]
    if len(sample_r)<int(VAL["blind_sampling_transfer"]["minimum_sample_centers_per_tile"]):
        raise RuntimeError(f"{tile}: too few blind sample centers {len(sample_r)}")

    flags=[]
    scale_summaries=[]
    for radius_m in RADII:
        rpx=max(1,int(round(radius_m/CELL_EXPECT)))
        size=2*rpx+1
        sigma=max(1.0,rpx/3.0)

        mean=ndimage.uniform_filter(zf,size=size,mode="nearest")
        mean2=ndimage.uniform_filter(zf*zf,size=size,mode="nearest")
        std=np.sqrt(np.maximum(0.0,mean2-mean*mean))
        relief=ndimage.maximum_filter(zf,size=size,mode="nearest")-ndimage.minimum_filter(zf,size=size,mode="nearest")
        smooth=ndimage.gaussian_filter(zf,sigma=sigma,mode="nearest",truncate=3.0)
        lap=np.abs(ndimage.laplace(smooth,mode="nearest"))/(CELL_EXPECT*CELL_EXPECT)
        g1=ndimage.gaussian_filter(zf,sigma=max(0.8,rpx/4.0),mode="nearest",truncate=3.0)
        g2=ndimage.gaussian_filter(zf,sigma=max(1.2,rpx/1.5),mode="nearest",truncate=3.0)
        dog=np.abs(g1-g2)
        gy,gx=np.gradient(smooth,CELL_EXPECT,CELL_EXPECT)
        grad=np.hypot(gx,gy)
        ts=max(1.0,rpx/2.0)
        jxx=ndimage.gaussian_filter(gx*gx,sigma=ts,mode="nearest",truncate=3.0)
        jyy=ndimage.gaussian_filter(gy*gy,sigma=ts,mode="nearest",truncate=3.0)
        jxy=ndimage.gaussian_filter(gx*gy,sigma=ts,mode="nearest",truncate=3.0)
        anis=np.sqrt((jxx-jyy)**2+4*jxy*jxy)/(jxx+jyy+1e-12)
        persist=np.abs(zf-mean)/(std+1e-6)

        grids={
          "local_relief_m":relief,
          "local_std_m":std,
          "abs_laplacian":lap,
          "abs_center_surround_dog":dog,
          "structure_tensor_anisotropy":anis,
          "gradient_magnitude":grad,
          "isolated_extremum_persistence":persist
        }
        vals={k:g[sample_r,sample_c].astype(np.float64) for k,g in grids.items()}
        thresholds={};robusts={}
        for k,v in vals.items():
            good=np.isfinite(v)
            if not np.all(good):
                if not np.any(good):
                    raise RuntimeError(f"{tile} radius {radius_m}: metric {k} all invalid")
                v=np.where(good,v,np.nanmedian(v[good]));vals[k]=v
            q=float(np.quantile(v,TAIL))
            rz,med,mad=robust_z(v)
            thresholds[k]={"q995":q,"median":med,"mad":mad}
            robusts[k]=rz

        zstack=np.column_stack([robusts[k] for k in sorted(robusts)])
        generic=0.5*np.max(zstack,axis=1)+0.5*np.mean(zstack,axis=1)
        generic_q=float(np.quantile(generic,TAIL))
        count_pre=0
        for i in range(len(sample_r)):
            votes=[];detail={}
            for fam,keys in FAMILIES.items():
                ex=[k for k in keys if vals[k][i]>=thresholds[k]["q995"]]
                if ex:
                    votes.append(fam);detail[fam]=ex
            if len(votes)>=2 or generic[i]>=generic_q:
                row=int(sample_r[i]);col=int(sample_c[i]);x,y=candidate_xy(transform,row,col)
                flags.append({
                  "tile":tile,"row":row,"col":col,"x_m":x,"y_m":y,
                  "radius_m":radius_m,"family_votes":votes,
                  "extreme_metrics_by_family":detail,
                  "generic_robust_score":float(generic[i]),
                  "metrics":{k:float(vals[k][i]) for k in vals}
                })
                count_pre+=1
        scale_summaries.append({
          "radius_m":radius_m,
          "sample_centers":int(len(sample_r)),
          "generic_q995":generic_q,
          "metric_thresholds":thresholds,
          "followup_count_pre_dedup":count_pre
        })

        del mean,mean2,std,relief,smooth,lap,g1,g2,dog,gy,gx,grad,jxx,jyy,jxy,anis,persist,grids,vals,robusts,zstack,generic
        gc.collect()

    flags.sort(key=lambda q:(-len(q["family_votes"]),-q["generic_robust_score"],q["radius_m"],q["row"],q["col"]))
    kept=[]
    for q in flags:
        if any(math.hypot(q["x_m"]-p["x_m"],q["y_m"]-p["y_m"])<DEDUP_M for p in kept):
            continue
        q=dict(q)
        q["candidate_id"]=f"CHELAN_{tile}_C{len(kept)+1:03d}"
        q["blind_rank"]=len(kept)+1
        kept.append(q)
        if len(kept)>=MAX_CAND:
            break
    all_candidates.extend(kept)

    tile_results[tile]={
      "input":{
        "url":spec["url"],"zip_bytes":len(blob),"zip_md5_verified":md5,"zip_sha256":sha256,
        "tif_name":rr["tif_name"],"tfw_name":rr["tfw_name"],"zip_members":rr["zip_members"],
        "shape":rr["shape"],"dtype":rr["dtype"],"nodata":rr["nodata"],
        "valid_cells":int(np.sum(valid)),"resolution_m":rr["resolution_m"],
        "georef_source":rr["georef_source"],"embedded_crs":rr["embedded_crs"],
        "embedded_transform":rr["embedded_transform"],"effective_transform":rr["effective_transform"],
        "world_file_values":rr["world_values"]
      },
      "sample_centers":int(len(sample_r)),
      "scale_summaries":scale_summaries,
      "pre_dedup_followup_total":len(flags),
      "frozen_candidate_count":len(kept),
      "frozen_candidates":kept
    }
    print(tile,"valid",int(np.sum(valid)),"sample",len(sample_r),"cand",len(kept),flush=True)
    del blob,z,valid,zf,edge,dist_valid,inner,R,C,sample_r,sample_c,rr
    gc.collect()

out={
 "artifact_id":"JANUS-KUSTO-LAKE-CHELAN-LC1-BLIND-BATHYMETRY-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "validity_addendum":VAL["artifact_id"],
 "source_receipt":SRC["artifact_id"],
 "truth_firewall":{
   "backscatter_zip_downloaded":False,
   "backscatter_values_read":False,
   "mass_transport_deposit_labels_read":False,
   "delta_labels_read":False
 },
 "tiles":tile_results,
 "frozen_candidate_total":len(all_candidates),
 "frozen_candidates":all_candidates,
 "next_gate":"LC2_FROZEN_CANDIDATE_BACKSCATTER_REPLAY",
 "claim_ceiling":"UNSEEN_BLIND_BATHYMETRIC_ANOMALY_CANDIDATE_GENERATION_ONLY"
}
p=OUT/"JANUS-KUSTO-LAKE-CHELAN-LC1-BLIND-BATHYMETRY-RUN-2026-09-23-v1.0.json"
raw=json.dumps(out,indent=2)
p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "candidate_total":out["frozen_candidate_total"],
 "tiles":{k:{
   "valid_cells":v["input"]["valid_cells"],
   "sample_centers":v["sample_centers"],
   "pre_dedup":v["pre_dedup_followup_total"],
   "candidates":v["frozen_candidate_count"],
   "zip_sha256":v["input"]["zip_sha256"]
 } for k,v in tile_results.items()},
 "backscatter_values_read":False,
 "mtd_labels_read":False,
 "candidate_json_sha256":hashlib.sha256(raw.encode()).hexdigest()
},indent=2))
