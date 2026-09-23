#!/usr/bin/env python3
from __future__ import annotations

import gc, hashlib, json, math, os, random, zipfile
from pathlib import Path

import numpy as np
import requests
import rasterio
from rasterio.windows import from_bounds
from rasterio.warp import reproject, Resampling, transform_bounds
from affine import Affine
from scipy import ndimage

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
TMP=OUT/"em1_tmp";TMP.mkdir(parents=True,exist_ok=True)

PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EM0-UNCHANGED-V2-PREREG-2026-09-24-v1.0.json").read_text())
DUAL=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EM0C-DUAL-METHOD-FREEZE-2026-09-24-v1.0.json").read_text())
MEM=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EM0B-ARCHIVE-MEMBER-RECEIPT-2026-09-24-v1.0.json").read_text())
DOMAIN=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EM0E-PAIRED-DOMAIN-REPROJECTION-FREEZE-2026-09-24-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-ElizabethMiddleton-EM1-blind/1.0"}
CELL=5.0
RADII=[20,40,80,160]
STRIDE=40.0
TAIL=.995
DEDUP=250.0
MAX_CAND=100
TARGET_CRS="EPSG:32757"
BUFFER_M=320.0
FAMILIES={
 "RELIEF_ROUGHNESS":["local_relief_m","local_std_m"],
 "CURVATURE":["abs_laplacian"],
 "RADIAL":["abs_center_surround_dog"],
 "LINEARITY":["structure_tensor_anisotropy"],
 "TERRACE_STEP":["gradient_magnitude"],
 "TOPOLOGIC":["isolated_extremum_persistence"]
}

def download_verified(spec,path):
    with requests.get(spec["url"],headers=UA,timeout=300,allow_redirects=True,stream=True) as r:
        r.raise_for_status()
        h=hashlib.sha256();n=0
        with open(path,"wb") as f:
            for chunk in r.iter_content(1024*1024):
                if chunk:
                    f.write(chunk);h.update(chunk);n+=len(chunk)
    if n!=int(spec["bytes"]):raise RuntimeError(f"archive size mismatch {n}")
    if h.hexdigest()!=spec["sha256"]:raise RuntimeError("archive sha mismatch")
    return h.hexdigest(),n

def snap_outward(bounds,cell):
    l,b,r,t=map(float,bounds)
    return (
      math.floor(l/cell)*cell,
      math.floor(b/cell)*cell,
      math.ceil(r/cell)*cell,
      math.ceil(t/cell)*cell
    )

def robust(v):
    med=float(np.median(v));mad=float(np.median(np.abs(v-med)))
    scale=max(1e-12,1.4826*mad)
    return np.abs((v-med)/scale),med,mad

bath_zip=TMP/"bathymetry.zip"
archive_sha,archive_bytes=download_verified(MEM["bathymetry_archive"],bath_zip)
member=MEM["bathymetry_archive"]["tif_members"][0]["name"]
spec=MEM["bathymetry_archive"]["tif_members"][0]
with zipfile.ZipFile(bath_zip) as zf:
    zi=zf.getinfo(member)
    if int(zi.file_size)!=int(spec["uncompressed_bytes"]):raise RuntimeError("member size mismatch")
    if int(zi.CRC)!=int(spec["crc32"]):raise RuntimeError("member crc mismatch")

uri=f"zip://{bath_zip}!{member}"
regions=[]
with rasterio.open(uri) as src:
    if str(src.crs)!="EPSG:4326":raise RuntimeError(f"unexpected bathy CRS {src.crs}")
    for rg in DOMAIN["paired_study_domain"]["regions"]:
        rid=rg["id"]; core=tuple(map(float,rg["core_bounds_m"]))
        expanded=(core[0]-BUFFER_M,core[1]-BUFFER_M,core[2]+BUFFER_M,core[3]+BUFFER_M)
        src_bounds=transform_bounds(TARGET_CRS,src.crs,*expanded,densify_pts=21)
        win=from_bounds(*src_bounds,transform=src.transform)
        win=win.round_offsets().round_lengths()
        # Clamp window to source dimensions.
        col0=max(0,int(win.col_off)); row0=max(0,int(win.row_off))
        col1=min(src.width,int(math.ceil(win.col_off+win.width)))
        row1=min(src.height,int(math.ceil(win.row_off+win.height)))
        win=rasterio.windows.Window(col0,row0,col1-col0,row1-row0)
        raw=src.read(1,window=win,masked=False)
        src_transform=src.window_transform(win)

        snapped=snap_outward(expanded,CELL)
        width=int(round((snapped[2]-snapped[0])/CELL))
        height=int(round((snapped[3]-snapped[1])/CELL))
        dst_transform=Affine(CELL,0.0,snapped[0],0.0,-CELL,snapped[3])
        dst=np.full((height,width),np.nan,dtype=np.float32)
        reproject(
          source=raw,destination=dst,
          src_transform=src_transform,src_crs=src.crs,src_nodata=src.nodata,
          dst_transform=dst_transform,dst_crs=TARGET_CRS,dst_nodata=np.nan,
          resampling=Resampling.nearest,num_threads=2
        )
        del raw;gc.collect()
        valid=np.isfinite(dst)
        if not np.any(valid):raise RuntimeError(f"{rid}: no valid bathy after projection")
        fill=float(np.median(dst[valid]))
        zf=dst.copy();zf[~valid]=fill
        maxpx=int(round(max(RADII)/CELL))
        dist=ndimage.distance_transform_edt(valid).astype(np.float32)
        # Sampling centers are restricted to the unexpanded paired footprint.
        step=max(1,int(round(STRIDE/CELL)))
        rows=np.arange(maxpx,height-maxpx,step,dtype=np.int32)
        cols=np.arange(maxpx,width-maxpx,step,dtype=np.int32)
        R,C=np.meshgrid(rows,cols,indexing="ij")
        sr=R.ravel();sc=C.ravel()
        xs=dst_transform.c+(sc+0.5)*CELL
        ys=dst_transform.f-(sr+0.5)*CELL
        in_core=(xs>=core[0])&(xs<=core[2])&(ys>=core[1])&(ys<=core[3])
        ok=in_core & valid[sr,sc] & (dist[sr,sc]>=maxpx)
        sr=sr[ok];sc=sc[ok]
        if len(sr)<1000:raise RuntimeError(f"{rid}: too few sample centers {len(sr)}")
        regions.append({
          "id":rid,"core_bounds_m":list(core),"expanded_bounds_m":list(expanded),
          "source_window":[int(win.col_off),int(win.row_off),int(win.width),int(win.height)],
          "projected_bounds_m":list(snapped),"transform":dst_transform,
          "z":dst,"zf":zf,"valid":valid,"sr":sr,"sc":sc
        })
        del dist,R,C,xs,ys,in_core,ok;gc.collect()

all_flags=[];scale_summaries=[]
for radius in RADII:
    print("EM1 radius",radius,flush=True)
    rpx=max(1,int(round(radius/CELL)));size=2*rpx+1;sigma=max(1.0,rpx/3.0)
    sampled={}
    offsets={}
    total=0
    for rg in regions:
        zf=rg["zf"]; sr=rg["sr"];sc=rg["sc"]
        vals={}
        mean=ndimage.uniform_filter(zf,size=size,mode="nearest")
        sq=zf*zf;mean2=ndimage.uniform_filter(sq,size=size,mode="nearest")
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
        lap=np.abs(ndimage.laplace(smooth,mode="nearest"))/(CELL*CELL)
        vals["abs_laplacian"]=lap[sr,sc].astype(np.float64);del lap;gc.collect()
        gy,gx=np.gradient(smooth,CELL,CELL)
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

        offsets[rg["id"]]=(total,total+len(sr))
        total+=len(sr)
        sampled[rg["id"]]=vals

    metrics=sorted(next(iter(sampled.values())).keys())
    pooled={m:np.concatenate([sampled[rg["id"]][m] for rg in regions]) for m in metrics}
    thresholds={}; rz={}
    for m,v in pooled.items():
        good=np.isfinite(v)
        if not np.all(good):
            if not np.any(good):raise RuntimeError(f"{radius} {m}: all invalid")
            v=np.where(good,v,np.nanmedian(v[good]));pooled[m]=v
        q=float(np.quantile(v,TAIL));zv,med,mad=robust(v)
        thresholds[m]={"q995":q,"median":med,"mad":mad};rz[m]=zv
    stack=np.column_stack([rz[m] for m in metrics])
    generic=.5*np.max(stack,axis=1)+.5*np.mean(stack,axis=1)
    gq=float(np.quantile(generic,TAIL))
    npre=0
    for rg in regions:
        lo,hi=offsets[rg["id"]]
        sr=rg["sr"];sc=rg["sc"];tr=rg["transform"]
        for ii in range(len(sr)):
            gi=lo+ii
            votes=[];detail={}
            for fam,keys in FAMILIES.items():
                ex=[m for m in keys if pooled[m][gi]>=thresholds[m]["q995"]]
                if ex:votes.append(fam);detail[fam]=ex
            if len(votes)>=2 or generic[gi]>=gq:
                row=int(sr[ii]);col=int(sc[ii])
                x=float(tr.c+(col+0.5)*CELL); y=float(tr.f-(row+0.5)*CELL)
                all_flags.append({
                  "region_id":rg["id"],"row":row,"col":col,"x_m":x,"y_m":y,"radius_m":radius,
                  "family_votes":votes,"extreme_metrics_by_family":detail,
                  "generic_robust_score":float(generic[gi]),
                  "metrics":{m:float(pooled[m][gi]) for m in metrics}
                });npre+=1
    scale_summaries.append({
      "radius_m":radius,
      "sample_centers_total":int(total),
      "sample_centers_by_region":{rg["id"]:int(len(rg["sr"])) for rg in regions},
      "generic_q995":gq,"metric_thresholds":thresholds,"followup_count_pre_dedup":npre
    })
    del sampled,pooled,rz,stack,generic;gc.collect()

all_flags.sort(key=lambda q:(-len(q["family_votes"]),-q["generic_robust_score"],q["radius_m"],q["region_id"],q["row"],q["col"]))
kept=[]
for q in all_flags:
    if any(math.hypot(q["x_m"]-p["x_m"],q["y_m"]-p["y_m"])<DEDUP for p in kept):continue
    q=dict(q);q["candidate_id"]=f"EM_BATHY_C{len(kept)+1:03d}";q["blind_rank"]=len(kept)+1
    kept.append(q)
    if len(kept)>=MAX_CAND:break

out={
 "artifact_id":"JANUS-KUSTO-ELIZABETH-MIDDLETON-EM1-BLIND-BATHYMETRY-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],
 "domain_freeze":DOMAIN["artifact_id"],
 "dual_method_freeze":DUAL["artifact_id"],
 "source":{
   "bathymetry_archive_sha256":archive_sha,"bathymetry_archive_bytes":archive_bytes,
   "bathymetry_member":member,"source_crs":"EPSG:4326","working_crs":TARGET_CRS,
   "reprojection":"NEAREST_TO_5M"
 },
 "regions":[{
   "id":rg["id"],"core_bounds_m":rg["core_bounds_m"],"expanded_bounds_m":rg["expanded_bounds_m"],
   "source_window":rg["source_window"],"projected_bounds_m":rg["projected_bounds_m"],
   "projected_shape":list(rg["z"].shape),"valid_cells":int(np.sum(rg["valid"])),
   "sample_centers":int(len(rg["sr"]))
 } for rg in regions],
 "truth_firewall":{
   "backscatter_values_read":False,"auv_imagery_read":False,"sediment_labels_read":False,
   "bruv_imagery_read":False,"candidate_coordinates_frozen":True,"candidate_radii_frozen":True
 },
 "blind_sampling":{"stride_m":STRIDE,"radii_m":RADII,"tail_quantile":TAIL,"threshold_pool":"BOTH_REGIONS"},
 "scale_summaries":scale_summaries,
 "pre_dedup_followup_total":len(all_flags),
 "frozen_candidate_count":len(kept),
 "frozen_candidates":kept,
 "next_gate":"EM2A_V2_AND_EM2B_V3_HELDOUT_BACKSCATTER_REPLAY",
 "claim_ceiling":"UNSEEN_INDEPENDENT_SENSOR_LINEAGE_BLIND_BATHYMETRIC_CANDIDATE_GENERATION_ONLY"
}
p=OUT/"JANUS-KUSTO-ELIZABETH-MIDDLETON-EM1-BLIND-BATHYMETRY-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "regions":[{"id":r["id"],"projected_shape":r["projected_shape"],"valid_cells":r["valid_cells"],"sample_centers":r["sample_centers"]} for r in out["regions"]],
 "pre_dedup":len(all_flags),"candidates":len(kept),
 "candidate_json_sha256":hashlib.sha256(raw.encode()).hexdigest(),
 "backscatter_values_read":False
},indent=2))
