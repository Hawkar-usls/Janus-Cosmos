#!/usr/bin/env python3
from __future__ import annotations

import gc, hashlib, io, json, math, zipfile
from pathlib import Path

import numpy as np
import requests
import rasterio
from rasterio.io import MemoryFile
from rasterio.windows import Window
from scipy import ndimage
from pyproj import Transformer

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-MONEY-SHOAL-UNSEEN-V2-PREREG-2026-09-23-v1.0.json").read_text())
COR=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-A0-PREREG-METHOD-INVARIANT-CORRECTION-2026-09-24-v1.0.json").read_text())
SRC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-A0E-EXACT-MEMBER-SELECTION-FREEZE-2026-09-24-v1.0.json").read_text())
BR=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-A0G-MONEY-SHOAL-REPRESENTATION-BRIDGE-RECEIPT-2026-09-24-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-Arafura-MoneyShoal-A1-blind-bathy/1.0"}
CELL=6.0
RADII=[24,48,96,192]
STRIDE=48.0
TAIL=.995
DEDUP=300.0
MAX_CAND=100
MAX_RADIUS=max(RADII)
FAMILIES={
 "RELIEF_ROUGHNESS":["local_relief_m","local_std_m"],
 "CURVATURE":["abs_laplacian"],
 "RADIAL":["abs_center_surround_dog"],
 "LINEARITY":["structure_tensor_anisotropy"],
 "TERRACE_STEP":["gradient_magnitude"],
 "TOPOLOGIC":["isolated_extremum_persistence"]
}
# Frozen Money Shoal held-out header rectangle from A0F/A0G. Header geometry only, no values.
BS_LEFT,BS_BOTTOM,BS_RIGHT,BS_TOP=241437.0,8850667.0,259581.0,8863317.0

def get(url):
    r=requests.get(url,headers=UA,timeout=300,allow_redirects=True); r.raise_for_status(); return r.content

def robust(v):
    med=float(np.median(v)); mad=float(np.median(np.abs(v-med)))
    scale=max(1e-12,1.4826*mad)
    return np.abs((v-med)/scale),med,mad

blob=get(SRC["bathymetry"]["archive_url"])
sha=hashlib.sha256(blob).hexdigest()
if sha.lower()!=SRC["bathymetry"]["archive_sha256"].lower(): raise RuntimeError("bathymetry archive SHA mismatch")
if len(blob)!=int(SRC["bathymetry"]["archive_bytes"]): raise RuntimeError("bathymetry archive byte length mismatch")
member=SRC["bathymetry"]["member"]
with zipfile.ZipFile(io.BytesIO(blob)) as zf:
    zi=zf.getinfo(member)
    if int(zi.file_size)!=int(SRC["bathymetry"]["member_uncompressed_bytes"]): raise RuntimeError("member size mismatch")
    if int(zi.CRC)!=int(SRC["bathymetry"]["member_crc32"]): raise RuntimeError("member CRC mismatch")
    tif=zf.read(member)

w=BR["conservative_bathymetry_source_window"]
win=Window(int(w["col_start"]),int(w["row_start"]),int(w["width"]),int(w["height"]))

with MemoryFile(tif) as mf:
  with mf.open() as ds:
    if str(ds.crs)!="EPSG:4326": raise RuntimeError(f"unexpected bathymetry CRS {ds.crs}")
    ma=ds.read(1,window=win,masked=True)
    z=np.asarray(ma.filled(np.nan),dtype=np.float32)
    valid=(~np.ma.getmaskarray(ma)) & np.isfinite(z)
    transform=ds.window_transform(win)
    nodata=ds.nodata
    dtype=str(ds.dtypes[0])
    source_shape=[ds.height,ds.width]
    source_transform=list(ds.transform)[:6]

if not np.any(valid): raise RuntimeError("no valid bathymetry cells in frozen window")
nrows,ncols=z.shape
fill=float(np.median(z[valid])); zf=z.copy(); zf[~valid]=fill
maxpx=int(round(MAX_RADIUS/CELL))

rr,cc=np.indices(z.shape,dtype=np.int32)
edge=np.minimum.reduce([rr,cc,nrows-1-rr,ncols-1-cc]); del rr,cc; gc.collect()
if np.all(valid): dist=np.full(z.shape,maxpx+1,dtype=np.float32)
else: dist=ndimage.distance_transform_edt(valid).astype(np.float32)
inner=valid & (edge>=maxpx) & (dist>=maxpx)
del edge,dist;gc.collect()

step=max(1,int(round(STRIDE/CELL)))
rows=np.arange(maxpx,nrows-maxpx,step,dtype=np.int32)
cols=np.arange(maxpx,ncols-maxpx,step,dtype=np.int32)
R,C=np.meshgrid(rows,cols,indexing="ij")
sr=R.ravel(); sc=C.ravel()
ok=inner[sr,sc]
sr=sr[ok]; sc=sc[ok]
del R,C,ok,inner;gc.collect()

# Convert only sample-center coordinates, never backscatter values.
lon,lat=rasterio.transform.xy(transform,sr,sc,offset="center")
lon=np.asarray(lon,dtype=np.float64); lat=np.asarray(lat,dtype=np.float64)
to_utm=Transformer.from_crs("EPSG:4326","EPSG:32753",always_xy=True)
ux,uy=to_utm.transform(lon,lat)
ux=np.asarray(ux); uy=np.asarray(uy)
domain=(ux>=BS_LEFT+MAX_RADIUS)&(ux<=BS_RIGHT-MAX_RADIUS)&(uy>=BS_BOTTOM+MAX_RADIUS)&(uy<=BS_TOP-MAX_RADIUS)
sr=sr[domain]; sc=sc[domain]; lon=lon[domain]; lat=lat[domain]; ux=ux[domain]; uy=uy[domain]
if len(sr)<1000: raise RuntimeError(f"too few frozen-domain sample centers {len(sr)}")

all_flags=[]; scale_summaries=[]
for radius in RADII:
    print("A1 radius",radius,flush=True)
    rpx=max(1,int(round(radius/CELL))); size=2*rpx+1; sigma=max(1.0,rpx/3.0)
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
    lap=np.abs(ndimage.laplace(smooth,mode="nearest"))/(CELL*CELL)
    vals["abs_laplacian"]=lap[sr,sc].astype(np.float64); del lap;gc.collect()

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

    thresholds={}; rz={}
    for k,v in vals.items():
        good=np.isfinite(v)
        if not np.all(good):
            if not np.any(good): raise RuntimeError(f"{radius} {k}: all invalid")
            v=np.where(good,v,np.nanmedian(v[good])); vals[k]=v
        q=float(np.quantile(v,TAIL)); zv,med,mad=robust(v)
        thresholds[k]={"q995":q,"median":med,"mad":mad}; rz[k]=zv

    stack=np.column_stack([rz[k] for k in sorted(rz)])
    generic=.5*np.max(stack,axis=1)+.5*np.mean(stack,axis=1)
    gq=float(np.quantile(generic,TAIL))
    npre=0
    for i in range(len(sr)):
        votes=[]; detail={}
        for fam,keys in FAMILIES.items():
            ex=[k for k in keys if vals[k][i]>=thresholds[k]["q995"]]
            if ex: votes.append(fam); detail[fam]=ex
        if len(votes)>=2 or generic[i]>=gq:
            local_row=int(sr[i]); local_col=int(sc[i])
            all_flags.append({
              "source_row":int(w["row_start"])+local_row,
              "source_col":int(w["col_start"])+local_col,
              "local_row":local_row,"local_col":local_col,
              "lon":float(lon[i]),"lat":float(lat[i]),
              "x_m":float(ux[i]),"y_m":float(uy[i]),
              "radius_m":radius,"family_votes":votes,"extreme_metrics_by_family":detail,
              "generic_robust_score":float(generic[i]),
              "metrics":{k:float(vals[k][i]) for k in vals}
            }); npre+=1
    scale_summaries.append({"radius_m":radius,"sample_centers":int(len(sr)),"generic_q995":gq,
                            "metric_thresholds":thresholds,"followup_count_pre_dedup":npre})
    del vals,rz,stack,generic;gc.collect()

all_flags.sort(key=lambda q:(-len(q["family_votes"]),-q["generic_robust_score"],q["radius_m"],q["source_row"],q["source_col"]))
kept=[]
for q in all_flags:
    if any(math.hypot(q["x_m"]-p["x_m"],q["y_m"]-p["y_m"])<DEDUP for p in kept): continue
    q=dict(q); q["candidate_id"]=f"ARAFURA_MONEY_C{len(kept)+1:03d}"; q["blind_rank"]=len(kept)+1
    kept.append(q)
    if len(kept)>=MAX_CAND: break

out={
 "artifact_id":"JANUS-KUSTO-ARAFURA-A1-MONEY-SHOAL-BLIND-BATHYMETRY-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],
 "method_invariant_correction":COR["artifact_id"],
 "source_freeze":SRC["artifact_id"],
 "representation_bridge":BR["artifact_id"],
 "input":{
   "archive_sha256":sha,"archive_bytes":len(blob),"member":member,
   "source_shape":source_shape,"source_transform":source_transform,
   "window_shape":list(z.shape),"window_transform":list(transform)[:6],
   "window_valid_cells":int(np.sum(valid)),"dtype":dtype,"nodata":nodata,
   "source_crs":"EPSG:4326","metric_crs":"EPSG:32753"
 },
 "blind_domain":{
   "money_shoal_header_rectangle_utm":[BS_LEFT,BS_BOTTOM,BS_RIGHT,BS_TOP],
   "boundary_shrink_m":MAX_RADIUS,
   "sample_centers":int(len(sr)),
   "nominal_cell_m":CELL,"stride_m":STRIDE,"radii_m":RADII,"tail_quantile":TAIL
 },
 "truth_firewall":{
   "backscatter_member_downloaded_or_opened":False,
   "backscatter_values_read":False,
   "seabed_samples_read":False,
   "underwater_imagery_read":False,
   "candidate_coordinates_changed_after_detection":False
 },
 "scale_summaries":scale_summaries,
 "pre_dedup_followup_total":len(all_flags),
 "spatial_dedup_m":DEDUP,
 "frozen_candidate_count":len(kept),
 "frozen_candidates":kept,
 "next_gate":"A2_FROZEN_CANDIDATE_BACKSCATTER_REPLAY_V2",
 "claim_ceiling":"UNSEEN_ARAFURA_MONEY_SHOAL_BLIND_BATHYMETRIC_CANDIDATE_GENERATION_ONLY"
}
p=OUT/"JANUS-KUSTO-ARAFURA-A1-MONEY-SHOAL-BLIND-BATHYMETRY-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2); p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"window_shape":out["input"]["window_shape"],
 "sample_centers":out["blind_domain"]["sample_centers"],
 "pre_dedup":out["pre_dedup_followup_total"],"candidates":len(kept),
 "candidate_json_sha256":hashlib.sha256(raw.encode()).hexdigest(),
 "backscatter_values_read":False
},indent=2))
