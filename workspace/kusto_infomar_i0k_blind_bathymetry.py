#!/usr/bin/env python3
from __future__ import annotations

import gc, hashlib, json, math, tempfile, zipfile
from pathlib import Path

import numpy as np
import requests
import rasterio
from rasterio.windows import Window, from_bounds
from scipy import ndimage
from pyproj import Transformer
from shapely.geometry import Polygon
from shapely import contains_xy

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I0K-BLIND-BATHYMETRY-CANDIDATE-GENERATION-PREREG-2026-09-24-v1.0.json").read_text())
FR=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I0I-PAIRED-RASTER-MEMBER-SELECTION-FREEZE-2026-09-24-v1.0.json").read_text())
GEO=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I0J-PAIRED-GRID-GEOMETRY-METADATA-RECEIPT-2026-09-24-v1.0.json").read_text())
BR=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I0J-NATIVE-CELL-REPRESENTATION-BRIDGE-2026-09-24-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-INFOMAR-I0K-blind-bathy/1.0","Accept-Encoding":"identity"}
CELL=float(BR["effective_native_cell_m"])
MULT=[4,8,16,32]
RADII=[CELL*k for k in MULT]
STRIDE_CELLS=8
STRIDE=CELL*STRIDE_CELLS
TAIL=.995
DEDUP=CELL*50
MAX_CAND=100
MAXPX=max(MULT)
MAX_RADIUS=max(RADII)
FAMILIES={
 "RELIEF_ROUGHNESS":["local_relief_m","local_std_m"],
 "CURVATURE":["abs_laplacian"],
 "RADIAL":["abs_center_surround_dog"],
 "LINEARITY":["structure_tensor_anisotropy"],
 "TERRACE_STEP":["gradient_magnitude"],
 "TOPOLOGIC":["isolated_extremum_persistence"]
}

def download(url,path):
    h=hashlib.sha256(); n=0
    with requests.get(url,headers=UA,timeout=300,allow_redirects=True,stream=True) as r:
        r.raise_for_status()
        with open(path,"wb") as f:
            for b in r.iter_content(1024*1024):
                if not b: continue
                f.write(b); h.update(b); n+=len(b)
    return h.hexdigest(),n

def robust(v):
    med=float(np.median(v)); mad=float(np.median(np.abs(v-med)))
    return np.abs((v-med)/max(1e-12,1.4826*mad)),med,mad

def dense_rect_polygon(bounds,src_crs,dst_crs,n=65):
    left,bottom,right,top=map(float,bounds)
    xs1=np.linspace(left,right,n); ys1=np.full(n,bottom)
    xs2=np.full(n,right); ys2=np.linspace(bottom,top,n)
    xs3=np.linspace(right,left,n); ys3=np.full(n,top)
    xs4=np.full(n,left); ys4=np.linspace(top,bottom,n)
    xs=np.concatenate([xs1,xs2[1:],xs3[1:],xs4[1:]])
    ys=np.concatenate([ys1,ys2[1:],ys3[1:],ys4[1:]])
    tr=Transformer.from_crs(src_crs,dst_crs,always_xy=True)
    x,y=tr.transform(xs,ys)
    poly=Polygon(np.column_stack([x,y]))
    if not poly.is_valid:
        poly=poly.buffer(0)
    if poly.is_empty or not poly.is_valid:
        raise RuntimeError("heldout header polygon transformation failed")
    return poly

bath=FR["selected"]["bathymetry"]
expected=GEO["bathymetry"]
back=GEO["heldout_backscatter"]

with tempfile.TemporaryDirectory(prefix="kusto_infomar_i0k_") as td:
    td=Path(td)
    zpath=td/"bathymetry.zip"
    archive_sha,archive_bytes=download(bath["archive_url"],zpath)
    if archive_sha.lower()!=expected["archive_sha256"].lower():
        raise RuntimeError(f"bathymetry archive SHA drift {archive_sha}")
    if archive_bytes!=int(expected["archive_bytes"]):
        raise RuntimeError("bathymetry archive size drift")

    tifpath=td/"bathymetry.tif"
    member_sha=hashlib.sha256(); member_bytes=0
    with zipfile.ZipFile(zpath) as zf:
        zi=zf.getinfo(bath["member"])
        if f"{zi.CRC:08x}".lower()!=str(bath["crc32"]).lower(): raise RuntimeError("bathymetry CRC drift")
        if int(zi.file_size)!=int(bath["uncompressed_size"]): raise RuntimeError("bathymetry member size drift")
        with zf.open(zi) as src, open(tifpath,"wb") as dst:
            while True:
                b=src.read(1024*1024)
                if not b: break
                dst.write(b); member_sha.update(b); member_bytes+=len(b)
    if member_sha.hexdigest().lower()!=expected["member_sha256"].lower():
        raise RuntimeError("bathymetry member SHA drift")

    with rasterio.open(tifpath) as ds:
        if str(ds.crs)!=expected["crs"]: raise RuntimeError(f"unexpected bathymetry CRS {ds.crs}")
        rx,ry=abs(float(ds.transform.a)),abs(float(ds.transform.e))
        if max(abs(rx-CELL),abs(ry-CELL))>1e-8:
            raise RuntimeError(f"native cell bridge drift {(rx,ry)} vs {CELL}")

        heldout_poly=dense_rect_polygon(back["bounds"],back["crs"],str(ds.crs))
        domain=heldout_poly.buffer(-MAX_RADIUS)
        if domain.is_empty: raise RuntimeError("heldout geometry erodes to empty domain")

        # Read only the bathymetry window required to evaluate the held-out geometry domain.
        l,b,r,t=heldout_poly.bounds
        fw=from_bounds(l,b,r,t,transform=ds.transform)
        c0=max(0,int(math.floor(fw.col_off))); r0=max(0,int(math.floor(fw.row_off)))
        c1=min(ds.width,int(math.ceil(fw.col_off+fw.width))); r1=min(ds.height,int(math.ceil(fw.row_off+fw.height)))
        if c1<=c0 or r1<=r0: raise RuntimeError("empty bathymetry overlap window")
        win=Window(c0,r0,c1-c0,r1-r0)
        ma=ds.read(1,window=win,masked=True)
        z=np.asarray(ma.filled(np.nan),dtype=np.float32)
        valid=(~np.ma.getmaskarray(ma)) & np.isfinite(z)
        transform=ds.window_transform(win)
        source_shape=[int(ds.height),int(ds.width)]
        nodata=None if ds.nodata is None else float(ds.nodata)
        dtype=str(ds.dtypes[0])

    if not np.any(valid): raise RuntimeError("no valid bathymetry cells in frozen overlap window")
    nrows,ncols=z.shape
    fill=float(np.median(z[valid])); zf=z.copy(); zf[~valid]=fill

    rr,cc=np.indices(z.shape,dtype=np.int32)
    edge=np.minimum.reduce([rr,cc,nrows-1-rr,ncols-1-cc]); del rr,cc; gc.collect()
    if np.all(valid): dist=np.full(z.shape,MAXPX+1,dtype=np.float32)
    else: dist=ndimage.distance_transform_edt(valid).astype(np.float32)
    inner=valid & (edge>=MAXPX) & (dist>=MAXPX)
    del edge,dist; gc.collect()

    rows=np.arange(MAXPX,nrows-MAXPX,STRIDE_CELLS,dtype=np.int32)
    cols=np.arange(MAXPX,ncols-MAXPX,STRIDE_CELLS,dtype=np.int32)
    R,C=np.meshgrid(rows,cols,indexing="ij")
    sr=R.ravel(); sc=C.ravel(); ok=inner[sr,sc]
    sr=sr[ok]; sc=sc[ok]
    del R,C,inner,ok; gc.collect()

    xx,yy=rasterio.transform.xy(transform,sr,sc,offset="center")
    xx=np.asarray(xx,dtype=np.float64); yy=np.asarray(yy,dtype=np.float64)
    ingeo=contains_xy(domain,xx,yy)
    sr=sr[ingeo]; sc=sc[ingeo]; xx=xx[ingeo]; yy=yy[ingeo]
    if len(sr)<1000: raise RuntimeError(f"too few frozen-domain sample centers {len(sr)}")

    to_wgs=Transformer.from_crs(expected["crs"],"EPSG:4326",always_xy=True)
    lon,lat=to_wgs.transform(xx,yy)
    lon=np.asarray(lon,dtype=np.float64); lat=np.asarray(lat,dtype=np.float64)

    all_flags=[]; scale_summaries=[]
    for mult,radius in zip(MULT,RADII):
        rpx=int(mult); size=2*rpx+1; sigma=max(1.0,rpx/3.0)
        print("I0K radius",radius,"px",rpx,flush=True)
        vals={}

        mean=ndimage.uniform_filter(zf,size=size,mode="nearest")
        sq=zf*zf
        mean2=ndimage.uniform_filter(sq,size=size,mode="nearest")
        std=np.sqrt(np.maximum(0.0,mean2-mean*mean))
        vals["local_std_m"]=std[sr,sc].astype(np.float64)
        persist=np.abs(zf-mean)/(std+1e-6)
        vals["isolated_extremum_persistence"]=persist[sr,sc].astype(np.float64)
        del sq,mean2,std,persist; gc.collect()

        zmax=ndimage.maximum_filter(zf,size=size,mode="nearest")
        zmin=ndimage.minimum_filter(zf,size=size,mode="nearest")
        vals["local_relief_m"]=(zmax[sr,sc]-zmin[sr,sc]).astype(np.float64)
        del zmax,zmin; gc.collect()

        smooth=ndimage.gaussian_filter(zf,sigma=sigma,mode="nearest",truncate=3.0)
        lap=np.abs(ndimage.laplace(smooth,mode="nearest"))/(CELL*CELL)
        vals["abs_laplacian"]=lap[sr,sc].astype(np.float64)
        del lap; gc.collect()

        gy,gx=np.gradient(smooth,CELL,CELL)
        grad=np.hypot(gx,gy)
        vals["gradient_magnitude"]=grad[sr,sc].astype(np.float64)
        del grad,smooth; gc.collect()

        ts=max(1.0,rpx/2.0)
        jxx=ndimage.gaussian_filter(gx*gx,sigma=ts,mode="nearest",truncate=3.0)
        jyy=ndimage.gaussian_filter(gy*gy,sigma=ts,mode="nearest",truncate=3.0)
        jxy=ndimage.gaussian_filter(gx*gy,sigma=ts,mode="nearest",truncate=3.0)
        anis=np.sqrt((jxx-jyy)**2+4*jxy*jxy)/(jxx+jyy+1e-12)
        vals["structure_tensor_anisotropy"]=anis[sr,sc].astype(np.float64)
        del gx,gy,jxx,jyy,jxy,anis; gc.collect()

        g1=ndimage.gaussian_filter(zf,sigma=max(.8,rpx/4.0),mode="nearest",truncate=3.0)
        g2=ndimage.gaussian_filter(zf,sigma=max(1.2,rpx/1.5),mode="nearest",truncate=3.0)
        vals["abs_center_surround_dog"]=np.abs(g1[sr,sc]-g2[sr,sc]).astype(np.float64)
        del g1,g2,mean; gc.collect()

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
        gq=float(np.quantile(generic,TAIL)); npre=0

        for i in range(len(sr)):
            votes=[]; detail={}
            for fam,keys in FAMILIES.items():
                ex=[k for k in keys if vals[k][i]>=thresholds[k]["q995"]]
                if ex: votes.append(fam); detail[fam]=ex
            if len(votes)>=2 or generic[i]>=gq:
                all_flags.append({
                  "source_row":r0+int(sr[i]),"source_col":c0+int(sc[i]),
                  "local_row":int(sr[i]),"local_col":int(sc[i]),
                  "x_m":float(xx[i]),"y_m":float(yy[i]),"lon":float(lon[i]),"lat":float(lat[i]),
                  "radius_m":float(radius),"radius_cells":int(mult),
                  "family_votes":votes,"extreme_metrics_by_family":detail,
                  "generic_robust_score":float(generic[i]),
                  "metrics":{k:float(vals[k][i]) for k in vals}
                }); npre+=1
        scale_summaries.append({
          "radius_cells":int(mult),"radius_m":float(radius),"sample_centers":int(len(sr)),
          "generic_q995":gq,"metric_thresholds":thresholds,"followup_count_pre_dedup":npre
        })
        del vals,rz,stack,generic; gc.collect()

    all_flags.sort(key=lambda q:(-len(q["family_votes"]),-q["generic_robust_score"],q["radius_cells"],q["source_row"],q["source_col"]))
    kept=[]
    for q in all_flags:
        if any(math.hypot(q["x_m"]-p["x_m"],q["y_m"]-p["y_m"])<DEDUP for p in kept): continue
        q=dict(q); q["candidate_id"]=f"INFOMAR_CB13_02_C{len(kept)+1:03d}"; q["blind_rank"]=len(kept)+1
        kept.append(q)
        if len(kept)>=MAX_CAND: break

    out={
      "artifact_id":"JANUS-KUSTO-INFOMAR-I0K-BLIND-BATHYMETRY-RUN-2026-09-24-v1.0",
      "prereg":PRE["artifact_id"],
      "representation_bridge":BR["artifact_id"],
      "source":{
        "survey_id":"CB13_02","bathymetry_archive_sha256":archive_sha,"bathymetry_archive_bytes":archive_bytes,
        "bathymetry_member":bath["member"],"bathymetry_member_sha256":member_sha.hexdigest(),
        "source_shape":source_shape,"window":[r0,c0,int(z.shape[0]),int(z.shape[1])],
        "window_transform":list(transform)[:6],"window_valid_cells":int(np.sum(valid)),
        "dtype":dtype,"nodata":nodata,"crs":expected["crs"]
      },
      "blind_domain":{
        "source":"heldout_backscatter_header_geometry_only",
        "backscatter_wgs84_bounds":back["bounds"],
        "domain_inward_buffer_m":MAX_RADIUS,
        "sample_centers":int(len(sr)),
        "native_cell_m":CELL,"radii_cells":MULT,"radii_m":RADII,
        "sampling_stride_cells":STRIDE_CELLS,"sampling_stride_m":STRIDE,
        "tail_quantile":TAIL,"spatial_dedup_m":DEDUP
      },
      "truth_firewall":{
        "backscatter_archive_downloaded":False,
        "backscatter_member_opened":False,
        "backscatter_values_read":False,
        "groundtruth_or_object_labels_read":False,
        "candidate_coordinates_changed_after_detection":False
      },
      "scale_summaries":scale_summaries,
      "pre_dedup_followup_total":len(all_flags),
      "spatial_dedup_m":DEDUP,
      "frozen_candidate_count":len(kept),
      "frozen_candidates":kept,
      "next_gate":"I0L_FROZEN_CANDIDATE_BACKSCATTER_REPLAY_V2",
      "claim_ceiling":"INFOMAR_BLIND_BATHYMETRIC_CANDIDATE_GENERATION_ONLY"
    }
    p=OUT/"JANUS-KUSTO-INFOMAR-I0K-BLIND-BATHYMETRY-RUN-2026-09-24-v1.0.json"
    raw=json.dumps(out,indent=2); p.write_text(raw)
    print(json.dumps({
      "artifact_id":out["artifact_id"],"sample_centers":out["blind_domain"]["sample_centers"],
      "pre_dedup":out["pre_dedup_followup_total"],"candidates":out["frozen_candidate_count"],
      "candidate_json_sha256":hashlib.sha256(raw.encode()).hexdigest(),
      "backscatter_archive_downloaded":False,"backscatter_values_read":False
    },indent=2))
