#!/usr/bin/env python3
from __future__ import annotations
import gc, hashlib, json, math, os, subprocess, tempfile
from pathlib import Path
import numpy as np, requests, rasterio
from rasterio.windows import Window, from_bounds
from scipy import ndimage
from pyproj import Transformer
from shapely.geometry import Polygon
from shapely import contains_xy

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
FREEZE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1E-PER-SURVEY-NATIVE-CELL-BLIND-FREEZE-2026-09-24-v1.0.json").read_text())
GEO=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1D-BATCH-PAIRED-RASTER-HEADER-GEOMETRY-RECEIPT-2026-09-24-v1.0.json").read_text())\nBATCH=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1B-PREOUTCOME-BATCH-METADATA-AND-TRANSPORT-RECEIPT-2026-09-24-v1.0.json").read_text())

SID="CB13_03"
SREC=next(x for x in GEO["results"] if x["SURVEY_ID"]==SID)
CFG=FREEZE["survey_freezes"][SID]
BATH=SREC["bathymetry"]\nBROW=next(x for x in BATCH["prospective"] if x["SURVEY_ID"]==SID)
CELL=float(CFG["native_cell_m"]); MULT=[4,8,16,32]; RADII=list(map(float,CFG["radii_m"]))
STRIDE_CELLS=8; DEDUP=float(CFG["spatial_dedup_m"]); TAIL=.995; MAX_CAND=100; MAXPX=32
FAMILIES=FREEZE["bathymetry_discovery"]["families"]
UA={"User-Agent":"JANUS-KUSTO-INFOMAR-I1F-CB13_03/1.0","Accept-Encoding":"identity"}
METRICS=["local_relief_m","local_std_m","abs_laplacian","abs_center_surround_dog","structure_tensor_anisotropy","gradient_magnitude","isolated_extremum_persistence"]
MI={k:i for i,k in enumerate(METRICS)}

def download(url,path,expected):
    h=hashlib.sha256(); n=0
    with requests.get(url,headers=UA,stream=True,timeout=900,allow_redirects=True) as r:
        r.raise_for_status()
        with open(path,"wb") as f:
            for b in r.iter_content(4*1024*1024):
                if not b: continue
                f.write(b); h.update(b); n+=len(b)
    if n!=int(expected): raise RuntimeError(f"archive size drift {n} != {expected}")
    return h.hexdigest(),n

def sha256_file(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()

def extract_7z(zpath,member,outdir):
    cp=subprocess.run(["7z","x","-y",f"-o{outdir}",str(zpath),member],capture_output=True,text=True)
    if cp.returncode: raise RuntimeError(cp.stdout+"\n"+cp.stderr)
    p=Path(outdir)/member
    if not p.exists(): raise RuntimeError("member extraction missing")
    return p

def dense_rect(bounds,src,dst,n=129):
    l,b,r,t=map(float,bounds)
    xs=np.concatenate([np.linspace(l,r,n),np.full(n-1,r),np.linspace(r,l,n)[1:],np.full(n-1,l)])
    ys=np.concatenate([np.full(n,b),np.linspace(b,t,n)[1:],np.full(n-1,t),np.linspace(t,b,n)[1:]])
    tr=Transformer.from_crs(src,dst,always_xy=True); x,y=tr.transform(xs,ys)
    p=Polygon(np.column_stack([x,y]))
    if not p.is_valid: p=p.buffer(0)
    if p.is_empty: raise RuntimeError("empty transformed overlap polygon")
    return p

def exact_median(mm_valid_values):
    n=len(mm_valid_values)
    if n<1: raise RuntimeError("no valid values")
    if n&1:
        k=n//2; mm_valid_values.partition(k); return float(mm_valid_values[k])
    a=n//2-1; b=n//2; mm_valid_values.partition((a,b)); return float((float(mm_valid_values[a])+float(mm_valid_values[b]))/2.0)

with tempfile.TemporaryDirectory(prefix="kusto_i1f_cb1303_") as td0:
    td=Path(td0); zpath=td/"bathy.zip"; ex=td/"ex"; ex.mkdir()
    asha,abytes=download(BROW["bathymetry_url"],zpath,BROW["bathymetry_archive_size_bytes"])
    if asha.lower()!=BATH["archive_sha256"].lower(): raise RuntimeError("archive SHA drift")
    tif=extract_7z(zpath,BATH["member"],ex)
    if tif.stat().st_size!=int(BATH["member_bytes"]): raise RuntimeError("member size drift")
    msha=sha256_file(tif)
    if msha.lower()!=BATH["member_sha256"].lower(): raise RuntimeError("member SHA drift")

    with rasterio.open(tif) as ds:
        if str(ds.crs)!=BATH["crs"]: raise RuntimeError(f"CRS drift {ds.crs}")
        overlap=dense_rect(CFG["wgs84_overlap_bounds"],"EPSG:4326",str(ds.crs))
        domain=overlap.buffer(-max(RADII))
        if domain.is_empty: raise RuntimeError("eroded overlap domain empty")
        l,b,r,t=overlap.bounds
        fw=from_bounds(l,b,r,t,transform=ds.transform)
        c0=max(0,int(math.floor(fw.col_off))); r0=max(0,int(math.floor(fw.row_off)))
        c1=min(ds.width,int(math.ceil(fw.col_off+fw.width))); r1=min(ds.height,int(math.ceil(fw.row_off+fw.height)))
        H,W=r1-r0,c1-c0
        if H<=0 or W<=0: raise RuntimeError("empty overlap window")
        transform=ds.window_transform(Window(c0,r0,W,H))

        z=np.memmap(td/"z.f32",dtype="float32",mode="w+",shape=(H,W))
        valid=np.memmap(td/"valid.u8",dtype="uint8",mode="w+",shape=(H,W))
        valid_count=0
        for a in range(0,H,512):
            n=min(512,H-a)
            ma=ds.read(1,window=Window(c0,r0+a,W,n),masked=True)
            arr=np.asarray(ma.filled(np.nan),dtype=np.float32)
            vm=(~np.ma.getmaskarray(ma)) & np.isfinite(arr)
            z[a:a+n]=arr; valid[a:a+n]=vm.astype(np.uint8); valid_count+=int(vm.sum())
        z.flush(); valid.flush()

    vv=np.memmap(td/"valid_values.f32",dtype="float32",mode="w+",shape=(valid_count,))
    pos=0
    for a in range(0,H,512):
        n=min(512,H-a); vm=np.asarray(valid[a:a+n],bool); arr=np.asarray(z[a:a+n])
        vals=arr[vm]; vv[pos:pos+len(vals)]=vals; pos+=len(vals)
    fill=exact_median(vv); vv.flush(); del vv; os.remove(td/"valid_values.f32")

    zf=np.memmap(td/"zf.f32",dtype="float32",mode="w+",shape=(H,W))
    for a in range(0,H,512):
        n=min(512,H-a); arr=np.asarray(z[a:a+n]).copy(); vm=np.asarray(valid[a:a+n],bool)
        arr[~vm]=fill; zf[a:a+n]=arr
    zf.flush(); del z; os.remove(td/"z.f32")

    rows_all=np.arange(MAXPX,H-MAXPX,STRIDE_CELLS,dtype=np.int32)
    cols_all=np.arange(MAXPX,W-MAXPX,STRIDE_CELLS,dtype=np.int32)
    sr_parts=[]; sc_parts=[]; x_parts=[]; y_parts=[]
    for q in range(0,len(rows_all),128):
        rg=rows_all[q:q+128]
        if len(rg)==0: continue
        wa=max(0,int(rg[0])-MAXPX); wb=min(H,int(rg[-1])+MAXPX+1)
        vm=np.asarray(valid[wa:wb],dtype=bool)
        dist=ndimage.distance_transform_edt(vm)
        good=dist[(rg-wa)[:,None],cols_all[None,:]]>=MAXPX
        ii,jj=np.nonzero(good)
        if len(ii)==0: continue
        rr=rg[ii]; cc=cols_all[jj]
        xx=transform.c+(cc.astype(float)+.5)*transform.a
        yy=transform.f+(rr.astype(float)+.5)*transform.e
        geo=contains_xy(domain,xx,yy)
        if np.any(geo):
            sr_parts.append(rr[geo]); sc_parts.append(cc[geo]); x_parts.append(xx[geo]); y_parts.append(yy[geo])
    if not sr_parts: raise RuntimeError("no valid sample centers")
    sr=np.concatenate(sr_parts).astype(np.int32); sc=np.concatenate(sc_parts).astype(np.int32)
    xx=np.concatenate(x_parts).astype(np.float64); yy=np.concatenate(y_parts).astype(np.float64)
    order=np.lexsort((sc,sr)); sr,sc,xx,yy=sr[order],sc[order],xx[order],yy[order]
    N=len(sr)
    if N<1000: raise RuntimeError(f"too few sample centers {N}")
    to_wgs=Transformer.from_crs(BATH["crs"],"EPSG:4326",always_xy=True)

    all_flags=[]; scale_summaries=[]
    unique_rows,row_starts=np.unique(sr,return_index=True)
    row_ends=np.r_[row_starts[1:],N]

    for mult,radius in zip(MULT,RADII):
        rpx=int(mult); size=2*rpx+1; sigma=max(1.0,rpx/3.0); ts=max(1.0,rpx/2.0)
        halo=max(4*rpx+8,96)
        metrics=np.memmap(td/f"metrics_{rpx}.f32",dtype="float32",mode="w+",shape=(N,len(METRICS)))
        for ug in range(0,len(unique_rows),96):
            urows=unique_rows[ug:ug+96]
            i0=int(row_starts[ug]); i1=int(row_ends[min(ug+95,len(row_ends)-1)])
            if ug+96 < len(unique_rows): i1=int(row_starts[ug+96])
            ra=max(0,int(urows[0])-halo); rb=min(H,int(urows[-1])+halo+1)
            tile=np.asarray(zf[ra:rb],dtype=np.float32).copy()
            lr=sr[i0:i1]-ra; lc=sc[i0:i1]

            mean=ndimage.uniform_filter(tile,size=size,mode="nearest")
            sq=tile*tile; mean2=ndimage.uniform_filter(sq,size=size,mode="nearest")
            std=np.sqrt(np.maximum(0.0,mean2-mean*mean))
            metrics[i0:i1,MI["local_std_m"]]=std[lr,lc]
            metrics[i0:i1,MI["isolated_extremum_persistence"]]=np.abs(tile[lr,lc]-mean[lr,lc])/(std[lr,lc]+1e-6)
            del sq,mean2,std,mean; gc.collect()

            zmax=ndimage.maximum_filter(tile,size=size,mode="nearest")
            zmin=ndimage.minimum_filter(tile,size=size,mode="nearest")
            metrics[i0:i1,MI["local_relief_m"]]=zmax[lr,lc]-zmin[lr,lc]
            del zmax,zmin; gc.collect()

            smooth=ndimage.gaussian_filter(tile,sigma=sigma,mode="nearest",truncate=3.0)
            lap=ndimage.laplace(smooth,mode="nearest")
            metrics[i0:i1,MI["abs_laplacian"]]=np.abs(lap[lr,lc])/(CELL*CELL)
            del lap; gc.collect()
            gy,gx=np.gradient(smooth,CELL,CELL)
            metrics[i0:i1,MI["gradient_magnitude"]]=np.hypot(gx[lr,lc],gy[lr,lc])
            del smooth; gc.collect()

            tmp=gx*gx; jj=ndimage.gaussian_filter(tmp,sigma=ts,mode="nearest",truncate=3.0); jxxs=jj[lr,lc].astype(np.float64); del tmp,jj
            tmp=gy*gy; jj=ndimage.gaussian_filter(tmp,sigma=ts,mode="nearest",truncate=3.0); jyys=jj[lr,lc].astype(np.float64); del tmp,jj
            tmp=gx*gy; jj=ndimage.gaussian_filter(tmp,sigma=ts,mode="nearest",truncate=3.0); jxys=jj[lr,lc].astype(np.float64); del tmp,jj,gx,gy
            metrics[i0:i1,MI["structure_tensor_anisotropy"]]=np.sqrt((jxxs-jyys)**2+4*jxys*jxys)/(jxxs+jyys+1e-12)
            del jxxs,jyys,jxys; gc.collect()

            g1=ndimage.gaussian_filter(tile,sigma=max(.8,rpx/4.0),mode="nearest",truncate=3.0)
            a1=g1[lr,lc].astype(np.float64); del g1
            g2=ndimage.gaussian_filter(tile,sigma=max(1.2,rpx/1.5),mode="nearest",truncate=3.0)
            metrics[i0:i1,MI["abs_center_surround_dog"]]=np.abs(a1-g2[lr,lc])
            del a1,g2,tile; gc.collect()
        metrics.flush()

        thresholds={}; meds=np.empty(len(METRICS)); scales=np.empty(len(METRICS))
        for k in METRICS:
            v=np.asarray(metrics[:,MI[k]],dtype=np.float64)
            q=float(np.quantile(v,TAIL)); med=float(np.median(v)); mad=float(np.median(np.abs(v-med)))
            thresholds[k]={"q995":q,"median":med,"mad":mad}; meds[MI[k]]=med; scales[MI[k]]=max(1e-12,1.4826*mad)
        generic=np.memmap(td/f"generic_{rpx}.f32",dtype="float32",mode="w+",shape=(N,))
        for a in range(0,N,250000):
            b=min(N,a+250000)
            m=np.asarray(metrics[a:b],dtype=np.float64)
            rz=np.abs((m-meds)/scales)
            generic[a:b]=.5*np.max(rz,axis=1)+.5*np.mean(rz,axis=1)
        generic.flush(); gq=float(np.quantile(np.asarray(generic,dtype=np.float64),TAIL))
        npre=0
        for a in range(0,N,250000):
            b=min(N,a+250000); m=np.asarray(metrics[a:b]); g=np.asarray(generic[a:b])
            fams=[
              (m[:,MI["local_relief_m"]]>=thresholds["local_relief_m"]["q995"]) | (m[:,MI["local_std_m"]]>=thresholds["local_std_m"]["q995"]),
              m[:,MI["abs_laplacian"]]>=thresholds["abs_laplacian"]["q995"],
              m[:,MI["abs_center_surround_dog"]]>=thresholds["abs_center_surround_dog"]["q995"],
              m[:,MI["structure_tensor_anisotropy"]]>=thresholds["structure_tensor_anisotropy"]["q995"],
              m[:,MI["gradient_magnitude"]]>=thresholds["gradient_magnitude"]["q995"],
              m[:,MI["isolated_extremum_persistence"]]>=thresholds["isolated_extremum_persistence"]["q995"]
            ]
            fc=np.sum(np.column_stack(fams),axis=1)
            take=np.flatnonzero((fc>=2)|(g>=gq))
            if len(take):
                ids=a+take; lon,lat=to_wgs.transform(xx[ids],yy[ids])
                for zidx,gi,lo,la in zip(take,ids,lon,lat):
                    votes=[]; detail={}
                    for fam,keys in FAMILIES.items():
                        ex=[k for k in keys if float(m[zidx,MI[k]])>=thresholds[k]["q995"]]
                        if ex: votes.append(fam); detail[fam]=ex
                    all_flags.append({
                      "source_row":r0+int(sr[gi]),"source_col":c0+int(sc[gi]),
                      "local_row":int(sr[gi]),"local_col":int(sc[gi]),
                      "x_m":float(xx[gi]),"y_m":float(yy[gi]),"lon":float(lo),"lat":float(la),
                      "radius_cells":rpx,"radius_m":float(radius),"family_votes":votes,"extreme_metrics_by_family":detail,
                      "generic_robust_score":float(g[zidx]),"metrics":{k:float(m[zidx,MI[k]]) for k in METRICS}
                    })
                npre+=len(take)
        scale_summaries.append({"radius_cells":rpx,"radius_m":float(radius),"sample_centers":N,"generic_q995":gq,
                                "metric_thresholds":thresholds,"followup_count_pre_dedup":npre})
        del metrics,generic; gc.collect()
        os.remove(td/f"metrics_{rpx}.f32"); os.remove(td/f"generic_{rpx}.f32")

    all_flags.sort(key=lambda q:(-len(q["family_votes"]),-q["generic_robust_score"],q["radius_cells"],q["source_row"],q["source_col"]))
    kept=[]
    for q in all_flags:
        if any(math.hypot(q["x_m"]-p["x_m"],q["y_m"]-p["y_m"])<DEDUP for p in kept): continue
        q=dict(q); q["candidate_id"]=f"INFOMAR_CB13_03_C{len(kept)+1:03d}"; q["blind_rank"]=len(kept)+1
        kept.append(q)
        if len(kept)>=MAX_CAND: break

    out={
      "artifact_id":"JANUS-KUSTO-INFOMAR-I1F-CB13_03-BLIND-BATHYMETRY-RUN-2026-09-24-v1.0",
      "freeze":FREEZE["artifact_id"],"SURVEY_ID":SID,"rank":2,
      "source":{"archive_sha256":asha,"archive_bytes":abytes,"member":BATH["member"],"member_sha256":msha,
                "window":[r0,c0,H,W],"window_valid_cells":valid_count,"fill_median":fill,"crs":BATH["crs"]},
      "blind_domain":{"wgs84_overlap_bounds":CFG["wgs84_overlap_bounds"],"native_cell_m":CELL,"radii_cells":MULT,
                      "radii_m":RADII,"sampling_stride_cells":8,"spatial_dedup_m":DEDUP,"sample_centers":N},
      "implementation":{"storage":"DISK_MEMMAP_PLUS_EXACT_ROW_CHUNK_FILTERS","scientific_math_changed":False},
      "scale_summaries":scale_summaries,"pre_dedup_followup_total":len(all_flags),
      "frozen_candidate_count":len(kept),"frozen_candidates":kept,
      "truth_firewall":{"backscatter_values_read":False,"groundtruth_read":False,"candidate_coordinates_frozen":True},
      "next_gate":"I1G_CB13_03_FROZEN_CANDIDATE_BACKSCATTER_REPLAY_V2",
      "claim_ceiling":"INFOMAR_CB13_03_BLIND_BATHYMETRIC_CANDIDATE_GENERATION_ONLY"
    }
    raw=json.dumps(out,indent=2); p=OUT/"JANUS-KUSTO-INFOMAR-I1F-CB13_03-BLIND-BATHYMETRY-RUN-2026-09-24-v1.0.json"; p.write_text(raw)
    print(json.dumps({"artifact_id":out["artifact_id"],"sample_centers":N,"pre_dedup":len(all_flags),"candidates":len(kept),
                      "candidate_json_sha256":hashlib.sha256(raw.encode()).hexdigest(),"backscatter_values_read":False},indent=2))
