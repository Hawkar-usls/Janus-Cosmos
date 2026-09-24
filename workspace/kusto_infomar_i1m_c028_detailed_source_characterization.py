#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json, math, shutil, subprocess, tempfile
from pathlib import Path

import numpy as np
import requests
import rasterio
from rasterio.windows import Window
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"/"i1m_c028"
OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1M-C028-DETAILED-CHARACTERIZATION-PREREG-2026-09-25-v1.0.json").read_text())
IMP=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1M-C028-DETAILED-SOURCE-IMPLEMENTATION-FREEZE-2026-09-25-v1.0.json").read_text())
T=IMP["target"]
UA={"User-Agent":"JANUS-KUSTO-I1K-C028/1.0","Accept-Encoding":"identity"}

BATH={
 "url":"https://gsi.geodata.gov.ie/downloads/Marine/Data/Downloads/2013/CB13_03/BY_CB13_03_DingleBay_2m_U29N_LAT_TIFF_Inshore_Ireland.zip",
 "archive_sha":"a498b9dce20752ed34af195a8a544176d5d25b40e506b1d28a77eaabba15c98a",
 "archive_bytes":122289603,
 "member":"BY_CB13_03_Dingle_2m_U29N.tif",
 "member_sha":"d2820c7ac62a910acb8d87fbfa9cd83d72576d475efec75572b4c6765e1f5c10",
 "member_bytes":2122964187
}
BACK={
 "url":"https://gsi.geodata.gov.ie/downloads/Marine/Data/Downloads/2013/CB13_03/BS_CB13_03_Dingle_2m_Inshore_IE_WGS84_LAT_TIFF.zip",
 "archive_sha":"a92894da8a5b7172433000a45b1aff531a1778eb56a8e3379a1d724daf6b8bf3",
 "archive_bytes":96624179,
 "member":"BS_CB13_03_Dingle_2m_Inshore_IE_WGS84_LAT.tif",
 "member_sha":"8c99066aa90fcdadf4f3dd69a86f0ec934a1e15a9b28de17a5f712102c654fce",
 "member_bytes":2775010886
}

def download(url,path):
    h=hashlib.sha256(); n=0
    with requests.get(url,headers=UA,timeout=600,stream=True,allow_redirects=True) as r:
        r.raise_for_status()
        with open(path,"wb") as f:
            for b in r.iter_content(1024*1024):
                if not b: continue
                f.write(b); h.update(b); n+=len(b)
    return h.hexdigest(),n

def extract(spec,label,tmp):
    z=tmp/f"{label}.zip"; od=tmp/f"{label}_x"; od.mkdir()
    sha,n=download(spec["url"],z)
    if sha.lower()!=spec["archive_sha"] or n!=spec["archive_bytes"]:
        raise RuntimeError(f"{label} archive identity drift {sha} {n}")
    seven=shutil.which("7z")
    if not seven: raise RuntimeError("7z missing")
    subprocess.run([seven,"e","-y",f"-o{od}",str(z),spec["member"]],check=True,stdout=subprocess.DEVNULL)
    p=od/Path(spec["member"]).name
    h=hashlib.sha256(); m=0
    with open(p,"rb") as f:
        while True:
            b=f.read(8*1024*1024)
            if not b: break
            h.update(b);m+=len(b)
    if h.hexdigest().lower()!=spec["member_sha"] or m!=spec["member_bytes"]:
        raise RuntimeError(f"{label} member identity drift")
    return p,{"archive_sha256":sha,"archive_bytes":n,"member_sha256":h.hexdigest(),"member_bytes":m}

def read_metric_window(ds,row,col,half_m,dx,dy):
    hx=int(math.ceil(half_m/dx)); hy=int(math.ceil(half_m/dy))
    win=Window(col-hx,row-hy,2*hx+1,2*hy+1)
    ma=ds.read(1,window=win,masked=True,boundless=True)
    a=np.asarray(ma.data,dtype=np.float64)
    valid=(~np.ma.getmaskarray(ma)) & np.isfinite(a)
    rr=(np.arange(a.shape[0])-hy)*dy
    cc=(np.arange(a.shape[1])-hx)*dx
    yy,xx=np.meshgrid(rr,cc,indexing="ij")
    d=np.hypot(xx,yy)
    return a,valid,d,xx,yy,win

def finite_stats(v):
    v=np.asarray(v,dtype=np.float64); v=v[np.isfinite(v)]
    return {"count":int(v.size),"min":float(np.min(v)),"q05":float(np.quantile(v,.05)),"median":float(np.median(v)),"q95":float(np.quantile(v,.95)),"max":float(np.max(v)),"mean":float(np.mean(v)),"std":float(np.std(v))}

def profile_line_artifact(a,valid,center_r,center_c):
    # gradient profiles across +/-100 pixels around target on 21 parallel lines.
    span=min(100,center_c-1,a.shape[1]-center_c-2)
    vspan=min(100,center_r-1,a.shape[0]-center_r-2)
    row_stats=[]; row_grads=[]
    for off in range(-10,11):
        r=center_r+off
        z=a[r,center_c-span:center_c+span+1].copy()
        ok=valid[r,center_c-span:center_c+span+1]
        g=np.abs(np.diff(z)); gok=ok[:-1]&ok[1:]
        gv=np.where(gok,g,np.nan)
        row_stats.append(float(np.nanmax(gv)))
        row_grads.append(gv)
    col_stats=[]; col_grads=[]
    for off in range(-10,11):
        c=center_c+off
        z=a[center_r-vspan:center_r+vspan+1,c].copy()
        ok=valid[center_r-vspan:center_r+vspan+1,c]
        g=np.abs(np.diff(z)); gok=ok[:-1]&ok[1:]
        gv=np.where(gok,g,np.nan)
        col_stats.append(float(np.nanmax(gv)))
        col_grads.append(gv)
    def pct(vals):
        target=vals[10]
        return float(100.0*np.mean(np.asarray(vals)<=target))
    def med_corr(gs):
        tgt=gs[10]
        cors=[]
        for i,g in enumerate(gs):
            if i==10: continue
            ok=np.isfinite(tgt)&np.isfinite(g)
            if np.sum(ok)>=5 and np.std(tgt[ok])>0 and np.std(g[ok])>0:
                cors.append(float(np.corrcoef(tgt[ok],g[ok])[0,1]))
        return None if not cors else float(np.median(cors))
    return {
      "target_row_max_abs_first_difference":row_stats[10],
      "target_row_percentile_vs_21_parallel_rows":pct(row_stats),
      "median_parallel_row_gradient_correlation":med_corr(row_grads),
      "target_col_max_abs_first_difference":col_stats[10],
      "target_col_percentile_vs_21_parallel_cols":pct(col_stats),
      "median_parallel_col_gradient_correlation":med_corr(col_grads)
    }

def render(arr,valid,extent,title,path,cmap,center_circle_m):
    x=np.array(arr,dtype=float); x[~valid]=np.nan
    fig,ax=plt.subplots(figsize=(8,8))
    lo,hi=np.nanquantile(x,[.02,.98])
    im=ax.imshow(x,origin="upper",extent=extent,cmap=cmap,vmin=lo,vmax=hi,interpolation="nearest")
    ax.add_patch(Circle((0,0),center_circle_m,fill=False,linewidth=1.5))
    ax.scatter([0],[0],marker="+",s=70)
    ax.set_xlabel("metres east of frozen target");ax.set_ylabel("metres north of frozen target")
    ax.set_title(title)
    fig.colorbar(im,ax=ax,shrink=.75)
    fig.tight_layout();fig.savefig(path,dpi=160);plt.close(fig)

def render_hillshade(a,valid,dx,dy,extent,path,r):
    z=a.copy()
    med=np.nanmedian(np.where(valid,z,np.nan))
    z[~valid]=med
    gy,gx=np.gradient(z,dy,dx)
    slope=np.pi/2-np.arctan(np.hypot(gx,gy))
    aspect=np.arctan2(-gx,gy)
    az=np.deg2rad(315.0); alt=np.deg2rad(45.0)
    hs=np.sin(alt)*np.sin(slope)+np.cos(alt)*np.cos(slope)*np.cos(az-aspect)
    fig,ax=plt.subplots(figsize=(8,8))
    ax.imshow(hs,origin="upper",extent=extent,cmap="gray",interpolation="nearest")
    ax.add_patch(Circle((0,0),r,fill=False,linewidth=1.5))
    ax.scatter([0],[0],marker="+",s=70)
    ax.set_xlabel("metres east of frozen target");ax.set_ylabel("metres north of frozen target")
    ax.set_title("C028 bathymetry hillshade, ±250 m")
    fig.tight_layout();fig.savefig(path,dpi=160);plt.close(fig)

with tempfile.TemporaryDirectory(prefix="kusto_i1k_") as td:
    tmp=Path(td)
    bp,bid=extract(BATH,"bath",tmp)
    kp,kid=extract(BACK,"back",tmp)

    with rasterio.open(bp) as ds:
        brow,bcol=rasterio.transform.rowcol(ds.transform,float(T["x_m"]),float(T["y_m"]))
        dx=abs(float(ds.transform.a)); dy=abs(float(ds.transform.e))
        bath_meta={"crs":str(ds.crs),"shape":[ds.height,ds.width],"resolution_m":[dx,dy],"nodata":None if ds.nodata is None else float(ds.nodata),"tags":ds.tags()}
        windows={}
        for h in (50,100,250):
            windows[h]=read_metric_window(ds,int(brow),int(bcol),h,dx,dy)
        a50,v50,d50,x50,y50,_=windows[50]
        a100,v100,d100,x100,y100,_=windows[100]
        a250,v250,d250,x250,y250,_=windows[250]
        r=float(T["radius_m"])
        disk=(d50<=r)&v50
        ann=(d50>1.25*r)&(d50<=2*r)&v50
        disk_med=float(np.median(a50[disk])); ann_med=float(np.median(a50[ann]))
        signed=disk_med-ann_med
        sign=1 if signed>0 else -1 if signed<0 else 0
        extremum_mask=(d50<=50)&v50
        if sign>=0:
            idx=np.nanargmax(np.where(extremum_mask,a50,np.nan))
        else:
            idx=np.nanargmin(np.where(extremum_mask,a50,np.nan))
        er,ec=np.unravel_index(idx,a50.shape)
        extremum=float(a50[er,ec])
        extremum_dist=float(d50[er,ec])
        threshold=ann_med+0.5*(extremum-ann_med)

        # Component in +/-100 m window containing local extremum mapped by metric offsets.
        if sign>=0: mask=(a100>=threshold)&v100
        else: mask=(a100<=threshold)&v100
        labels,nlab=ndimage.label(mask,structure=np.ones((3,3),dtype=np.uint8))
        c100r=a100.shape[0]//2; c100c=a100.shape[1]//2
        # Convert extremum displacement to nearest +/-100m pixel.
        de=float(x50[er,ec]); dn=float(y50[er,ec])
        cer=int(round(c100r+dn/dy)); cec=int(round(c100c+de/dx))
        lab=int(labels[cer,cec])
        comp=(labels==lab) if lab>0 else np.zeros_like(labels,dtype=bool)
        area=float(np.sum(comp)*dx*dy)
        eqd=float(math.sqrt(4*area/math.pi)) if area>0 else 0.0
        touch=bool(np.any(comp[0,:]) or np.any(comp[-1,:]) or np.any(comp[:,0]) or np.any(comp[:,-1]))
        pts=np.column_stack([x100[comp],y100[comp]]) if np.any(comp) else np.empty((0,2))
        if len(pts)>=2:
            cen=pts.mean(axis=0); q=pts-cen
            cov=np.cov(q,rowvar=False)
            vals,vecs=np.linalg.eigh(cov); order=np.argsort(vals)[::-1]; vecs=vecs[:,order]
            proj=q@vecs
            spans=np.ptp(proj,axis=0)
            major=float(max(spans)); minor=float(min(spans)); aspect=None if minor<=0 else float(major/minor)
        else:
            major=minor=0.0;aspect=None
        radial=[]
        for lo in range(0,100,5):
            sel=(d100>=lo)&(d100<lo+5)&v100
            radial.append({"r0_m":lo,"r1_m":lo+5,"median":None if not np.any(sel) else float(np.median(a100[sel])),"count":int(np.sum(sel))})
        linecheck=profile_line_artifact(a250,v250,a250.shape[0]//2,a250.shape[1]//2)
        possible_line=bool(
          ((linecheck["target_row_percentile_vs_21_parallel_rows"]>=99 and (linecheck["median_parallel_row_gradient_correlation"] or -1)>=.90) or
           (linecheck["target_col_percentile_vs_21_parallel_cols"]>=99 and (linecheck["median_parallel_col_gradient_correlation"] or -1)>=.90))
        )
        if possible_line:
            morph="POSSIBLE_RASTER_LINE_ARTEFACT"
        elif touch:
            morph="BROAD_GEOLOGIC_RELIEF"
        elif area>0:
            morph="COMPACT_RELIEF"
        else:
            morph="UNRESOLVED"

        bath={
          "target_pixel":[int(brow),int(bcol)],
          "center_value":float(a50[a50.shape[0]//2,a50.shape[1]//2]),
          "disk_stats":finite_stats(a50[disk]),"annulus_stats":finite_stats(a50[ann]),
          "disk_median":disk_med,"annulus_median":ann_med,"signed_disk_minus_annulus_median":signed,
          "native_value_relief_sign":"POSITIVE" if sign>0 else "NEGATIVE" if sign<0 else "ZERO",
          "local_extremum_value_within_50m":extremum,"local_extremum_distance_from_target_m":extremum_dist,
          "half_relief_threshold":threshold,"half_relief_component_area_m2":area,"half_relief_equivalent_diameter_m":eqd,
          "half_relief_principal_axis_span_m":[major,minor],"half_relief_aspect_ratio":aspect,
          "half_relief_component_touches_100m_window":touch,
          "valid_fraction":{"50m":float(np.mean(v50)),"100m":float(np.mean(v100)),"250m":float(np.mean(v250))},
          "context_stats":{"50m":finite_stats(a50[v50]),"100m":finite_stats(a100[v100]),"250m":finite_stats(a250[v250])},
          "radial_median_profile_5m_bins":radial,
          "raster_line_artefact_check":linecheck,
          "preliminary_morphology_class":morph
        }
        extent=[float(np.min(x250)),float(np.max(x250)),float(np.min(y250)),float(np.max(y250))]
        render(a250,v250,extent,"C028 native bathymetry, ±250 m",OUT/"c028_bathymetry_native_250m.png","viridis",r)
        extent50=[float(np.min(x50)),float(np.max(x50)),float(np.min(y50)),float(np.max(y50))]
        render(a50,v50,extent50,"C028 native bathymetry, ±50 m",OUT/"c028_bathymetry_native_50m.png","viridis",r)
        render_hillshade(a250,v250,dx,dy,extent,OUT/"c028_bathymetry_hillshade_250m.png",r)

    with rasterio.open(kp) as ds:
        krow,kcol=rasterio.transform.rowcol(ds.transform,float(T["lon"]),float(T["lat"]))
        # Frozen physical axis metres at header midpoint.
        kdx=1.6349820375381001; kdy=2.775327722669724
        ka,kv,kd,kx,ky,_=read_metric_window(ds,int(krow),int(kcol),250,kdx,kdy)
        r=float(T["radius_m"]); kdisk=(kd<=r)&kv; kann=(kd>1.25*r)&(kd<=2*r)&kv
        back={
          "target_pixel":[int(krow),int(kcol)],"crs":str(ds.crs),"nodata":None if ds.nodata is None else float(ds.nodata),
          "disk_stats":finite_stats(ka[kdisk]),"annulus_stats":finite_stats(ka[kann]),
          "disk_median":float(np.median(ka[kdisk])),"annulus_median":float(np.median(ka[kann])),
          "signed_disk_minus_annulus_median":float(np.median(ka[kdisk])-np.median(ka[kann])),
          "context_250m_stats":finite_stats(ka[kv]),"valid_fraction_250m":float(np.mean(kv)),
          "raster_line_artefact_check":profile_line_artifact(ka,kv,ka.shape[0]//2,ka.shape[1]//2)
        }
        kextent=[float(np.min(kx)),float(np.max(kx)),float(np.min(ky)),float(np.max(ky))]
        render(ka,kv,kextent,"C028 native backscatter, ±250 m",OUT/"c028_backscatter_native_250m.png","gray",r)

    out={
      "artifact_id":"JANUS-KUSTO-INFOMAR-I1M-C028-DETAILED-SOURCE-CHARACTERIZATION-RUN-2026-09-25-v1.0",
      "prereg":PRE["artifact_id"],"implementation_freeze":IMP["artifact_id"],
      "target":T,"source_identity":{"bathymetry":bid,"backscatter":kid},
      "bathymetry_metadata":bath_meta,"bathymetry":bath,"backscatter":back,
      "context_blind":False,
      "target_or_radius_changed":False,
      "external_identity_claim":None,
      "claim_ceiling":"POSTRESULT_C028_DETAILED_SOURCE_CHARACTERIZATION_ONLY"
    }
    raw=json.dumps(out,indent=2,ensure_ascii=False)
    (OUT/"JANUS-KUSTO-INFOMAR-I1K-C028-TARGET-CENTERED-SOURCE-VALUE-CHARACTERIZATION-RUN-2026-09-25-v1.0.json").write_text(raw)
    print(json.dumps({
      "artifact_id":out["artifact_id"],
      "target":T,
      "bathymetry_center_value":bath["center_value"],
      "bathymetry_disk_median":bath["disk_median"],
      "bathymetry_annulus_median":bath["annulus_median"],
      "signed_disk_minus_annulus_median":bath["signed_disk_minus_annulus_median"],
      "local_extremum":bath["local_extremum_value_within_50m"],
      "local_extremum_distance_m":bath["local_extremum_distance_from_target_m"],
      "equivalent_diameter_m":bath["half_relief_equivalent_diameter_m"],
      "principal_axis_span_m":bath["half_relief_principal_axis_span_m"],
      "aspect_ratio":bath["half_relief_aspect_ratio"],
      "component_touches_100m_window":bath["half_relief_component_touches_100m_window"],
      "preliminary_morphology_class":bath["preliminary_morphology_class"],
      "backscatter_disk_median":back["disk_median"],
      "backscatter_annulus_median":back["annulus_median"],
      "output_json_sha256":hashlib.sha256(raw.encode()).hexdigest()
    },indent=2))
