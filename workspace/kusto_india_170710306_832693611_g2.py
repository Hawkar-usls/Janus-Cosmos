#!/usr/bin/env python3
import json, math, hashlib, re
from pathlib import Path
import numpy as np
import requests
import rasterio
from rasterio.windows import from_bounds
from pyproj import Geod
import matplotlib.pyplot as plt

LAT=17.0710306
LON=83.2693611
S=16.9710306
N=17.1710306
W=83.1693611
E=83.3693611
COG="https://data.source.coop/ausantarctic/gebco/GEBCO_2026.tif"
WMS="https://wms.gebco.net/mapserv?"
OUT=Path("workspace/kusto_global_groundtruth_out/india_g2")
OUT.mkdir(parents=True,exist_ok=True)
GEOD=Geod(ellps="WGS84")

def finite_or_none(x):
    x=float(x)
    return x if math.isfinite(x) else None

def get_tid():
    # WMS 1.3.0 EPSG:4326 axis order is lat,lon.
    half=.01
    params={
      "SERVICE":"WMS","VERSION":"1.3.0","REQUEST":"GetFeatureInfo",
      "LAYERS":"gebco_latest_tid_2","QUERY_LAYERS":"gebco_latest_tid_2",
      "CRS":"EPSG:4326",
      "BBOX":f"{LAT-half},{LON-half},{LAT+half},{LON+half}",
      "WIDTH":"101","HEIGHT":"101","I":"50","J":"50",
      "INFO_FORMAT":"text/plain"
    }
    rec={"params":params}
    try:
        r=requests.get(WMS,params=params,timeout=60,headers={"User-Agent":"JANUS-KUSTO-G2/1.0"})
        rec.update({"status":r.status_code,"url":r.url,"content_type":r.headers.get("content-type"),
                    "text":r.text[:2000]})
        m=re.search(r"(?i)(?:value|pixel|tid)[^0-9-]*(-?\d+)",r.text)
        rec["parsed_tid"]=int(m.group(1)) if m else None
    except Exception as e:
        rec["error"]=repr(e); rec["parsed_tid"]=None
    return rec

with rasterio.Env(
    GDAL_HTTP_MULTIRANGE="YES",
    GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
    CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif"
):
    with rasterio.open(COG) as ds:
        if ds.crs is None:
            raise RuntimeError("GEBCO COG has no CRS")
        win=from_bounds(W,S,E,N,transform=ds.transform).round_offsets().round_lengths()
        ma=ds.read(1,window=win,masked=True)
        arr=np.asarray(ma.data,dtype=np.float64)
        valid=(~np.ma.getmaskarray(ma)) & np.isfinite(arr)
        tr=ds.window_transform(win)
        if not valid.any():
            raise RuntimeError("No valid GEBCO cells in frozen window")

        # Target cell in full raster and local window
        full_row,full_col=ds.index(LON,LAT)
        local_row=int(full_row-int(win.row_off))
        local_col=int(full_col-int(win.col_off))
        if not (0<=local_row<arr.shape[0] and 0<=local_col<arr.shape[1]):
            raise RuntimeError("Target not inside extracted window")
        target_valid=bool(valid[local_row,local_col])
        target=float(arr[local_row,local_col]) if target_valid else None

        # Local metric pixel axes
        x0,y0=ds.xy(full_row,full_col)
        x1,y1=ds.xy(full_row,full_col+1)
        x2,y2=ds.xy(full_row+1,full_col)
        _,_,dx=GEOD.inv(x0,y0,x1,y1)
        _,_,dy=GEOD.inv(x0,y0,x2,y2)
        dx,dy=abs(dx),abs(dy)

        rr=(np.arange(arr.shape[0])-local_row)[:,None]*dy
        cc=(np.arange(arr.shape[1])-local_col)[None,:]*dx
        d=np.hypot(rr,cc)

        radii=[]
        for radius in [500,1000,2000,5000,10000]:
            disk=d<=radius
            ann=(d>1.25*radius)&(d<=2.0*radius)
            dv=arr[disk&valid]
            av=arr[ann&valid]
            rec={"radius_m":radius,"disk_valid_n":int(dv.size),"annulus_valid_n":int(av.size)}
            if dv.size:
                rec.update({
                  "local_min_m":float(np.min(dv)),
                  "local_max_m":float(np.max(dv)),
                  "local_relief_m":float(np.max(dv)-np.min(dv)),
                  "disk_median_m":float(np.median(dv)),
                  "target_local_percentile":float(np.mean(dv<=target)) if target_valid else None
                })
            if av.size:
                rec["annulus_median_m"]=float(np.median(av))
                if target_valid:
                    rec["target_minus_annulus_median_m"]=float(target-np.median(av))
            radii.append(rec)

        rings=[]
        edges=[0,1000,2000,3000,5000,7500,10000]
        for a,b in zip(edges[:-1],edges[1:]):
            m=(d>a)&(d<=b)&valid
            vals=arr[m]
            rings.append({"inner_m":a,"outer_m":b,"n":int(vals.size),
                          "median_m":float(np.median(vals)) if vals.size else None,
                          "min_m":float(np.min(vals)) if vals.size else None,
                          "max_m":float(np.max(vals)) if vals.size else None})

        # Shape proxies around 5 km: covariance of gradient magnitude locations in high-gradient tail
        gy,gx=np.gradient(arr,dy,dx)
        grad=np.hypot(gx,gy)
        m5=(d<=5000)&valid&np.isfinite(grad)
        gv=grad[m5]
        shape={}
        if gv.size:
            q95=float(np.quantile(gv,.95))
            yy,xx=np.nonzero(m5 & (grad>=q95))
            if len(xx)>=3:
                X=np.column_stack([(xx-local_col)*dx,(yy-local_row)*dy]).astype(float)
                cov=np.cov(X,rowvar=False)
                ev=np.linalg.eigvalsh(cov)
                anis=float((ev[-1]-ev[0])/(ev[-1]+ev[0]+1e-12))
            else: anis=None
            shape={"gradient_q95":q95,"high_gradient_count":int(len(xx)),"high_gradient_spatial_anisotropy":anis}

        tid=get_tid()

        out={
          "artifact_id":"JANUS-KUSTO-INDIA-170710306N-832693611E-G2-GEBCO2026-RUN-2026-09-25-v1.0",
          "target":{"lat":LAT,"lon":LON},
          "source":{
            "cog":COG,"crs":str(ds.crs),"shape":[int(ds.height),int(ds.width)],
            "native_resolution_deg":[abs(float(ds.transform.a)),abs(float(ds.transform.e))],
            "window":[int(win.row_off),int(win.col_off),int(win.height),int(win.width)]
          },
          "grid":{
            "pixel_axis_m_at_target":[dx,dy],
            "target_row_col_full":[full_row,full_col],
            "target_row_col_window":[local_row,local_col],
            "target_valid":target_valid,
            "target_depth_m":target,
            "window_valid_fraction":float(valid.mean()),
            "window_min_m":float(np.min(arr[valid])),
            "window_max_m":float(np.max(arr[valid])),
            "radii":radii,
            "radial_profile":rings,
            "shape_proxy_5km":shape
          },
          "tid_query":tid,
          "claim_ceiling":"GEBCO2026_15_ARCSECOND_REGIONAL_MORPHOLOGY_ONLY"
        }
        raw=json.dumps(out,indent=2,allow_nan=False)
        (OUT/"JANUS-KUSTO-INDIA-170710306N-832693611E-G2-GEBCO2026-RUN-2026-09-25-v1.0.json").write_text(raw)

        # Image: native data only, no interpolation
        vals=arr[valid]
        lo=float(np.quantile(vals,.02)); hi=float(np.quantile(vals,.98))
        fig,ax=plt.subplots(figsize=(7,6))
        im=ax.imshow(np.where(valid,arr,np.nan),origin="upper",cmap="terrain",vmin=lo,vmax=hi,interpolation="nearest")
        ax.plot(local_col,local_row,"r+",ms=15,mew=2)
        ax.set_title(f"GEBCO 2026 | {LAT:.7f}, {LON:.7f}\n15 arc-sec native cells")
        ax.set_axis_off()
        fig.colorbar(im,ax=ax,label="elevation (m)")
        fig.tight_layout()
        fig.savefig(OUT/"gebco2026_target_window.png",dpi=180,bbox_inches="tight")
        plt.close(fig)

        print(json.dumps({
          "target_depth_m":target,
          "pixel_axis_m":[dx,dy],
          "window_min_m":out["grid"]["window_min_m"],
          "window_max_m":out["grid"]["window_max_m"],
          "radii":radii,
          "radial_profile":rings,
          "shape_proxy_5km":shape,
          "tid_query":tid,
          "output_sha256":hashlib.sha256(raw.encode()).hexdigest()
        },indent=2))

