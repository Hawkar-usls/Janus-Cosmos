#!/usr/bin/env python3
import json, math, hashlib, time, tempfile
from pathlib import Path
import numpy as np, requests, rasterio
from pyproj import Geod
import matplotlib.pyplot as plt

LAT=17.0710306; LON=83.2693611
S=16.9710306; N=17.1710306; W=83.1693611; E=83.3693611
OUT=Path("workspace/kusto_global_groundtruth_out/india_g1"); OUT.mkdir(parents=True,exist_ok=True)
URL="https://www.gmrt.org/services/GridServer"
UA={"User-Agent":"JANUS-KUSTO-INDIA-G1/1.0","Accept-Encoding":"identity"}
GEOD=Geod(ellps="WGS84")
params={"north":N,"south":S,"west":W,"east":E,"layer":"topo-mask","format":"geotiff","resolution":"max"}

attempts=[]; raw=None
for i in range(30):
    r=requests.get(URL,params=params,headers=UA,timeout=900,allow_redirects=True)
    b=r.content
    attempts.append({"attempt":i+1,"status":r.status_code,"bytes":len(b),"url":r.url,
                     "magic":b[:8].hex(),"content_type":r.headers.get("content-type")})
    if r.status_code==200 and b[:4] in (b"II*"+bytes([0]), b"MM"+bytes([0])+b"*"):
        raw=b; break
    time.sleep(1)
if raw is None:
    raise RuntimeError("No valid GMRT topo-mask TIFF after frozen 30 attempts: "+json.dumps(attempts,indent=2))

tif=OUT/"gmrt_topo_mask_target.tif"; tif.write_bytes(raw)
sha=hashlib.sha256(raw).hexdigest()

with rasterio.open(tif) as ds:
    ma=ds.read(1,masked=True)
    arr=np.asarray(ma.data,dtype=float)
    valid=(~np.ma.getmaskarray(ma)) & np.isfinite(arr)
    row,col=ds.index(LON,LAT)
    x0,y0=ds.xy(row,col)
    x1,y1=ds.xy(row,col+1); x2,y2=ds.xy(row+1,col)
    _,_,dx=GEOD.inv(x0,y0,x1,y1); _,_,dy=GEOD.inv(x0,y0,x2,y2)
    dx,dy=abs(dx),abs(dy)
    target_valid=bool(0<=row<ds.height and 0<=col<ds.width and valid[row,col])
    target=float(arr[row,col]) if target_valid else None
    radii=[]
    for radius in [250,500,1000,2000,5000]:
        rx=max(1,int(math.ceil(radius/dx))); ry=max(1,int(math.ceil(radius/dy)))
        r0=max(0,row-2*ry); r1=min(ds.height,row+2*ry+1)
        c0=max(0,col-2*rx); c1=min(ds.width,col+2*rx+1)
        sub=arr[r0:r1,c0:c1]; vm=valid[r0:r1,c0:c1]
        rr=(np.arange(r0,r1)-row)[:,None]*dy
        cc=(np.arange(c0,c1)-col)[None,:]*dx
        d=np.hypot(rr,cc)
        disk=d<=radius; ann=(d>1.25*radius)&(d<=2*radius)
        dv=sub[disk&vm]; av=sub[ann&vm]
        rec={"radius_m":radius,"disk_valid_n":int(len(dv)),"annulus_valid_n":int(len(av))}
        if len(dv):
            rec.update({"local_min_m":float(np.min(dv)),"local_max_m":float(np.max(dv)),
                        "local_relief_m":float(np.max(dv)-np.min(dv)),
                        "disk_median_m":float(np.median(dv)),
                        "target_local_percentile":float(np.mean(dv<=target)) if target_valid else None})
        if len(av):
            rec["annulus_median_m"]=float(np.median(av))
            if target_valid: rec["target_minus_annulus_median_m"]=float(target-np.median(av))
        radii.append(rec)

    # nearest valid cell if target itself is masked
    nearest=None
    if not target_valid and valid.any():
        rr,cc=np.nonzero(valid)
        dist2=((rr-row)*dy)**2+((cc-col)*dx)**2
        k=int(np.argmin(dist2))
        nr,nc=int(rr[k]),int(cc[k])
        nlon,nlat=ds.xy(nr,nc)
        nearest={"row":nr,"col":nc,"distance_m":float(np.sqrt(dist2[k])),
                 "lon":float(nlon),"lat":float(nlat),"depth_m":float(arr[nr,nc])}

    # quick morphology image, native grid
    v=arr[valid]
    lo=float(np.quantile(v,.02)); hi=float(np.quantile(v,.98))
    fig,ax=plt.subplots(figsize=(7,6))
    im=ax.imshow(np.where(valid,arr,np.nan),origin="upper",vmin=lo,vmax=hi,cmap="gray")
    ax.plot(col,row,"r+",ms=14,mew=2)
    ax.set_title(f"GMRT topo-mask | {LAT:.7f}, {LON:.7f}")
    ax.set_axis_off()
    fig.colorbar(im,ax=ax,label="elevation/depth (m)")
    fig.tight_layout()
    png=OUT/"gmrt_topo_mask_target.png"; fig.savefig(png,dpi=180,bbox_inches="tight"); plt.close(fig)

    out={"artifact_id":"JANUS-KUSTO-INDIA-170710306N-832693611E-G1-HIGHRES-RUN-2026-09-25-v1.0",
         "target":{"lat":LAT,"lon":LON},"source":{"url":URL,"params":params,"attempts":attempts,
         "tiff_sha256":sha,"tiff_bytes":len(raw)},
         "grid":{"shape":[ds.height,ds.width],"crs":str(ds.crs),"bounds":list(ds.bounds),
                 "resolution_deg":[abs(ds.transform.a),abs(ds.transform.e)],
                 "pixel_axis_m":[dx,dy],"valid_fraction":float(valid.mean()),
                 "target_row_col":[row,col],"target_valid":target_valid,
                 "target_depth_m":target,"nearest_valid_if_masked":nearest,
                 "radii":radii},
         "claim_ceiling":"HIGH_RES_BATHYMETRY_MORPHOLOGY_ONLY"}
    (OUT/"JANUS-KUSTO-INDIA-170710306N-832693611E-G1-HIGHRES-RUN-2026-09-25-v1.0.json").write_text(json.dumps(out,indent=2,allow_nan=False))
    print(json.dumps({"target_valid":target_valid,"target_depth_m":target,"nearest_valid_if_masked":nearest,
                      "pixel_axis_m":[dx,dy],"valid_fraction":float(valid.mean()),
                      "radii":radii,"tiff_sha256":sha,"successful_attempt":len(attempts)},indent=2))
