#!/usr/bin/env python3
import json, math, hashlib, tempfile
from pathlib import Path
import numpy as np, requests, rasterio
from pyproj import Geod

LAT=17.0710306; LON=83.2693611
S=16.9710306; N=17.1710306; W=83.1693611; E=83.3693611
OUT=Path("workspace/kusto_global_groundtruth_out"); OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-KUSTO-INDIA-G0/1.0","Accept-Encoding":"identity"}
BASES=["https://www.gmrt.org/services/GridServer","http://www.gmrt.org/services/GridServer"]
GEOD=Geod(ellps="WGS84")

def dl(layer,path):
    attempts=[]
    variants=[
      {"resolution":"max"},
      {"mresolution":"100"}
    ]
    for base in BASES:
      for extra in variants:
        params={"north":N,"south":S,"west":W,"east":E,"layer":layer,"format":"geotiff",**extra}
        try:
          r=requests.get(base,params=params,headers=UA,timeout=900,allow_redirects=True)
          head=r.content[:160]
          rec={"request_url":r.url,"status":r.status_code,"bytes":len(r.content),
               "content_type":r.headers.get("content-type"),"first_bytes_hex":head[:32].hex(),
               "first_text":head.decode("utf-8","replace")}
          attempts.append(rec)
          r.raise_for_status()
          b=r.content
          # TIFF little/big endian magic
          if b[:4] in (b"II*\\x00",b"MM\\x00*"):
            path.write_bytes(b)
            return {"url":r.url,"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest(),
                    "content_type":r.headers.get("content-type"),"attempts":attempts}
          # Some services return ZIP containing a tif
          if b[:2]==b"PK":
            import io,zipfile
            with zipfile.ZipFile(io.BytesIO(b)) as zf:
              tif=[n for n in zf.namelist() if n.lower().endswith((".tif",".tiff"))]
              if len(tif)==1:
                raw=zf.read(tif[0]); path.write_bytes(raw)
                return {"url":r.url,"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest(),
                        "content_type":r.headers.get("content-type"),"archive_member":tif[0],"attempts":attempts}
        except Exception as e:
          attempts.append({"error":repr(e)})
    if layer=="topo-mask" and attempts and all(int(x.get("bytes",-1))==0 for x in attempts if "bytes" in x):\n        return {"empty_response":True,"attempts":attempts}\n    raise RuntimeError("GMRT did not return a GeoTIFF: "+json.dumps(attempts,indent=2))

def radius_stats(ds,arr,valid,radius):
    row,col=ds.index(LON,LAT)
    # metres per pixel at target from geodesic transform
    x0,y0=ds.xy(row,col)
    x1,y1=ds.xy(row,col+1)
    x2,y2=ds.xy(row+1,col)
    _,_,dx=GEOD.inv(x0,y0,x1,y1); _,_,dy=GEOD.inv(x0,y0,x2,y2)
    dx=abs(dx); dy=abs(dy)
    rx=max(1,int(math.ceil(radius/dx))); ry=max(1,int(math.ceil(radius/dy)))
    r0=max(0,row-2*ry); r1=min(ds.height,row+2*ry+1); c0=max(0,col-2*rx); c1=min(ds.width,col+2*rx+1)
    sub=arr[r0:r1,c0:c1]; vm=valid[r0:r1,c0:c1]
    rr=(np.arange(r0,r1)-row)[:,None]*dy; cc=(np.arange(c0,c1)-col)[None,:]*dx
    d=np.hypot(rr,cc)
    disk=d<=radius; ann=(d>1.25*radius)&(d<=2*radius)
    dv=sub[disk&vm]; av=sub[ann&vm]
    if len(dv)<1:
        return {"radius_m":radius,"valid":False}
    target=float(arr[row,col]) if valid[row,col] else None
    relief=float(np.max(dv)-np.min(dv))
    med=float(np.median(dv))
    amed=float(np.median(av)) if len(av) else None
    contrast=(target-amed) if target is not None and amed is not None else None
    # percentile of target in local disk: 0 low, 1 high
    rank=float(np.mean(dv<=target)) if target is not None else None
    return {"radius_m":radius,"valid":True,"disk_valid_n":int(len(dv)),"annulus_valid_n":int(len(av)),
            "local_relief_m":relief,"disk_median_m":med,"annulus_median_m":amed,
            "target_minus_annulus_median_m":contrast,"target_local_percentile":rank,
            "pixel_axis_m":[dx,dy]}

def inspect(path,label):
    with rasterio.open(path) as ds:
        ma=ds.read(1,masked=True)
        arr=np.asarray(ma.data,dtype=float)
        valid=(~np.ma.getmaskarray(ma)) & np.isfinite(arr)
        row,col=ds.index(LON,LAT)
        inb=0<=row<ds.height and 0<=col<ds.width
        target_valid=bool(inb and valid[row,col])
        result={"label":label,"shape":[ds.height,ds.width],"crs":str(ds.crs),"bounds":list(ds.bounds),
                "resolution_deg":[abs(ds.transform.a),abs(ds.transform.e)],"valid_fraction":float(valid.mean()),
                "target_row_col":[row,col],"target_valid":target_valid,
                "target_depth_m":float(arr[row,col]) if target_valid else None,
                "radii":[radius_stats(ds,arr,valid,r) for r in [250,500,1000,2000,5000]]}
        # central 5x5 values for sanity
        if inb:
            a=max(0,row-2); b=min(ds.height,row+3); c=max(0,col-2); d=min(ds.width,col+3)
            vals=np.where(valid[a:b,c:d],arr[a:b,c:d],np.nan)
            result["target_5x5_m"]=vals.tolist()
        return result

with tempfile.TemporaryDirectory() as td:
    td=Path(td)
    outputs={}
    for layer in ["topo-mask","topo"]:
        p=td/(layer.replace("-","_")+".tif")
        tr=dl(layer,p)\n        outputs[layer]={"transport":tr,"grid":None if tr.get("empty_response") else inspect(p,layer)}
    masked=outputs["topo-mask"]["grid"]
    unmasked=outputs["topo"]["grid"]
    if masked["target_valid"]:
        coverage="HIGH_RES_GMRT_AT_TARGET"
        ceiling="HIGH_RES_BATHYMETRY_MORPHOLOGY_ONLY"
    elif masked["valid_fraction"]>0:
        coverage="HIGH_RES_GMRT_IN_BOX_BUT_NOT_AT_TARGET"
        ceiling="LOW_RES_TARGET_CONTEXT_PLUS_NEARBY_HIGH_RES_COVERAGE"
    else:
        coverage="NO_HIGH_RES_GMRT_IN_BOX"
        ceiling="LOW_RES_GLOBAL_CONTEXT_ONLY"
    out={"artifact_id":"JANUS-KUSTO-INDIA-170710306N-832693611E-G0-RUN-2026-09-25-v1.0",
         "target":{"lat":LAT,"lon":LON},"box":{"south":S,"north":N,"west":W,"east":E},
         "gmrt_version_context":"v4.5 service queried 2026-09-25",
         "outputs":outputs,"coverage_verdict":coverage,"claim_ceiling":ceiling}
    raw=json.dumps(out,indent=2,allow_nan=False)
    (OUT/"JANUS-KUSTO-INDIA-170710306N-832693611E-G0-RUN-2026-09-25-v1.0.json").write_text(raw)
    print(json.dumps({
      "coverage_verdict":coverage,
      "masked_target_valid":masked["target_valid"],
      "masked_valid_fraction":masked["valid_fraction"],
      "unmasked_target_depth_m":unmasked["target_depth_m"],
      "unmasked_radii":unmasked["radii"],
      "output_sha256":hashlib.sha256(raw.encode()).hexdigest()
    },indent=2))
