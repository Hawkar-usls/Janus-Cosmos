#!/usr/bin/env python3
from __future__ import annotations
import hashlib,io,json,tempfile,zipfile
from pathlib import Path
import numpy as np
import requests,rasterio
from rasterio.enums import Resampling
from scipy import ndimage

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
SRC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EMR0C-EXACT-SOURCE-SELECTION-FREEZE-2026-09-24-v1.0.json").read_text())
REC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EMR0D-RASTER-HEADER-RECEIPT-2026-09-24-v1.0.json").read_text())
D=SRC["discovery_source"]
UA={"User-Agent":"JANUS-KUSTO-ElizabethMiddleton-EMR0E-support-mask/1.0"}

r=requests.get(D["url"],headers=UA,timeout=360,allow_redirects=True);r.raise_for_status()
blob=r.content
if hashlib.sha256(blob).hexdigest()!=D["archive_sha256"]:raise RuntimeError("archive SHA mismatch")

with tempfile.TemporaryDirectory() as td:
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        zi=zf.getinfo(D["member"])
        if zi.file_size!=D["member_uncompressed_bytes"] or zi.CRC!=D["member_crc32"]:raise RuntimeError("member mismatch")
        outp=Path(td)/"bathy.tif"
        with zf.open(D["member"]) as src, open(outp,"wb") as dst:
            while True:
                b=src.read(8*1024*1024)
                if not b:break
                dst.write(b)
    with rasterio.open(outp) as ds:
        overviews=ds.overviews(1)
        target_h=2048
        target_w=max(1,int(round(ds.width*(target_h/ds.height))))
        mask=ds.read_masks(1,out_shape=(target_h,target_w),resampling=Resampling.nearest)
        valid=mask>0
        labels,n=ndimage.label(valid,structure=np.ones((3,3),dtype=np.uint8))
        objs=ndimage.find_objects(labels)
        comps=[]
        sy=ds.height/target_h;sx=ds.width/target_w
        for lab,sl in enumerate(objs,start=1):
            if sl is None:continue
            ys,xs=sl
            count=int(np.sum(labels[sl]==lab))
            # map coarse half-open box to full-res pixel window
            r0=max(0,int(np.floor(ys.start*sy)));r1=min(ds.height,int(np.ceil(ys.stop*sy)))
            c0=max(0,int(np.floor(xs.start*sx)));c1=min(ds.width,int(np.ceil(xs.stop*sx)))
            # pixel corners -> geographic bounds
            x0,y0=ds.transform*(c0,r0);x1,y1=ds.transform*(c1,r1)
            comps.append({
              "component_id":lab,
              "coarse_valid_pixel_count":count,
              "coarse_bbox":[int(xs.start),int(ys.start),int(xs.stop),int(ys.stop)],
              "source_pixel_window":[r0,r1,c0,c1],
              "source_crs_bounds":[float(min(x0,x1)),float(min(y0,y1)),float(max(x0,x1)),float(max(y0,y1))]
            })
        comps.sort(key=lambda q:q["coarse_valid_pixel_count"],reverse=True)
        summary={
          "source_shape":[ds.height,ds.width],
          "source_crs":str(ds.crs),
          "source_bounds":[float(ds.bounds.left),float(ds.bounds.bottom),float(ds.bounds.right),float(ds.bounds.top)],
          "block_shapes":[list(x) for x in ds.block_shapes],
          "overviews":overviews,
          "coarse_mask_shape":[target_h,target_w],
          "coarse_valid_pixels":int(np.sum(valid)),
          "coarse_valid_fraction":float(np.mean(valid)),
          "component_count":len(comps)
        }

out={
 "artifact_id":"JANUS-KUSTO-ELIZABETH-MIDDLETON-EMR0E-DISCOVERY-SUPPORT-MASK-AUDIT-2026-09-24-v1.0",
 "source_freeze":SRC["artifact_id"],
 "header_receipt":REC["artifact_id"],
 "support_mask_summary":summary,
 "components":comps,
 "numeric_depth_values_read":False,
 "heldout_backscatter_geometry_used_for_component_detection":False,
 "backscatter_values_read":False,
 "groundtruth_assets_read":False,
 "claim_ceiling":"DISCOVERY_CHANNEL_SUPPORT_MASK_GEOMETRY_ONLY"
}
p=OUT/"JANUS-KUSTO-ELIZABETH-MIDDLETON-EMR0E-DISCOVERY-SUPPORT-MASK-AUDIT-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"summary":summary,
 "top_components":comps[:20],
 "numeric_depth_values_read":False,
 "heldout_backscatter_geometry_used":False
},indent=2))
