#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json, shutil, tempfile, zipfile
from pathlib import Path
import requests, rasterio
from rasterio.warp import transform_bounds

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PR=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I0J-PAIRED-GRID-GEOMETRY-METADATA-PREREG-2026-09-24-v1.0.json").read_text())
FR=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I0I-PAIRED-RASTER-MEMBER-SELECTION-FREEZE-2026-09-24-v1.0.json").read_text())
S=requests.Session(); S.headers.update({"User-Agent":"JANUS-KUSTO-INFOMAR-I0J/1.0","Accept-Encoding":"identity"})

def download(url,path):
    h=hashlib.sha256(); n=0
    with S.get(url,stream=True,timeout=300,allow_redirects=True) as r:
        r.raise_for_status()
        with open(path,"wb") as f:
            for chunk in r.iter_content(1024*1024):
                if not chunk: continue
                f.write(chunk); h.update(chunk); n+=len(chunk)
    return h.hexdigest(),n

def inspect(label,spec,tmp):
    zpath=tmp/f"{label}.zip"
    archive_sha,archive_bytes=download(spec["archive_url"],zpath)
    member=spec["member"]
    outpath=tmp/f"{label}.tif"
    member_sha=hashlib.sha256(); member_bytes=0
    with zipfile.ZipFile(zpath) as zf:
        zi=zf.getinfo(member)
        if f"{zi.CRC:08x}".lower()!=str(spec["crc32"]).lower():
            raise RuntimeError(f"{label}: CRC drift")
        if int(zi.compress_size)!=int(spec["compressed_size"]):
            raise RuntimeError(f"{label}: compressed size drift")
        if int(zi.file_size)!=int(spec["uncompressed_size"]):
            raise RuntimeError(f"{label}: uncompressed size drift")
        with zf.open(zi,"r") as src, open(outpath,"wb") as dst:
            while True:
                b=src.read(1024*1024)
                if not b: break
                dst.write(b); member_sha.update(b); member_bytes+=len(b)
    if member_bytes!=int(spec["uncompressed_size"]):
        raise RuntimeError(f"{label}: extracted size mismatch")
    with rasterio.open(outpath) as ds:
        # HEADER/METADATA ONLY. Deliberately no ds.read(), ds.sample(), statistics or histogram calls.
        crs=str(ds.crs) if ds.crs else None
        epsg=ds.crs.to_epsg() if ds.crs else None
        tr=ds.transform
        b=ds.bounds
        meta={
          "width":int(ds.width),"height":int(ds.height),"count":int(ds.count),
          "dtypes":list(ds.dtypes),"nodatavals":[None if x is None else float(x) for x in ds.nodatavals],
          "crs":crs,"epsg":epsg,
          "transform":[float(tr.a),float(tr.b),float(tr.c),float(tr.d),float(tr.e),float(tr.f)],
          "bounds":[float(b.left),float(b.bottom),float(b.right),float(b.top)],
          "resolution":[float(abs(tr.a)),float(abs(tr.e))],
          "block_shapes":[list(map(int,x)) for x in ds.block_shapes]
        }
    if not meta["crs"]:
        raise RuntimeError(f"{label}: missing CRS")
    wb=transform_bounds(meta["crs"],"EPSG:4326",*meta["bounds"],densify_pts=21)
    return {
      "archive_url":spec["archive_url"],"archive_sha256":archive_sha,"archive_bytes":archive_bytes,
      "member":member,"member_sha256":member_sha.hexdigest(),"member_bytes":member_bytes,
      "member_crc32_verified":str(spec["crc32"]).lower(),
      "header":meta,"wgs84_bounds":[float(x) for x in wb],
      "raster_pixel_values_read":False
    }

with tempfile.TemporaryDirectory(prefix="kusto_infomar_i0j_") as td:
    tmp=Path(td)
    bathy=inspect("bathymetry",FR["selected"]["bathymetry"],tmp)
    back=inspect("heldout_backscatter",FR["selected"]["heldout_backscatter"],tmp)

l=max(bathy["wgs84_bounds"][0],back["wgs84_bounds"][0])
b=max(bathy["wgs84_bounds"][1],back["wgs84_bounds"][1])
r=min(bathy["wgs84_bounds"][2],back["wgs84_bounds"][2])
t=min(bathy["wgs84_bounds"][3],back["wgs84_bounds"][3])
overlap=bool(r>l and t>b)
out={
 "artifact_id":"JANUS-KUSTO-INFOMAR-I0J-PAIRED-GRID-GEOMETRY-METADATA-RECEIPT-2026-09-24-v1.0",
 "prereg":PR["artifact_id"],"member_freeze":FR["artifact_id"],"survey_id":FR["survey_id"],
 "bathymetry":bathy,"heldout_backscatter":back,
 "wgs84_overlap_bounds":[l,b,r,t] if overlap else None,
 "nonempty_overlap":overlap,
 "archive_bodies_consumed":True,
 "raster_members_materialized_for_header_open":True,
 "raster_pixel_values_read":False,
 "result":"PASS" if overlap else "FAIL",
 "next_gate":PR["next_gate_on_pass"] if overlap else PR["next_gate_on_fail"],
 "claim_ceiling":PR["claim_ceiling"]
}
p=OUT/"JANUS-KUSTO-INFOMAR-I0J-PAIRED-GRID-GEOMETRY-METADATA-RECEIPT-2026-09-24-v1.0.json"
p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps(out,indent=2,ensure_ascii=False))
