#!/usr/bin/env python3
from __future__ import annotations
import io, json, zipfile, hashlib, math
from pathlib import Path
import requests
import rasterio
from rasterio.io import MemoryFile
from pyproj import Transformer

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
F=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-A0E-EXACT-MEMBER-SELECTION-FREEZE-2026-09-24-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-Arafura-A0G-domain-bridge/1.0"}

def fetch_member(spec, member_key="member"):
    r=requests.get(spec["archive_url"],headers=UA,timeout=300,allow_redirects=True); r.raise_for_status()
    blob=r.content
    sha=hashlib.sha256(blob).hexdigest()
    if sha.lower()!=spec["archive_sha256"].lower(): raise RuntimeError("archive sha mismatch")
    if len(blob)!=int(spec["archive_bytes"]): raise RuntimeError("archive size mismatch")
    member=spec[member_key]
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        zi=zf.getinfo(member)
        raw=zf.read(member)
    return raw,zi

braw,_=fetch_member(F["bathymetry"])
hraw,_=fetch_member(F["heldout_backscatter"])

with MemoryFile(braw) as mf:
  with mf.open() as ds:
    bh={"crs":str(ds.crs),"transform":list(ds.transform)[:6],"bounds":[ds.bounds.left,ds.bounds.bottom,ds.bounds.right,ds.bounds.top],"width":ds.width,"height":ds.height,"nodata":ds.nodata}
with MemoryFile(hraw) as mf:
  with mf.open() as ds:
    hh={"crs":str(ds.crs),"transform":list(ds.transform)[:6],"bounds":[ds.bounds.left,ds.bounds.bottom,ds.bounds.right,ds.bounds.top],"width":ds.width,"height":ds.height,"nodata":ds.nodata}

if bh["crs"]!="EPSG:4326": raise RuntimeError(f"unexpected bathy CRS {bh['crs']}")
if hh["crs"]!="EPSG:32753": raise RuntimeError(f"unexpected heldout CRS {hh['crs']}")

left,bottom,right,top=hh["bounds"]
# Densify each edge so the transformed footprint is frozen without reading raster values.
n=64
ring=[]
for i in range(n+1):
    t=i/n; ring.append((left+t*(right-left), top))
for i in range(1,n+1):
    t=i/n; ring.append((right, top-t*(top-bottom)))
for i in range(1,n+1):
    t=i/n; ring.append((right-t*(right-left), bottom))
for i in range(1,n):
    t=i/n; ring.append((left, bottom+t*(top-bottom)))
ring.append(ring[0])

to_wgs=Transformer.from_crs("EPSG:32753","EPSG:4326",always_xy=True)
lon,lat=to_wgs.transform([p[0] for p in ring],[p[1] for p in ring])
poly=[[float(x),float(y)] for x,y in zip(lon,lat)]
minlon=min(x for x,y in poly); maxlon=max(x for x,y in poly)
minlat=min(y for x,y in poly); maxlat=max(y for x,y in poly)

# Compute a conservative source-grid window from header geometry only.
a,b,c,d,e,f=bh["transform"]
if abs(b)>1e-15 or abs(d)>1e-15: raise RuntimeError("rotated bathy transform not supported")
def rc(x,y):
    col=(x-c)/a
    row=(y-f)/e
    return row,col
r0,c0=rc(minlon,maxlat); r1,c1=rc(maxlon,minlat)
row0=max(0,int(math.floor(min(r0,r1)))-2); row1=min(bh["height"],int(math.ceil(max(r0,r1)))+2)
col0=max(0,int(math.floor(min(c0,c1)))-2); col1=min(bh["width"],int(math.ceil(max(c0,c1)))+2)

out={
 "artifact_id":"JANUS-KUSTO-ARAFURA-A0G-MONEY-SHOAL-REPRESENTATION-BRIDGE-RUN-2026-09-24-v1.0",
 "source_freeze":F["artifact_id"],
 "selected_subarea":"Money Shoal",
 "bathymetry_header":bh,
 "heldout_backscatter_header":hh,
 "money_shoal_header_footprint_wgs84_polygon":poly,
 "money_shoal_header_footprint_wgs84_bbox":[minlon,minlat,maxlon,maxlat],
 "conservative_bathymetry_source_window":{"row_start":row0,"row_stop":row1,"col_start":col0,"col_stop":col1,"height":row1-row0,"width":col1-col0},
 "representation_rule":{
   "morphology_grid":"ORIGINAL_BATHYMETRY_SOURCE_GRID_NO_RESAMPLING",
   "nominal_cell_m":6,
   "candidate_domain":"BATHYMETRY_CELL_CENTER_MUST_LIE_INSIDE_HEADER_ONLY_MONEY_SHOAL_FOOTPRINT_POLYGON",
   "metric_coordinate_bridge":"EPSG:4326_CELL_CENTER_TO_EPSG:32753_FOR_DEDUP_AND_LATER_PAIRED_SCORING_ONLY",
   "heldout_backscatter_values_used":False
 },
 "pixel_values_read":False,
 "groundtruth_assets_read":False,
 "next_gate":"A1_MONEY_SHOAL_BATHYMETRY_ONLY_BLIND_CANDIDATE_FREEZE",
 "claim_ceiling":"HEADER_ONLY_SUBAREA_AND_COORDINATE_REPRESENTATION_BRIDGE"
}
p=OUT/"JANUS-KUSTO-ARAFURA-A0G-MONEY-SHOAL-REPRESENTATION-BRIDGE-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2); p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"bbox":out["money_shoal_header_footprint_wgs84_bbox"],
 "window":out["conservative_bathymetry_source_window"],"pixel_values_read":False,
 "heldout_backscatter_values_used":False
},indent=2))
