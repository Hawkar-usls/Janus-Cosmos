#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json, os, shutil, subprocess, tempfile
from pathlib import Path
import requests, rasterio
from rasterio.warp import transform_bounds

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
BATCH=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1B-PREOUTCOME-BATCH-METADATA-AND-TRANSPORT-RECEIPT-2026-09-24-v1.0.json").read_text())
INV=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1C-BATCH-ARCHIVE-MEMBER-INVENTORY-RECEIPT-2026-09-24-v1.0.json").read_text())
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1D-BATCH-PAIRED-RASTER-HEADER-GEOMETRY-PREREG-2026-09-24-v1.0.json").read_text())
S=requests.Session(); S.headers.update({"User-Agent":"JANUS-KUSTO-INFOMAR-I1D/1.0","Accept-Encoding":"identity"})

rows={x["SURVEY_ID"]:x for x in BATCH["prospective"]}
inv={x["SURVEY_ID"]:x for x in INV["batch"]}

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()

def download(url,path,expected_size):
    h=hashlib.sha256(); n=0
    with S.get(url,stream=True,timeout=600,allow_redirects=True) as r:
        r.raise_for_status()
        with open(path,"wb") as f:
            for b in r.iter_content(4*1024*1024):
                if not b: continue
                f.write(b); h.update(b); n+=len(b)
    if n!=int(expected_size): raise RuntimeError(f"archive size drift {n} vs {expected_size}")
    return h.hexdigest(),n

def extract_7z(zip_path,member,outdir):
    cmd=["7z","x","-y",f"-o{outdir}",str(zip_path),member]
    cp=subprocess.run(cmd,text=True,capture_output=True)
    if cp.returncode!=0: raise RuntimeError(f"7z failed rc={cp.returncode}: {cp.stdout}\n{cp.stderr}")
    p=Path(outdir)/member
    if not p.exists(): raise RuntimeError(f"7z did not produce {member}")
    return p

def header(p):
    with rasterio.open(p) as ds:
        tr=ds.transform; b=ds.bounds
        return {
          "width":int(ds.width),"height":int(ds.height),"count":int(ds.count),
          "dtypes":list(ds.dtypes),"nodatavals":[None if x is None else float(x) for x in ds.nodatavals],
          "crs":str(ds.crs) if ds.crs else None,"epsg":ds.crs.to_epsg() if ds.crs else None,
          "transform":[float(tr.a),float(tr.b),float(tr.c),float(tr.d),float(tr.e),float(tr.f)],
          "bounds":[float(b.left),float(b.bottom),float(b.right),float(b.top)],
          "resolution":[float(abs(tr.a)),float(abs(tr.e))],
          "block_shapes":[list(map(int,x)) for x in ds.block_shapes]
        }

results=[]
with tempfile.TemporaryDirectory(prefix="kusto_infomar_i1d_") as td:
    td=Path(td)
    for sid in ["CB13_03","CB13_05"]:
        row=rows[sid]; ii=inv[sid]
        rec={"SURVEY_ID":sid,"rank":row["rank"],"state":"PENDING","channels":{}}
        for channel,urlkey,sizekey in [
          ("bathymetry","bathymetry_url","bathymetry_archive_size_bytes"),
          ("backscatter","backscatter_url","backscatter_archive_size_bytes")
        ]:
            spec=ii[channel]; z=td/f"{sid}_{channel}.zip"; ex=td/f"{sid}_{channel}_ex"; ex.mkdir()
            asha,abytes=download(row[urlkey],z,row[sizekey])
            p=extract_7z(z,spec["filename"],ex)
            st=p.stat().st_size
            if st!=int(spec["uncompressed_size"]): raise RuntimeError(f"{sid} {channel} extracted size drift {st}")
            msha=sha256_file(p); h=header(p)
            if not h["crs"]: raise RuntimeError(f"{sid} {channel} missing CRS")
            wb=transform_bounds(h["crs"],"EPSG:4326",*h["bounds"],densify_pts=21)
            rec["channels"][channel]={
              "archive_url":row[urlkey],"archive_sha256":asha,"archive_bytes":abytes,
              "member":spec["filename"],"member_crc32":spec["crc32"],"member_sha256":msha,"member_bytes":st,
              "header":h,"wgs84_bounds":[float(x) for x in wb],
              "raster_pixel_values_read":False
            }
            shutil.rmtree(ex); z.unlink()
        a=rec["channels"]["bathymetry"]["wgs84_bounds"]; b=rec["channels"]["backscatter"]["wgs84_bounds"]
        overlap=[max(a[0],b[0]),max(a[1],b[1]),min(a[2],b[2]),min(a[3],b[3])]
        ok=overlap[2]>overlap[0] and overlap[3]>overlap[1]
        rec["wgs84_overlap_bounds"]=overlap if ok else None
        rec["nonempty_overlap"]=ok
        rec["state"]="PASS_HEADER_GEOMETRY" if ok else "GEOMETRY_BLOCKED__NOT_REPLACED"
        results.append(rec)

out={
 "artifact_id":"JANUS-KUSTO-INFOMAR-I1D-BATCH-PAIRED-RASTER-HEADER-GEOMETRY-RECEIPT-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],"results":results,
 "transport_blocks_preserved":["CB14_01","CB14_02"],
 "header_geometry_pass_count":sum(r["state"]=="PASS_HEADER_GEOMETRY" for r in results),
 "raster_pixel_values_read":False,"groundtruth_read":False,
 "next_gate":"I1E_PER_SURVEY_NATIVE_CELL_BLIND_FREEZE",
 "claim_ceiling":PRE["claim_ceiling"]
}
p=OUT/"JANUS-KUSTO-INFOMAR-I1D-BATCH-PAIRED-RASTER-HEADER-GEOMETRY-RECEIPT-2026-09-24-v1.0.json"
p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "results":[{"SURVEY_ID":r["SURVEY_ID"],"state":r["state"],
             "bathymetry":r["channels"]["bathymetry"]["header"],
             "backscatter":r["channels"]["backscatter"]["header"],
             "wgs84_overlap_bounds":r["wgs84_overlap_bounds"]} for r in results],
 "raster_pixel_values_read":False
},indent=2,ensure_ascii=False))
