#!/usr/bin/env python3
import ftplib, hashlib, json, math, statistics
from pathlib import Path

HOST="livftp.noc.ac.uk"
ROOT="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/EM12"
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)

C={"lat":-4.015075679897318,"lon":-12.29915403590432}
M=111320.0; LAT0=-4.0; COS0=math.cos(math.radians(LAT0)); CELL=100.0
PRODUCTS={
 "LINE31_PRODUCT":{"name":"B1-81-1_Acceptl28-33.xyz.ascii","sha256":"0e8cac049b53dd6ba0e5013e8aaea2bfb48347db887edc3c4e3e9f697e96dcf0"},
 "LINE18_PRODUCT":{"name":"B1-81-1_Acceptl15-22.xyz.ascii","sha256":"ad2e4d60401f7ffaa5e14b3e43aa0e45da9d320c6b34fa7f3b873f4024be6007"}
}
CELLS=[
 {"id":"CELL_A","lat":-4.0165853,"lon":-12.2985657,"z":4025.08},
 {"id":"CELL_B","lat":-4.0152462,"lon":-12.299318,"z":4036.31}
]

def ftp():
    f=ftplib.FTP(timeout=120); f.connect(HOST,21); f.login("anonymous","janus-kusto@example.invalid"); f.voidcmd("TYPE I"); return f

def gxy(lon,lat): return lon*M*COS0,lat*M
CX,CY=gxy(C["lon"],C["lat"])

def parse_xyz(raw):
    s=raw.decode("ascii","ignore").strip()
    if not s or s[0] in "#;!": return None
    p=s.replace(","," ").split()
    if len(p)<3:return None
    try:lon=float(p[0]);lat=float(p[1]);z=float(p[2])
    except:return None
    if not all(math.isfinite(v) for v in (lon,lat,z)):return None
    x,y=gxy(lon,lat)
    dx=x-CX;dy=y-CY
    ix=math.floor((dx+1000.0)/CELL);iy=math.floor((dy+1000.0)/CELL)
    return {"lon":lon,"lat":lat,"z":z,"x":x,"y":y,"dx":dx,"dy":dy,"ix":ix,"iy":iy}

def load(name,expected):
    f=ftp();h=hashlib.sha256();rows=[]
    try:
        sock=f.transfercmd("RETR "+ROOT+"/"+name);stream=sock.makefile("rb")
        try:
            for raw in stream:
                h.update(raw)
                q=parse_xyz(raw)
                if q is not None: rows.append(q)
        finally: stream.close();sock.close()
    finally:
        try:f.quit()
        except:f.close()
    sha=h.hexdigest()
    if sha!=expected:raise RuntimeError(f"{name} SHA mismatch {sha}")
    return rows,sha

def median_abs_dev(vals):
    if not vals:return None
    med=statistics.median(vals)
    return statistics.median([abs(v-med) for v in vals])

def pct(vals,p):
    s=sorted(vals)
    if not s:return None
    q=(len(s)-1)*p;i=int(math.floor(q));f=q-i
    return s[i] if i==len(s)-1 else s[i]*(1-f)+s[i+1]*f

def cell_center(ix,iy):
    return -1000+(ix+.5)*CELL,-1000+(iy+.5)*CELL

def stats(rows,ix,iy):
    rr=[r for r in rows if r["ix"]==ix and r["iy"]==iy]
    if not rr:return {"n":0}
    z=[r["z"] for r in rr]
    xc,yc=cell_center(ix,iy)
    dcenter=[math.hypot(r["dx"]-xc,r["dy"]-yc) for r in rr]
    xs=[r["dx"] for r in rr];ys=[r["dy"] for r in rr]
    return {
      "n":len(rr),
      "depth_m":{
        "median":statistics.median(z),"mean":statistics.fmean(z),"min":min(z),"max":max(z),
        "p10":pct(z,.10),"p90":pct(z,.90),"iqr":pct(z,.75)-pct(z,.25),"mad":median_abs_dev(z),
        "range":max(z)-min(z)
      },
      "geometry":{
        "cell_center_offset_m":[xc,yc],
        "min_distance_to_cell_center_m":min(dcenter),
        "max_distance_to_cell_center_m":max(dcenter),
        "east_span_m":max(xs)-min(xs),"north_span_m":max(ys)-min(ys)
      },
      "soundings":[{"lat":r["lat"],"lon":r["lon"],"depth_m":r["z"],"distance_to_cell_center_m":d}
                   for r,d in sorted(zip(rr,dcenter),key=lambda q:q[1])]
    }

loaded={}
manifest={}
for label,p in PRODUCTS.items():
    rows,sha=load(p["name"],p["sha256"]);loaded[label]=rows;manifest[label]={"file":p["name"],"sha256":sha,"parsed_xyz":len(rows)}

results={}
for cell in CELLS:
    x,y=gxy(cell["lon"],cell["lat"]); dx=x-CX;dy=y-CY
    ix=math.floor((dx+1000)/CELL);iy=math.floor((dy+1000)/CELL)
    xc,yc=cell_center(ix,iy)
    primary=stats(loaded["LINE31_PRODUCT"],ix,iy)
    if primary["n"]==0: raise RuntimeError(f"{cell['id']} exact line31 cell unexpectedly empty")
    rep=min(primary["soundings"],key=lambda s:s["distance_to_cell_center_m"])
    rep_minus_med=cell["z"]-primary["depth_m"]["median"]
    if abs(rep_minus_med)<=10:label="REPRESENTATIVE_CONSISTENT"
    elif abs(rep_minus_med)>30:label="REPRESENTATIVE_OUTLIER_LIKE"
    else:label="INTERMEDIATE"
    adjacent={}
    for oy in [-1,0,1]:
        for ox in [-1,0,1]:
            if ox==0 and oy==0:continue
            adjacent[f"{ox:+d},{oy:+d}"]=stats(loaded["LINE31_PRODUCT"],ix+ox,iy+oy)
    line18=stats(loaded["LINE18_PRODUCT"],ix,iy)
    results[cell["id"]]={
      "frozen_representative":{"lat":cell["lat"],"lon":cell["lon"],"depth_m":cell["z"],"offset_m":[dx,dy],
                               "distance_to_cell_center_m":math.hypot(dx-xc,dy-yc)},
      "grid":{"ix":ix,"iy":iy,"center_offset_m":[xc,yc]},
      "line31_exact_cell":primary,
      "recovered_nearest_to_center_sounding":rep,
      "representative_minus_cell_median_m":rep_minus_med,
      "representative_consistency_label":label,
      "line31_adjacent_8_cells":adjacent,
      "line18_same_exact_cell":line18,
      "same_cell_cross_product_support":line18.get("n",0)>0
    }

out={
 "artifact_id":"JANUS-KUSTO-CAND003-LINE31-DOMINANT-CELLS-INTERNAL-CONSISTENCY-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CAND003-LINE31-DOMINANT-CELLS-INTERNAL-CONSISTENCY-PREREG-2026-09-23-v1.0.json",
 "candidate":C,
 "manifest":manifest,
 "results":results,
 "strict_parent_fail":"IMMUTABLE",
 "cells_removed":False,
 "claim_ceiling":"PROCESSED_CELL_INTERNAL_CONSISTENCY_DIAGNOSTIC_ONLY"
}
p=OUT/"JANUS-KUSTO-CAND003-LINE31-DOMINANT-CELLS-INTERNAL-CONSISTENCY-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "manifest":manifest,
 "results":{k:{
   "frozen_representative":v["frozen_representative"],
   "grid":v["grid"],
   "line31_exact_cell":{kk:vv for kk,vv in v["line31_exact_cell"].items() if kk!="soundings"},
   "recovered_nearest_to_center_sounding":v["recovered_nearest_to_center_sounding"],
   "representative_minus_cell_median_m":v["representative_minus_cell_median_m"],
   "representative_consistency_label":v["representative_consistency_label"],
   "line18_same_exact_cell":{kk:vv for kk,vv in v["line18_same_exact_cell"].items() if kk!="soundings"},
   "same_cell_cross_product_support":v["same_cell_cross_product_support"],
   "adjacent_summary":{a:{"n":s.get("n"),"median":s.get("depth_m",{}).get("median")} for a,s in v["line31_adjacent_8_cells"].items()}
 } for k,v in results.items()},
 "strict_parent_fail":"IMMUTABLE"
},indent=2))
