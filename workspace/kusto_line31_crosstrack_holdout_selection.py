#!/usr/bin/env python3
from __future__ import annotations
import concurrent.futures, ftplib, hashlib, json, math, re, time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import requests
import shapely
from shapely.geometry import LineString, Point
from shapely.strtree import STRtree
from scipy.spatial import cKDTree

HOST="livftp.noc.ac.uk"
ROOT="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281"
NAV_PATH=ROOT+"/Nav/cd169leg1_1min.listit"
PRODUCT_PATH=ROOT+"/EM12/B1-81-1_Acceptl28-33.xyz.ascii"
PRODUCT_SHA="0e8cac049b53dd6ba0e5013e8aaea2bfb48347db887edc3c4e3e9f697e96dcf0"
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
M=111320.0; LAT0=-4.0; COS0=math.cos(math.radians(LAT0))
UA={"User-Agent":"JANUS-KUSTO-line31-holdout-selection/1.0"}

SURVEYS={
 "KN192-07":"https://data.ngdc.noaa.gov/platforms/ocean/ships/knorr/KN192-07/multibeam/data/version1/MB/generated/",
 "KNOX15RR":"https://data.ngdc.noaa.gov/platforms/ocean/ships/roger_revelle/KNOX15RR/multibeam/data/version1/MB/generated/"
}
WINDOWS={
 28:("2005-02-25T07:30:00Z","2005-02-25T21:10:00Z"),
 29:("2005-02-25T21:10:00Z","2005-02-26T07:34:00Z"),
 30:("2005-02-26T08:22:00Z","2005-02-26T09:35:00Z"),
 31:("2005-02-26T19:11:00Z","2005-02-27T14:40:00Z"),
 32:("2005-02-27T14:40:00Z","2005-02-28T07:31:00Z"),
 33:("2005-02-28T07:31:00Z","2005-02-28T14:59:00Z")
}
EXCLUDE=[
 ("CAND003_FAIL",-4.015075679897318,-12.29915403590432,2000.0),
 ("MAR4DEG02S_PASS",-4.03,-12.25,2000.0)
]
BANDS={
 "CENTER":(0.20,0.40,0.30),
 "EDGE":(0.95,0.99,0.975)
}

def xy(lon,lat): return lon*M*COS0,lat*M
def ll(x,y): return y/M,x/(M*COS0)
def parse_iso(s): return datetime.fromisoformat(s.replace("Z","+00:00")).timestamp()

def ftp_fetch(path):
    f=ftplib.FTP(timeout=120); f.connect(HOST,21); f.login("anonymous","janus-kusto@example.invalid"); f.voidcmd("TYPE I"); chunks=[]
    try:f.retrbinary("RETR "+path,chunks.append)
    finally:
        try:f.quit()
        except:f.close()
    return b"".join(chunks)

def fetch(url,timeout=120):
    last=None
    for attempt in range(7):
        try:
            r=requests.get(url,headers=UA,timeout=timeout)
            if r.status_code==429:
                last=RuntimeError(f"429 {url}"); time.sleep(min(2**attempt,32)); continue
            r.raise_for_status(); return r.content
        except requests.exceptions.RequestException as e:
            last=e
            if attempt==6: raise
            time.sleep(min(2**attempt,32))
    raise last

def list_fnv(base):
    html=fetch(base,60).decode("latin1","replace")
    return sorted(set(re.findall(r'href="([^"]+\.fnv)"',html,re.I)))

def parse_fnv_segments(data):
    segs=[]
    for line in data.decode("utf-8","replace").splitlines():
        p=line.split()
        if len(p)<19: continue
        try:
            plon,plat,slon,slat=map(float,(p[15],p[16],p[17],p[18]))
        except: continue
        ax,ay=xy(plon,plat); bx,by=xy(slon,slat)
        if all(map(math.isfinite,(ax,ay,bx,by))):
            segs.append(LineString([(ax,ay),(bx,by)]))
    return segs

# CD169 navigation geometry.
navb=ftp_fetch(NAV_PATH)
base=datetime(2005,1,1,tzinfo=timezone.utc); nav=[]
for raw in navb.decode("ascii","ignore").splitlines():
    p=raw.split()
    if len(p)<7:continue
    try:
        jd=int(p[1]); hh,mm,ss=map(int,p[2].split(":")); lat=float(p[3]); lon=float(p[5])
    except:continue
    dt=base+timedelta(days=jd-1,hours=hh,minutes=mm,seconds=ss)
    x,y=xy(lon,lat); nav.append((dt.timestamp(),x,y))
line_geom={}
for line,(a,b) in WINDOWS.items():
    ta,tb=parse_iso(a),parse_iso(b)
    coords=[(x,y) for t,x,y in nav if ta<=t<=tb]
    if len(coords)<2: raise RuntimeError(f"line {line} insufficient nav")
    line_geom[line]=LineString(coords)

# Parse only XY from processed product. Depth column is intentionally ignored.
prod=ftp_fetch(PRODUCT_PATH)
sha=hashlib.sha256(prod).hexdigest()
if sha!=PRODUCT_SHA: raise RuntimeError(f"product SHA mismatch {sha}")
lons=[]; lats=[]; xs=[]; ys=[]
for raw in prod.decode("ascii","ignore").splitlines():
    s=raw.strip()
    if not s or s[0] in "#;!":continue
    p=s.replace(","," ").split()
    if len(p)<2:continue
    try:lon=float(p[0]);lat=float(p[1])
    except:continue
    if not(math.isfinite(lon) and math.isfinite(lat)):continue
    x,y=xy(lon,lat); lons.append(lon);lats.append(lat);xs.append(x);ys.append(y)
x=np.asarray(xs,float); y=np.asarray(ys,float); lon=np.asarray(lons,float); lat=np.asarray(lats,float)
pts=shapely.points(x,y)
lines=sorted(WINDOWS)
dist=np.empty((len(lines),len(x)),float)
for j,line in enumerate(lines):
    dist[j,:]=shapely.distance(pts,line_geom[line])
arg=np.argmin(dist,axis=0)
mind=dist[arg,np.arange(len(x))]
assigned=np.asarray([lines[i] for i in arg],dtype=int)
mask31=assigned==31
idx31=np.where(mask31)[0]
ref=mind[mask31]; ref_sorted=np.sort(ref)
pctl31=np.searchsorted(ref_sorted,ref,side="right")/len(ref_sorted)
coords31=np.column_stack([x[mask31],y[mask31]])
support_tree=cKDTree(coords31)

# Load 2008 FNV cross-swath segments only.
coverage={}
fnv_manifest={}
for survey,baseurl in SURVEYS.items():
    names=list_fnv(baseurl); segs=[]; manifest=[]
    def load_one(name):
        b=fetch(urljoin(baseurl,name),90)
        return name,hashlib.sha256(b).hexdigest(),parse_fnv_segments(b)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        futs=[ex.submit(load_one,n) for n in names]
        for fut in concurrent.futures.as_completed(futs):
            name,h,ss=fut.result(); segs.extend(ss); manifest.append({"file":name,"sha256":h,"segments":len(ss)})
    if not segs: raise RuntimeError(f"{survey}: no FNV segments")
    coverage[survey]={"segments":segs,"tree":STRtree(segs)}
    fnv_manifest[survey]=sorted(manifest,key=lambda z:z["file"])

def nearest_coverage_m(survey,px,py):
    obj=Point(px,py); rec=coverage[survey]
    k=rec["tree"].nearest(obj)
    return float(obj.distance(rec["segments"][int(k)]))

def excluded(px,py):
    for _,ela,elo,r in EXCLUDE:
        ex,ey=xy(elo,ela)
        if math.hypot(px-ex,py-ey)<r:return True
    return False

# Candidate records are only line31-attributed points.
records=[]
for j,orig in enumerate(idx31):
    records.append({
      "orig":int(orig),"x":float(x[orig]),"y":float(y[orig]),"lat":float(lat[orig]),"lon":float(lon[orig]),
      "cross_track_m":float(mind[orig]),"percentile":float(pctl31[j])
    })

selected={band:[] for band in BANDS}
selection_audit={band:{"band_candidates":0,"after_exclusion":0,"local_support_pass":0,"dual_coverage_pass":0} for band in BANDS}
for band,(lo,hi,targetp) in BANDS.items():
    cand=[r for r in records if lo<=r["percentile"]<=hi]
    selection_audit[band]["band_candidates"]=len(cand)
    cand.sort(key=lambda r:(abs(r["percentile"]-targetp),r["lat"],r["lon"]))
    for r in cand:
        if len(selected[band])>=3: break
        if excluded(r["x"],r["y"]): continue
        selection_audit[band]["after_exclusion"]+=1
        local_n=len(support_tree.query_ball_point([r["x"],r["y"]],1000.0))
        if local_n<100: continue
        selection_audit[band]["local_support_pass"]+=1
        cov={s:nearest_coverage_m(s,r["x"],r["y"]) for s in SURVEYS}
        if not all(v<=100.0 for v in cov.values()): continue
        selection_audit[band]["dual_coverage_pass"]+=1
        if any(math.hypot(r["x"]-q["x"],r["y"]-q["y"])<3000.0 for q in selected[band]): continue
        rec=dict(r);rec["local_line31_support_within_1000m"]=int(local_n);rec["coverage_distance_m"]=cov
        selected[band].append(rec)

for band in selected:
    for r in selected[band]:
        r.pop("x",None);r.pop("y",None);r.pop("orig",None)

out={
 "artifact_id":"JANUS-KUSTO-CD169-EM12-LINE31-CROSSTRACK-HOLDOUT-SELECTION-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CD169-EM12-LINE31-CROSSTRACK-HOLDOUT-SELECTION-PREREG-2026-09-23-v1.0.json",
 "depth_values_parsed":False,
 "navigation":{"sha256":hashlib.sha256(navb).hexdigest(),"parsed_rows":len(nav)},
 "processed_product":{"sha256":sha,"parsed_xy":len(x),"line31_attributed_soundings":int(np.sum(mask31))},
 "fnv_manifest":fnv_manifest,
 "selection_audit":selection_audit,
 "selected":selected,
 "selected_counts":{k:len(v) for k,v in selected.items()},
 "thresholds_relaxed":False,
 "claim_ceiling":"GEOMETRY_ONLY_HOLDOUT_SELECTION"
}
p=OUT/"JANUS-KUSTO-CD169-EM12-LINE31-CROSSTRACK-HOLDOUT-SELECTION-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({"selection_audit":selection_audit,"selected":selected,"selected_counts":out["selected_counts"]},indent=2))
