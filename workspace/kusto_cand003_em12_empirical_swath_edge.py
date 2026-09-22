#!/usr/bin/env python3
from __future__ import annotations
import ftplib, hashlib, json, math
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import LineString

C={"id":"KN19207_CAND_003","lat":-4.015075679897318,"lon":-12.29915403590432}
M=111320.0; LAT0=-4.0; COS0=math.cos(math.radians(LAT0))
HOST="livftp.noc.ac.uk"
ROOT="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281"
EM12=ROOT+"/EM12"
NAV=ROOT+"/Nav/cd169leg1_1min.listit"
OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)

PRODUCTS={
 "B1-81-1_Acceptl15-22.xyz.ascii":{"sha":"ad2e4d60401f7ffaa5e14b3e43aa0e45da9d320c6b34fa7f3b873f4024be6007","lines":list(range(15,23)),"target_line":18},
 "B1-81-1_Acceptl28-33.xyz.ascii":{"sha":"0e8cac049b53dd6ba0e5013e8aaea2bfb48347db887edc3c4e3e9f697e96dcf0","lines":list(range(28,34)),"target_line":31}
}
WINDOWS={
15:("2005-02-23T17:37:00Z","2005-02-23T19:37:00Z"),16:("2005-02-23T19:37:00Z","2005-02-23T21:30:00Z"),
17:("2005-02-23T21:30:00Z","2005-02-24T00:02:00Z"),18:("2005-02-24T00:02:00Z","2005-02-24T01:59:00Z"),
19:("2005-02-24T01:59:00Z","2005-02-24T04:11:00Z"),20:("2005-02-24T04:11:00Z","2005-02-24T06:00:00Z"),
21:("2005-02-24T06:00:00Z","2005-02-24T06:25:00Z"),22:("2005-02-24T06:25:00Z","2005-02-24T06:35:00Z"),
28:("2005-02-25T07:30:00Z","2005-02-25T21:10:00Z"),29:("2005-02-25T21:10:00Z","2005-02-26T07:34:00Z"),
30:("2005-02-26T08:22:00Z","2005-02-26T09:35:00Z"),31:("2005-02-26T19:11:00Z","2005-02-27T14:40:00Z"),
32:("2005-02-27T14:40:00Z","2005-02-28T07:31:00Z"),33:("2005-02-28T07:31:00Z","2005-02-28T14:59:00Z")
}

def xy(lon,lat):return lon*M*COS0,lat*M
CX,CY=xy(C["lon"],C["lat"])

def ftp_fetch(path):
    f=ftplib.FTP(timeout=120);f.connect(HOST,21);f.login("anonymous","janus-kusto@example.invalid");f.voidcmd("TYPE I");chunks=[]
    try:f.retrbinary("RETR "+path,chunks.append)
    finally:
        try:f.quit()
        except:f.close()
    return b"".join(chunks)

def parse_iso(s):return datetime.fromisoformat(s.replace("Z","+00:00")).timestamp()

navb=ftp_fetch(NAV)
nav=[]
base=datetime(2005,1,1,tzinfo=timezone.utc)
for line in navb.decode("ascii","ignore").splitlines():
    p=line.split()
    if len(p)<7:continue
    try:
        yy=int(p[0]);jd=int(p[1]);hh,mm,ss=map(int,p[2].split(":"));lat=float(p[3]);lon=float(p[5])
    except:continue
    dt=base+timedelta(days=jd-1,hours=hh,minutes=mm,seconds=ss)
    x,y=xy(lon,lat);nav.append((dt.timestamp(),x,y))

line_geom={}
for line,(a,b) in WINDOWS.items():
    ta,tb=parse_iso(a),parse_iso(b)
    coords=[(x,y) for t,x,y in nav if ta<=t<=tb]
    if len(coords)<2:raise RuntimeError(f"line {line} has insufficient nav")
    line_geom[line]=LineString(coords)

def parse_xy_only(data):
    xs=[];ys=[]
    for raw in data.decode("ascii","ignore").splitlines():
        s=raw.strip()
        if not s or s[0] in "#;!":continue
        p=s.replace(","," ").split()
        if len(p)<2:continue
        try:lon=float(p[0]);lat=float(p[1])
        except:continue
        if not (math.isfinite(lon) and math.isfinite(lat)):continue
        x,y=xy(lon,lat);xs.append(x);ys.append(y)
    return np.asarray(xs,float),np.asarray(ys,float)

def summarize_product(name,meta):
    b=ftp_fetch(EM12+"/"+name)
    sha=hashlib.sha256(b).hexdigest()
    if sha!=meta["sha"]:raise RuntimeError(f"{name} SHA mismatch")
    x,y=parse_xy_only(b);n=len(x)
    # Shapely 2 vectorized exact point-to-LineString distances.
    pts=shapely.points(x,y)
    lines=meta["lines"]
    dist=np.empty((len(lines),n),dtype=np.float64)
    for j,line in enumerate(lines):
        dist[j,:]=shapely.distance(pts,line_geom[line])
    arg=np.argmin(dist,axis=0)
    mind=dist[arg,np.arange(n)]
    assigned=np.asarray([lines[i] for i in arg],dtype=int)
    tl=meta["target_line"]
    ref=mind[assigned==tl]
    tmask=(assigned==tl)&(np.hypot(x-CX,y-CY)<=1000.0)
    td=mind[tmask]
    ref_sorted=np.sort(ref)
    pctl=np.searchsorted(ref_sorted,td,side="right")/len(ref_sorted) if len(ref_sorted) and len(td) else np.array([])
    out={
      "sha256":sha,"parsed_xy":n,"target_line":tl,
      "all_soundings_attributed_to_target_line":int(len(ref)),
      "target_region_soundings_attributed_to_target_line":int(len(td)),
      "target_region_median_cross_track_distance_m":None if not len(td) else float(np.median(td)),
      "target_region_median_empirical_cross_track_percentile":None if not len(pctl) else float(np.median(pctl)),
      "target_region_fraction_percentile_ge_0p80":None if not len(pctl) else float(np.mean(pctl>=.80)),
      "target_region_fraction_percentile_ge_0p90":None if not len(pctl) else float(np.mean(pctl>=.90)),
      "target_region_fraction_percentile_ge_0p95":None if not len(pctl) else float(np.mean(pctl>=.95)),
      "target_region_percentile_p10_p50_p90":None if not len(pctl) else [float(v) for v in np.percentile(pctl,[10,50,90])],
      "reference_cross_track_m_p50_p80_p90_p95_p99":None if not len(ref) else [float(v) for v in np.percentile(ref,[50,80,90,95,99])]
    }
    return out

res={name:summarize_product(name,meta) for name,meta in PRODUCTS.items()}
a=res["B1-81-1_Acceptl15-22.xyz.ascii"]["target_region_median_empirical_cross_track_percentile"]
b=res["B1-81-1_Acceptl28-33.xyz.ascii"]["target_region_median_empirical_cross_track_percentile"]
positive=bool(a is not None and b is not None and b>=.80 and a<.80 and (b-a)>=.15)
out={
 "artifact_id":"JANUS-KUSTO-CAND003-CD169-EM12-EMPIRICAL-SWATH-EDGE-DIAGNOSTIC-RUN-2026-09-22-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CAND003-CD169-EM12-EMPIRICAL-SWATH-EDGE-DIAGNOSTIC-PREREG-2026-09-22-v1.0.json",
 "target":C,
 "nav_sha256":hashlib.sha256(navb).hexdigest(),
 "distance_engine":"Shapely vectorized exact Euclidean point-to-LineString distance in local equirectangular metric coordinates",
 "products":res,
 "median_percentile_gap_line31_minus_line18":None if a is None or b is None else float(b-a),
 "diagnostic_signature_pass":positive,
 "verdict":"EMPIRICAL_SWATH_EDGE_MECHANISM_SUPPORTED" if positive else "NO_EMPIRICAL_SWATH_EDGE_MECHANISM",
 "warning":"Empirical cross-track percentile is not an exact beam number or raw ping provenance.",
 "strict_parent_fail_immutable":True,
 "claim_ceiling":"GEOMETRIC_MECHANISM_DIAGNOSTIC_ONLY"
}
p=OUT/"JANUS-KUSTO-CAND003-CD169-EM12-EMPIRICAL-SWATH-EDGE-DIAGNOSTIC-RUN-2026-09-22-v1.0.json";p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
