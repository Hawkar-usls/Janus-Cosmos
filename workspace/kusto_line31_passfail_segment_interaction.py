#!/usr/bin/env python3
from __future__ import annotations
import ftplib, hashlib, json, math, statistics
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import LineString

HOST="livftp.noc.ac.uk"
ROOT="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281"
NAV_PATH=ROOT+"/Nav/cd169leg1_1min.listit"
PRODUCT_PATH=ROOT+"/EM12/B1-81-1_Acceptl28-33.xyz.ascii"
PRODUCT_SHA="0e8cac049b53dd6ba0e5013e8aaea2bfb48347db887edc3c4e3e9f697e96dcf0"
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
M=111320.0; LAT0=-4.0; COS0=math.cos(math.radians(LAT0)); R=6371008.8

TARGETS={
 "CAND003_FAIL":{"lat":-4.015075679897318,"lon":-12.29915403590432,"inherited_temporal_status":"FAIL_LINE31_VS_BOTH_2008_LINEAGES"},
 "MAR4DEG02S_PASS":{"lat":-4.03,"lon":-12.25,"inherited_temporal_status":"PASS_LINE31_VS_BOTH_2008_LINEAGES"}
}
WINDOWS={
 28:("2005-02-25T07:30:00Z","2005-02-25T21:10:00Z"),
 29:("2005-02-25T21:10:00Z","2005-02-26T07:34:00Z"),
 30:("2005-02-26T08:22:00Z","2005-02-26T09:35:00Z"),
 31:("2005-02-26T19:11:00Z","2005-02-27T14:40:00Z"),
 32:("2005-02-27T14:40:00Z","2005-02-28T07:31:00Z"),
 33:("2005-02-28T07:31:00Z","2005-02-28T14:59:00Z")
}

def xy(lon,lat): return lon*M*COS0,lat*M
def parse_iso(s): return datetime.fromisoformat(s.replace("Z","+00:00")).timestamp()
def iso(t): return datetime.fromtimestamp(t,tz=timezone.utc).isoformat().replace("+00:00","Z")

def ftp_fetch(path):
    f=ftplib.FTP(timeout=120); f.connect(HOST,21); f.login("anonymous","janus-kusto@example.invalid"); f.voidcmd("TYPE I"); chunks=[]
    try: f.retrbinary("RETR "+path,chunks.append)
    finally:
        try:f.quit()
        except:f.close()
    return b"".join(chunks)

def hav(lat1,lon1,lat2,lon2):
    p1,p2=math.radians(lat1),math.radians(lat2)
    dp=p2-p1; dl=math.radians(lon2-lon1)
    a=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R*math.asin(min(1.0,math.sqrt(a)))

def bearing(lat1,lon1,lat2,lon2):
    p1,p2=math.radians(lat1),math.radians(lat2); dl=math.radians(lon2-lon1)
    y=math.sin(dl)*math.cos(p2)
    x=math.cos(p1)*math.sin(p2)-math.sin(p1)*math.cos(p2)*math.cos(dl)
    return (math.degrees(math.atan2(y,x))+360)%360

def angle_diff(a,b): return ((a-b+180)%360)-180

def pct(xs,p):
    if not xs:return None
    return float(np.percentile(np.asarray(xs,float),p*100))

def circ_std(xs):
    if not xs:return None
    sx=sum(math.cos(math.radians(x)) for x in xs)/len(xs)
    sy=sum(math.sin(math.radians(x)) for x in xs)/len(xs)
    rho=max(1e-15,min(1.0,math.hypot(sx,sy)))
    return math.degrees(math.sqrt(max(0.0,-2*math.log(rho))))

navb=ftp_fetch(NAV_PATH)
nav=[]; base=datetime(2005,1,1,tzinfo=timezone.utc)
for raw in navb.decode("ascii","ignore").splitlines():
    p=raw.split()
    if len(p)<7: continue
    try:
        jd=int(p[1]); hh,mm,ss=map(int,p[2].split(":")); lat=float(p[3]); lon=float(p[5])
    except: continue
    dt=base+timedelta(days=jd-1,hours=hh,minutes=mm,seconds=ss)
    x,y=xy(lon,lat)
    nav.append({"t":dt.timestamp(),"lat":lat,"lon":lon,"x":x,"y":y})
nav.sort(key=lambda r:r["t"])

line_rows={}
line_geom={}
for line,(a,b) in WINDOWS.items():
    ta,tb=parse_iso(a),parse_iso(b)
    rows=[r for r in nav if ta<=r["t"]<=tb]
    line_rows[line]=rows
    coords=[(r["x"],r["y"]) for r in rows]
    if len(coords)<2: raise RuntimeError(f"line {line} insufficient nav")
    line_geom[line]=LineString(coords)

prod=ftp_fetch(PRODUCT_PATH)
sha=hashlib.sha256(prod).hexdigest()
if sha!=PRODUCT_SHA: raise RuntimeError(f"product SHA mismatch {sha}")
xs=[];ys=[]
for raw in prod.decode("ascii","ignore").splitlines():
    s=raw.strip()
    if not s or s[0] in "#;!": continue
    p=s.replace(","," ").split()
    if len(p)<2: continue
    try: lon=float(p[0]); lat=float(p[1])
    except: continue
    if not(math.isfinite(lon) and math.isfinite(lat)): continue
    x,y=xy(lon,lat); xs.append(x); ys.append(y)
x=np.asarray(xs,float); y=np.asarray(ys,float)
pts=shapely.points(x,y)
lines=sorted(WINDOWS)
dist=np.empty((len(lines),len(x)),float)
for j,line in enumerate(lines):
    dist[j,:]=shapely.distance(pts,line_geom[line])
arg=np.argmin(dist,axis=0)
mind=dist[arg,np.arange(len(x))]
assigned=np.asarray([lines[i] for i in arg],dtype=int)
ref=mind[assigned==31]
ref_sorted=np.sort(ref)

def navstate(target):
    rows=line_rows[31]
    nearest=min(rows,key=lambda r:hav(target["lat"],target["lon"],r["lat"],r["lon"]))
    nd=hav(target["lat"],target["lon"],nearest["lat"],nearest["lon"])
    a=max(parse_iso(WINDOWS[31][0]),nearest["t"]-1800)
    b=min(parse_iso(WINDOWS[31][1]),nearest["t"]+1800)
    win=[r for r in rows if a<=r["t"]<=b]
    speeds=[]; heads=[]; dts=[]
    for u,v in zip(win,win[1:]):
        dt=v["t"]-u["t"]
        if dt<=0:continue
        d=hav(u["lat"],u["lon"],v["lat"],v["lon"])
        speeds.append(d/dt*1.9438444924406)
        heads.append(bearing(u["lat"],u["lon"],v["lat"],v["lon"]))
        dts.append(dt)
    rates=[]
    for i in range(len(heads)-1):
        dm=max(1e-9,(dts[i]+dts[i+1])/120.0)
        rates.append(abs(angle_diff(heads[i+1],heads[i]))/dm)
    return {
      "nearest_nav_utc":iso(nearest["t"]),
      "nearest_ship_track_distance_m":nd,
      "analysis_window_utc":[iso(a),iso(b)],
      "analysis_window_points":len(win),
      "median_speed_knots":None if not speeds else float(statistics.median(speeds)),
      "speed_p10_knots":pct(speeds,.10),
      "speed_p90_knots":pct(speeds,.90),
      "circular_heading_std_deg":circ_std(heads),
      "median_abs_heading_change_deg_per_min":None if not rates else float(statistics.median(rates)),
      "p95_abs_heading_change_deg_per_min":pct(rates,.95)
    }

def geometry(target):
    cx,cy=xy(target["lon"],target["lat"])
    mask=(assigned==31)&(np.hypot(x-cx,y-cy)<=1000.0)
    td=mind[mask]
    pctl=np.searchsorted(ref_sorted,td,side="right")/len(ref_sorted) if len(td) and len(ref_sorted) else np.asarray([])
    return {
      "line31_all_attributed_soundings":int(len(ref)),
      "target_region_line31_soundings":int(len(td)),
      "median_cross_track_distance_to_line31_m":None if not len(td) else float(np.median(td)),
      "median_empirical_cross_track_percentile":None if not len(pctl) else float(np.median(pctl)),
      "fraction_empirical_percentile_ge_0p95":None if not len(pctl) else float(np.mean(pctl>=.95)),
      "empirical_percentile_p10_p50_p90":None if not len(pctl) else [float(v) for v in np.percentile(pctl,[10,50,90])],
      "line31_reference_cross_track_m_p50_p80_p90_p95_p99":[float(v) for v in np.percentile(ref,[50,80,90,95,99])]
    }

results={}
for name,t in TARGETS.items():
    results[name]={
      "inherited_temporal_status":t["inherited_temporal_status"],
      "geometry":geometry(t),
      "navstate":navstate(t)
    }

a=results["CAND003_FAIL"]; b=results["MAR4DEG02S_PASS"]
def ratio(v1,v2):
    return None if v1 is None or v2 in (None,0) else float(v1/v2)
comparison={
 "cross_track_distance_ratio_fail_over_pass":ratio(a["geometry"]["median_cross_track_distance_to_line31_m"],b["geometry"]["median_cross_track_distance_to_line31_m"]),
 "empirical_percentile_difference_fail_minus_pass":None if a["geometry"]["median_empirical_cross_track_percentile"] is None or b["geometry"]["median_empirical_cross_track_percentile"] is None else float(a["geometry"]["median_empirical_cross_track_percentile"]-b["geometry"]["median_empirical_cross_track_percentile"]),
 "nearest_ship_track_distance_ratio_fail_over_pass":ratio(a["navstate"]["nearest_ship_track_distance_m"],b["navstate"]["nearest_ship_track_distance_m"]),
 "heading_std_ratio_fail_over_pass":ratio(a["navstate"]["circular_heading_std_deg"],b["navstate"]["circular_heading_std_deg"]),
 "p95_turn_rate_ratio_fail_over_pass":ratio(a["navstate"]["p95_abs_heading_change_deg_per_min"],b["navstate"]["p95_abs_heading_change_deg_per_min"]),
 "median_speed_ratio_fail_over_pass":ratio(a["navstate"]["median_speed_knots"],b["navstate"]["median_speed_knots"])
}
out={
 "artifact_id":"JANUS-KUSTO-CD169-EM12-LINE31-PASSFAIL-SEGMENT-INTERACTION-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CD169-EM12-LINE31-PASSFAIL-SEGMENT-INTERACTION-PREREG-2026-09-23-v1.0.json",
 "navigation":{"path":NAV_PATH,"sha256":hashlib.sha256(navb).hexdigest(),"parsed_rows":len(nav)},
 "processed_product":{"path":PRODUCT_PATH,"sha256":sha,"parsed_xy":len(x)},
 "results":results,
 "comparison":comparison,
 "depth_used":False,
 "interpretation":"Descriptive same-line31 PASS-vs-FAIL segment comparison only; differences do not establish causality.",
 "claim_ceiling":"SAME_LINE31_LOCAL_SEGMENT_INTERACTION_DIAGNOSTIC_ONLY"
}
p=OUT/"JANUS-KUSTO-CD169-EM12-LINE31-PASSFAIL-SEGMENT-INTERACTION-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
