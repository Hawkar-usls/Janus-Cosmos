#!/usr/bin/env python3
from __future__ import annotations
import ftplib, hashlib, json, math, statistics
from datetime import datetime, timedelta, timezone
from pathlib import Path

HOST="livftp.noc.ac.uk"
NAV_PATH="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/Nav/cd169leg1_1min.listit"
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
LINE={"start":"2005-02-26T19:11:00Z","end":"2005-02-27T14:40:00Z"}
TARGETS=[
 {"id":"CAND003_FAIL","lat":-4.015075679897318,"lon":-12.29915403590432,"role":"temporal_fail"},
 {"id":"MAR4DEG02S_KNOWNPOS_PASS","lat":-4.03,"lon":-12.25,"role":"temporal_pass"}
]
R=6371008.8

def ftp():
    f=ftplib.FTP(timeout=120); f.connect(HOST,21); f.login("anonymous","janus-kusto@example.invalid"); f.voidcmd("TYPE I"); return f

def iso_ts(s): return datetime.fromisoformat(s.replace("Z","+00:00")).timestamp()
def iso(t): return datetime.fromtimestamp(t,tz=timezone.utc).isoformat().replace("+00:00","Z")

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

def percentile(xs,p):
    s=sorted(xs)
    if not s:return None
    q=(len(s)-1)*p; i=int(math.floor(q)); f=q-i
    return s[i] if i==len(s)-1 else s[i]*(1-f)+s[i+1]*f

def circular_mean_deg(xs):
    if not xs:return None
    sx=sum(math.cos(math.radians(x)) for x in xs); sy=sum(math.sin(math.radians(x)) for x in xs)
    return (math.degrees(math.atan2(sy,sx))+360)%360

def circular_std_deg(xs):
    if not xs:return None
    sx=sum(math.cos(math.radians(x)) for x in xs)/len(xs)
    sy=sum(math.sin(math.radians(x)) for x in xs)/len(xs)
    rho=max(1e-15,min(1.0,math.hypot(sx,sy)))
    return math.degrees(math.sqrt(max(0.0,-2*math.log(rho))))

def load_nav():
    f=ftp(); chunks=[]
    try:f.retrbinary("RETR "+NAV_PATH,chunks.append)
    finally:
        try:f.quit()
        except:f.close()
    b=b"".join(chunks)
    rows=[]; base=datetime(2005,1,1,tzinfo=timezone.utc)
    for raw in b.decode("ascii","ignore").splitlines():
        p=raw.split()
        if len(p)<7: continue
        try:
            jd=int(p[1]); hh,mm,ss=map(int,p[2].split(":"))
            lat=float(p[3]); lon=float(p[5])
        except: continue
        dt=base+timedelta(days=jd-1,hours=hh,minutes=mm,seconds=ss)
        rows.append({"t":dt.timestamp(),"lat":lat,"lon":lon})
    rows.sort(key=lambda r:r["t"])
    return rows,hashlib.sha256(b).hexdigest()

def summarize_target(line_rows,target):
    nearest=min(line_rows,key=lambda r:hav(target["lat"],target["lon"],r["lat"],r["lon"]))
    nd=hav(target["lat"],target["lon"],nearest["lat"],nearest["lon"])
    a,b=iso_ts(LINE["start"]),iso_ts(LINE["end"])
    w0=max(a,nearest["t"]-1800); w1=min(b,nearest["t"]+1800)
    win=[r for r in line_rows if w0<=r["t"]<=w1]
    speeds=[]; heads=[]; seg_lengths=[]; seg_dt=[]
    for x,y in zip(win,win[1:]):
        dt=y["t"]-x["t"]
        if dt<=0: continue
        d=hav(x["lat"],x["lon"],y["lat"],y["lon"])
        speeds.append(d/dt*1.9438444924406)
        heads.append(bearing(x["lat"],x["lon"],y["lat"],y["lon"]))
        seg_lengths.append(d); seg_dt.append(dt)
    rates=[]
    for (h0,h1),(dt0,dt1) in zip(zip(heads,heads[1:]),zip(seg_dt,seg_dt[1:])):
        dm=max(1e-9,(dt0+dt1)/120.0)
        rates.append(abs(angle_diff(h1,h0))/dm)
    return {
      "target":target,
      "nearest_track_point":{"utc":iso(nearest["t"]),"lat":nearest["lat"],"lon":nearest["lon"],"distance_m":nd},
      "analysis_window_utc":[iso(w0),iso(w1)],
      "analysis_window_epoch":[w0,w1],
      "analysis_window_points":len(win),
      "derived_speed_knots":{
        "median":statistics.median(speeds) if speeds else None,
        "p10":percentile(speeds,.10),"p90":percentile(speeds,.90)
      },
      "heading":{
        "circular_mean_deg":circular_mean_deg(heads),
        "circular_std_deg":circular_std_deg(heads),
        "median_abs_change_deg_per_min":statistics.median(rates) if rates else None,
        "p95_abs_change_deg_per_min":percentile(rates,.95)
      },
      "track_path_length_m":sum(seg_lengths),
      "window_duration_s": (win[-1]["t"]-win[0]["t"]) if len(win)>=2 else 0
    }

nav,sha=load_nav()
a,b=iso_ts(LINE["start"]),iso_ts(LINE["end"])
line=[r for r in nav if a<=r["t"]<=b]
res={t["id"]:summarize_target(line,t) for t in TARGETS}
A=res["CAND003_FAIL"]; B=res["MAR4DEG02S_KNOWNPOS_PASS"]

a0,a1=A["analysis_window_epoch"]; b0,b1=B["analysis_window_epoch"]
overlap=max(0.0,min(a1,b1)-max(a0,b0))
union=max(a1,b1)-min(a0,b0)
overlap_fraction=overlap/union if union>0 else None

def ratio(x,y):
    return None if y in (None,0) or x is None else x/y

comparison={
  "nearest_approach_time_difference_s":abs(
      datetime.fromisoformat(A["nearest_track_point"]["utc"].replace("Z","+00:00")).timestamp()
      - datetime.fromisoformat(B["nearest_track_point"]["utc"].replace("Z","+00:00")).timestamp()
  ),
  "analysis_window_overlap_s":overlap,
  "analysis_window_overlap_fraction_of_union":overlap_fraction,
  "cand003_over_knownpos_heading_std_ratio":ratio(A["heading"]["circular_std_deg"],B["heading"]["circular_std_deg"]),
  "cand003_over_knownpos_p95_turn_rate_ratio":ratio(A["heading"]["p95_abs_change_deg_per_min"],B["heading"]["p95_abs_change_deg_per_min"]),
  "cand003_over_knownpos_median_speed_ratio":ratio(A["derived_speed_knots"]["median"],B["derived_speed_knots"]["median"]),
  "cand003_minus_knownpos_nearest_track_distance_m":A["nearest_track_point"]["distance_m"]-B["nearest_track_point"]["distance_m"]
}
out={
 "artifact_id":"JANUS-KUSTO-LINE31-CAND003-VS-KNOWNPOS-NAVSEGMENT-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-LINE31-CAND003-VS-KNOWNPOS-NAVSEGMENT-PREREG-2026-09-23-v1.0.json",
 "navigation":{"path":NAV_PATH,"sha256":sha,"parsed_rows":len(nav),"line31_rows":len(line)},
 "line31_window_utc":[LINE["start"],LINE["end"]],
 "targets":res,
 "comparison":comparison,
 "depth_used":False,
 "processed_xyz_used":False,
 "interpretation":"Same-line navigation-state comparison only. Association does not establish causality.",
 "claim_ceiling":"SAME_LINE_LOCATION_SPECIFIC_NAVSTATE_MECHANISM_DIAGNOSTIC_ONLY"
}
p=OUT/"JANUS-KUSTO-LINE31-CAND003-VS-KNOWNPOS-NAVSEGMENT-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
