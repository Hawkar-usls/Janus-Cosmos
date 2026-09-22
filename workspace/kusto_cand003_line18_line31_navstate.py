#!/usr/bin/env python3
from __future__ import annotations
import ftplib, hashlib, json, math, statistics
from datetime import datetime, timedelta, timezone
from pathlib import Path

HOST="livftp.noc.ac.uk"
NAV_PATH="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/Nav/cd169leg1_1min.listit"
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
TARGET={"id":"KN19207_CAND_003","lat":-4.015075679897318,"lon":-12.29915403590432}
LINES={
  "18":{"start":"2005-02-24T00:02:00Z","end":"2005-02-24T01:59:00Z","context":"planned swath survey / WP4-era line"},
  "31":{"start":"2005-02-26T19:11:00Z","end":"2005-02-27T14:40:00Z","context":"transit to TOBI 02 re-deployment"}
}
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

def circular_median_deg(xs):
    if not xs:return None
    mu=circular_mean_deg(xs)
    unwrapped=[mu+angle_diff(x,mu) for x in xs]
    return statistics.median(unwrapped)%360

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

def summarize_line(allrows,spec):
    a,b=iso_ts(spec["start"]),iso_ts(spec["end"])
    line=[r for r in allrows if a<=r["t"]<=b]
    nearest=min(line,key=lambda r:hav(TARGET["lat"],TARGET["lon"],r["lat"],r["lon"]))
    nd=hav(TARGET["lat"],TARGET["lon"],nearest["lat"],nearest["lon"])
    w0=max(a,nearest["t"]-1800); w1=min(b,nearest["t"]+1800)
    win=[r for r in line if w0<=r["t"]<=w1]
    speeds=[]; heads=[]; seg_lengths=[]; seg_dt=[]
    for x,y in zip(win,win[1:]):
        dt=y["t"]-x["t"]
        if dt<=0: continue
        d=hav(x["lat"],x["lon"],y["lat"],y["lon"])
        speeds.append(d/dt*1.9438444924406)
        heads.append(bearing(x["lat"],x["lon"],y["lat"],y["lon"]))
        seg_lengths.append(d); seg_dt.append(dt)
    changes=[]
    change_rates=[]
    for (h0,h1),(dt0,dt1) in zip(zip(heads,heads[1:]),zip(seg_dt,seg_dt[1:])):
        # Heading samples represent adjacent segments; normalize by average segment interval in minutes.
        dm=max(1e-9,(dt0+dt1)/120.0)
        ch=abs(angle_diff(h1,h0))
        changes.append(ch); change_rates.append(ch/dm)
    return {
      "line_window_utc":[spec["start"],spec["end"]],
      "source_context":spec["context"],
      "line_nav_points":len(line),
      "nearest_track_point":{
        "utc":iso(nearest["t"]),"lat":nearest["lat"],"lon":nearest["lon"],"distance_m":nd
      },
      "analysis_window_utc":[iso(w0),iso(w1)],
      "analysis_window_points":len(win),
      "window_is_clipped_by_line_boundary": bool(w0>nearest["t"]-1800 or w1<nearest["t"]+1800),
      "derived_speed_knots":{
        "n":len(speeds),
        "median":statistics.median(speeds) if speeds else None,
        "p10":percentile(speeds,.10),"p90":percentile(speeds,.90)
      },
      "heading":{
        "n":len(heads),
        "circular_median_deg":circular_median_deg(heads),
        "circular_mean_deg":circular_mean_deg(heads),
        "circular_std_deg":circular_std_deg(heads),
        "median_abs_change_deg_per_min":statistics.median(change_rates) if change_rates else None,
        "p95_abs_change_deg_per_min":percentile(change_rates,.95)
      },
      "track_path_length_m_in_window":sum(seg_lengths),
      "window_duration_s": (win[-1]["t"]-win[0]["t"]) if len(win)>=2 else 0
    }

nav,sha=load_nav()
res={line:summarize_line(nav,spec) for line,spec in LINES.items()}
l18,l31=res["18"],res["31"]
comparison={
 "nearest_distance_difference_line31_minus_line18_m":l31["nearest_track_point"]["distance_m"]-l18["nearest_track_point"]["distance_m"],
 "median_speed_ratio_line31_over_line18":l31["derived_speed_knots"]["median"]/l18["derived_speed_knots"]["median"] if l18["derived_speed_knots"]["median"] else None,
 "heading_std_ratio_line31_over_line18":l31["heading"]["circular_std_deg"]/l18["heading"]["circular_std_deg"] if l18["heading"]["circular_std_deg"] else None,
 "median_turn_rate_ratio_line31_over_line18":l31["heading"]["median_abs_change_deg_per_min"]/l18["heading"]["median_abs_change_deg_per_min"] if l18["heading"]["median_abs_change_deg_per_min"] else None
}
out={
 "artifact_id":"JANUS-KUSTO-CAND003-CD169-EM12-LINE18-VS-LINE31-NAVSTATE-RUN-2026-09-22-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CAND003-CD169-EM12-LINE18-VS-LINE31-NAVSTATE-PREREG-2026-09-22-v1.0.json",
 "target":TARGET,
 "navigation":{"path":NAV_PATH,"sha256":sha,"parsed_rows":len(nav)},
 "lines":res,
 "comparison":comparison,
 "interpretation_rule":"Operational/nav-state differences may identify a mechanism candidate but do not prove bathymetric causality.",
 "depth_used":False,
 "processed_xyz_used":False,
 "claim_ceiling":"NAVIGATION_STATE_MECHANISM_DIAGNOSTIC_ONLY"
}
p=OUT/"JANUS-KUSTO-CAND003-CD169-EM12-LINE18-VS-LINE31-NAVSTATE-RUN-2026-09-22-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
