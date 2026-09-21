#!/usr/bin/env python3
"""Reproduce KUSTO Hannah/CD169 real-world navigation calibration.

Inputs are the original externally supplied cd169.veh_nav and commands.cfg.
The script never reads sonar intensity and emits calibration metrics only.
CD169 becomes calibration-only after use and cannot validate Synthetic-v2.
"""
from __future__ import annotations

import argparse, hashlib, json, math, re, statistics
from datetime import datetime, timezone
from pathlib import Path

R = 6371008.8
TARGET = (-3.8654180644718967, -12.142441475)
S0 = (-3.869, -12.145)
S1 = (-3.855, -12.135)
SYNTHETIC_CLOSEST = datetime.fromisoformat("2005-02-28T01:07:25.275+00:00")
SYNTHETIC_MC_P995_M = 60.852179
OLD_V2 = (-3.8649594146460964, -12.14180870234793)
ANCHORS = [
    ("00:00",-3.907,-12.153),("00:15",-3.895,-12.167),("00:29",-3.886,-12.151),
    ("01:00",-3.869,-12.145),("01:29",-3.855,-12.135),("02:00",-3.835,-12.128),
    ("02:30",-3.814,-12.128),("03:00",-3.795,-12.131),("03:30",-3.774,-12.134),
    ("04:00",-3.752,-12.136),("04:30",-3.734,-12.140),("05:00",-3.715,-12.146),
    ("05:30",-3.683,-12.181),("06:00",-3.670,-12.151),("06:31",-3.649,-12.154)
]

def sha256(p: Path):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda:f.read(1<<20),b""): h.update(c)
    return h.hexdigest()

def hav(a,b):
    p1,l1=map(math.radians,a);p2,l2=map(math.radians,b)
    dp=p2-p1;dl=l2-l1
    q=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R*math.asin(min(1.0,math.sqrt(q)))

def bearing(a,b):
    p1,l1=map(math.radians,a);p2,l2=map(math.radians,b);dl=l2-l1
    y=math.sin(dl)*math.cos(p2)
    x=math.cos(p1)*math.sin(p2)-math.sin(p1)*math.cos(p2)*math.cos(dl)
    return (math.degrees(math.atan2(y,x))+360)%360

def en(origin,point):
    p0,l0=map(math.radians,origin);p,l=map(math.radians,point)
    return ((l-l0)*math.cos((p+p0)/2)*R,(p-p0)*R)

def pct(xs,p):
    s=sorted(xs);q=(len(s)-1)*p;i=int(math.floor(q));f=q-i
    return s[i] if i==len(s)-1 else s[i]*(1-f)+s[i+1]*f

def stat(xs):
    return {"n":len(xs),"min_m":min(xs),"median_m":statistics.median(xs),
            "mean_m":statistics.mean(xs),"p90_m":pct(xs,.9),"p95_m":pct(xs,.95),
            "max_m":max(xs)}

def parse(path: Path):
    out=[]
    for line in path.read_text(encoding="utf-8",errors="replace").splitlines():
        p=line.split()
        if len(p)!=9: continue
        try:
            dt=datetime(2000+int(p[1][:2]),int(p[1][2:4]),int(p[1][4:6]),
                        int(p[2][:2]),int(p[2][2:4]),tzinfo=timezone.utc)
            vals=list(map(float,p[3:]))
        except Exception: continue
        out.append({"utc":dt,"ship":(vals[0],vals[1]),"wire":vals[2],
                    "aux":vals[3],"veh":(vals[4],vals[5])})
    return sorted(out,key=lambda x:x["utc"])

def iso(x): return x.isoformat().replace("+00:00","Z")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--veh-nav",required=True,type=Path)
    ap.add_argument("--commands",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    a=ap.parse_args()

    rows=parse(a.veh_nav)
    bytime={r["utc"].strftime("%H:%M"):r for r in rows if r["utc"].date()==datetime(2005,2,28).date()}
    line=[r for r in rows if datetime(2005,2,28,1,0,tzinfo=timezone.utc)<=r["utc"]<=datetime(2005,2,28,1,29,tzinfo=timezone.utc)]
    if len(line)!=30: raise SystemExit(f"expected 30 line19 minute rows, got {len(line)}")

    cmd=a.commands.read_text(encoding="utf-8",errors="replace")
    m=re.search(r"^\s*mrgnav_inertia\b[^\n]*?\s-u\s+(\d+(?:\.\d+)?)\s+-n\s+(\S+)",cmd,re.M)

    h=bearing(S0,S1);hr=math.radians(h)
    residual=[];along=[];cross=[]
    for k,r in enumerate(line):
        f=k/29
        s=(S0[0]+f*(S1[0]-S0[0]),S0[1]+f*(S1[1]-S0[1]))
        residual.append(hav(s,r["veh"]))
        e,n=en(s,r["veh"])
        along.append(e*math.sin(hr)+n*math.cos(hr))
        cross.append(e*math.cos(hr)-n*math.sin(hr))

    anchor_res=[]
    for hhmm,lat,lon in ANCHORS:
        r=bytime[hhmm]
        anchor_res.append(hav((lat,lon),r["veh"]))

    best=None
    for i in range(len(line)-1):
        for sec in range(60):
            f=sec/60
            pt=(line[i]["veh"][0]+f*(line[i+1]["veh"][0]-line[i]["veh"][0]),
                line[i]["veh"][1]+f*(line[i+1]["veh"][1]-line[i]["veh"][1]))
            d=hav(TARGET,pt)
            dt=line[i]["utc"].replace(second=sec)
            if best is None or d<best[0]: best=(d,dt,pt,i)

    r7,r8=line[7],line[8]
    native=(r7["veh"][0]+25/60*(r8["veh"][0]-r7["veh"][0]),
            r7["veh"][1]+25/60*(r8["veh"][1]-r7["veh"][1]))
    i=best[3]
    aa=line[max(0,i-1)]["veh"];bb=line[min(len(line)-1,i+1)]["veh"]
    heading=bearing(aa,bb);tb=bearing(best[2],TARGET)
    rel=((tb-heading+180)%360)-180

    out={
      "status":"CALIBRATION_CORPUS_READY",
      "source":{"veh_nav_sha256":sha256(a.veh_nav),"veh_nav_rows":len(rows),
                "commands_cfg_sha256":sha256(a.commands),
                "mrgnav_u_parameter_m":float(m.group(1)) if m else None},
      "line19":{"position_residual":stat(residual),"along_error":stat(along),"cross_error":stat(cross)},
      "science_log_anchors":{"position_residual":stat(anchor_res)},
      "target":{"native_closest_utc":iso(best[1]),"native_closest_distance_m":best[0],
                "time_shift_vs_synthetic_s":(best[1]-SYNTHETIC_CLOSEST).total_seconds(),
                "native_track_heading_deg":heading,"bearing_to_target_deg":tb,
                "relative_bearing_deg":rel,"side":"PORT" if rel<0 else "STARBOARD",
                "native_distance_at_synthetic_time_m":hav(TARGET,native),
                "old_v2_to_native_at_synthetic_time_m":hav(OLD_V2,native),
                "synthetic_rounding_only_mc_p995_m":SYNTHETIC_MC_P995_M},
      "firewall":{"role":"CALIBRATION_ONLY","cd169_can_validate_v2":False,
                  "sonar_intensity_used":False,"morphology_claim":False}
    }
    a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2))

if __name__=="__main__": main()
