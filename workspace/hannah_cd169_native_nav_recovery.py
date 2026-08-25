#!/usr/bin/env python3
"""Recover CD169 TOBI PRISM navigation/cable inputs around the frozen target.

This stage intentionally runs before sample-level sonar inspection. It hashes
and parses only README.txt, cd169.nav and cd169.cable from the archived TOBI
PRISM package, and reports time-bracketing/interpolation without recentering.
"""
from __future__ import annotations

import argparse
import ftplib
import hashlib
import json
import math
import re
import statistics
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

HOST = "livftp.noc.ac.uk"
BASE = "/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281"
README = BASE + "/README.txt"
NAV = BASE + "/TOBI/cd169.nav"
CABLE = BASE + "/TOBI/cd169.cable"
TARGET = datetime(2005, 2, 28, 1, 7, 25, tzinfo=timezone.utc)
FROZEN = (-3.8654180644718967, -12.142441475)


def ftp() -> ftplib.FTP:
    f = ftplib.FTP(timeout=60)
    f.connect(HOST, 21)
    f.login("anonymous", "janus-probe@example.invalid")
    f.voidcmd("TYPE I")
    return f


def get(path: str) -> bytes:
    f = ftp(); out: list[bytes] = []
    try:
        f.retrbinary("RETR " + path, out.append, blocksize=262144)
    finally:
        try: f.quit()
        except Exception:
            try: f.close()
            except Exception: pass
    return b"".join(out)


def dt_from_fields(date: str, hm: str) -> datetime | None:
    try:
        if not re.fullmatch(r"\d{6}", date) or not re.fullmatch(r"\d{4}", hm): return None
        yy, mo, dd = int(date[:2]), int(date[2:4]), int(date[4:6])
        hh, mm = int(hm[:2]), int(hm[2:4])
        return datetime(2000 + yy, mo, dd, hh, mm, tzinfo=timezone.utc)
    except Exception: return None


def parse_nav(raw: bytes) -> list[dict[str, Any]]:
    out=[]
    for no, line in enumerate(raw.decode("utf-8",errors="replace").splitlines(),1):
        p=line.split()
        if len(p) < 5: continue
        dt=dt_from_fields(p[1],p[2])
        if dt is None: continue
        try: lat=float(p[3]); lon=float(p[4])
        except ValueError: continue
        if not (-90<=lat<=90 and -180<=lon<=180): continue
        out.append({"line":no,"datetime_utc":dt,"latitude":lat,"longitude":lon,"raw":line[:300]})
    return out


def parse_cable(raw: bytes) -> list[dict[str, Any]]:
    out=[]
    for no,line in enumerate(raw.decode("utf-8",errors="replace").splitlines(),1):
        p=line.split()
        if len(p)<4: continue
        dt=dt_from_fields(p[1],p[2])
        if dt is None: continue
        nums=[]
        for x in p[3:]:
            try: nums.append(float(x))
            except ValueError: pass
        if not nums: continue
        # Archive preview shows the final numeric field as the cable series.
        # Retain every numeric field so this choice remains auditable.
        out.append({"line":no,"datetime_utc":dt,"numeric_fields":nums,"selected_last_numeric":nums[-1],"raw":line[:300]})
    return out


def epoch(x: dict[str,Any]) -> float: return x["datetime_utc"].timestamp()


def bracket(rows: list[dict[str,Any]], target: datetime) -> tuple[dict|None,dict|None]:
    before=None; after=None; t=target.timestamp()
    for r in rows:
        e=epoch(r)
        if e<=t and (before is None or e>epoch(before)): before=r
        if e>=t and (after is None or e<epoch(after)): after=r
    return before,after


def lerp(a: float,b: float,f: float)->float: return a+(b-a)*f


def interp_nav(a: dict|None,b: dict|None,target:datetime)->dict|None:
    if a is None or b is None:return None
    ta,tb=epoch(a),epoch(b); t=target.timestamp()
    if tb==ta:f=0.0
    else:f=(t-ta)/(tb-ta)
    return {"fraction":f,"latitude":lerp(a["latitude"],b["latitude"],f),"longitude":lerp(a["longitude"],b["longitude"],f)}


def interp_cable(a:dict|None,b:dict|None,target:datetime)->dict|None:
    if a is None or b is None:return None
    ta,tb=epoch(a),epoch(b); t=target.timestamp(); f=0.0 if tb==ta else (t-ta)/(tb-ta)
    return {"fraction":f,"selected_last_numeric":lerp(a["selected_last_numeric"],b["selected_last_numeric"],f)}


def hav_km(a:tuple[float,float],b:tuple[float,float])->float:
    r=6371.0088; p1,l1=map(math.radians,a); p2,l2=map(math.radians,b)
    dp=p2-p1; dl=l2-l1
    h=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*r*math.asin(min(1,math.sqrt(h)))


def serial(r:dict|None)->dict|None:
    if r is None:return None
    o=dict(r); o["datetime_utc"]=r["datetime_utc"].isoformat().replace("+00:00","Z"); return o


def cadence(rows:list[dict[str,Any]])->dict:
    ds=[epoch(b)-epoch(a) for a,b in zip(rows,rows[1:]) if 0<epoch(b)-epoch(a)<86400]
    if not ds:return {"median_seconds":None,"mad_seconds":None}
    m=statistics.median(ds); return {"median_seconds":float(m),"mad_seconds":float(statistics.median(abs(x-m) for x in ds))}


def around(rows:list[dict[str,Any]], minutes:int=15)->list[dict]:
    lo=TARGET-timedelta(minutes=minutes); hi=TARGET+timedelta(minutes=minutes)
    return [serial(r) for r in rows if lo<=r["datetime_utc"]<=hi]


def main()->int:
    ap=argparse.ArgumentParser(); ap.add_argument("--output",required=True,type=Path); a=ap.parse_args()
    result={"schema":"janus.cosmos.cousteau.cd169_native_nav_recovery.v1","target_utc":TARGET.isoformat().replace("+00:00","Z"),"frozen_coordinate":{"latitude":FROZEN[0],"longitude":FROZEN[1]},"status":"STARTED","sonar_samples_inspected":False,"no_recenter":True,"scientific_claim":False}
    try:
        rr=get(README); nr=get(NAV); cr=get(CABLE)
        readme=rr.decode("utf-8",errors="replace")
        nav=parse_nav(nr); cable=parse_cable(cr)
        na,nb=bracket(nav,TARGET); ca,cb=bracket(cable,TARGET)
        ni=interp_nav(na,nb,TARGET); ci=interp_cable(ca,cb,TARGET)
        result["source_integrity"]={
          "readme":{"path":"sd11281/README.txt","size_bytes":len(rr),"sha256":hashlib.sha256(rr).hexdigest()},
          "nav":{"path":"sd11281/TOBI/cd169.nav","size_bytes":len(nr),"sha256":hashlib.sha256(nr).hexdigest(),"parsed_rows":len(nav),"cadence":cadence(nav)},
          "cable":{"path":"sd11281/TOBI/cd169.cable","size_bytes":len(cr),"sha256":hashlib.sha256(cr).hexdigest(),"parsed_rows":len(cable),"cadence":cadence(cable)}
        }
        result["archive_readme_text"]=readme[:4000]
        result["archive_readme_prism_statement_present"]=(".nav" in readme and ".cable" in readme and "PRISM" in readme.upper())
        result["nav_target"]={"before":serial(na),"after":serial(nb),"interpolated":ni,"rows_plus_minus_15min":around(nav)}
        result["cable_target"]={"before":serial(ca),"after":serial(cb),"interpolated":ci,"rows_plus_minus_15min":around(cable)}
        if ni:
            result["nav_interpolated_to_frozen_distance_km"]=hav_km((ni["latitude"],ni["longitude"]),FROZEN)
        result["classification"]={
          "nav_archive_role":"PRISM_TOBI_NAV_INPUT",
          "basis":"README states CDF files generated from raw TOBI files are accompanied by relevant .nav/.cable files required for PRISM processing; files are stored under TOBI/.",
          "coordinate_semantics":"TO_BE_VALIDATED_BY_TARGET_ALIGNMENT_AND_PRISM/CDF_METADATA",
          "may_pass_G5_yet":False
        }
        result["status"]="NATIVE_PRISM_NAV_AND_CABLE_TARGET_SLICE_READY" if ni and ci else "NAV_OR_CABLE_TARGET_BRACKET_INCOMPLETE"
    except Exception as e:
        result["status"]="NATIVE_NAV_RECOVERY_FAILED"; result["error_type"]=type(e).__name__; result["error"]=str(e)
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"nav_target":(result.get("nav_target") or {}).get("interpolated"),"cable_target":(result.get("cable_target") or {}).get("interpolated"),"distance_km":result.get("nav_interpolated_to_frozen_distance_km")},indent=2))
    return 0 if result["status"] in {"NATIVE_PRISM_NAV_AND_CABLE_TARGET_SLICE_READY","NAV_OR_CABLE_TARGET_BRACKET_INCOMPLETE"} else 2

if __name__=="__main__": raise SystemExit(main())
