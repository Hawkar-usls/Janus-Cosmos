#!/usr/bin/env python3
import bisect, ftplib, hashlib, json, math, struct, zipfile
from datetime import datetime, timezone
from pathlib import Path

HOST="livftp.noc.ac.uk"
ZIP_PATH="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11280/EM12/EM12raw.zip"
ZIP_SHA="060b65865312330e866f842959a89674b744ff44569dcb1fa13daa7b5b0cefb0"
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
LOCAL=Path("/tmp/CD169_EM12raw.zip")
TARGETS=[
 {"id":"CELL_A","lat":-4.0165853,"lon":-12.2985657,"processed_2005_depth_m":4025.08,"KN192_2008_depth_m":3903.7895286232233,"KNOX15RR_2008_depth_m":3925.6959423840044},
 {"id":"CELL_B","lat":-4.0152462,"lon":-12.299318,"processed_2005_depth_m":4036.31,"KN192_2008_depth_m":3915.9501554071903,"KNOX15RR_2008_depth_m":3937.9770237505436}
]
BANDS=[50,100,250,500,1000]
BATH_LABELS={0x0294:"EM12DS_BATH",0x0295:"EM12DP_BATH",0x0296:"EM12S_BATH"}
POS_LABEL=0x0293
POS_SIZE=93
BATH_SIZE=926
R=6371008.8

# Exact ZIP acquisition
h=hashlib.sha256()
f=ftplib.FTP(timeout=180); f.connect(HOST,21); f.login("anonymous","janus-kusto@example.invalid"); f.voidcmd("TYPE I")
with LOCAL.open("wb") as out:
    def cb(b): out.write(b); h.update(b)
    f.retrbinary("RETR "+ZIP_PATH,cb,blocksize=1024*1024)
try:f.quit()
except:f.close()
if h.hexdigest()!=ZIP_SHA: raise RuntimeError("ZIP SHA mismatch")

def digits(b,a,z):
    s=b[a:z]
    if len(s)!=(z-a) or not all(48<=x<=57 for x in s): return None
    return int(s.decode("ascii"))

def parse_time(payload, pos_mode):
    try:
        dd=digits(payload,0,2); mm=digits(payload,2,4); yy=digits(payload,4,6)
        if pos_mode:
            hh=digits(payload,7,9); mi=digits(payload,9,11); ss=digits(payload,11,13); cs=digits(payload,13,15)
        else:
            hh=digits(payload,6,8); mi=digits(payload,8,10); ss=digits(payload,10,12); cs=digits(payload,12,14)
        if None in (dd,mm,yy,hh,mi,ss,cs): return None
        year=2000+yy if yy<80 else 1900+yy
        dt=datetime(year,mm,dd,hh,mi,ss,cs*10000,tzinfo=timezone.utc)
        return dt.timestamp(),dt.isoformat().replace("+00:00","Z")
    except Exception:
        return None

def parse_float_ascii(b,a,z):
    try:return float(b[a:z].decode("ascii","ignore").strip())
    except:return None

def parse_pos(payload,offset):
    tt=parse_time(payload,True)
    if tt is None:return None
    try:
        deg=digits(payload,16,18); minute=parse_float_ascii(payload,18,25); hemi=chr(payload[25])
        ldeg=digits(payload,27,30); lmin=parse_float_ascii(payload,30,37); lhemi=chr(payload[37])
        if None in (deg,minute,ldeg,lmin):return None
        lat=deg+minute/60.0
        if hemi.lower()=="s":lat=-lat
        lon=ldeg+lmin/60.0
        if lhemi.lower()=="w":lon=-lon
        speed=parse_float_ascii(payload,80,84) or 0.0
        line_heading=parse_float_ascii(payload,85,90) or 0.0
        if not (-90<=lat<=90 and -180<=lon<=180):return None
        return {"offset":offset,"epoch":tt[0],"utc":tt[1],"lat":lat,"lon":lon,"speed_mps":speed,"line_heading_deg":line_heading}
    except:return None

def parse_bath_geometry(payload,offset,tname):
    tt=parse_time(payload,False)
    if tt is None:return None
    if len(payload)<BATH_SIZE:return None
    try:
        ping_number=struct.unpack_from("<H",payload,14)[0]
        bath_res=payload[16]
        heading_raw=struct.unpack_from("<H",payload,20)[0]
        heading_deg=0.1*heading_raw
        if bath_res==1:
            offset_scale=0.2
            depth_scale=0.1
        elif bath_res==2:
            offset_scale=0.5
            depth_scale=0.2
        else:
            return None
        beams=[]
        for i in range(81):
            base=32+11*i
            stored_bath=struct.unpack_from("<h",payload,base)[0]
            across=struct.unpack_from("<h",payload,base+2)[0]*offset_scale
            along=struct.unpack_from("<h",payload,base+4)[0]*offset_scale
            quality=payload[base+9]
            depth_m=None if stored_bath==0 else depth_scale*stored_bath
            beams.append((i,across,along,quality,stored_bath,depth_m))
        return {"offset":offset,"epoch":tt[0],"utc":tt[1],"type":tname,"ping_number":ping_number,
                "bath_res":bath_res,"heading_deg":heading_deg,"beams":beams}
    except:return None

def coor_scale(lat):
    C1=111412.84; C2=-93.5; C3=0.118
    C4=111132.92; C5=-559.82; C6=1.175; C7=0.0023
    r=math.radians(lat)
    mtodeglon=1.0/abs(C1*math.cos(r)+C2*math.cos(3*r)+C3*math.cos(5*r))
    mtodeglat=1.0/abs(C4+C5*math.cos(2*r)+C6*math.cos(4*r)+C7*math.cos(6*r))
    return mtodeglon,mtodeglat

def nav_interp(fixes,t,current_speed,current_line_heading):
    if not fixes:return None
    times=[x["epoch"] for x in fixes]
    mtodeglon,mtodeglat=coor_scale(fixes[-1]["lat"])
    if len(fixes)>1 and fixes[0]["epoch"]<=t<=fixes[-1]["epoch"]:
        j=bisect.bisect_left(times,t)
        if j<=0:j=1
        a,b=fixes[j-1],fixes[j]
        den=b["epoch"]-a["epoch"]
        fac=0.0 if den==0 else (t-a["epoch"])/den
        return {"lon":a["lon"]+fac*(b["lon"]-a["lon"]),"lat":a["lat"]+fac*(b["lat"]-a["lat"]),
                "mode":"INTERPOLATED","nav_fix_before_utc":a["utc"],"nav_fix_after_utc":b["utc"]}
    anchor=fixes[-1] if t>fixes[-1]["epoch"] else fixes[0]
    dt=t-anchor["epoch"]
    dd=dt*current_speed
    hx=math.sin(math.radians(current_line_heading)); hy=math.cos(math.radians(current_line_heading))
    return {"lon":anchor["lon"]+hx*mtodeglon*dd,"lat":anchor["lat"]+hy*mtodeglat*dd,
            "mode":"EXTRAPOLATED","anchor_utc":anchor["utc"],"dt_s":dt,"raw_speed_mps":current_speed,
            "line_heading_deg":current_line_heading}

def beam_lonlat(nav,heading,across,along):
    mlon,mlat=coor_scale(nav["lat"])
    hx=math.sin(math.radians(heading)); hy=math.cos(math.radians(heading))
    lon=nav["lon"] + hy*mlon*across + hx*mlon*along
    lat=nav["lat"] - hx*mlat*across + hy*mlat*along
    return lon,lat

def hav(lat1,lon1,lat2,lon2):
    p1,p2=math.radians(lat1),math.radians(lat2)
    dp=p2-p1; dl=math.radians(lon2-lon1)
    a=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R*math.asin(min(1.0,math.sqrt(a)))

def scan_records(data):
    recs=[]
    n=len(data)
    for i in range(n-18):
        if data[i]!=0x02:continue
        typ=(data[i]<<8)|data[i+1]
        if typ==POS_LABEL and i+2+POS_SIZE<=n:
            p=parse_pos(data[i+2:i+2+POS_SIZE],i)
            if p and datetime(2005,2,1,tzinfo=timezone.utc).timestamp()<=p["epoch"]<=datetime(2005,4,1,tzinfo=timezone.utc).timestamp():
                recs.append((i,"POS",p))
        elif typ in BATH_LABELS and i+2+BATH_SIZE<=n:
            b=parse_bath_geometry(data[i+2:i+2+BATH_SIZE],i,BATH_LABELS[typ])
            if b and datetime(2005,2,1,tzinfo=timezone.utc).timestamp()<=b["epoch"]<=datetime(2005,4,1,tzinfo=timezone.utc).timestamp():
                recs.append((i,"BATH",b))
    # de-dup exact label offsets
    unique={}
    for r in recs: unique[(r[0],r[1])]=r
    return sorted(unique.values(),key=lambda r:r[0])

state={t["id"]:{"target":t,"nearest":None,"counts":{str(b):0 for b in BANDS},"selected_le100m":[]} for t in TARGETS}
entry_summaries=[]
total_bath=0; total_beams=0; navless_bath=0

with zipfile.ZipFile(LOCAL,"r") as z:
    for info in z.infolist():
        data=z.read(info.filename)
        records=scan_records(data)
        fixes=[]; current_speed=0.0; current_line_heading=0.0
        npos=0; nbath=0; nbeam=0
        for _,kind,obj in records:
            if kind=="POS":
                fixes.append(obj); current_speed=obj["speed_mps"]; current_line_heading=obj["line_heading_deg"]; npos+=1
            else:
                nbath+=1; total_bath+=1
                nav=nav_interp(fixes,obj["epoch"],current_speed,current_line_heading)
                if nav is None:
                    navless_bath+=1; continue
                for beam_index,across,along,quality,stored_bath,depth_m in obj["beams"]:
                    lon,lat=beam_lonlat(nav,obj["heading_deg"],across,along)
                    nbeam+=1; total_beams+=1
                    for t in TARGETS:
                        d=hav(t["lat"],t["lon"],lat,lon)
                        s=state[t["id"]]
                        for b in BANDS:
                            if d<=b:s["counts"][str(b)]+=1
                        cand={
                          "distance_m":d,"raw_file":info.filename,"ping_utc":obj["utc"],"ping_number":obj["ping_number"],
                          "bath_type":obj["type"],"beam_index":beam_index,"bath_res":obj["bath_res"],
                          "ping_heading_deg":obj["heading_deg"],"nav_lon":nav["lon"],"nav_lat":nav["lat"],
                          "nav_mode":nav["mode"],"raw_across_m":across,"raw_along_m":along,
                          "beam_lon":lon,"beam_lat":lat,"raw_quality_byte":quality,
                          "stored_bath":stored_bath,"raw_depth_m":depth_m
                        }
                        # Geometry-only selection was frozen before raw depth read.
                        if d<=100 and stored_bath!=0:
                            s["selected_le100m"].append(cand)
                        if s["nearest"] is None or d<s["nearest"]["distance_m"]:
                            s["nearest"]=cand
        entry_summaries.append({"filename":info.filename,"position_records":npos,"bath_records":nbath,"beam_geometries":nbeam})

def median(xs):
    xs=sorted(xs)
    n=len(xs)
    if n==0:return None
    return xs[n//2] if n%2 else 0.5*(xs[n//2-1]+xs[n//2])

def summarize_cell(s):
    beams=sorted(s["selected_le100m"],key=lambda x:(x["distance_m"],x["raw_file"],x["ping_utc"],x["beam_index"]))
    depths=[x["raw_depth_m"] for x in beams if x["raw_depth_m"] is not None]
    by_file={}
    for x in beams:
        by_file.setdefault(x["raw_file"],[]).append(x)
    groups={}
    for fn,bb in sorted(by_file.items()):
        dd=[x["raw_depth_m"] for x in bb if x["raw_depth_m"] is not None]
        groups[fn]={
          "n":len(bb),
          "median_depth_m":median(dd),
          "min_depth_m":min(dd) if dd else None,
          "max_depth_m":max(dd) if dd else None,
          "ping_times_utc":sorted(set(x["ping_utc"] for x in bb)),
          "beam_indices":sorted(set(x["beam_index"] for x in bb))
        }
    t=s["target"]
    rawmed=median(depths)
    proc=t["processed_2005_depth_m"]; kn=t["KN192_2008_depth_m"]; kx=t["KNOX15RR_2008_depth_m"]
    if rawmed is None:
        verdict="NO_NON_NULL_RAW_DEPTH_IN_GEOMETRY_SELECTED_BEAMS"
    else:
        dp=abs(rawmed-proc); dkn=abs(rawmed-kn); dkx=abs(rawmed-kx)
        if dkn < dp and dkx < dp:
            verdict="PROCESSING_SHIFT_SUPPORTED"
        elif dp < dkn and dp < dkx:
            verdict="RAW_ACQUISITION_SUPPORTS_PROCESSED"
        else:
            verdict="MIXED_OR_AMBIGUOUS"
    nearest_nonnull=next((x for x in beams if x["raw_depth_m"] is not None),None)
    return {
      "target":t,
      "geometry_selected_nonnull_beam_count":len(beams),
      "beams":beams,
      "raw_depth_median_m":rawmed,
      "raw_depth_min_m":min(depths) if depths else None,
      "raw_depth_max_m":max(depths) if depths else None,
      "by_raw_file":groups,
      "nearest_geometry_nonnull_beam":nearest_nonnull,
      "comparisons":None if rawmed is None else {
        "processed_minus_nearest_raw_depth_m":proc-nearest_nonnull["raw_depth_m"] if nearest_nonnull else None,
        "processed_minus_raw_median_depth_m":proc-rawmed,
        "KN192_minus_raw_median_depth_m":kn-rawmed,
        "KNOX15RR_minus_raw_median_depth_m":kx-rawmed,
        "abs_raw_median_to_processed_m":abs(rawmed-proc),
        "abs_raw_median_to_KN192_m":abs(rawmed-kn),
        "abs_raw_median_to_KNOX15RR_m":abs(rawmed-kx)
      },
      "diagnostic_verdict":verdict
    }

summaries={k:summarize_cell(v) for k,v in state.items()}
out={
 "artifact_id":"JANUS-KUSTO-CAND003-DOMINANT-CELLS-RAW-DEPTH-STAGEB-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CAND003-DOMINANT-CELLS-RAW-DEPTH-STAGEB-PREREG-2026-09-23-v1.0.json",
 "source":{"path":ZIP_PATH,"sha256":ZIP_SHA,"entries":66},
 "source_authority":{
   "raw_reader":"MB-System MBF_EMOLDRAW",
   "depth_semantics":"EM12 bath_res=1 => 0.1 m/unit; bath_res=2 => 0.2 m/unit",
   "null_semantics":"stored bath == 0 => MB_FLAG_NULL; nonzero => MB_FLAG_NONE",
   "quality_byte":"reported but not filtered"
 },
 "selection":{
   "geometry_radius_m":100,
   "depth_used_for_selection":False,
   "quality_used_for_selection":False,
   "scan_all_66_files":True
 },
 "global_counts":{"bath_records":total_bath,"beam_geometries":total_beams,"bath_records_without_prior_nav":navless_bath},
 "cells":summaries,
 "strict_parent_temporal_fail":"IMMUTABLE",
 "claim_ceiling":"RAW_EM12_DEPTH_PROVENANCE_DIAGNOSTIC_FOR_TWO_FROZEN_CELLS"
}
p=OUT/"JANUS-KUSTO-CAND003-DOMINANT-CELLS-RAW-DEPTH-STAGEB-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
  "global_counts":out["global_counts"],
  "cells":{k:{
    "geometry_selected_nonnull_beam_count":v["geometry_selected_nonnull_beam_count"],
    "raw_depth_median_m":v["raw_depth_median_m"],
    "raw_depth_min_m":v["raw_depth_min_m"],
    "raw_depth_max_m":v["raw_depth_max_m"],
    "by_raw_file":v["by_raw_file"],
    "nearest_geometry_nonnull_beam":v["nearest_geometry_nonnull_beam"],
    "comparisons":v["comparisons"],
    "diagnostic_verdict":v["diagnostic_verdict"]
  } for k,v in summaries.items()}
},indent=2))
