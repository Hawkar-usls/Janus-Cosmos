#!/usr/bin/env python3
from __future__ import annotations
import bisect, ftplib, hashlib, json, math, struct, zipfile
from datetime import datetime, timezone
from pathlib import Path
import numpy as np

C={"lat":-4.015075679897318,"lon":-12.29915403590432}
FIXED=[
 {"id":"CELL_A","lat":-4.0165853,"lon":-12.2985657},
 {"id":"CELL_B","lat":-4.0152462,"lon":-12.299318}
]
ANALYSIS_R=5000.0
CELL=100.0
LAT0=C["lat"]; M=111320.0; COS0=math.cos(math.radians(LAT0))
HOST="livftp.noc.ac.uk"
ZIP_PATH="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11280/EM12/EM12raw.zip"
ZIP_SHA="060b65865312330e866f842959a89674b744ff44569dcb1fa13daa7b5b0cefb0"
LOCAL=Path("/tmp/CD169_EM12raw.zip")
RAW_FILES={"LINE18_ERA":"0018_240205_000335_raw.all","LINE32_ERA":"0032_270205_144041_raw.all"}
POS_LABEL=0x0293
BATH_LABELS={0x0294:"EM12DS_BATH",0x0295:"EM12DP_BATH",0x0296:"EM12S_BATH"}
POS_SIZE=93; BATH_SIZE=926
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)

def get_zip():
    h=hashlib.sha256()
    f=ftplib.FTP(timeout=180); f.connect(HOST,21); f.login("anonymous","janus-kusto@example.invalid"); f.voidcmd("TYPE I")
    with LOCAL.open("wb") as o:
        def cb(b): o.write(b); h.update(b)
        f.retrbinary("RETR "+ZIP_PATH,cb,blocksize=1024*1024)
    try:f.quit()
    except:f.close()
    if h.hexdigest()!=ZIP_SHA: raise RuntimeError("ZIP SHA mismatch")

def digits(b,a,z):
    s=b[a:z]
    if len(s)!=(z-a) or not all(48<=x<=57 for x in s): return None
    return int(s.decode("ascii"))

def parse_time(payload,pos_mode=False):
    try:
        dd=digits(payload,0,2); mm=digits(payload,2,4); yy=digits(payload,4,6)
        if pos_mode:
            hh=digits(payload,7,9); mi=digits(payload,9,11); ss=digits(payload,11,13); cs=digits(payload,13,15)
        else:
            hh=digits(payload,6,8); mi=digits(payload,8,10); ss=digits(payload,10,12); cs=digits(payload,12,14)
        if None in (dd,mm,yy,hh,mi,ss,cs): return None
        dt=datetime(2000+yy if yy<80 else 1900+yy,mm,dd,hh,mi,ss,cs*10000,tzinfo=timezone.utc)
        return dt.timestamp()
    except:return None

def af(b,a,z):
    try:return float(b[a:z].decode("ascii","ignore").strip())
    except:return None

def parse_pos(payload,offset):
    t=parse_time(payload,True)
    if t is None:return None
    try:
        deg=digits(payload,16,18); minute=af(payload,18,25); hemi=chr(payload[25])
        ldeg=digits(payload,27,30); lmin=af(payload,30,37); lhemi=chr(payload[37])
        if None in (deg,minute,ldeg,lmin):return None
        lat=deg+minute/60.0; lon=ldeg+lmin/60.0
        if hemi.lower()=="s":lat=-lat
        if lhemi.lower()=="w":lon=-lon
        speed=af(payload,80,84) or 0.0; head=af(payload,85,90) or 0.0
        return {"offset":offset,"epoch":t,"lat":lat,"lon":lon,"speed_mps":speed,"line_heading_deg":head}
    except:return None

def parse_bath(payload,offset,tname):
    t=parse_time(payload,False)
    if t is None or len(payload)<BATH_SIZE:return None
    try:
        res=payload[16]; heading=.1*struct.unpack_from("<H",payload,20)[0]
        if res==1: ds=.1; xs=.2
        elif res==2: ds=.2; xs=.5
        else:return None
        beams=[]
        for i in range(81):
            k=32+11*i; bath=struct.unpack_from("<H",payload,k)[0]
            if bath==0:continue
            across=struct.unpack_from("<h",payload,k+2)[0]*xs
            along=struct.unpack_from("<h",payload,k+4)[0]*xs
            beams.append((i,ds*bath,across,along))
        return {"offset":offset,"epoch":t,"heading_deg":heading,"beams":beams}
    except:return None

def coor_scale(lat):
    C1=111412.84; C2=-93.5; C3=.118; C4=111132.92; C5=-559.82; C6=1.175; C7=.0023
    r=math.radians(lat)
    return (1.0/abs(C1*math.cos(r)+C2*math.cos(3*r)+C3*math.cos(5*r)),
            1.0/abs(C4+C5*math.cos(2*r)+C6*math.cos(4*r)+C7*math.cos(6*r)))

def nav_interp(fixes,t,speed,heading):
    if not fixes:return None
    times=[x["epoch"] for x in fixes]; mlon,mlat=coor_scale(fixes[-1]["lat"])
    if len(fixes)>1 and fixes[0]["epoch"]<=t<=fixes[-1]["epoch"]:
        j=max(1,bisect.bisect_left(times,t)); a,b=fixes[j-1],fixes[j]
        den=b["epoch"]-a["epoch"]; fac=0 if den==0 else (t-a["epoch"])/den
        return {"lon":a["lon"]+fac*(b["lon"]-a["lon"]),"lat":a["lat"]+fac*(b["lat"]-a["lat"])}
    anchor=fixes[-1] if t>fixes[-1]["epoch"] else fixes[0]; dt=t-anchor["epoch"]; dd=dt*speed
    hx=math.sin(math.radians(heading)); hy=math.cos(math.radians(heading))
    return {"lon":anchor["lon"]+hx*mlon*dd,"lat":anchor["lat"]+hy*mlat*dd}

def beam_lonlat(nav,heading,across,along):
    mlon,mlat=coor_scale(nav["lat"]); hx=math.sin(math.radians(heading)); hy=math.cos(math.radians(heading))
    return (nav["lon"]+hy*mlon*across+hx*mlon*along,
            nav["lat"]-hx*mlat*across+hy*mlat*along)

def xy(lon,lat):
    return ((lon-C["lon"])*M*COS0,(lat-C["lat"])*M)

def scan_raw(data):
    records=[]; n=len(data)
    for i in range(n-18):
        if data[i]!=0x02:continue
        typ=(data[i]<<8)|data[i+1]
        if typ==POS_LABEL and i+2+POS_SIZE<=n:
            p=parse_pos(data[i+2:i+2+POS_SIZE],i)
            if p:records.append((i,"POS",p))
        elif typ in BATH_LABELS and i+2+BATH_SIZE<=n:
            b=parse_bath(data[i+2:i+2+BATH_SIZE],i,BATH_LABELS[typ])
            if b:records.append((i,"BATH",b))
    records=sorted({(r[0],r[1]):r for r in records}.values(),key=lambda r:r[0])
    npos=sum(1 for _,k,_ in records if k=="POS"); nbath=sum(1 for _,k,_ in records if k=="BATH")
    fixes=[]; speed=0.; linehead=0.; pts=[]
    for _,kind,o in records:
        if kind=="POS":
            fixes.append(o); speed=o["speed_mps"]; linehead=o["line_heading_deg"]
        else:
            nav=nav_interp(fixes,o["epoch"],speed,linehead)
            if nav is None:continue
            for _,depth,across,along in o["beams"]:
                lon,lat=beam_lonlat(nav,o["heading_deg"],across,along)
                x,y=xy(lon,lat)
                if math.hypot(x,y)<=ANALYSIS_R:
                    pts.append((x,y,depth,across))
    return np.asarray(pts,float),{"position_records":npos,"bath_records":nbath}

def grid(pts):
    d={}
    for x,y,z,across in pts:
        ix=math.floor((x+ANALYSIS_R)/CELL); iy=math.floor((y+ANALYSIS_R)/CELL)
        xc=-ANALYSIS_R+(ix+.5)*CELL; yc=-ANALYSIS_R+(iy+.5)*CELL
        if math.hypot(xc,yc)>ANALYSIS_R:continue
        d.setdefault((ix,iy),[]).append((float(z),float(across)))
    out={}
    for k,v in d.items():
        zz=np.array([q[0] for q in v]); aa=np.array([q[1] for q in v])
        out[k]={
          "x":-ANALYSIS_R+(k[0]+.5)*CELL,
          "y":-ANALYSIS_R+(k[1]+.5)*CELL,
          "z":float(np.median(zz)),
          "signed_across_m":float(np.median(aa)),
          "abs_across_m":float(np.median(np.abs(aa))),
          "n":len(v)
        }
    return out

def midrank_percentile(values,val):
    a=np.asarray(values,float)
    less=np.sum(a<val); equal=np.sum(a==val)
    return float((less+0.5*equal)/len(a))

def target_key(t):
    x,y=xy(t["lon"],t["lat"])
    return (math.floor((x+ANALYSIS_R)/CELL),math.floor((y+ANALYSIS_R)/CELL)),x,y

def components(keys):
    keys=set(keys); comps=[]
    while keys:
        seed=keys.pop(); stack=[seed]; n=0
        while stack:
            q=stack.pop(); n+=1
            for dx in (-1,0,1):
                for dy in (-1,0,1):
                    if dx==0 and dy==0:continue
                    nb=(q[0]+dx,q[1]+dy)
                    if nb in keys:
                        keys.remove(nb); stack.append(nb)
        comps.append(n)
    return sorted(comps,reverse=True)

get_zip()
with zipfile.ZipFile(LOCAL,"r") as z:
    raw={}; selftest={}
    expected={"LINE18_ERA":{"position_records":6867,"bath_records":533},"LINE32_ERA":{"position_records":60275,"bath_records":4784}}
    for lab,name in RAW_FILES.items():
        pts,cnt=scan_raw(z.read(name)); raw[lab]=pts; selftest[lab]=cnt
        if cnt!=expected[lab]: raise RuntimeError(f"EMOLDRAW_SELFTEST_FAIL {lab}: {cnt} != {expected[lab]}")
g={k:grid(v) for k,v in raw.items()}
common=sorted(set(g["LINE18_ERA"])&set(g["LINE32_ERA"]))
if len(common)!=2096: raise RuntimeError(f"PARENT_COMMON_CELL_COUNT_MISMATCH {len(common)} != 2096")

rows=[]
for k in common:
    a=g["LINE18_ERA"][k]; b=g["LINE32_ERA"][k]
    diff=b["z"]-a["z"]
    rows.append({
      "key":[k[0],k[1]],"x":a["x"],"y":a["y"],
      "line18_depth_m":a["z"],"line32_depth_m":b["z"],
      "signed_difference_m":diff,"abs_difference_m":abs(diff),
      "line18_signed_across_m":a["signed_across_m"],
      "line32_signed_across_m":b["signed_across_m"],
      "line18_abs_across_m":a["abs_across_m"],
      "line32_abs_across_m":b["abs_across_m"]
    })
absvals=[r["abs_difference_m"] for r in rows]
signed=[r["signed_difference_m"] for r in rows]
ge100=[tuple(r["key"]) for r in rows if r["abs_difference_m"]>=100]
ge150=[tuple(r["key"]) for r in rows if r["abs_difference_m"]>=150]
comp100=components(ge100)

fixed={}
for t in FIXED:
    k,x,y=target_key(t)
    r=next((q for q in rows if tuple(q["key"])==k),None)
    if r is None: raise RuntimeError(f"{t['id']} not in parent common cells: {k}")
    others100=[q for q in rows if q["abs_difference_m"]>=100 and tuple(q["key"])!=k]
    nearest100=min((math.hypot(q["x"]-r["x"],q["y"]-r["y"]) for q in others100),default=None)
    local={}
    for rad in [250,500,1000]:
        vals=[q["abs_difference_m"] for q in rows if math.hypot(q["x"]-r["x"],q["y"]-r["y"])<=rad]
        local[str(rad)]={
          "n":len(vals),"median_abs_difference_m":float(np.median(vals)) if vals else None,
          "p90_abs_difference_m":float(np.percentile(vals,90)) if vals else None,
          "max_abs_difference_m":max(vals) if vals else None,
          "fraction_ge100":float(np.mean(np.asarray(vals)>=100)) if vals else None
        }
    fixed[t["id"]]={
      "coordinate":{"lat":t["lat"],"lon":t["lon"]},
      "grid_key":[k[0],k[1]],"grid_center_xy_m":[r["x"],r["y"]],
      "line18_depth_m":r["line18_depth_m"],"line32_depth_m":r["line32_depth_m"],
      "signed_difference_m":r["signed_difference_m"],"abs_difference_m":r["abs_difference_m"],
      "abs_difference_midrank_percentile":midrank_percentile(absvals,r["abs_difference_m"]),
      "signed_difference_midrank_percentile":midrank_percentile(signed,r["signed_difference_m"]),
      "line18_signed_across_m":r["line18_signed_across_m"],
      "line32_signed_across_m":r["line32_signed_across_m"],
      "nearest_other_ge100_cell_distance_m":nearest100,
      "local_neighborhoods":local
    }

opp=[r for r in rows if r["line18_signed_across_m"]>=5500 and r["line32_signed_across_m"]<=-5500]
contrast=[r for r in rows if r not in opp]
def grp(rr):
    v=np.array([q["abs_difference_m"] for q in rr],float)
    return {
      "n":len(rr),
      "median_abs_difference_m":float(np.median(v)) if len(v) else None,
      "p90_abs_difference_m":float(np.percentile(v,90)) if len(v) else None,
      "fraction_abs_ge100":float(np.mean(v>=100)) if len(v) else None,
      "fraction_abs_ge150":float(np.mean(v>=150)) if len(v) else None
    }
og=grp(opp); cg=grp(contrast)

if og["n"]<20:
    verdict="INSUFFICIENT_OPPOSITE_EXTREME_SUPPORT"
elif any(fixed[x]["abs_difference_midrank_percentile"]<0.95 for x in fixed):
    verdict="NOT_EXTREME_IN_PARENT_FIELD"
elif all(fixed[x]["abs_difference_midrank_percentile"]>=0.99 for x in fixed) and og["median_abs_difference_m"]>cg["median_abs_difference_m"] and og["p90_abs_difference_m"]>cg["p90_abs_difference_m"]:
    verdict="LOCAL_OUTER_BEAM_CROSSOVER_SUPPORTED"
else:
    verdict="ISOLATED_LOCAL_OUTLIERS"

out={
 "artifact_id":"JANUS-KUSTO-CAND003-LOCAL-EXTREME-OUTER-BEAM-TAIL-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CAND003-LOCAL-EXTREME-OUTER-BEAM-TAIL-PREREG-2026-09-23-v1.0.json",
 "parser_selftest":selftest,
 "parent_common_cell_count":len(rows),
 "global_tail":{
   "abs_difference_median_m":float(np.median(absvals)),
   "abs_difference_p95_m":float(np.percentile(absvals,95)),
   "abs_difference_p99_m":float(np.percentile(absvals,99)),
   "abs_difference_max_m":max(absvals),
   "count_ge100":len(ge100),"fraction_ge100":len(ge100)/len(rows),
   "count_ge150":len(ge150),"fraction_ge150":len(ge150)/len(rows),
   "ge100_connected_component_sizes_8neighbor":comp100
 },
 "fixed_cells":fixed,
 "opposite_extreme_group":og,
 "contrast_group":cg,
 "opposite_extreme_definition":"line18_signed_across>=+5500m AND line32_signed_across<=-5500m",
 "authoritative_verdict":verdict,
 "claim_ceiling":"LOCAL_OUTER_BEAM_TAIL_DIAGNOSTIC_ONLY__NO_REFRACTION_CAUSALITY"
}
p=OUT/"JANUS-KUSTO-CAND003-LOCAL-EXTREME-OUTER-BEAM-TAIL-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
