#!/usr/bin/env python3
from __future__ import annotations
import bisect, concurrent.futures, ftplib, hashlib, json, math, re, struct, time, zipfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin
import numpy as np
import requests
from scipy.stats import spearmanr

C={"lat":-4.015075679897318,"lon":-12.29915403590432}
R_EARTH=6371008.8
ANALYSIS_R=5000.0
FETCH_R=6500.0
CELL=100.0
LAT0=C["lat"]; M=111320.0; COS0=math.cos(math.radians(LAT0))
HOST="livftp.noc.ac.uk"
ZIP_PATH="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11280/EM12/EM12raw.zip"
ZIP_SHA="060b65865312330e866f842959a89674b744ff44569dcb1fa13daa7b5b0cefb0"
LOCAL=Path("/tmp/CD169_EM12raw.zip")
RAW_FILES={"LINE18_ERA":"0018_240205_000335_raw.all","LINE32_ERA":"0032_270205_144041_raw.all"}
SURVEYS={
 "KN192-07_2008":"https://data.ngdc.noaa.gov/platforms/ocean/ships/knorr/KN192-07/multibeam/data/version1/MB/generated/",
 "KNOX15RR_2008":"https://data.ngdc.noaa.gov/platforms/ocean/ships/roger_revelle/KNOX15RR/multibeam/data/version1/MB/generated/"
}
UA={"User-Agent":"JANUS-KUSTO-raw-line18-line32-overlap/1.0"}
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)

# EMOLDRAW constants
POS_LABEL=0x0293; BATH_LABELS={0x0296:"EM12S_BATH",0x0295:"EM12DP_BATH",0x0294:"EM12DS_BATH"}
POS_SIZE=93; BATH_SIZE=926

def get_zip():
    h=hashlib.sha256()
    f=ftplib.FTP(timeout=180); f.connect(HOST,21); f.login("anonymous","janus-kusto@example.invalid"); f.voidcmd("TYPE I")
    with LOCAL.open("wb") as o:
        def cb(b):o.write(b);h.update(b)
        f.retrbinary("RETR "+ZIP_PATH,cb,blocksize=1024*1024)
    try:f.quit()
    except:f.close()
    if h.hexdigest()!=ZIP_SHA: raise RuntimeError("ZIP SHA mismatch")

def digits(b,a,z):
    x=b[a:z]
    if len(x)!=(z-a) or not all(48<=q<=57 for q in x):return None
    return int(x.decode("ascii"))

def parse_time(payload,pos_mode=False):
    try:
        dd=digits(payload,0,2);mm=digits(payload,2,4);yy=digits(payload,4,6)
        if pos_mode:
            hh=digits(payload,7,9);mi=digits(payload,9,11);ss=digits(payload,11,13);cs=digits(payload,13,15)
        else:
            hh=digits(payload,6,8);mi=digits(payload,8,10);ss=digits(payload,10,12);cs=digits(payload,12,14)
        if None in (dd,mm,yy,hh,mi,ss,cs):return None
        year=2000+yy if yy<80 else 1900+yy
        dt=datetime(year,mm,dd,hh,mi,ss,cs*10000,tzinfo=timezone.utc)
        return dt.timestamp(),dt.isoformat().replace("+00:00","Z")
    except:return None

def af(b,a,z):
    try:return float(b[a:z].decode("ascii","ignore").strip())
    except:return None

def parse_pos(payload,offset):
    tt=parse_time(payload,True)
    if tt is None:return None
    try:
        deg=digits(payload,16,18);minute=af(payload,18,25);hemi=chr(payload[25])
        ldeg=digits(payload,27,30);lmin=af(payload,30,37);lhemi=chr(payload[37])
        if None in (deg,minute,ldeg,lmin):return None
        lat=deg+minute/60.0
        if hemi.lower()=="s":lat=-lat
        lon=ldeg+lmin/60.0
        if lhemi.lower()=="w":lon=-lon
        speed=af(payload,80,84) or 0.0
        head=af(payload,85,90) or 0.0
        return {"offset":offset,"epoch":tt[0],"utc":tt[1],"lat":lat,"lon":lon,"speed_mps":speed,"line_heading_deg":head}
    except:return None

def parse_bath(payload,offset,tname):
    tt=parse_time(payload,False)
    if tt is None or len(payload)<BATH_SIZE:return None
    try:
        ping=struct.unpack_from("<H",payload,14)[0];res=payload[16]
        heading=.1*struct.unpack_from("<H",payload,20)[0]
        if res==1:ds=.1;xs=.2
        elif res==2:ds=.2;xs=.5
        else:return None
        beams=[]
        for i in range(81):
            k=32+11*i
            bath=struct.unpack_from("<H",payload,k)[0]
            if bath==0:continue
            across=struct.unpack_from("<h",payload,k+2)[0]*xs
            along=struct.unpack_from("<h",payload,k+4)[0]*xs
            beams.append((i,ds*bath,across,along))
        return {"offset":offset,"epoch":tt[0],"utc":tt[1],"type":tname,"ping_number":ping,
                "bath_res":res,"heading_deg":heading,"beams":beams}
    except:return None

def coor_scale(lat):
    C1=111412.84;C2=-93.5;C3=.118;C4=111132.92;C5=-559.82;C6=1.175;C7=.0023
    r=math.radians(lat)
    return (1.0/abs(C1*math.cos(r)+C2*math.cos(3*r)+C3*math.cos(5*r)),
            1.0/abs(C4+C5*math.cos(2*r)+C6*math.cos(4*r)+C7*math.cos(6*r)))

def nav_interp(fixes,t,speed,heading):
    if not fixes:return None
    times=[x["epoch"] for x in fixes]
    mlon,mlat=coor_scale(fixes[-1]["lat"])
    if len(fixes)>1 and fixes[0]["epoch"]<=t<=fixes[-1]["epoch"]:
        j=bisect.bisect_left(times,t);j=max(1,j)
        a,b=fixes[j-1],fixes[j];den=b["epoch"]-a["epoch"];fac=0 if den==0 else (t-a["epoch"])/den
        return {"lon":a["lon"]+fac*(b["lon"]-a["lon"]),"lat":a["lat"]+fac*(b["lat"]-a["lat"])}
    anchor=fixes[-1] if t>fixes[-1]["epoch"] else fixes[0]
    dt=t-anchor["epoch"];dd=dt*speed;hx=math.sin(math.radians(heading));hy=math.cos(math.radians(heading))
    return {"lon":anchor["lon"]+hx*mlon*dd,"lat":anchor["lat"]+hy*mlat*dd}

def beam_lonlat(nav,heading,across,along):
    mlon,mlat=coor_scale(nav["lat"]);hx=math.sin(math.radians(heading));hy=math.cos(math.radians(heading))
    return (nav["lon"]+hy*mlon*across+hx*mlon*along,
            nav["lat"]-hx*mlat*across+hy*mlat*along)

def hav(lat1,lon1,lat2,lon2):
    p1,p2=math.radians(lat1),math.radians(lat2);dp=p2-p1;dl=math.radians(lon2-lon1)
    a=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R_EARTH*math.asin(min(1,math.sqrt(a)))

def xy(lon,lat):
    return ((lon-C["lon"])*M*COS0,(lat-C["lat"])*M)

def scan_raw_file(data):
    records=[];n=len(data)
    for i in range(n-18):
        if data[i]!=0x02:continue
        typ=(data[i]<<8)|data[i+1]
        if typ==POS_LABEL and i+2+POS_SIZE<=n:
            p=parse_pos(data[i+2:i+2+POS_SIZE],i)
            if p:records.append((i,"POS",p))
        elif typ in BATH_LABELS and i+2+BATH_SIZE<=n:
            b=parse_bath(data[i+2:i+2+BATH_SIZE],i,BATH_LABELS[typ])
            if b:records.append((i,"BATH",b))
    u={(r[0],r[1]):r for r in records}
    records=sorted(u.values(),key=lambda r:r[0])
    fixes=[];speed=0.;linehead=0.;pts=[]
    for _,kind,o in records:
        if kind=="POS":
            fixes.append(o);speed=o["speed_mps"];linehead=o["line_heading_deg"]
        else:
            nav=nav_interp(fixes,o["epoch"],speed,linehead)
            if nav is None:continue
            for bi,depth,across,along in o["beams"]:
                lon,lat=beam_lonlat(nav,o["heading_deg"],across,along)
                xx,yy=xy(lon,lat)
                if math.hypot(xx,yy)<=ANALYSIS_R:
                    pts.append((xx,yy,depth,across,bi,o["epoch"],lon,lat))
    return np.asarray(pts,float)

def grid_raw(pts):
    d={}
    for row in pts:
        x,y,z,across=row[:4]
        ix=math.floor((x+ANALYSIS_R)/CELL);iy=math.floor((y+ANALYSIS_R)/CELL)
        xc=-ANALYSIS_R+(ix+.5)*CELL;yc=-ANALYSIS_R+(iy+.5)*CELL
        if math.hypot(xc,yc)>ANALYSIS_R:continue
        d.setdefault((ix,iy),[]).append((z,across))
    out={}
    for k,v in d.items():
        z=np.asarray([q[0] for q in v]);ac=np.asarray([q[1] for q in v])
        ix,iy=k;xc=-ANALYSIS_R+(ix+.5)*CELL;yc=-ANALYSIS_R+(iy+.5)*CELL
        out[k]={"x":xc,"y":yc,"z":float(np.median(z)),"n":len(v),"abs_across_m":float(np.median(np.abs(ac))),"signed_across_m":float(np.median(ac))}
    return out

def fit(x,y,z):
    A=np.column_stack([np.ones(len(x)),x,y]);coef=np.linalg.lstsq(A,z,rcond=None)[0];res=z-A@coef
    return coef,res

def mad(v):
    v=np.asarray(v);m=np.median(v);return float(np.median(np.abs(v-m)))

# FBT
def fetch(url,timeout=120):
    last=None
    for attempt in range(7):
        try:
            r=requests.get(url,headers=UA,timeout=timeout)
            if r.status_code==429:
                time.sleep(min(2**attempt,32));continue
            r.raise_for_status();return r.content
        except Exception as e:
            last=e;time.sleep(min(2**attempt,32))
    raise last

def list_files(base,ext):
    html=fetch(base,60).decode("latin1","replace")
    return sorted(set(h for h in re.findall(r'href="([^"]+)"',html,re.I) if h.lower().endswith("."+ext.lower())))

def parse_fnv(text):
    rows=[]
    for line in text.splitlines():
        p=line.split()
        if len(p)<19:continue
        try:rows.append((float(p[15]),float(p[16]),float(p[17]),float(p[18])))
        except:pass
    return rows

def segdist(px,py,ax,ay,bx,by):
    vx,vy=bx-ax,by-ay;wx,wy=px-ax,py-ay;den=vx*vx+vy*vy;t=0 if den==0 else max(0,min(1,(wx*vx+wy*vy)/den))
    return math.hypot(px-(ax+t*vx),py-(ay+t*vy))

def min_cross(rows):
    best=float("inf")
    for plon,plat,slon,slat in rows:
        ax,ay=xy(plon,plat);bx,by=xy(slon,slat);best=min(best,segdist(0,0,ax,ay,bx,by))
    return best

def parse_fbt(data):
    off=0;chunks=[];ids={b"V4":(4,90),b"V5":(5,98),b"cc":("c",2),b"##":("c",30)}
    while off+2<=len(data):
        tag=data[off:off+2]
        if tag not in ids:raise RuntimeError(f"FBT bad tag {tag!r} @{off}")
        ver,hs=ids[tag]
        if isinstance(ver,str):off+=hs+128;continue
        h=data[off:off+hs]
        t,lon,lat,sd,alt=struct.unpack_from(">5d",h,2)
        heading=struct.unpack_from(">f",h,42)[0]
        if ver==4:
            nb,na,ns,head=struct.unpack_from(">4h",h,70);ds,xs=struct.unpack_from(">2f",h,78)
        else:
            nb,na,ns,head=struct.unpack_from(">4i",h,70);ds,xs=struct.unpack_from(">2f",h,86)
        p=off+hs;flags=np.frombuffer(data,dtype=np.uint8,count=nb,offset=p).copy();p+=nb
        def arr(n):
            nonlocal p
            a=np.frombuffer(data,dtype=">i2",count=n,offset=p).astype(float);p+=2*n;return a
        bath=arr(nb);across=arr(nb);along=arr(nb);p+=2*na+6*ns
        navx,navy=xy(lon,lat);hr=math.radians(heading);sh=math.sin(hr);ch=math.cos(hr)
        xt=xs*across;lt=xs*along
        xx=navx+lt*sh+xt*ch;yy=navy+lt*ch-xt*sh;zz=ds*bath+sd
        good=np.where(flags==0)[0]
        if len(good):chunks.append(np.column_stack([xx[good],yy[good],zz[good]]))
        off=p
    return np.vstack(chunks) if chunks else np.empty((0,3))

def load_2008(label,base):
    fnvs=list_files(base,"fnv");selected=[];manifest=[]
    def one(n):
        b=fetch(urljoin(base,n),90);rows=parse_fnv(b.decode("utf-8","replace"));return n,min_cross(rows)
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:
        fs=[ex.submit(one,n) for n in fnvs]
        for f in concurrent.futures.as_completed(fs):
            n,d=f.result()
            if d<=FETCH_R:selected.append(n[:-4]+".fbt")
    selected=sorted(set(selected))
    parts=[]
    for n in selected:
        b=fetch(urljoin(base,n),180);p=parse_fbt(b)
        if len(p):
            rr=np.hypot(p[:,0],p[:,1]);parts.append(p[rr<=ANALYSIS_R])
        manifest.append(n)
    pts=np.vstack(parts) if parts else np.empty((0,3))
    return pts,manifest

def grid_xyz(pts):
    d={}
    for x,y,z in pts:
        ix=math.floor((x+ANALYSIS_R)/CELL);iy=math.floor((y+ANALYSIS_R)/CELL)
        xc=-ANALYSIS_R+(ix+.5)*CELL;yc=-ANALYSIS_R+(iy+.5)*CELL
        if math.hypot(xc,yc)>ANALYSIS_R:continue
        d.setdefault((ix,iy),[]).append(float(z))
    return {k:{"x":-ANALYSIS_R+(k[0]+.5)*CELL,"y":-ANALYSIS_R+(k[1]+.5)*CELL,"z":float(np.median(v)),"n":len(v)} for k,v in d.items()}

def pair_stats(ga,gb):
    common=sorted(set(ga)&set(gb));out={"common_cells":len(common)}
    if len(common)<30:return {**out,"support_pass":False}
    x=np.array([ga[k]["x"] for k in common]);y=np.array([ga[k]["y"] for k in common]);a=np.array([ga[k]["z"] for k in common]);b=np.array([gb[k]["z"] for k in common])
    _,ra=fit(x,y,a);_,rb=fit(x,y,b)
    diff=b-a
    out.update({
      "support_pass":True,
      "median_signed_B_minus_A_m":float(np.median(diff)),
      "median_absolute_residual_m":float(np.median(np.abs(diff))),
      "pearson_after_plane":float(np.corrcoef(ra,rb)[0,1]),
      "spearman_after_plane":float(spearmanr(ra,rb).statistic)
    })
    return out

get_zip()
with zipfile.ZipFile(LOCAL,"r") as z:
    raw={}
    for lab,name in RAW_FILES.items():
        raw[lab]=scan_raw_file(z.read(name))
        print(lab,"raw points",len(raw[lab]))
gr={k:grid_raw(v) for k,v in raw.items()}
common=sorted(set(gr["LINE18_ERA"])&set(gr["LINE32_ERA"]))
if len(common)>=30:
    x=np.array([gr["LINE18_ERA"][k]["x"] for k in common]);y=np.array([gr["LINE18_ERA"][k]["y"] for k in common])
    z18=np.array([gr["LINE18_ERA"][k]["z"] for k in common]);z32=np.array([gr["LINE32_ERA"][k]["z"] for k in common])
    d=z32-z18;_,r18=fit(x,y,z18);_,r32=fit(x,y,z32)
    ac18=np.array([gr["LINE18_ERA"][k]["abs_across_m"] for k in common]);ac32=np.array([gr["LINE32_ERA"][k]["abs_across_m"] for k in common])
    overlap={
      "support_pass":True,"common_cell_count":len(common),
      "median_line32_minus_line18_depth_m":float(np.median(d)),
      "p10_p90_line32_minus_line18_depth_m":[float(np.percentile(d,10)),float(np.percentile(d,90))],
      "MAD_depth_difference_m":mad(d),
      "pearson_after_independent_plane":float(np.corrcoef(r18,r32)[0,1]),
      "spearman_after_independent_plane":float(spearmanr(r18,r32).statistic),
      "fraction_line32_deeper_than_line18":float(np.mean(d>0)),
      "line18_median_abs_across_track_m":float(np.median(ac18)),
      "line32_median_abs_across_track_m":float(np.median(ac32)),
      "cells":[
        {"ix":k[0],"iy":k[1],"x":gr["LINE18_ERA"][k]["x"],"y":gr["LINE18_ERA"][k]["y"],
         "line18_depth_m":gr["LINE18_ERA"][k]["z"],"line32_depth_m":gr["LINE32_ERA"][k]["z"],
         "line32_minus_line18_m":gr["LINE32_ERA"][k]["z"]-gr["LINE18_ERA"][k]["z"],
         "line18_abs_across_m":gr["LINE18_ERA"][k]["abs_across_m"],"line32_abs_across_m":gr["LINE32_ERA"][k]["abs_across_m"]}
        for k in common
      ]
    }
else:overlap={"support_pass":False,"common_cell_count":len(common)}

surv={};g2008={}
for lab,base in SURVEYS.items():
    print("loading",lab)
    pts,manifest=load_2008(lab,base);surv[lab]={"points":len(pts),"selected_fbt":manifest};g2008[lab]=grid_xyz(pts)

tie={}
for rawlab in ["LINE18_ERA","LINE32_ERA"]:
    tie[rawlab]={}
    for lab in SURVEYS:
        tie[rawlab][lab]=pair_stats(gr[rawlab],g2008[lab])

verdict="MIXED_OR_INSUFFICIENT"
if overlap.get("support_pass"):
    med=overlap["median_line32_minus_line18_depth_m"];frac=overlap["fraction_line32_deeper_than_line18"]
    l18closer=all(tie["LINE18_ERA"][s].get("support_pass") and tie["LINE32_ERA"][s].get("support_pass") and tie["LINE18_ERA"][s]["median_absolute_residual_m"] < tie["LINE32_ERA"][s]["median_absolute_residual_m"] for s in SURVEYS)
    l32closer=all(tie["LINE18_ERA"][s].get("support_pass") and tie["LINE32_ERA"][s].get("support_pass") and tie["LINE32_ERA"][s]["median_absolute_residual_m"] < tie["LINE18_ERA"][s]["median_absolute_residual_m"] for s in SURVEYS)
    if med>=100 and frac>=.80 and l18closer: verdict="LINE32_SYSTEMATIC_DEEP_BIAS_SUPPORTED"
    elif med>=100 and frac>=.80 and l32closer: verdict="RAW_LINEAGE_SYSTEMATIC_DISCORDANCE_SUPPORTED"
    elif med<=-100 and frac<=.20 and l32closer: verdict="LINE18_SYSTEMATIC_SHALLOW_BIAS_SUPPORTED"
    elif abs(med)>=100: verdict="RAW_LINEAGE_SYSTEMATIC_DISCORDANCE_SUPPORTED"
    elif abs(med)<50: verdict="LOCAL_CELL_ONLY"

out={
 "artifact_id":"JANUS-KUSTO-CAND003-RAW-LINE18-VS-LINE32-SPATIAL-OVERLAP-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CAND003-RAW-LINE18-VS-LINE32-SPATIAL-OVERLAP-PREREG-2026-09-23-v1.0.json",
 "target_center":C,"analysis_radius_m":ANALYSIS_R,"cell_m":CELL,
 "raw_point_counts":{k:int(len(v)) for k,v in raw.items()},
 "raw_grid_cell_counts":{k:len(v) for k,v in gr.items()},
 "raw_overlap":overlap,
 "independent_2008":surv,
 "raw_vs_2008":tie,
 "authoritative_verdict":verdict,
 "strict_temporal_fail":"IMMUTABLE",
 "claim_ceiling":"RAW_EM12_LINEAGE_LEVEL_SPATIAL_SYSTEMATIC_DIAGNOSTIC_ONLY"
}
p=OUT/"JANUS-KUSTO-CAND003-RAW-LINE18-VS-LINE32-SPATIAL-OVERLAP-RUN-2026-09-23-v1.0.json";p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "raw_point_counts":out["raw_point_counts"],"raw_grid_cell_counts":out["raw_grid_cell_counts"],
 "raw_overlap":{k:v for k,v in overlap.items() if k!="cells"},
 "independent_2008":surv,
 "raw_vs_2008":tie,
 "authoritative_verdict":verdict
},indent=2))
