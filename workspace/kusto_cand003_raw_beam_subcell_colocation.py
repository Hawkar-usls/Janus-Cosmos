#!/usr/bin/env python3
import bisect, ftplib, hashlib, json, math, struct, zipfile
import concurrent.futures, re, time
from urllib.parse import urljoin
import numpy as np
import requests
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


C={"id":"KN19207_CAND_003","lat":-4.015075679897318,"lon":-12.29915403590432}
M=111320.0; LAT0=-4.0; COS0=math.cos(math.radians(LAT0)); SELECT_R=1600.0
UA={"User-Agent":"JANUS-KUSTO-CAND003-subcell-colocation/1.0"}
SURVEYS={
 "KN192-07_2008":"https://data.ngdc.noaa.gov/platforms/ocean/ships/knorr/KN192-07/multibeam/data/version1/MB/generated/",
 "KNOX15RR_2008":"https://data.ngdc.noaa.gov/platforms/ocean/ships/roger_revelle/KNOX15RR/multibeam/data/version1/MB/generated/"
}
LINE_FILES=["0018_240205_000335_raw.all","0032_270205_144041_raw.all"]
RADII=[25,50,100]
PRIMARY=50

def gxy(lon,lat): return lon*M*COS0,lat*M
def cxy(): return gxy(C["lon"],C["lat"])

def fetch(url,timeout=120):
    import time
    last=None
    for attempt in range(7):
        try:
            r=requests.get(url,headers=UA,timeout=timeout)
            if r.status_code==429:
                time.sleep(min(30,2**attempt)); continue
            r.raise_for_status(); return r.content
        except Exception as e:
            last=e
            if attempt==6: raise
            time.sleep(min(30,2**attempt))
    raise last

def list_files(base,ext):
    html=fetch(base,60).decode("latin1","replace")
    hrefs=re.findall(r'href="([^"]+)"',html,re.I)
    suffix='.'+ext.lower()
    return sorted(set(h for h in hrefs if h.lower().endswith(suffix)))

def parse_fnv(text):
    rows=[]
    for line in text.splitlines():
        p=line.split()
        if len(p)<19:continue
        try:rows.append({"portlon":float(p[15]),"portlat":float(p[16]),"stbdlon":float(p[17]),"stbdlat":float(p[18])})
        except:pass
    return rows

def segdist(px,py,ax,ay,bx,by):
    vx,vy=bx-ax,by-ay; wx,wy=px-ax,py-ay; den=vx*vx+vy*vy
    t=0 if den==0 else max(0,min(1,(wx*vx+wy*vy)/den))
    return math.hypot(px-(ax+t*vx),py-(ay+t*vy))

def min_fnv(rows):
    px,py=cxy();best=float('inf')
    for r in rows:
        ax,ay=gxy(r["portlon"],r["portlat"]);bx,by=gxy(r["stbdlon"],r["stbdlat"])
        best=min(best,segdist(px,py,ax,ay,bx,by))
    return best

def parse_fbt(data,name):
    off=0; chunks=[]; nrec=0
    ids={b"V4":(4,90),b"V5":(5,98),b"cc":("comment",2),b"##":("oldcomment",30)}
    while off+2<=len(data):
        tag=data[off:off+2]
        if tag not in ids:raise RuntimeError(f"{name} bad tag {tag!r} @{off}")
        ver,hs=ids[tag]
        if isinstance(ver,str):off+=hs+128;continue
        h=data[off:off+hs]
        t,lon,lat,sd,alt=struct.unpack_from(">5d",h,2)
        heading,speed,roll,pitch,heave,bxw,blw=struct.unpack_from(">7f",h,42)
        if ver==4:
            nb,na,ns,head=struct.unpack_from(">4h",h,70);dscale,xscale=struct.unpack_from(">2f",h,78)
        else:
            nb,na,ns,head=struct.unpack_from(">4i",h,70);dscale,xscale=struct.unpack_from(">2f",h,86)
        p=off+hs
        flags=np.frombuffer(data,dtype=np.uint8,count=nb,offset=p).copy();p+=nb
        def arr(n):
            nonlocal p
            a=np.frombuffer(data,dtype=">i2",count=n,offset=p).astype(float);p+=2*n;return a
        bath=arr(nb);across=arr(nb);along=arr(nb);p+=2*na+6*ns
        gx,gy=gxy(lon,lat);hr=math.radians(heading);sh=math.sin(hr);ch=math.cos(hr)
        xt=xscale*across;lt=xscale*along
        x=gx+lt*sh+xt*ch;y=gy+lt*ch-xt*sh;z=dscale*bath+sd
        good=np.where(flags==0)[0]
        if len(good):chunks.append(np.column_stack([x[good],y[good],z[good]]))
        off=p;nrec+=1
    return (np.vstack(chunks) if chunks else np.empty((0,3))),nrec

def load_2008(label,base):
    names=list_files(base,"fnv"); selected=[]; fnv_manifest=[]
    def getfnv(n):
        b=fetch(urljoin(base,n),90);return n,hashlib.sha256(b).hexdigest(),parse_fnv(b.decode("utf-8","replace"))
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
        for n,sha,rows in (f.result() for f in concurrent.futures.as_completed([ex.submit(getfnv,n) for n in names])):
            d=min_fnv(rows);fnv_manifest.append({"file":n,"sha256":sha,"min_cross_m":d})
            if d<=SELECT_R:selected.append(n[:-4]+".fbt")
    selected=sorted(set(selected));fnv_manifest.sort(key=lambda r:r["file"])
    parts=[];fbt_manifest=[]
    def getfbt(n):
        b=fetch(urljoin(base,n),180);pts,nrec=parse_fbt(b,n);return n,hashlib.sha256(b).hexdigest(),len(b),pts,nrec
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        for n,sha,size,pts,nrec in (f.result() for f in concurrent.futures.as_completed([ex.submit(getfbt,n) for n in selected])):
            parts.append(pts);fbt_manifest.append({"file":n,"sha256":sha,"bytes":size,"records":nrec,"good_beams":int(len(pts))})
    allpts=np.vstack(parts) if parts else np.empty((0,3))
    cx,cy=cxy()
    rr=np.hypot(allpts[:,0]-cx,allpts[:,1]-cy) if len(allpts) else np.array([])
    allpts=allpts[rr<=SELECT_R] if len(allpts) else allpts
    return allpts,{"selected_fbt":selected,"fbt_manifest":sorted(fbt_manifest,key=lambda r:r["file"]),"fnv_manifest":fnv_manifest}



def collect_raw_selected():
    selected={t["id"]:[] for t in TARGETS}
    with zipfile.ZipFile(LOCAL,"r") as z:
        for fn in LINE_FILES:
            data=z.read(fn)
            records=scan_records(data)
            fixes=[]; current_speed=0.0; current_line_heading=0.0
            for _,kind,obj in records:
                if kind=="POS":
                    fixes.append(obj); current_speed=obj["speed_mps"]; current_line_heading=obj["line_heading_deg"]
                    continue
                nav=nav_interp(fixes,obj["epoch"],current_speed,current_line_heading)
                if nav is None: continue
                for beam_index,across,along,quality,stored_bath,depth_m in obj["beams"]:
                    if stored_bath==0: continue
                    lon,lat=beam_lonlat(nav,obj["heading_deg"],across,along)
                    for t in TARGETS:
                        d=hav(t["lat"],t["lon"],lat,lon)
                        if d<=100:
                            selected[t["id"]].append({
                              "cell_distance_m":d,"raw_file":fn,"lineage":"LINE18_ERA" if fn.startswith("0018_") else "LINE32_ERA",
                              "ping_utc":obj["utc"],"ping_number":obj["ping_number"],"beam_index":beam_index,
                              "raw_depth_m":depth_m,"raw_across_m":across,"raw_along_m":along,
                              "beam_lon":lon,"beam_lat":lat,"raw_quality_byte":quality
                            })
    for k in selected:
        selected[k]=sorted(selected[k],key=lambda x:(x["lineage"],x["ping_utc"],x["beam_index"]))
    return selected

def estimate_at_beam(pts,beam,radius):
    if len(pts)==0:return {"n":0,"median_depth_m":None}
    bx,by=gxy(beam["beam_lon"],beam["beam_lat"])
    rr=np.hypot(pts[:,0]-bx,pts[:,1]-by)
    z=pts[rr<=radius,2]
    return {"n":int(len(z)),"median_depth_m":float(np.median(z)) if len(z) else None}

def pair_distance(a,b):
    return hav(a["beam_lat"],a["beam_lon"],b["beam_lat"],b["beam_lon"])

rawsel=collect_raw_selected()
loaded={}; manifest={}
for lab,base in SURVEYS.items():
    print("Loading",lab,flush=True)
    pts,man=load_2008(lab,base)
    loaded[lab]=pts;manifest[lab]=man
    print(lab,"points",len(pts),flush=True)

cells={}
for cid,beams in rawsel.items():
    enriched=[]
    for b in beams:
        e=dict(b); e["2008"]={}
        for lab,pts in loaded.items():
            e["2008"][lab]={str(rad):estimate_at_beam(pts,b,rad) for rad in RADII}
        enriched.append(e)
    line18=[b for b in enriched if b["lineage"]=="LINE18_ERA"]
    line32=[b for b in enriched if b["lineage"]=="LINE32_ERA"]
    closest=None
    for a in line18:
        for b in line32:
            d=pair_distance(a,b)
            if closest is None or d<closest["horizontal_separation_m"]:
                closest={"horizontal_separation_m":d,"line18":a,"line32":b,
                         "raw_line32_minus_line18_depth_m":b["raw_depth_m"]-a["raw_depth_m"]}
    if closest:
        closest["2008_depth_difference_at_own_beam_coordinates"]={}
        for lab in loaded:
            x=closest["line18"]["2008"][lab][str(PRIMARY)]["median_depth_m"]
            y=closest["line32"]["2008"][lab][str(PRIMARY)]["median_depth_m"]
            closest["2008_depth_difference_at_own_beam_coordinates"][lab]=None if x is None or y is None else y-x

    lineage_summary={}
    for lin,bb in [("LINE18_ERA",line18),("LINE32_ERA",line32)]:
        lineage_summary[lin]={}
        for lab in loaded:
            lineage_summary[lin][lab]={}
            for rad in RADII:
                residuals=[]
                for b in bb:
                    z=b["2008"][lab][str(rad)]["median_depth_m"]
                    if z is not None: residuals.append(abs(b["raw_depth_m"]-z))
                lineage_summary[lin][lab][str(rad)]={
                  "supported_beams":len(residuals),
                  "median_absolute_raw_minus_2008_m":float(np.median(residuals)) if residuals else None,
                  "min_absolute_residual_m":float(min(residuals)) if residuals else None,
                  "max_absolute_residual_m":float(max(residuals)) if residuals else None
                }
    cells[cid]={
      "raw_beam_count":len(enriched),
      "raw_beams":enriched,
      "closest_opposite_lineage_pair":closest,
      "lineage_summary":lineage_summary
    }

out={
 "artifact_id":"JANUS-KUSTO-CAND003-RAW-BEAM-SUBCELL-COLOCATION-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CAND003-RAW-BEAM-SUBCELL-COLOCATION-PREREG-2026-09-23-v1.0.json",
 "target":C,
 "raw_source":{"path":ZIP_PATH,"sha256":ZIP_SHA},
 "line_files":LINE_FILES,
 "radii_m":RADII,
 "primary_radius_m":PRIMARY,
 "cells":cells,
 "ncei_manifest":manifest,
 "claim_ceiling":"RAW_BEAM_SUBCELL_COLOCATION_DIAGNOSTIC_ONLY"
}
p=OUT/"JANUS-KUSTO-CAND003-RAW-BEAM-SUBCELL-COLOCATION-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "cells":{cid:{
   "raw_beam_count":v["raw_beam_count"],
   "closest_opposite_lineage_pair":v["closest_opposite_lineage_pair"],
   "lineage_summary":v["lineage_summary"]
 } for cid,v in cells.items()}
},indent=2))
