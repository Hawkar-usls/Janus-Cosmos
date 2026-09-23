#!/usr/bin/env python3
from __future__ import annotations

import concurrent.futures
import json
import math
import re
import struct
import time
from collections import defaultdict
from pathlib import Path
from urllib.parse import urljoin, urlparse

import numpy as np
import requests
from scipy import ndimage
from shapely.geometry import LineString, box

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHARK-ROCK-POSTHOC-G2B-ELIGIBILITY-DIAGNOSTIC-PREREG-2026-09-23-v1.0.json").read_text())
G3=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHARK-ROCK-G3-TRUTH-UNLOCK-SCORING-RECEIPT-2026-09-23-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-SharkRock-eligibility-diagnostic/1.0"}
LAYER="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
T=PRE["frozen_g2b_geometry"]["tile"]; TILE=box(T["lon_min"],T["lat_min"],T["lon_max"],T["lat_max"])
GRID=float(PRE["frozen_g2b_geometry"]["grid_m"])
STRIDE=float(PRE["frozen_g2b_geometry"]["sampling_stride_m"])
RADII=[int(x) for x in PRE["frozen_g2b_geometry"]["radii_m"]]
TRUTH=PRE["truth_coordinate"]
UNC=17.67766952966369
M_PER_DEG=111320.0
LON0=(T["lon_min"]+T["lon_max"])/2.0
LAT0=(T["lat_min"]+T["lat_max"])/2.0
COS0=math.cos(math.radians(LAT0))
LOCAL_HALF_M=800.0

def get(url,params=None):
    last=None
    for a in range(8):
        try:
            r=requests.get(url,params=params,headers=UA,timeout=180)
            if r.status_code==429:
                time.sleep(min(45,2**a)); continue
            r.raise_for_status(); return r
        except Exception as e:
            last=e
            if a==7: raise
            time.sleep(min(45,2**a))
    raise last

def hrefs(u):
    return sorted(set(urljoin(u,h) for h in re.findall(r'href=["\']([^"\']+)["\']',get(u).text,re.I)))

def survey_features(sid):
    d=get(LAYER,params={
      "where":f"SURVEY_ID='{sid}'",
      "outFields":"SURVEY_ID,PLATFORM,INSTRUMENT,SURVEY_YEAR,DOWNLOAD_URL",
      "returnGeometry":"false","f":"json"
    }).json()
    return [f.get("attributes") or {} for f in d.get("features",[])]

def derive_base(download_url):
    p=urlparse(download_url); parts=[x for x in p.path.split("/") if x]
    ship=parts[-2]; survey=parts[-1]
    if survey.endswith("_mb.html"): survey=survey[:-8]
    else: survey=survey.split(".")[0]
    return f"https://data.ngdc.noaa.gov/platforms/ocean/ships/{ship}/{survey}/multibeam/data/"

def list_fnv(base,maxdepth=4):
    seen=set(); out=[]
    def walk(u,d):
        if u in seen or d>maxdepth:return
        seen.add(u)
        for x in hrefs(u):
            if not x.startswith(base) or x.rstrip("/")==u.rstrip("/"):continue
            if x.endswith("/"):walk(x,d+1)
            elif x.lower().endswith(".fnv"):out.append(x)
    walk(base,0)
    return sorted(set(out))

def parse_fnv(text):
    rows=[]
    for line in text.splitlines():
        p=line.split()
        if len(p)<19:continue
        try:rows.append((float(p[15]),float(p[16]),float(p[17]),float(p[18])))
        except:pass
    return rows

def fnv_hit(u):
    try:
        for plon,plat,slon,slat in parse_fnv(get(u).text):
            if LineString([(plon,plat),(slon,slat)]).intersects(TILE):return u
    except Exception:
        return None
    return None

def fbt_url(u):return u[:-4]+".fbt"

def xy(lon,lat):
    return (lon-LON0)*M_PER_DEG*COS0,(lat-LAT0)*M_PER_DEG

XMIN,YMIN=xy(T["lon_min"],T["lat_min"])
XMAX,YMAX=xy(T["lon_max"],T["lat_max"])
NX=int(math.ceil((XMAX-XMIN)/GRID))
NY=int(math.ceil((YMAX-YMIN)/GRID))
TX,TY=xy(TRUTH["lon"],TRUTH["lat"])
TC=int(math.floor((TX-XMIN)/GRID))
TR=int(math.floor((TY-YMIN)/GRID))
MARGIN_C=int(math.ceil(LOCAL_HALF_M/GRID))
C0=max(0,TC-MARGIN_C); C1=min(NX-1,TC+MARGIN_C)
R0=max(0,TR-MARGIN_C); R1=min(NY-1,TR+MARGIN_C)
LW=C1-C0+1; LH=R1-R0+1

def parse_fbt_local(data,counts):
    off=0; min_d=float("inf"); used=0; recs=0
    ids={b"V4":(4,90),b"V5":(5,98),b"cc":("comment",2),b"##":("oldcomment",30)}
    while off+2<=len(data):
        tag=data[off:off+2]
        if tag not in ids:break
        ver,hs=ids[tag]
        if isinstance(ver,str):
            off+=hs+128;continue
        h=data[off:off+hs]
        if len(h)<hs:break
        t,lon,lat,sd,alt=struct.unpack_from(">5d",h,2)
        heading,speed,roll,pitch,heave,bx,bl=struct.unpack_from(">7f",h,42)
        if ver==4:
            nb,na,ns,shd=struct.unpack_from(">4h",h,70);dscale,xscale=struct.unpack_from(">2f",h,78)
        else:
            nb,na,ns,shd=struct.unpack_from(">4i",h,70);dscale,xscale=struct.unpack_from(">2f",h,86)
        if min(nb,na,ns)<0 or nb>100000:break
        p=off+hs
        flags=np.frombuffer(data,dtype=np.uint8,count=nb,offset=p);p+=nb
        bath=np.frombuffer(data,dtype=">i2",count=nb,offset=p);p+=2*nb
        across=np.frombuffer(data,dtype=">i2",count=nb,offset=p).astype(np.float64);p+=2*nb
        along=np.frombuffer(data,dtype=">i2",count=nb,offset=p).astype(np.float64);p+=2*nb
        p+=2*na+6*ns
        if p>len(data):break
        gi=np.where(flags==0)[0]
        if len(gi):
            navx,navy=xy(lon,lat)
            hr=math.radians(heading);sh=math.sin(hr);ch=math.cos(hr)
            xt=xscale*across[gi];lt=xscale*along[gi]
            east=navx+lt*sh+xt*ch
            north=navy+lt*ch-xt*sh
            d=np.hypot(east-TX,north-TY)
            if len(d):min_d=min(min_d,float(np.min(d)))
            col=np.floor((east-XMIN)/GRID).astype(np.int64)
            row=np.floor((north-YMIN)/GRID).astype(np.int64)
            m=(row>=R0)&(row<=R1)&(col>=C0)&(col<=C1)
            if np.any(m):
                rr=row[m]-R0;cc=col[m]-C0
                flat=rr*LW+cc
                np.add.at(counts.ravel(),flat,1)
                used+=int(np.sum(m))
        recs+=1;off=p
    return min_d,used,recs

g3_support={sid:bool(v["truth_supported"]) for sid,v in G3["per_survey"].items()}
survey_out={}

for sid in sorted(G3["per_survey"]):
    feats=survey_features(sid)
    bases=sorted(set(derive_base(f["DOWNLOAD_URL"]) for f in feats if f.get("DOWNLOAD_URL")))
    fnvs=[]
    for b in bases:fnvs.extend(list_fnv(b))
    fnvs=sorted(set(fnvs))
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        hits=[x for x in ex.map(fnv_hit,fnvs) if x]
    counts=np.zeros((LH,LW),dtype=np.int64)
    min_good=float("inf");file_errors=0;local_beams=0
    for k,u in enumerate(sorted(hits),1):
        try:
            b=get(fbt_url(u)).content
            md,used,recs=parse_fbt_local(b,counts)
            min_good=min(min_good,md);local_beams+=used
        except Exception:
            file_errors+=1
        if k%25==0:print(sid,k,"/",len(hits),flush=True)
    valid=counts>0
    dist=ndimage.distance_transform_edt(valid)
    truth_lr=TR-R0;truth_lc=TC-C0
    truth_cell_occupied=bool(0<=truth_lr<LH and 0<=truth_lc<LW and valid[truth_lr,truth_lc])
    per_radius={}
    any_chance=False
    step=max(1,int(round(STRIDE/GRID)))
    for rad in RADII:
        rpx=max(1,int(round(rad/GRID)))
        global_rows=np.arange(rpx,NY-rpx,step,dtype=int)
        global_cols=np.arange(rpx,NX-rpx,step,dtype=int)
        rows=global_rows[(global_rows>=R0)&(global_rows<=R1)]
        cols=global_cols[(global_cols>=C0)&(global_cols<=C1)]
        nearest_sample=float("inf");nearest_eligible=float("inf");eligible_count=0
        for rr in rows:
            for cc in cols:
                cx=XMIN+(cc+.5)*GRID;cy=YMIN+(rr+.5)*GRID
                dd=math.hypot(cx-TX,cy-TY)
                nearest_sample=min(nearest_sample,dd)
                lr=rr-R0;lc=cc-C0
                if valid[lr,lc] and dist[lr,lc]>=rpx:
                    eligible_count+=1
                    nearest_eligible=min(nearest_eligible,dd)
        centers_y=(np.arange(R0,R1+1)+.5)*GRID+YMIN
        centers_x=(np.arange(C0,C1+1)+.5)*GRID+XMIN
        YY,XX=np.meshgrid(centers_y,centers_x,indexing="ij")
        mask=np.hypot(XX-TX,YY-TY)<=rad
        occupied_within=int(np.sum(valid & mask))
        threshold=rad+UNC
        opportunity=bool(math.isfinite(nearest_eligible) and nearest_eligible<=threshold)
        any_chance=any_chance or opportunity
        per_radius[str(rad)]={
          "radius_m":rad,
          "occupied_grid_cells_within_radius":occupied_within,
          "nearest_G2B_sample_center_distance_m":None if not math.isfinite(nearest_sample) else nearest_sample,
          "nearest_eligible_G2B_sample_center_distance_m":None if not math.isfinite(nearest_eligible) else nearest_eligible,
          "eligible_sample_center_count_in_local_window":eligible_count,
          "truth_support_opportunity_under_frozen_geometry":opportunity,
          "support_threshold_m":threshold
        }
    if g3_support[sid]:
        cls="RECOVERED"
    elif any_chance:
        cls="DETECTOR_SELECTION_MISS"
    else:
        cls="COVERAGE_OR_VALIDITY_LIMITED"
    survey_out[sid]={
      "tile_hit_fnv_count":len(hits),
      "fbt_file_error_count":file_errors,
      "minimum_good_sounding_distance_to_truth_m":None if not math.isfinite(min_good) else min_good,
      "local_good_sounding_count":local_beams,
      "truth_25m_cell_occupied":truth_cell_occupied,
      "g3_truth_supported":g3_support[sid],
      "per_radius":per_radius,
      "diagnostic_classification":cls
    }
    print(sid,cls,"min_beam",survey_out[sid]["minimum_good_sounding_distance_to_truth_m"],flush=True)

out={
 "artifact_id":"JANUS-KUSTO-SHARK-ROCK-POSTHOC-G2B-ELIGIBILITY-DIAGNOSTIC-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "parent_g3_receipt":G3["artifact_id"],
 "truth_coordinate":TRUTH,
 "local_window_halfwidth_m":LOCAL_HALF_M,
 "surveys":survey_out,
 "classification_counts":{k:sum(1 for v in survey_out.values() if v["diagnostic_classification"]==k) for k in ["RECOVERED","DETECTOR_SELECTION_MISS","COVERAGE_OR_VALIDITY_LIMITED"]},
 "scientific_promotion":False,
 "interpretation_rule":"Diagnostic only. Any algorithm changes motivated by this result require unseen validation.",
 "claim_ceiling":PRE["claim_ceiling"]
}
p=OUT/"JANUS-KUSTO-SHARK-ROCK-POSTHOC-G2B-ELIGIBILITY-DIAGNOSTIC-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "classifications":{k:v["diagnostic_classification"] for k,v in survey_out.items()},
 "minimum_good_sounding_distance_m":{k:v["minimum_good_sounding_distance_to_truth_m"] for k,v in survey_out.items()},
 "truth_cell_occupied":{k:v["truth_25m_cell_occupied"] for k,v in survey_out.items()},
 "nearest_eligible_by_radius":{k:{r:x["nearest_eligible_G2B_sample_center_distance_m"] for r,x in v["per_radius"].items()} for k,v in survey_out.items()},
 "scientific_promotion":False
},indent=2))
