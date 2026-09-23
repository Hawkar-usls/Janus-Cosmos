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
from scipy.spatial import cKDTree
from shapely.geometry import LineString, box

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"
OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHARK-ROCK-G2B-BLIND-MULTISURVEY-MORPHOLOGY-PREREG-2026-09-23-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-shark-rock-G2B/1.0"}
LAYER="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
T=PRE["blind_tile"]
TILE=box(T["lon_min"],T["lat_min"],T["lon_max"],T["lat_max"])
GRID=float(PRE["grid_m"])
RADII=[int(x) for x in PRE["radii_m"]]
STRIDE=float(PRE["sampling_stride_m"])
TAIL=float(PRE["selection"]["per_metric_upper_tail_quantile"])
DEDUP=float(PRE["selection"]["spatial_dedup_min_m"])
MAX_CAND=int(PRE["selection"]["maximum_candidates_per_survey"])
MATCH=float(PRE["cross_survey_persistence"]["match_tolerance_m"])
MIN_CLUSTER_SURVEYS=int(PRE["cross_survey_persistence"]["minimum_distinct_surveys_for_cluster"])
M_PER_DEG=111320.0
LON0=(T["lon_min"]+T["lon_max"])/2.0
LAT0=(T["lat_min"]+T["lat_max"])/2.0
COS0=math.cos(math.radians(LAT0))

FAMILIES=PRE["morphology_families"]

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
    d=get(LAYER,params={"where":f"SURVEY_ID='{sid}'","outFields":"SURVEY_ID,PLATFORM,INSTRUMENT,SURVEY_YEAR,DOWNLOAD_URL","returnGeometry":"false","f":"json"}).json()
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

def fbt_url(u): return u[:-4]+".fbt"

def xy(lon,lat):
    return (lon-LON0)*M_PER_DEG*COS0,(lat-LAT0)*M_PER_DEG

def ll(x,y):
    return LON0+x/(M_PER_DEG*COS0),LAT0+y/M_PER_DEG

XMIN,YMIN=xy(T["lon_min"],T["lat_min"])
XMAX,YMAX=xy(T["lon_max"],T["lat_max"])
NX=int(math.ceil((XMAX-XMIN)/GRID))
NY=int(math.ceil((YMAX-YMIN)/GRID))

def parse_fbt_into_grid(data,sum_grid,count_grid):
    off=0; recs=0; good=0; used=0
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
        bath=np.frombuffer(data,dtype=">i2",count=nb,offset=p).astype(np.float64);p+=2*nb
        across=np.frombuffer(data,dtype=">i2",count=nb,offset=p).astype(np.float64);p+=2*nb
        along=np.frombuffer(data,dtype=">i2",count=nb,offset=p).astype(np.float64);p+=2*nb
        p+=2*na+6*ns
        if p>len(data):break
        good_idx=np.where(flags==0)[0]
        good+=len(good_idx)
        if len(good_idx):
            navx,navy=xy(lon,lat)
            hr=math.radians(heading);sh=math.sin(hr);ch=math.cos(hr)
            xt=xscale*across[good_idx];lt=xscale*along[good_idx]
            east=navx+lt*sh+xt*ch
            north=navy+lt*ch-xt*sh
            depth=dscale*bath[good_idx]+sd
            col=np.floor((east-XMIN)/GRID).astype(np.int64)
            row=np.floor((north-YMIN)/GRID).astype(np.int64)
            m=(row>=0)&(row<NY)&(col>=0)&(col<NX)&np.isfinite(depth)
            if np.any(m):
                rr=row[m];cc=col[m];zz=depth[m]
                flat=rr*NX+cc
                np.add.at(sum_grid.ravel(),flat,zz)
                np.add.at(count_grid.ravel(),flat,1)
                used+=int(np.sum(m))
        recs+=1;off=p
    return recs,good,used

def robust_abs_z(v):
    med=float(np.median(v));mad=float(np.median(np.abs(v-med)))
    scale=max(1e-12,1.4826*mad)
    return np.abs((v-med)/scale),med,mad

def discover_on_grid(sid,z,valid):
    if np.sum(valid)==0:raise RuntimeError(f"{sid}: no valid rasterized cells")
    fill=float(np.median(z[valid]));zf=z.copy();zf[~valid]=fill
    if np.all(valid):
        dist=np.full(z.shape,max(RADII)/GRID+2,dtype=np.float32)
    else:
        dist=ndimage.distance_transform_edt(valid).astype(np.float32)
    step=max(1,int(round(STRIDE/GRID)))
    records=[]
    scale_summary=[]
    for rad in RADII:
        rpx=max(1,int(round(rad/GRID)))
        inner=valid & (dist>=rpx)
        rows=np.arange(rpx,NY-rpx,step,dtype=int)
        cols=np.arange(rpx,NX-rpx,step,dtype=int)
        R,C=np.meshgrid(rows,cols,indexing="ij")
        sr=R.ravel();sc=C.ravel();m=inner[sr,sc]
        sr=sr[m];sc=sc[m]
        if len(sr)<int(PRE["validity"]["minimum_blind_sample_centers_per_survey"]):
            scale_summary.append({"radius_m":rad,"sample_centers":int(len(sr)),"status":"INSUFFICIENT_VALID_CENTERS"})
            continue
        size=2*rpx+1;sigma=max(1.0,rpx/3.0)
        mean=ndimage.uniform_filter(zf,size=size,mode="nearest")
        mean2=ndimage.uniform_filter(zf*zf,size=size,mode="nearest")
        std=np.sqrt(np.maximum(0.0,mean2-mean*mean))
        rel=ndimage.maximum_filter(zf,size=size,mode="nearest")-ndimage.minimum_filter(zf,size=size,mode="nearest")
        smooth=ndimage.gaussian_filter(zf,sigma=sigma,mode="nearest",truncate=3.0)
        lap=np.abs(ndimage.laplace(smooth,mode="nearest"))/(GRID*GRID)
        g1=ndimage.gaussian_filter(zf,sigma=max(.8,rpx/4.0),mode="nearest",truncate=3.0)
        g2=ndimage.gaussian_filter(zf,sigma=max(1.2,rpx/1.5),mode="nearest",truncate=3.0)
        dog=np.abs(g1-g2)
        gy,gx=np.gradient(smooth,GRID,GRID);grad=np.hypot(gx,gy)
        ts=max(1.0,rpx/2.0)
        jxx=ndimage.gaussian_filter(gx*gx,sigma=ts,mode="nearest",truncate=3.0)
        jyy=ndimage.gaussian_filter(gy*gy,sigma=ts,mode="nearest",truncate=3.0)
        jxy=ndimage.gaussian_filter(gx*gy,sigma=ts,mode="nearest",truncate=3.0)
        anis=np.sqrt((jxx-jyy)**2+4*jxy*jxy)/(jxx+jyy+1e-12)
        persist=np.abs(zf-mean)/(std+1e-6)
        grids={
          "local_relief_m":rel,"local_std_m":std,"abs_laplacian":lap,
          "abs_center_surround_dog":dog,"structure_tensor_anisotropy":anis,
          "gradient_magnitude":grad,"isolated_extremum_persistence":persist
        }
        vals={k:g[sr,sc].astype(np.float64) for k,g in grids.items()}
        th={};rz={}
        for k,v in vals.items():
            if not np.all(np.isfinite(v)):
                good=np.isfinite(v);v=np.where(good,v,np.nanmedian(v[good]));vals[k]=v
            th[k]=float(np.quantile(v,TAIL))
            rz[k],med,mad=robust_abs_z(v)
        zstack=np.column_stack([rz[k] for k in sorted(rz)])
        generic=0.5*np.max(zstack,axis=1)+0.5*np.mean(zstack,axis=1)
        generic_q=float(np.quantile(generic,TAIL))
        nfollow=0
        for i in range(len(sr)):
            votes=[];detail={}
            for fam,keys in FAMILIES.items():
                ex=[k for k in keys if vals[k][i]>=th[k]]
                if ex:votes.append(fam);detail[fam]=ex
            if len(votes)>=2 or generic[i]>=generic_q:
                cx=XMIN+(int(sc[i])+.5)*GRID;cy=YMIN+(int(sr[i])+.5)*GRID
                lon,lat=ll(cx,cy)
                records.append({
                  "survey_id":sid,"lon":lon,"lat":lat,"x_m":cx,"y_m":cy,"radius_m":rad,
                  "family_votes":votes,"extreme_metrics_by_family":detail,
                  "generic_robust_score":float(generic[i]),
                  "metrics":{k:float(vals[k][i]) for k in vals}
                });nfollow+=1
        scale_summary.append({"radius_m":rad,"sample_centers":int(len(sr)),"followup_pre_dedup":nfollow,"generic_q995":generic_q})
    records.sort(key=lambda q:(-len(q["family_votes"]),-q["generic_robust_score"],q["radius_m"],q["lat"],q["lon"]))
    kept=[]
    for q in records:
        if any(math.hypot(q["x_m"]-p["x_m"],q["y_m"]-p["y_m"])<DEDUP for p in kept):continue
        zq=dict(q);zq["candidate_id"]=f"{sid}_C{len(kept)+1:03d}";zq["blind_rank"]=len(kept)+1
        kept.append(zq)
        if len(kept)>=MAX_CAND:break
    return kept,scale_summary,len(records)

survey_out={}
all_candidates=[]
for sid in PRE["survey_ids"]:
    feats=survey_features(sid)
    bases=sorted(set(derive_base(f["DOWNLOAD_URL"]) for f in feats if f.get("DOWNLOAD_URL")))
    fnvs=[]
    for b in bases:
        fnvs.extend(list_fnv(b))
    fnvs=sorted(set(fnvs))
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        hits=[x for x in ex.map(fnv_hit,fnvs) if x]
    hits=sorted(hits)
    sum_grid=np.zeros((NY,NX),dtype=np.float64)
    count_grid=np.zeros((NY,NX),dtype=np.int64)
    files=[];total_records=0;total_good=0;total_used=0
    for k,fnv in enumerate(hits,1):
        fu=fbt_url(fnv)
        try:
            b=get(fu).content
            recs,good,used=parse_fbt_into_grid(b,sum_grid,count_grid)
            files.append({"fbt":fu,"bytes":len(b),"records":recs,"good_beams":good,"beams_in_tile":used})
            total_records+=recs;total_good+=good;total_used+=used
        except Exception as e:
            files.append({"fbt":fu,"error":type(e).__name__+": "+str(e)})
        if k%25==0:print(sid,"FBT",k,"/",len(hits),flush=True)
    valid=count_grid>0
    z=np.full((NY,NX),np.nan,dtype=np.float64)
    z[valid]=sum_grid[valid]/count_grid[valid]
    cand,scales,pre=discover_on_grid(sid,z,valid)
    all_candidates.extend(cand)
    survey_out[sid]={
      "metadata_features":feats,"tile_hit_fnv_count":len(hits),"files":files,
      "fbt_success_count":sum(1 for f in files if "error" not in f),
      "fbt_error_count":sum(1 for f in files if "error" in f),
      "records":total_records,"good_beams_total":total_good,"good_beams_in_tile":total_used,
      "occupied_grid_cells":int(np.sum(valid)),"scale_summaries":scales,
      "pre_dedup_followup_total":pre,"frozen_candidate_count":len(cand),"frozen_candidates":cand
    }
    print(sid,"occupied",int(np.sum(valid)),"candidates",len(cand),flush=True)

# Blind same-radius cross-survey persistence graph.
clusters=[]
for rad in RADII:
    pts=[q for q in all_candidates if int(q["radius_m"])==rad]
    if len(pts)<2:continue
    arr=np.array([[q["x_m"],q["y_m"]] for q in pts],dtype=float)
    tree=cKDTree(arr)
    pairs=tree.query_pairs(MATCH)
    parent=list(range(len(pts)))
    def find(a):
        while parent[a]!=a:
            parent[a]=parent[parent[a]];a=parent[a]
        return a
    def union(a,b):
        ra,rb=find(a),find(b)
        if ra!=rb:parent[rb]=ra
    for a,b in pairs:
        if pts[a]["survey_id"]!=pts[b]["survey_id"]:union(a,b)
    comp=defaultdict(list)
    for i in range(len(pts)):comp[find(i)].append(i)
    rr=[]
    for inds in comp.values():
        members=[pts[i] for i in inds]
        surveys=sorted(set(q["survey_id"] for q in members))
        if len(surveys)<MIN_CLUSTER_SURVEYS:continue
        x=float(np.mean([q["x_m"] for q in members]));y=float(np.mean([q["y_m"] for q in members]))
        lon,lat=ll(x,y)
        rr.append({
          "radius_m":rad,"lon":lon,"lat":lat,"x_m":x,"y_m":y,
          "distinct_survey_count":len(surveys),"surveys":surveys,
          "member_count":len(members),
          "members":[{"survey_id":q["survey_id"],"candidate_id":q["candidate_id"],"blind_rank":q["blind_rank"],"lon":q["lon"],"lat":q["lat"]} for q in members]
        })
    rr.sort(key=lambda q:(-q["distinct_survey_count"],min(m["blind_rank"] for m in q["members"]),q["lat"],q["lon"]))
    for q in rr[:int(PRE["cross_survey_persistence"]["maximum_clusters_per_radius"])]:
        q["cluster_id"]=f"R{rad}_P{len([x for x in clusters if x['radius_m']==rad])+1:03d}"
        clusters.append(q)

out={
 "artifact_id":"JANUS-KUSTO-SHARK-ROCK-G2B-BLIND-MULTISURVEY-MORPHOLOGY-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "truth_firewall":{"shark_rock_coordinate_read":False,"rov_interpretation_read":False,"video_stills_read":False},
 "grid":{"cell_m":GRID,"nx":NX,"ny":NY,"tile":T},
 "surveys":survey_out,
 "candidate_total":len(all_candidates),
 "persistent_cluster_count":len(clusters),
 "persistent_clusters":clusters,
 "persistence_support_histogram":{str(k):sum(1 for q in clusters if q["distinct_survey_count"]==k) for k in range(2,len(PRE["survey_ids"])+1)},
 "next_gate":"G3_UNLOCK_ROV_TRUTH_AND_SCORE_FROZEN_CANDIDATES_AND_PERSISTENT_CLUSTERS",
 "claim_ceiling":"BLIND_MULTISURVEY_MORPHOLOGY_CANDIDATE_AND_PERSISTENCE_GENERATION_ONLY"
}
p=OUT/"JANUS-KUSTO-SHARK-ROCK-G2B-BLIND-MULTISURVEY-MORPHOLOGY-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "survey_candidates":{k:v["frozen_candidate_count"] for k,v in survey_out.items()},
 "occupied_cells":{k:v["occupied_grid_cells"] for k,v in survey_out.items()},
 "candidate_total":out["candidate_total"],
 "persistent_cluster_count":out["persistent_cluster_count"],
 "persistence_support_histogram":out["persistence_support_histogram"],
 "truth_coordinate_read":False
},indent=2))
