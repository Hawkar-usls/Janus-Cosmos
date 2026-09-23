#!/usr/bin/env python3
from __future__ import annotations

import hashlib, io, json, math, os, random, zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import requests
import rasterio
from affine import Affine
from rasterio.io import MemoryFile
from scipy.stats import mannwhitneyu

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)

PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-CROSS-SOUND-ALASKA-UNSEEN-DUALCHANNEL-PREREG-2026-09-23-v1.0.json").read_text())
SRC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-CROSS-SOUND-CS0-SOURCE-BINDING-RECEIPT-2026-09-23-v1.0.json").read_text())
CAND_AUTH=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-CROSS-SOUND-CS1-BLIND-BATHYMETRY-RECEIPT-2026-09-23-v1.0.json").read_text())
IMP=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-CROSS-SOUND-CS2-BACKSCATTER-IMPLEMENTATION-ADDENDUM-2026-09-23-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-CrossSound-CS2-unseen-backscatter/1.0"}
ARTIFACT_ID=int(CAND_AUTH["run"]["artifact_id"])
EXPECTED_ARTIFACT_SHA=CAND_AUTH["run"]["artifact_zip_sha256"]
EXPECTED_JSON_SHA=CAND_AUTH["run"]["candidate_json_sha256"]
EXPECTED_CANDIDATES=int(CAND_AUTH["blind_execution"]["frozen_candidate_count"])
NOMINAL_CELL=float(IMP["payload_gate"]["georef_nominal_cell_m"])
REL_TOL=float(IMP["payload_gate"]["georef_relative_tolerance"])
MIN_VALID=float(IMP["method"]["minimum_patch_valid_fraction"])
RNG=random.Random(int(IMP["method"]["control_rng_seed"]))
METRICS=list(IMP["method"]["metrics"])
CONTROLS_PER=int(IMP["method"]["controls_per_candidate"])
MIN_CONTROLS=int(IMP["method"]["minimum_controls_per_radius"])
MIN_DIST=float(IMP["method"]["minimum_distance_from_any_frozen_candidate_m"])
MIN_ACTIVE=int(IMP["method"]["minimum_active_metrics_per_stratum"])
ITERS=int(IMP["method"]["permutation_iterations"])
ALPHA=float(IMP["method"]["alpha"])

def get(url,headers=None,timeout=180):
    h=dict(UA)
    if headers:h.update(headers)
    r=requests.get(url,headers=h,timeout=timeout,allow_redirects=True)
    r.raise_for_status()
    return r

def get_candidate_json():
    token=os.environ.get("GITHUB_TOKEN")
    repo=os.environ.get("GITHUB_REPOSITORY","hawkar-usls/Janus-Cosmos")
    if not token:raise RuntimeError("GITHUB_TOKEN required")
    u=f"https://api.github.com/repos/{repo}/actions/artifacts/{ARTIFACT_ID}/zip"
    blob=get(u,headers={
      "Authorization":f"Bearer {token}",
      "Accept":"application/vnd.github+json",
      "X-GitHub-Api-Version":"2022-11-28"
    }).content
    zsha=hashlib.sha256(blob).hexdigest()
    if zsha!=EXPECTED_ARTIFACT_SHA:
        raise RuntimeError(f"candidate artifact SHA mismatch {zsha} != {EXPECTED_ARTIFACT_SHA}")
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        js=[n for n in zf.namelist() if n.lower().endswith(".json") and not n.startswith("__MACOSX/")]
        if len(js)!=1:raise RuntimeError(f"expected one candidate JSON, found {js}")
        raw=zf.read(js[0])
    jsha=hashlib.sha256(raw).hexdigest()
    if jsha!=EXPECTED_JSON_SHA:
        raise RuntimeError(f"candidate JSON SHA mismatch {jsha} != {EXPECTED_JSON_SHA}")
    j=json.loads(raw)
    if int(j["frozen_candidate_count"])!=EXPECTED_CANDIDATES:
        raise RuntimeError("frozen candidate count mismatch")
    return j,zsha,jsha

def world_affine(text):
    vals=[float(x.strip()) for x in text.splitlines() if x.strip()]
    if len(vals)!=6:raise RuntimeError(f"world file requires 6 values, got {vals}")
    A,D,B,E,C,F=vals
    return vals,Affine(A,B,C-.5*A-.5*B,D,E,F-.5*D-.5*E)

def load_backscatter():
    url=IMP["backscatter_source"]["item_download_url"]
    blob=get(url).content
    if not blob:
        raise RuntimeError("empty backscatter payload")
    sha=hashlib.sha256(blob).hexdigest()
    tif=None; tname=None; members=[]; container=None
    if blob[:2]==b"PK":
        container="ZIP"
        with zipfile.ZipFile(io.BytesIO(blob)) as zf:
            members=zf.namelist()
            tifs=[n for n in members
                  if n.lower().endswith((".tif",".tiff"))
                  and not n.startswith("__MACOSX/")
                  and "/._" not in n
                  and "_qv" not in n.lower()
                  and "quick" not in n.lower()]
            if len(tifs)!=1:
                raise RuntimeError(f"expected exactly one primary TIFF in backscatter bundle, found {tifs}")
            tname=tifs[0]; tif=zf.read(tname)
    elif blob[:4] in (b"II*\\x00",b"MM\\x00*"):
        container="RAW_TIFF"; tif=blob; tname="ITEM_LEVEL_RAW_TIFF"
    else:
        raise RuntimeError(f"unsupported backscatter payload magic {blob[:16]!r}")
    with MemoryFile(tif) as mf:
        with mf.open() as src:
            ma=src.read(1,masked=True)
            data=np.asarray(ma.filled(np.nan),dtype=np.float32)
            valid=(~np.ma.getmaskarray(ma)) & np.isfinite(data)
            emb=src.transform
            crs=str(src.crs) if src.crs is not None else None
            nodata=src.nodata; dtype=str(src.dtypes[0])
    if emb.is_identity:
        raise RuntimeError("Cross Sound backscatter raster has identity geotransform")
    nominal=float(IMP["payload_gate"]["georef_nominal_cell_m"])
    tol=float(IMP["payload_gate"]["georef_relative_tolerance"])
    res=(abs(float(emb.a)),abs(float(emb.e)))
    if max(abs(res[0]-nominal),abs(res[1]-nominal))/nominal>tol:
        raise RuntimeError(f"effective backscatter resolution outside frozen tolerance: {res}")
    return {
      "data":data,"valid":valid,"transform":emb,"blob":blob,
      "sha256":sha,"tif_name":tname,"members":members,"payload_container":container,
      "embedded_transform":list(emb)[:6],"effective_transform":list(emb)[:6],
      "georef_source":"EMBEDDED_TRANSFORM","crs":crs,"nodata":nodata,"dtype":dtype,
      "resolution_m":list(res),"cell_eff_m":float(sum(res)/2.0),"url":url
    }

def patch_metrics(arr,valid,row,col,rpx):
    outer=2*rpx
    r0=row-outer;r1=row+outer+1;c0=col-outer;c1=col+outer+1
    if r0<0 or c0<0 or r1>arr.shape[0] or c1>arr.shape[1]:return None
    a=arr[r0:r1,c0:c1];v=valid[r0:r1,c0:c1]
    yy,xx=np.ogrid[-outer:outer+1,-outer:outer+1]
    d2=xx*xx+yy*yy
    disk=d2<=rpx*rpx
    ann=(d2>(1.25*rpx)**2)&(d2<=(2*rpx)**2)
    need=disk|ann
    if np.sum(need)==0 or float(np.mean(v[need]))<MIN_VALID:return None
    dv=a[disk&v].astype(np.float64);av=a[ann&v].astype(np.float64)
    if len(dv)<10 or len(av)<10:return None
    dstd=float(np.std(dv))
    contrast=float(abs(np.mean(dv)-np.mean(av))/(np.std(av)+1e-6))
    drange=float(np.quantile(dv,.95)-np.quantile(dv,.05))
    core=a.astype(np.float64)
    gx=np.abs(np.diff(core,axis=1));gy=np.abs(np.diff(core,axis=0))
    vx=v[:,:-1]&v[:,1:];vy=v[:-1,:]&v[1:,:]
    maskx=disk[:,:-1]&disk[:,1:]&vx
    masky=disk[:-1,:]&disk[1:,:]&vy
    edges=np.concatenate([gx[maskx],gy[masky]])
    edge=float(np.median(edges)) if len(edges) else float("nan")
    vals=[dstd,contrast,drange,edge]
    if not all(np.isfinite(vals)):return None
    return {
      "disk_standard_deviation":dstd,
      "absolute_center_vs_annulus_contrast_over_annulus_sd":contrast,
      "disk_q95_minus_q05_range":drange,
      "median_first_difference_edge_energy":edge
    }

candj,artifact_sha,candidate_sha=get_candidate_json()
candidates=list(candj["frozen_candidates"])
if len(candidates)!=EXPECTED_CANDIDATES:raise RuntimeError("candidate list length mismatch")
bs=load_backscatter()
arr=bs["data"];valid=bs["valid"];transform=bs["transform"];cell=bs["cell_eff_m"]
frozen_xy=np.array([[float(c["x_m"]),float(c["y_m"])] for c in candidates],dtype=float)

candidate_results=[];by_radius=defaultdict(list)
for c in candidates:
    row,col=rasterio.transform.rowcol(transform,float(c["x_m"]),float(c["y_m"]))
    radius=int(c["radius_m"]);rpx=max(1,int(round(radius/cell)))
    m=patch_metrics(arr,valid,int(row),int(col),rpx)
    rec={
      "candidate_id":c["candidate_id"],"blind_rank":c["blind_rank"],
      "x_m":c["x_m"],"y_m":c["y_m"],"radius_m":radius,
      "row":int(row),"col":int(col),"metrics":m
    }
    candidate_results.append(rec);by_radius[radius].append(rec)

def far_from_candidates(x,y):
    d2=(frozen_xy[:,0]-x)**2+(frozen_xy[:,1]-y)**2
    return bool(np.min(d2)>=MIN_DIST**2)

strata={};blocked=[];all_candidate_scores=[];all_control_scores=[]
for radius in sorted(by_radius):
    key=str(radius);frozen=by_radius[radius];scored=[q for q in frozen if q["metrics"] is not None]
    if not scored:
        strata[key]={
          "radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":0,
          "status":"COVERAGE_BLOCKED","active_metrics":[],"excluded_zero_mad_metrics":[],
          "candidate_support_count":0
        }
        blocked.append({"stratum":key,"reason":"COVERAGE_BLOCKED"});continue
    target=max(MIN_CONTROLS,len(scored)*CONTROLS_PER)
    controls=[];tries=0;rpx=max(1,int(round(radius/cell)));outer=2*rpx
    h,w=arr.shape
    while len(controls)<target and tries<target*1000:
        tries+=1
        if h<=2*outer+1 or w<=2*outer+1:break
        row=RNG.randint(outer,h-outer-1);col=RNG.randint(outer,w-outer-1)
        if not valid[row,col]:continue
        x,y=rasterio.transform.xy(transform,row,col,offset="center")
        if not far_from_candidates(float(x),float(y)):continue
        m=patch_metrics(arr,valid,row,col,rpx)
        if m is None:continue
        controls.append({"x_m":float(x),"y_m":float(y),"metrics":m})
    if len(controls)<target:raise RuntimeError(f"{key}: controls {len(controls)}/{target}")

    control_mat=np.array([[q["metrics"][m] for m in METRICS] for q in controls],dtype=float)
    med=np.median(control_mat,axis=0);mad=np.median(np.abs(control_mat-med),axis=0)
    active_idx=[i for i,x in enumerate(mad) if float(x)>0.0]
    excluded_idx=[i for i,x in enumerate(mad) if float(x)==0.0]
    active=[METRICS[i] for i in active_idx];excluded=[METRICS[i] for i in excluded_idx]
    if len(active)<MIN_ACTIVE:
        strata[key]={
          "radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":len(scored),
          "control_count":len(controls),"status":"DEGENERATE_METRICS_BLOCKED",
          "active_metrics":active,"excluded_zero_mad_metrics":excluded,
          "control_metric_median":{m:float(med[i]) for i,m in enumerate(METRICS)},
          "control_metric_mad":{m:float(mad[i]) for i,m in enumerate(METRICS)},
          "candidate_support_count":0
        }
        blocked.append({"stratum":key,"reason":"DEGENERATE_METRICS_BLOCKED"});continue

    meda=med[active_idx];mada=mad[active_idx];scale=1.4826*mada
    ca=control_mat[:,active_idx]
    cz=np.maximum(0.0,(ca-meda)/scale)
    cscores=.5*np.max(cz,axis=1)+.5*np.mean(cz,axis=1)
    q95=float(np.quantile(cscores,.95))
    cand_mat=np.array([[q["metrics"][m] for m in METRICS] for q in scored],dtype=float)[:,active_idx]
    kz=np.maximum(0.0,(cand_mat-meda)/scale)
    kscores=.5*np.max(kz,axis=1)+.5*np.mean(kz,axis=1)
    support=0
    for q,zv,score in zip(scored,kz,kscores):
        q["active_metrics"]=active;q["excluded_zero_mad_metrics"]=excluded
        q["metric_positive_robust_z"]={m:float(zv[i]) for i,m in enumerate(active)}
        q["aggregate_score"]=float(score);q["same_stratum_control_q95"]=q95
        q["cross_channel_support"]=bool(score>=q95);support+=int(score>=q95)
    all_candidate_scores.extend(float(x) for x in kscores)
    all_control_scores.extend(float(x) for x in cscores)
    strata[key]={
      "radius_m":radius,"frozen_candidate_count":len(frozen),
      "scored_candidate_count":len(scored),"unscored_candidate_count":len(frozen)-len(scored),
      "control_count":len(controls),"status":"SCORABLE",
      "control_metric_median":{m:float(med[i]) for i,m in enumerate(METRICS)},
      "control_metric_mad":{m:float(mad[i]) for i,m in enumerate(METRICS)},
      "active_metrics":active,"excluded_zero_mad_metrics":excluded,
      "control_score_q95":q95,"candidate_support_count":support,
      "_control_scores":[float(x) for x in cscores]
    }

scorable=[k for k,v in strata.items() if v["status"]=="SCORABLE"]
observed=float(np.mean(all_candidate_scores)) if all_candidate_scores else float("nan")
control_mean=float(np.mean(all_control_scores)) if all_control_scores else float("nan")
if blocked or not scorable:
    perm_p=None
else:
    ge=0
    for _ in range(ITERS):
        sim=[]
        for key in scorable:
            s=strata[key];n=int(s["scored_candidate_count"]);pool=s["_control_scores"]
            if n>len(pool):raise RuntimeError(f"{key}: permutation n > controls")
            sim.extend(RNG.sample(pool,n))
        if float(np.mean(sim))>=observed:ge+=1
    perm_p=(ge+1)/(ITERS+1)

if all_candidate_scores and all_control_scores:
    mw=mannwhitneyu(all_candidate_scores,all_control_scores,alternative="greater",method="asymptotic")
    mw_u=float(mw.statistic);mw_p=float(mw.pvalue)
else:
    mw_u=None;mw_p=None

support_count=sum(int(v.get("candidate_support_count",0)) for v in strata.values())
scored_count=sum(int(v.get("scored_candidate_count",0)) for v in strata.values() if v["status"]=="SCORABLE")
support_fraction=(support_count/scored_count) if scored_count else None
enrichment=(support_fraction/.05) if support_fraction is not None else None
primary_pass=bool(
  not blocked and perm_p is not None and perm_p<ALPHA and
  np.isfinite(observed) and np.isfinite(control_mean) and observed>control_mean and
  enrichment is not None and enrichment>1.0
)
for v in strata.values():v.pop("_control_scores",None)

if blocked:interpretation="PRIMARY_VALIDATION_BLOCKED_BY_COVERAGE_OR_DEGENERATE_STRATUM"
elif primary_pass:interpretation="PASS_PREREGISTERED_UNSEEN_SAME_SURVEY_CROSS_CHANNEL_VALIDATION"
else:interpretation="NO_PREREGISTERED_UNSEEN_SAME_SURVEY_CROSS_CHANNEL_VALIDATION"

out={
 "artifact_id":"JANUS-KUSTO-CROSS-SOUND-CS2-UNSEEN-BACKSCATTER-REPLAY-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],"implementation_addendum":IMP["artifact_id"],
 "candidate_authority":{
   "receipt":CAND_AUTH["artifact_id"],"artifact_id":ARTIFACT_ID,
   "artifact_zip_sha256_verified":artifact_sha,"candidate_json_sha256_verified":candidate_sha,
   "frozen_candidate_count":len(candidates)
 },
 "backscatter_input":{
   "url":bs["url"],"zip_bytes":len(bs["blob"]),
   "zip_sha256":bs["sha256"],"tif_name":bs["tif_name"],
   "zip_members":bs["members"],"shape":list(arr.shape),"dtype":bs["dtype"],"nodata":bs["nodata"],
   "valid_pixels":int(np.sum(valid)),"resolution_m":bs["resolution_m"],"effective_cell_m":cell,
   "georef_source":bs["georef_source"],"crs":bs["crs"],
   "embedded_transform":bs["embedded_transform"],"effective_transform":bs["effective_transform"],
   "payload_container":bs["payload_container"]
 },
 "strata":strata,"candidate_results":candidate_results,
 "coverage":{"frozen_candidate_count":len(candidates),"scorable_candidate_count":scored_count,"blocked_strata":blocked},
 "population":{
   "candidate_mean_aggregate_score":observed if np.isfinite(observed) else None,
   "control_mean_aggregate_score":control_mean if np.isfinite(control_mean) else None,
   "candidate_support_count":support_count,"candidate_support_fraction":support_fraction,
   "nominal_control_tail_fraction":0.05,"support_enrichment_ratio":enrichment,
   "stratified_permutation_iterations":ITERS,"stratified_permutation_p":perm_p,
   "mann_whitney_u":mw_u,"mann_whitney_one_sided_p":mw_p,
   "primary_pass":primary_pass,"scientific_promotion":primary_pass
 },
 "truth_firewall":{
   "landslide_or_fault_labels_read":False,"candidate_coordinates_changed":False,
   "candidate_radii_changed":False
 },
 "interpretation":interpretation,
 "claim_ceiling":"PREREGISTERED_UNSEEN_OPEN_OCEAN_SAME_SURVEY_CROSS_CHANNEL_VALIDATION_ONLY__NO_INDEPENDENT_SURVEY_REPLICATION"
}
p=OUT/"JANUS-KUSTO-CROSS-SOUND-CS2-UNSEEN-BACKSCATTER-REPLAY-RUN-2026-09-23-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"candidate_hash_verified":candidate_sha,
 "frozen_candidates":len(candidates),"scorable_candidates":scored_count,
 "blocked_strata":blocked,"support_count":support_count,"support_fraction":support_fraction,
 "support_enrichment_ratio":enrichment,"candidate_mean_score":out["population"]["candidate_mean_aggregate_score"],
 "control_mean_score":out["population"]["control_mean_aggregate_score"],
 "permutation_p":perm_p,"mann_whitney_p":mw_p,"primary_pass":primary_pass,
 "scientific_promotion":primary_pass,"interpretation":interpretation,
 "output_json_sha256":hashlib.sha256(raw.encode()).hexdigest()
},indent=2))
