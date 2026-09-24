#!/usr/bin/env python3
from __future__ import annotations

import hashlib, io, json, math, os, random, zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import requests
import rasterio
from rasterio.io import MemoryFile
from rasterio.windows import Window
from pyproj import Transformer
from scipy.stats import mannwhitneyu

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)

PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-MONEY-SHOAL-UNSEEN-V2-PREREG-2026-09-23-v1.0.json").read_text())
COR=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-A0-PREREG-METHOD-INVARIANT-CORRECTION-2026-09-24-v1.0.json").read_text())
SRC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-A0E-EXACT-MEMBER-SELECTION-FREEZE-2026-09-24-v1.0.json").read_text())
BR=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-A0G-MONEY-SHOAL-REPRESENTATION-BRIDGE-RECEIPT-2026-09-24-v1.0.json").read_text())
AUTH=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-A1-MONEY-SHOAL-BLIND-BATHYMETRY-RECEIPT-2026-09-24-v1.0.json").read_text())
IMP=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-A2-UNCHANGED-V2-IMPLEMENTATION-FREEZE-2026-09-24-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-Arafura-MoneyShoal-A2-v2/1.0"}
ARTIFACT_ID=int(AUTH["run"]["artifact_id"])
EXPECTED_ARTIFACT_SHA=AUTH["run"]["artifact_zip_sha256"]
EXPECTED_JSON_SHA=AUTH["run"]["candidate_json_sha256"]
EXPECTED_CANDIDATES=int(AUTH["blind_execution"]["frozen_candidate_count"])

BATHY_CELL=6.0
BACK_CELL=2.0
BATHY_SUPPORT_MIN=.80
BACK_SUPPORT_MIN=.80
CONTROLS_PER=20
MIN_CONTROLS=200
MIN_DIST=100.0
RNG=random.Random(260923)
METRICS=[
 "disk_standard_deviation",
 "absolute_center_vs_annulus_contrast_over_annulus_sd",
 "disk_q95_minus_q05_range",
 "median_first_difference_edge_energy"
]
MIN_ACTIVE=2
ITERS=5000
ALPHA=.05
to_wgs=Transformer.from_crs("EPSG:32753","EPSG:4326",always_xy=True)

def get(url,headers=None,timeout=300):
    h=dict(UA)
    if headers:h.update(headers)
    r=requests.get(url,headers=h,timeout=timeout,allow_redirects=True)
    r.raise_for_status()
    return r

def get_candidate_json():
    token=os.environ.get("GITHUB_TOKEN")
    repo=os.environ.get("GITHUB_REPOSITORY","hawkar-usls/Janus-Cosmos")
    if not token: raise RuntimeError("GITHUB_TOKEN required")
    u=f"https://api.github.com/repos/{repo}/actions/artifacts/{ARTIFACT_ID}/zip"
    blob=get(u,headers={
      "Authorization":f"Bearer {token}",
      "Accept":"application/vnd.github+json",
      "X-GitHub-Api-Version":"2022-11-28"
    }).content
    zsha=hashlib.sha256(blob).hexdigest()
    if zsha!=EXPECTED_ARTIFACT_SHA: raise RuntimeError(f"candidate artifact SHA mismatch {zsha}")
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        js=[n for n in zf.namelist() if n.lower().endswith(".json") and not n.startswith("__MACOSX/")]
        if len(js)!=1: raise RuntimeError(f"expected one candidate JSON, found {js}")
        raw=zf.read(js[0])
    jsha=hashlib.sha256(raw).hexdigest()
    if jsha!=EXPECTED_JSON_SHA: raise RuntimeError(f"candidate JSON SHA mismatch {jsha}")
    j=json.loads(raw)
    if int(j["frozen_candidate_count"])!=EXPECTED_CANDIDATES: raise RuntimeError("candidate count mismatch")
    return j,zsha,jsha

def fetch_member(spec):
    blob=get(spec["archive_url"]).content
    sha=hashlib.sha256(blob).hexdigest()
    if sha.lower()!=spec["archive_sha256"].lower(): raise RuntimeError("source archive SHA mismatch")
    if len(blob)!=int(spec["archive_bytes"]): raise RuntimeError("source archive byte length mismatch")
    member=spec["member"]
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        zi=zf.getinfo(member)
        if int(zi.file_size)!=int(spec["member_uncompressed_bytes"]): raise RuntimeError("member size mismatch")
        if int(zi.CRC)!=int(spec["member_crc32"]): raise RuntimeError("member CRC mismatch")
        raw=zf.read(member)
    return blob,raw,sha

def load_bathy():
    blob,raw,sha=fetch_member(SRC["bathymetry"])
    w=BR["conservative_bathymetry_source_window"]
    win=Window(int(w["col_start"]),int(w["row_start"]),int(w["width"]),int(w["height"]))
    with MemoryFile(raw) as mf:
      with mf.open() as ds:
        if str(ds.crs)!="EPSG:4326": raise RuntimeError(f"unexpected bathy CRS {ds.crs}")
        ma=ds.read(1,window=win,masked=True)
        data=np.asarray(ma.filled(np.nan),dtype=np.float32)
        valid=(~np.ma.getmaskarray(ma)) & np.isfinite(data)
        transform=ds.window_transform(win)
        meta={"source_url":SRC["bathymetry"]["archive_url"],"archive_sha256":sha,"archive_bytes":len(blob),
              "member":SRC["bathymetry"]["member"],"shape":list(data.shape),"valid_pixels":int(np.sum(valid)),
              "crs":str(ds.crs),"transform":list(transform)[:6],"nominal_cell_m":BATHY_CELL}
    return data,valid,transform,meta

def load_back():
    blob,raw,sha=fetch_member(SRC["heldout_backscatter"])
    with MemoryFile(raw) as mf:
      with mf.open() as ds:
        if str(ds.crs)!="EPSG:32753": raise RuntimeError(f"unexpected backscatter CRS {ds.crs}")
        ma=ds.read(1,masked=True)
        data=np.asarray(ma.filled(np.nan),dtype=np.float32)
        valid=(~np.ma.getmaskarray(ma)) & np.isfinite(data)
        transform=ds.transform
        res=(abs(float(transform.a)),abs(float(transform.e)))
        if max(abs(res[0]-BACK_CELL),abs(res[1]-BACK_CELL))/BACK_CELL>.001:
            raise RuntimeError(f"unexpected backscatter resolution {res}")
        meta={"source_url":SRC["heldout_backscatter"]["archive_url"],"archive_sha256":sha,"archive_bytes":len(blob),
              "member":SRC["heldout_backscatter"]["member"],"shape":list(data.shape),"valid_pixels":int(np.sum(valid)),
              "crs":str(ds.crs),"transform":list(transform)[:6],"resolution_m":list(res)}
    return data,valid,transform,meta

def patch_geometry(radius_m,cell_m):
    rpx=max(1,int(round(radius_m/cell_m))); outer=2*rpx
    yy,xx=np.ogrid[-outer:outer+1,-outer:outer+1]
    d2=xx*xx+yy*yy
    disk=d2<=rpx*rpx
    ann=(d2>(1.25*rpx)**2)&(d2<=(2*rpx)**2)
    return rpx,outer,disk,ann,(disk|ann)

def bathy_support(valid,transform,x,y,radius):
    lon,lat=to_wgs.transform(float(x),float(y))
    row,col=rasterio.transform.rowcol(transform,lon,lat)
    _,outer,disk,ann,need=patch_geometry(radius,BATHY_CELL)
    r0,r1=row-outer,row+outer+1; c0,c1=col-outer,col+outer+1
    if r0<0 or c0<0 or r1>valid.shape[0] or c1>valid.shape[1]: return 0.0,int(row),int(col)
    v=valid[r0:r1,c0:c1]
    return float(np.mean(v[need])) if np.sum(need) else 0.0,int(row),int(col)

def back_metrics(arr,valid,transform,x,y,radius):
    row,col=rasterio.transform.rowcol(transform,float(x),float(y))
    _,outer,disk,ann,need=patch_geometry(radius,BACK_CELL)
    r0,r1=row-outer,row+outer+1; c0,c1=col-outer,col+outer+1
    if r0<0 or c0<0 or r1>arr.shape[0] or c1>arr.shape[1]: return None,0.0,int(row),int(col)
    a=arr[r0:r1,c0:c1]; v=valid[r0:r1,c0:c1]
    frac=float(np.mean(v[need])) if np.sum(need) else 0.0
    if frac<BACK_SUPPORT_MIN: return None,frac,int(row),int(col)
    dv=a[disk&v].astype(np.float64); av=a[ann&v].astype(np.float64)
    if len(dv)<10 or len(av)<10:return None,frac,int(row),int(col)
    dstd=float(np.std(dv))
    contrast=float(abs(np.mean(dv)-np.mean(av))/(np.std(av)+1e-6))
    drange=float(np.quantile(dv,.95)-np.quantile(dv,.05))
    core=a.astype(np.float64)
    gx=np.abs(np.diff(core,axis=1)); gy=np.abs(np.diff(core,axis=0))
    vx=v[:,:-1]&v[:,1:]; vy=v[:-1,:]&v[1:,:]
    maskx=disk[:,:-1]&disk[:,1:]&vx; masky=disk[:-1,:]&disk[1:,:]&vy
    edges=np.concatenate([gx[maskx],gy[masky]])
    edge=float(np.median(edges)) if len(edges) else float("nan")
    vals=[dstd,contrast,drange,edge]
    if not all(np.isfinite(vals)):return None,frac,int(row),int(col)
    return {
      "disk_standard_deviation":dstd,
      "absolute_center_vs_annulus_contrast_over_annulus_sd":contrast,
      "disk_q95_minus_q05_range":drange,
      "median_first_difference_edge_energy":edge
    },frac,int(row),int(col)

candj,artifact_sha,candidate_sha=get_candidate_json()
candidates=list(candj["frozen_candidates"])
bathy,bvalid,btransform,bmeta=load_bathy()
back,backvalid,backtransform,backmeta=load_back()
frozen_xy=np.array([[float(c["x_m"]),float(c["y_m"])] for c in candidates],dtype=float)

def far(x,y):
    d2=(frozen_xy[:,0]-x)**2+(frozen_xy[:,1]-y)**2
    return bool(np.min(d2)>=MIN_DIST**2)

candidate_results=[]; by_radius=defaultdict(list)
for c in candidates:
    x,y=float(c["x_m"]),float(c["y_m"]); radius=int(c["radius_m"])
    bsup,brow,bcol=bathy_support(bvalid,btransform,x,y,radius)
    m=None; ksup=0.0; krow=kcol=None
    if bsup>=BATHY_SUPPORT_MIN:
        m,ksup,krow,kcol=back_metrics(back,backvalid,backtransform,x,y,radius)
    rec={"candidate_id":c["candidate_id"],"blind_rank":c["blind_rank"],"x_m":x,"y_m":y,
         "lon":c["lon"],"lat":c["lat"],"radius_m":radius,
         "bathymetry_support_fraction":bsup,"backscatter_native_support_fraction":ksup,
         "bathymetry_local_row":brow,"bathymetry_local_col":bcol,
         "backscatter_row":krow,"backscatter_col":kcol,"metrics":m}
    candidate_results.append(rec);by_radius[radius].append(rec)

strata={}; blocked=[]; all_candidate_scores=[];all_control_scores=[]

for radius in sorted(by_radius):
    key=str(radius); frozen=by_radius[radius]
    scored=[q for q in frozen if q["metrics"] is not None and q["bathymetry_support_fraction"]>=BATHY_SUPPORT_MIN and q["backscatter_native_support_fraction"]>=BACK_SUPPORT_MIN]
    if not scored:
        strata[key]={"radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":0,
                     "status":"COVERAGE_BLOCKED","candidate_support_count":0,"active_metrics":[],"excluded_zero_mad_metrics":[]}
        blocked.append({"stratum":key,"reason":"COVERAGE_BLOCKED"});continue

    target=max(MIN_CONTROLS,len(scored)*CONTROLS_PER)
    controls=[];tries=0
    _,outer,_,_,_=patch_geometry(radius,BACK_CELL)
    h,w=backvalid.shape
    while len(controls)<target and tries<target*2000:
        tries+=1
        row=RNG.randint(outer,h-outer-1); col=RNG.randint(outer,w-outer-1)
        x,y=rasterio.transform.xy(backtransform,row,col,offset="center");x=float(x);y=float(y)
        if not far(x,y):continue
        bsup,_,_=bathy_support(bvalid,btransform,x,y,radius)
        if bsup<BATHY_SUPPORT_MIN:continue
        m,ksup,_,_=back_metrics(back,backvalid,backtransform,x,y,radius)
        if m is None or ksup<BACK_SUPPORT_MIN:continue
        controls.append({"x_m":x,"y_m":y,"bathymetry_support_fraction":bsup,
                         "backscatter_native_support_fraction":ksup,"metrics":m})

    if len(controls)<target:
        strata[key]={"radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":len(scored),
                     "control_count":len(controls),"required_control_count":target,"status":"CONTROL_SUPPORT_BLOCKED",
                     "candidate_support_count":0,"active_metrics":[],"excluded_zero_mad_metrics":[]}
        blocked.append({"stratum":key,"reason":"CONTROL_SUPPORT_BLOCKED"});continue

    cm=np.array([[q["metrics"][m] for m in METRICS] for q in controls],dtype=float)
    med=np.median(cm,axis=0); mad=np.median(np.abs(cm-med),axis=0)
    active_idx=[i for i,x in enumerate(mad) if float(x)>0.0]
    excluded_idx=[i for i,x in enumerate(mad) if float(x)==0.0]
    active=[METRICS[i] for i in active_idx]; excluded=[METRICS[i] for i in excluded_idx]
    if len(active)<MIN_ACTIVE:
        strata[key]={"radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":len(scored),
                     "control_count":len(controls),"status":"DEGENERATE_METRICS_BLOCKED","active_metrics":active,
                     "excluded_zero_mad_metrics":excluded,"candidate_support_count":0,
                     "control_metric_median":{m:float(med[i]) for i,m in enumerate(METRICS)},
                     "control_metric_mad":{m:float(mad[i]) for i,m in enumerate(METRICS)}}
        blocked.append({"stratum":key,"reason":"DEGENERATE_METRICS_BLOCKED"});continue

    meda=med[active_idx]; mada=mad[active_idx]; scale=1.4826*mada
    ca=cm[:,active_idx]; cz=np.maximum(0.0,(ca-meda)/scale)
    cscores=.5*np.max(cz,axis=1)+.5*np.mean(cz,axis=1); q95=float(np.quantile(cscores,.95))
    km=np.array([[q["metrics"][m] for m in METRICS] for q in scored],dtype=float)[:,active_idx]
    kz=np.maximum(0.0,(km-meda)/scale)
    kscores=.5*np.max(kz,axis=1)+.5*np.mean(kz,axis=1)

    support=0
    for q,zv,score in zip(scored,kz,kscores):
        q["active_metrics"]=active;q["excluded_zero_mad_metrics"]=excluded
        q["metric_positive_robust_z"]={m:float(zv[i]) for i,m in enumerate(active)}
        q["aggregate_score"]=float(score);q["same_stratum_control_q95"]=q95
        q["cross_channel_support"]=bool(score>=q95);support+=int(score>=q95)

    all_candidate_scores.extend(float(x) for x in kscores);all_control_scores.extend(float(x) for x in cscores)
    strata[key]={"radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":len(scored),
                 "unscored_candidate_count":len(frozen)-len(scored),"control_count":len(controls),"status":"SCORABLE",
                 "control_metric_median":{m:float(med[i]) for i,m in enumerate(METRICS)},
                 "control_metric_mad":{m:float(mad[i]) for i,m in enumerate(METRICS)},
                 "active_metrics":active,"excluded_zero_mad_metrics":excluded,"control_score_q95":q95,
                 "candidate_support_count":support,"_control_scores":[float(x) for x in cscores]}

scorable=[k for k,v in strata.items() if v["status"]=="SCORABLE"]
observed=float(np.mean(all_candidate_scores)) if all_candidate_scores else float("nan")
control_mean=float(np.mean(all_control_scores)) if all_control_scores else float("nan")

if blocked or not scorable: perm_p=None
else:
    ge=0
    for _ in range(ITERS):
        sim=[]
        for key in scorable:
            s=strata[key];n=int(s["scored_candidate_count"]);pool=s["_control_scores"]
            sim.extend(RNG.sample(pool,n))
        if float(np.mean(sim))>=observed:ge+=1
    perm_p=(ge+1)/(ITERS+1)

if all_candidate_scores and all_control_scores:
    mw=mannwhitneyu(all_candidate_scores,all_control_scores,alternative="greater",method="asymptotic")
    mw_u=float(mw.statistic);mw_p=float(mw.pvalue)
else: mw_u=mw_p=None

support_count=sum(int(v.get("candidate_support_count",0)) for v in strata.values())
scored_count=sum(int(v.get("scored_candidate_count",0)) for v in strata.values() if v["status"]=="SCORABLE")
support_fraction=support_count/scored_count if scored_count else None
enrichment=support_fraction/.05 if support_fraction is not None else None
primary_pass=bool(not blocked and perm_p is not None and perm_p<ALPHA and
                  np.isfinite(observed) and np.isfinite(control_mean) and observed>control_mean and
                  enrichment is not None and enrichment>1.0)
for v in strata.values():v.pop("_control_scores",None)

interpretation=("PRIMARY_VALIDATION_BLOCKED_BY_SUPPORT_COVERAGE_OR_DEGENERATE_STRATUM" if blocked else
                "PASS_PREREGISTERED_UNSEEN_ARAFURA_MONEY_SHOAL_PAIRED_SUPPORT_DOMAIN_V2_VALIDATION" if primary_pass else
                "NO_PREREGISTERED_UNSEEN_ARAFURA_MONEY_SHOAL_PAIRED_SUPPORT_DOMAIN_V2_VALIDATION")

out={
 "artifact_id":"JANUS-KUSTO-ARAFURA-A2-MONEY-SHOAL-UNSEEN-BACKSCATTER-V2-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],"method_invariant_correction":COR["artifact_id"],"implementation_freeze":IMP["artifact_id"],
 "candidate_authority":{"receipt":AUTH["artifact_id"],"artifact_id":ARTIFACT_ID,
                        "artifact_zip_sha256_verified":artifact_sha,"candidate_json_sha256_verified":candidate_sha,
                        "frozen_candidate_count":len(candidates)},
 "paired_inputs":{"bathymetry":bmeta,"backscatter":backmeta},
 "support_domain_v2":{"bathymetry_support_min_fraction":BATHY_SUPPORT_MIN,
                      "backscatter_native_mask_min_fraction":BACK_SUPPORT_MIN,
                      "candidate_and_control_rule_identical":True,
                      "backscatter_value_based_support_inference_used":False},
 "strata":strata,"candidate_results":candidate_results,
 "coverage":{"frozen_candidate_count":len(candidates),"scorable_candidate_count":scored_count,"blocked_strata":blocked},
 "population":{"candidate_mean_aggregate_score":observed if np.isfinite(observed) else None,
               "control_mean_aggregate_score":control_mean if np.isfinite(control_mean) else None,
               "candidate_support_count":support_count,"candidate_support_fraction":support_fraction,
               "nominal_control_tail_fraction":.05,"support_enrichment_ratio":enrichment,
               "stratified_permutation_iterations":ITERS,"stratified_permutation_p":perm_p,
               "mann_whitney_u":mw_u,"mann_whitney_one_sided_p":mw_p,
               "primary_pass":primary_pass,"scientific_promotion":primary_pass},
 "truth_firewall":{"seabed_samples_read":False,"underwater_imagery_read":False,
                   "candidate_coordinates_changed":False,"candidate_radii_changed":False,
                   "pillar_bank_backscatter_used":False},
 "interpretation":interpretation,
 "claim_ceiling":"PREREGISTERED_UNSEEN_ARAFURA_MONEY_SHOAL_SAME_SURVEY_CROSSCHANNEL_V2_VALIDATION_ONLY__NO_INDEPENDENT_SURVEY_REPLICATION"
}
p=OUT/"JANUS-KUSTO-ARAFURA-A2-MONEY-SHOAL-UNSEEN-BACKSCATTER-V2-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"candidate_hash_verified":candidate_sha,"frozen_candidates":len(candidates),
 "scorable_candidates":scored_count,"blocked_strata":blocked,"support_count":support_count,
 "support_fraction":support_fraction,"support_enrichment_ratio":enrichment,
 "candidate_mean_score":out["population"]["candidate_mean_aggregate_score"],
 "control_mean_score":out["population"]["control_mean_aggregate_score"],
 "permutation_p":perm_p,"mann_whitney_p":mw_p,"primary_pass":primary_pass,
 "interpretation":interpretation,"output_json_sha256":hashlib.sha256(raw.encode()).hexdigest()
},indent=2))
