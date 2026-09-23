#!/usr/bin/env python3
from __future__ import annotations

import hashlib, io, json, math, os, random, zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import requests
import rasterio
from affine import Affine
from pyproj import CRS, Geod, Transformer
from rasterio.io import MemoryFile
from scipy.stats import mannwhitneyu

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)

PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOUTHWEST-CORNER-SWC0-UNCHANGED-V2-PREREG-2026-09-24-v1.0.json").read_text())
SRC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOUTHWEST-CORNER-SWC0D-EXACT-MEMBER-BINDING-RECEIPT-2026-09-24-v1.0.json").read_text())
CAND_AUTH=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOUTHWEST-CORNER-SWC1-BLIND-BATHYMETRY-RECEIPT-2026-09-24-v1.0.json").read_text())
IMP=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOUTHWEST-CORNER-SWC2-UNCHANGED-V2-IMPLEMENTATION-FREEZE-2026-09-24-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-SouthwestCorner-SWC2-support-domain-v2/1.0"}
ARTIFACT_ID=int(CAND_AUTH["run"]["artifact_id"])
EXPECTED_ARTIFACT_SHA=CAND_AUTH["run"]["artifact_zip_sha256"]
EXPECTED_JSON_SHA=CAND_AUTH["run"]["candidate_json_sha256"]
EXPECTED_CANDIDATES=int(CAND_AUTH["blind_execution"]["frozen_candidate_count"])

BATHY_CELL=5.0
BACK_CELL=4.0
BATHY_SUPPORT_MIN=float(IMP["support_domain_v2"]["bathymetry_support_min_fraction"])
BACK_MASK_MIN=float(IMP["support_domain_v2"]["backscatter_native_mask_min_fraction"])
RNG=random.Random(int(IMP["scoring"]["control_rng_seed"]))
METRICS=list(IMP["scoring"]["metrics"])
CONTROLS_PER=int(IMP["scoring"]["controls_per_candidate"])
MIN_CONTROLS=int(IMP["scoring"]["minimum_controls_per_radius"])
MIN_DIST=float(IMP["scoring"]["minimum_distance_from_any_frozen_candidate_m"])
MIN_ACTIVE=int(IMP["scoring"]["minimum_active_metrics_per_stratum"])
ITERS=int(IMP["scoring"]["permutation_iterations"])
ALPHA=float(IMP["scoring"]["alpha"])
METRIC_CRS=CRS.from_user_input("EPSG:32750")
WGS84=CRS.from_user_input("EPSG:4326")
TO_METRIC=Transformer.from_crs(WGS84,METRIC_CRS,always_xy=True)
GEOD=Geod(ellps="WGS84")

def get(url,headers=None,timeout=300):
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
    if zsha!=EXPECTED_ARTIFACT_SHA:raise RuntimeError("candidate artifact SHA mismatch")
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        js=[n for n in zf.namelist() if n.lower().endswith(".json") and not n.startswith("__MACOSX/")]
        if len(js)!=1:raise RuntimeError(f"expected one candidate JSON, found {js}")
        raw=zf.read(js[0])
    jsha=hashlib.sha256(raw).hexdigest()
    if jsha!=EXPECTED_JSON_SHA:raise RuntimeError("candidate JSON SHA mismatch")
    j=json.loads(raw)
    if int(j["frozen_candidate_count"])!=EXPECTED_CANDIDATES:raise RuntimeError("candidate count mismatch")
    for c in j["frozen_candidates"]:
        if "lon" not in c or "lat" not in c or "x_m" not in c or "y_m" not in c:
            raise RuntimeError("candidate georef fields missing")
    return j,zsha,jsha

def world_affine(text):
    vals=[float(x.strip()) for x in text.splitlines() if x.strip()]
    if len(vals)!=6:raise RuntimeError(f"world file requires 6 values, got {vals}")
    A,D,B,E,C,F=vals
    return vals,Affine(A,B,C-.5*A-.5*B,D,E,F-.5*D-.5*E)

def load_archive_raster(kind):
    arc=SRC[f"{kind}_archive"]
    blob=get(arc["url"]).content
    sha=hashlib.sha256(blob).hexdigest()
    if sha.lower()!=arc["sha256"].lower():raise RuntimeError(f"{kind} archive SHA mismatch")
    if len(blob)!=int(arc["bytes"]):raise RuntimeError(f"{kind} archive byte length mismatch")
    member=arc["member"]["name"]
    nominal=float(arc["member"]["published_resolution_m"])
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        names=zf.namelist()
        if member not in names:raise RuntimeError(f"{kind} member missing")
        zi=zf.getinfo(member)
        if int(zi.file_size)!=int(arc["member"]["uncompressed_bytes"]):raise RuntimeError(f"{kind} member size mismatch")
        if int(zi.CRC)!=int(arc["member"]["crc32"]):raise RuntimeError(f"{kind} member CRC mismatch")
        tif=zf.read(member)
        stem=member.rsplit(".",1)[0].lower()
        tfw_names=[n for n in names if n.lower()==stem+".tfw" or n.lower()==stem+".tifw" or n.lower()==stem+".wld"]
        if len(tfw_names)>1:raise RuntimeError(f"{kind}: multiple exact world files")
        tfw=zf.read(tfw_names[0]).decode("utf-8","replace") if tfw_names else None
    with MemoryFile(tif) as mf:
        with mf.open() as src:
            ma=src.read(1,masked=True)
            data=np.asarray(ma.filled(np.nan),dtype=np.float32)
            valid=(~np.ma.getmaskarray(ma)) & np.isfinite(data)
            emb=src.transform
            crs=str(src.crs) if src.crs is not None else None
            nodata=src.nodata
            dtype=str(src.dtypes[0])
    if crs is None:raise RuntimeError(f"{kind}: CRS missing")
    src_crs=CRS.from_user_input(crs)
    wvals=None
    if not emb.is_identity:
        transform=emb;geo="EMBEDDED_TRANSFORM"
    else:
        if tfw is None:raise RuntimeError(f"{kind}: identity transform and no exact world file")
        wvals,transform=world_affine(tfw);geo="WORLD_FILE"
    to_wgs=Transformer.from_crs(src_crs,WGS84,always_xy=True)
    from_wgs=Transformer.from_crs(WGS84,src_crs,always_xy=True)

    # Diagnostic physical dimensions only; never retune the frozen nominal cell.
    h,w=data.shape; rr=h//2;cc=w//2
    sx0,sy0=rasterio.transform.xy(transform,rr,cc,offset="center")
    sxe,sye=rasterio.transform.xy(transform,rr,cc+1,offset="center")
    sxn,syn=rasterio.transform.xy(transform,rr+1,cc,offset="center")
    lon0,lat0=to_wgs.transform(float(sx0),float(sy0))
    lone,late=to_wgs.transform(float(sxe),float(sye))
    lonn,latn=to_wgs.transform(float(sxn),float(syn))
    _,_,east_m=GEOD.inv(lon0,lat0,lone,late)
    _,_,north_m=GEOD.inv(lon0,lat0,lonn,latn)
    diag={"east_west":abs(float(east_m)),"north_south":abs(float(north_m))}
    if not (0.45*nominal <= diag["north_south"] <= 1.55*nominal and 0.45*nominal <= diag["east_west"] <= 1.55*nominal):
        raise RuntimeError(f"{kind}: diagnostic physical cell inconsistent with nominal {nominal} m: {diag}")

    return {
      "data":data,"valid":valid,"transform":transform,
      "crs":src_crs,"to_wgs":to_wgs,"from_wgs":from_wgs,
      "archive_sha256":sha,"archive_bytes":len(blob),"members":names,
      "tif_name":member,"tfw_name":tfw_names[0] if tfw_names else None,
      "dtype":dtype,"nodata":nodata,"embedded_transform":list(emb)[:6],
      "effective_transform":list(transform)[:6],"world_values":wvals,
      "georef_source":geo,"nominal_cell_m":nominal,
      "diagnostic_center_pixel_geodesic_m":diag
    }

def patch_geometry(radius_m,cell_m):
    rpx=max(1,int(round(radius_m/cell_m)))
    outer=2*rpx
    yy,xx=np.ogrid[-outer:outer+1,-outer:outer+1]
    d2=xx*xx+yy*yy
    disk=d2<=rpx*rpx
    ann=(d2>(1.25*rpx)**2)&(d2<=(2*rpx)**2)
    return rpx,outer,disk,ann,(disk|ann)

def rowcol_from_lonlat(r,lon,lat):
    x,y=r["from_wgs"].transform(float(lon),float(lat))
    row,col=rasterio.transform.rowcol(r["transform"],float(x),float(y))
    return int(row),int(col)

def lonlat_from_rowcol(r,row,col):
    x,y=rasterio.transform.xy(r["transform"],int(row),int(col),offset="center")
    lon,lat=r["to_wgs"].transform(float(x),float(y))
    return float(lon),float(lat)

def support_fraction(r,lon,lat,radius_m):
    row,col=rowcol_from_lonlat(r,lon,lat)
    _,outer,_,_,need=patch_geometry(radius_m,r["nominal_cell_m"])
    r0,r1=row-outer,row+outer+1;c0,c1=col-outer,col+outer+1
    if r0<0 or c0<0 or r1>r["valid"].shape[0] or c1>r["valid"].shape[1]:
        return 0.0,row,col
    v=r["valid"][r0:r1,c0:c1]
    return float(np.mean(v[need])) if np.sum(need) else 0.0,row,col

def backscatter_metrics(r,lon,lat,radius_m):
    row,col=rowcol_from_lonlat(r,lon,lat)
    _,outer,disk,ann,need=patch_geometry(radius_m,r["nominal_cell_m"])
    r0,r1=row-outer,row+outer+1;c0,c1=col-outer,col+outer+1
    if r0<0 or c0<0 or r1>r["data"].shape[0] or c1>r["data"].shape[1]:
        return None,0.0,row,col
    a=r["data"][r0:r1,c0:c1]
    v=r["valid"][r0:r1,c0:c1]
    mask_fraction=float(np.mean(v[need])) if np.sum(need) else 0.0
    if mask_fraction<BACK_MASK_MIN:return None,mask_fraction,row,col
    dv=a[disk&v].astype(np.float64);av=a[ann&v].astype(np.float64)
    if len(dv)<10 or len(av)<10:return None,mask_fraction,row,col
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
    if not all(np.isfinite(vals)):return None,mask_fraction,row,col
    return {
      "disk_standard_deviation":dstd,
      "absolute_center_vs_annulus_contrast_over_annulus_sd":contrast,
      "disk_q95_minus_q05_range":drange,
      "median_first_difference_edge_energy":edge
    },mask_fraction,row,col

# Candidate authority is verified before held-out backscatter is opened.
candj,artifact_sha,candidate_sha=get_candidate_json()
candidates=list(candj["frozen_candidates"])
if len(candidates)!=EXPECTED_CANDIDATES:raise RuntimeError("candidate list length mismatch")

bathy=load_archive_raster("bathymetry")
back=load_archive_raster("backscatter")

frozen_xy=np.array([[float(c["x_m"]),float(c["y_m"])] for c in candidates],dtype=float)
def far_from_candidates(mx,my):
    d2=(frozen_xy[:,0]-mx)**2+(frozen_xy[:,1]-my)**2
    return bool(np.min(d2)>=MIN_DIST**2)

candidate_results=[];by_radius=defaultdict(list)
for c in candidates:
    lon,lat=float(c["lon"]),float(c["lat"])
    radius=int(c["radius_m"])
    bsup,brow,bcol=support_fraction(bathy,lon,lat,radius)
    metrics=None;backsup=0.0;backrow=backcol=None
    if bsup>=BATHY_SUPPORT_MIN:
        metrics,backsup,backrow,backcol=backscatter_metrics(back,lon,lat,radius)
    rec={
      "candidate_id":c["candidate_id"],"blind_rank":c["blind_rank"],
      "x_m":float(c["x_m"]),"y_m":float(c["y_m"]),"lon":lon,"lat":lat,
      "radius_m":radius,
      "bathymetry_support_fraction":bsup,
      "backscatter_native_support_fraction":backsup,
      "bathymetry_row":brow,"bathymetry_col":bcol,
      "backscatter_row":backrow,"backscatter_col":backcol,
      "metrics":metrics
    }
    candidate_results.append(rec);by_radius[radius].append(rec)

strata={};blocked=[];all_candidate_scores=[];all_control_scores=[]
for radius in sorted(by_radius):
    key=str(radius);frozen=by_radius[radius]
    scored=[q for q in frozen if q["metrics"] is not None and q["bathymetry_support_fraction"]>=BATHY_SUPPORT_MIN and q["backscatter_native_support_fraction"]>=BACK_MASK_MIN]
    if not scored:
        strata[key]={"radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":0,
                     "status":"COVERAGE_BLOCKED","candidate_support_count":0,"active_metrics":[],"excluded_zero_mad_metrics":[]}
        blocked.append({"stratum":key,"reason":"COVERAGE_BLOCKED"});continue

    target=max(MIN_CONTROLS,len(scored)*CONTROLS_PER)
    controls=[];tries=0;h,w=bathy["valid"].shape
    margin=int(math.ceil(2*radius/bathy["nominal_cell_m"]))+1
    while len(controls)<target and tries<target*1500:
        tries+=1
        if h<=2*margin+1 or w<=2*margin+1:break
        brow=RNG.randint(margin,h-margin-1);bcol=RNG.randint(margin,w-margin-1)
        lon,lat=lonlat_from_rowcol(bathy,brow,bcol)
        mx,my=TO_METRIC.transform(lon,lat)
        if not far_from_candidates(float(mx),float(my)):continue
        bsup,_,_=support_fraction(bathy,lon,lat,radius)
        if bsup<BATHY_SUPPORT_MIN:continue
        m,backsup,_,_=backscatter_metrics(back,lon,lat,radius)
        if m is None or backsup<BACK_MASK_MIN:continue
        controls.append({"x_m":float(mx),"y_m":float(my),"lon":lon,"lat":lat,
                         "bathymetry_support_fraction":bsup,"backscatter_native_support_fraction":backsup,"metrics":m})

    if len(controls)<target:
        strata[key]={"radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":len(scored),
                     "control_count":len(controls),"required_control_count":target,"status":"CONTROL_SUPPORT_BLOCKED",
                     "candidate_support_count":0,"active_metrics":[],"excluded_zero_mad_metrics":[]}
        blocked.append({"stratum":key,"reason":"CONTROL_SUPPORT_BLOCKED"});continue

    control_mat=np.array([[q["metrics"][m] for m in METRICS] for q in controls],dtype=float)
    med=np.median(control_mat,axis=0);mad=np.median(np.abs(control_mat-med),axis=0)
    active_idx=[i for i,x in enumerate(mad) if float(x)>0.0]
    excluded_idx=[i for i,x in enumerate(mad) if float(x)==0.0]
    active=[METRICS[i] for i in active_idx];excluded=[METRICS[i] for i in excluded_idx]
    if len(active)<MIN_ACTIVE:
        strata[key]={"radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":len(scored),
          "control_count":len(controls),"status":"DEGENERATE_METRICS_BLOCKED","active_metrics":active,
          "excluded_zero_mad_metrics":excluded,
          "control_metric_median":{m:float(med[i]) for i,m in enumerate(METRICS)},
          "control_metric_mad":{m:float(mad[i]) for i,m in enumerate(METRICS)},"candidate_support_count":0}
        blocked.append({"stratum":key,"reason":"DEGENERATE_METRICS_BLOCKED"});continue

    meda=med[active_idx];mada=mad[active_idx];scale=1.4826*mada
    ca=control_mat[:,active_idx];cz=np.maximum(0.0,(ca-meda)/scale)
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
    all_candidate_scores.extend(float(x) for x in kscores);all_control_scores.extend(float(x) for x in cscores)
    strata[key]={"radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":len(scored),
      "unscored_candidate_count":len(frozen)-len(scored),"control_count":len(controls),"status":"SCORABLE",
      "control_metric_median":{m:float(med[i]) for i,m in enumerate(METRICS)},
      "control_metric_mad":{m:float(mad[i]) for i,m in enumerate(METRICS)},
      "active_metrics":active,"excluded_zero_mad_metrics":excluded,
      "control_score_q95":q95,"candidate_support_count":support,"_control_scores":[float(x) for x in cscores]}

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
            st=strata[key];n=int(st["scored_candidate_count"]);pool=st["_control_scores"]
            if n>len(pool):raise RuntimeError(f"{key}: permutation n > controls")
            sim.extend(RNG.sample(pool,n))
        if float(np.mean(sim))>=observed:ge+=1
    perm_p=(ge+1)/(ITERS+1)

if all_candidate_scores and all_control_scores:
    mw=mannwhitneyu(all_candidate_scores,all_control_scores,alternative="greater",method="asymptotic")
    mw_u=float(mw.statistic);mw_p=float(mw.pvalue)
else:mw_u=mw_p=None

support_count=sum(int(v.get("candidate_support_count",0)) for v in strata.values())
scored_count=sum(int(v.get("scored_candidate_count",0)) for v in strata.values() if v["status"]=="SCORABLE")
support_fraction=(support_count/scored_count) if scored_count else None
enrichment=(support_fraction/.05) if support_fraction is not None else None
primary_pass=bool(not blocked and perm_p is not None and perm_p<ALPHA and np.isfinite(observed) and np.isfinite(control_mean)
                  and observed>control_mean and enrichment is not None and enrichment>1.0)
for v in strata.values():v.pop("_control_scores",None)

if blocked:interpretation="PRIMARY_VALIDATION_BLOCKED_BY_SUPPORT_COVERAGE_OR_DEGENERATE_STRATUM"
elif primary_pass:interpretation="PASS_PREREGISTERED_UNSEEN_WESTERN_AUSTRALIAN_MARINE_PAIRED_SUPPORT_DOMAIN_V2_VALIDATION"
else:interpretation="NO_PREREGISTERED_UNSEEN_WESTERN_AUSTRALIAN_MARINE_PAIRED_SUPPORT_DOMAIN_V2_VALIDATION"

def meta(r):
    return {
      "archive_sha256":r["archive_sha256"],"archive_bytes":r["archive_bytes"],
      "tif_name":r["tif_name"],"tfw_name":r["tfw_name"],"shape":list(r["data"].shape),
      "valid_pixels":int(np.sum(r["valid"])),"dtype":r["dtype"],"nodata":r["nodata"],
      "nominal_cell_m":r["nominal_cell_m"],"diagnostic_center_pixel_geodesic_m":r["diagnostic_center_pixel_geodesic_m"],
      "georef_source":r["georef_source"],"crs":r["crs"].to_string(),
      "embedded_transform":r["embedded_transform"],"effective_transform":r["effective_transform"],
      "world_file_values":r["world_values"]
    }

out={
 "artifact_id":"JANUS-KUSTO-SOUTHWEST-CORNER-SWC2-UNSEEN-BACKSCATTER-V2-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],"implementation_freeze":IMP["artifact_id"],
 "candidate_authority":{"receipt":CAND_AUTH["artifact_id"],"artifact_id":ARTIFACT_ID,
   "artifact_zip_sha256_verified":artifact_sha,"candidate_json_sha256_verified":candidate_sha,
   "frozen_candidate_count":len(candidates)},
 "paired_inputs":{"bathymetry":meta(bathy),"backscatter":meta(back)},
 "support_domain_v2":{"bathymetry_support_min_fraction":BATHY_SUPPORT_MIN,
   "backscatter_native_mask_min_fraction":BACK_MASK_MIN,"candidate_and_control_rule_identical":True,
   "backscatter_value_based_support_inference_used":False},
 "georeference":{"metric_candidate_crs":METRIC_CRS.to_string(),
   "candidate_native_lookup":"frozen WGS84 lon/lat transformed independently into each raster CRS"},
 "strata":strata,"candidate_results":candidate_results,
 "coverage":{"frozen_candidate_count":len(candidates),"scorable_candidate_count":scored_count,"blocked_strata":blocked},
 "population":{"candidate_mean_aggregate_score":observed if np.isfinite(observed) else None,
   "control_mean_aggregate_score":control_mean if np.isfinite(control_mean) else None,
   "candidate_support_count":support_count,"candidate_support_fraction":support_fraction,
   "nominal_control_tail_fraction":0.05,"support_enrichment_ratio":enrichment,
   "stratified_permutation_iterations":ITERS,"stratified_permutation_p":perm_p,
   "mann_whitney_u":mw_u,"mann_whitney_one_sided_p":mw_p,
   "primary_pass":primary_pass,"scientific_promotion":primary_pass},
 "truth_firewall":{"underwater_imagery_read":False,"benthic_interpretation_read":False,
   "candidate_coordinates_changed":False,"candidate_radii_changed":False},
 "interpretation":interpretation,
 "claim_ceiling":"PREREGISTERED_UNSEEN_WESTERN_AUSTRALIAN_MARINE_SAME_SURVEY_CROSSCHANNEL_V2_VALIDATION_ONLY__NO_VISUAL_OR_INDEPENDENT_SURVEY_REPLICATION"
}
p=OUT/"JANUS-KUSTO-SOUTHWEST-CORNER-SWC2-UNSEEN-BACKSCATTER-V2-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"candidate_hash_verified":candidate_sha,
 "frozen_candidates":len(candidates),"scorable_candidates":scored_count,"blocked_strata":blocked,
 "support_count":support_count,"support_fraction":support_fraction,"support_enrichment_ratio":enrichment,
 "candidate_mean_score":out["population"]["candidate_mean_aggregate_score"],
 "control_mean_score":out["population"]["control_mean_aggregate_score"],
 "permutation_p":perm_p,"mann_whitney_p":mw_p,
 "primary_pass":primary_pass,"scientific_promotion":primary_pass,
 "interpretation":interpretation,"output_json_sha256":hashlib.sha256(raw.encode()).hexdigest()
},indent=2))
