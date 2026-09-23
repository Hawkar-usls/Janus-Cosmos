#!/usr/bin/env python3
from __future__ import annotations

import gc, hashlib, io, json, math, os, random, zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import requests
import rasterio
from rasterio.io import MemoryFile
from rasterio.windows import from_bounds
from rasterio.warp import reproject, Resampling, transform_bounds
from affine import Affine
from scipy import ndimage
from scipy.stats import mannwhitneyu

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
TMP=OUT/"em2_tmp";TMP.mkdir(parents=True,exist_ok=True)

DUAL=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EM0C-DUAL-METHOD-FREEZE-2026-09-24-v1.0.json").read_text())
MEM=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EM0B-ARCHIVE-MEMBER-RECEIPT-2026-09-24-v1.0.json").read_text())
DOMAIN=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EM0E-PAIRED-DOMAIN-REPROJECTION-FREEZE-2026-09-24-v1.0.json").read_text())
CAND_AUTH=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EM1-BLIND-BATHYMETRY-RECEIPT-2026-09-24-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-ElizabethMiddleton-EM2-dual-v2-v3/1.0"}
ARTIFACT_ID=int(CAND_AUTH["run"]["artifact_id"])
EXPECTED_ARTIFACT_SHA=CAND_AUTH["run"]["artifact_zip_sha256"]
EXPECTED_JSON_SHA=CAND_AUTH["run"]["candidate_json_sha256"]
EXPECTED_CANDIDATES=int(CAND_AUTH["blind_execution"]["frozen_candidate_count"])

BATH_CELL=5.0
BACK_CELL=4.0
RADII=[20,40,80,160]
BUFFER_M=320.0
TARGET_CRS="EPSG:32757"
BATHY_SUPPORT_MIN=0.80
BACK_MASK_MIN=0.80
CONTROLS_PER=20
MIN_CONTROLS=200
MIN_DIST=100.0
MIN_ACTIVE=2
ITERS=5000
INTERNAL_ALPHA=0.05
PROMOTION_ALPHA=0.025
CONTROL_SEED=260923
RELATIVE_MAD_FLOOR=0.01
METRICS=[
 "disk_standard_deviation",
 "absolute_center_vs_annulus_contrast_over_annulus_sd",
 "disk_q95_minus_q05_range",
 "median_first_difference_edge_energy"
]

def get(url,headers=None,timeout=300,stream=False):
    h=dict(UA)
    if headers:h.update(headers)
    r=requests.get(url,headers=h,timeout=timeout,allow_redirects=True,stream=stream)
    r.raise_for_status();return r

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
        js=[n for n in zf.namelist() if n.lower().endswith(".json")]
        if len(js)!=1:raise RuntimeError(f"expected one candidate JSON {js}")
        raw=zf.read(js[0])
    jsha=hashlib.sha256(raw).hexdigest()
    if jsha!=EXPECTED_JSON_SHA:raise RuntimeError("candidate JSON SHA mismatch")
    j=json.loads(raw)
    if int(j["frozen_candidate_count"])!=EXPECTED_CANDIDATES:raise RuntimeError("candidate count mismatch")
    return j,zsha,jsha

def download_verified(spec,path):
    with get(spec["url"],stream=True) as r:
        h=hashlib.sha256();n=0
        with open(path,"wb") as f:
            for chunk in r.iter_content(1024*1024):
                if chunk:
                    f.write(chunk);h.update(chunk);n+=len(chunk)
    if n!=int(spec["bytes"]):raise RuntimeError(f"{path.name}: size mismatch")
    if h.hexdigest()!=spec["sha256"]:raise RuntimeError(f"{path.name}: sha mismatch")
    return h.hexdigest(),n

def snap_outward(bounds,cell):
    l,b,r,t=map(float,bounds)
    return math.floor(l/cell)*cell,math.floor(b/cell)*cell,math.ceil(r/cell)*cell,math.ceil(t/cell)*cell

def load_bathymetry_regions():
    zpath=TMP/"bathymetry.zip"
    sha,n=download_verified(MEM["bathymetry_archive"],zpath)
    spec=MEM["bathymetry_archive"]["tif_members"][0]
    member=spec["name"]
    with zipfile.ZipFile(zpath) as zf:
        zi=zf.getinfo(member)
        if int(zi.file_size)!=int(spec["uncompressed_bytes"]) or int(zi.CRC)!=int(spec["crc32"]):
            raise RuntimeError("bathymetry member authority mismatch")
    uri=f"zip://{zpath}!{member}"
    out={}
    with rasterio.open(uri) as src:
        for rg in DOMAIN["paired_study_domain"]["regions"]:
            rid=rg["id"];core=tuple(map(float,rg["core_bounds_m"]))
            expanded=(core[0]-BUFFER_M,core[1]-BUFFER_M,core[2]+BUFFER_M,core[3]+BUFFER_M)
            sb=transform_bounds(TARGET_CRS,src.crs,*expanded,densify_pts=21)
            win=from_bounds(*sb,transform=src.transform).round_offsets().round_lengths()
            col0=max(0,int(win.col_off));row0=max(0,int(win.row_off))
            col1=min(src.width,int(math.ceil(win.col_off+win.width)))
            row1=min(src.height,int(math.ceil(win.row_off+win.height)))
            win=rasterio.windows.Window(col0,row0,col1-col0,row1-row0)
            raw=src.read(1,window=win,masked=False)
            st=src.window_transform(win)
            snapped=snap_outward(expanded,BATH_CELL)
            width=int(round((snapped[2]-snapped[0])/BATH_CELL))
            height=int(round((snapped[3]-snapped[1])/BATH_CELL))
            tr=Affine(BATH_CELL,0,snapped[0],0,-BATH_CELL,snapped[3])
            dst=np.full((height,width),np.nan,dtype=np.float32)
            reproject(
              source=raw,destination=dst,src_transform=st,src_crs=src.crs,src_nodata=src.nodata,
              dst_transform=tr,dst_crs=TARGET_CRS,dst_nodata=np.nan,resampling=Resampling.nearest,num_threads=2
            )
            del raw;gc.collect()
            valid=np.isfinite(dst)
            out[rid]={"data":dst,"valid":valid,"transform":tr,"core_bounds":core}
    return out,{"url":MEM["bathymetry_archive"]["url"],"sha256":sha,"bytes":n,"member":member}

def decode_backscatter_member(zf,spec,rid):
    member=spec["name"];zi=zf.getinfo(member)
    if int(zi.file_size)!=int(spec["uncompressed_bytes"]) or int(zi.CRC)!=int(spec["crc32"]):
        raise RuntimeError(f"{rid}: backscatter member authority mismatch")
    tif=zf.read(member)
    with MemoryFile(tif) as mf:
        with mf.open() as src:
            ma=src.read(1,masked=True)
            data=np.asarray(ma.filled(np.nan),dtype=np.float32)
            valid=(~np.ma.getmaskarray(ma))&np.isfinite(data)
            tr=src.transform;crs=str(src.crs) if src.crs else None
            res=(abs(float(tr.a)),abs(float(tr.e)))
            if crs!=TARGET_CRS:raise RuntimeError(f"{rid}: unexpected backscatter CRS {crs}")
            if max(abs(res[0]-BACK_CELL),abs(res[1]-BACK_CELL))/BACK_CELL>0.001:
                raise RuntimeError(f"{rid}: backscatter resolution {res}")
    return {"data":data,"valid":valid,"transform":tr,"member":member,"resolution_m":list(res)}

def load_backscatter():
    zpath=TMP/"backscatter.zip"
    sha,n=download_verified(MEM["backscatter_archive"],zpath)
    specs={x["name"]:x for x in MEM["backscatter_archive"]["tif_members"]}
    names={
      "ELIZABETH":"32Bit_Geotifs/bksct_elizabeth_clip_20200331.tif",
      "MIDDLETON":"32Bit_Geotifs/bksct_middleton_clip_20200331.tif"
    }
    out={}
    with zipfile.ZipFile(zpath) as zf:
        for rid,name in names.items():
            out[rid]=decode_backscatter_member(zf,specs[name],rid)
    return out,{"url":MEM["backscatter_archive"]["url"],"sha256":sha,"bytes":n,"members":names}

def patch_geometry(radius_m,cell_m):
    rpx=max(1,int(round(radius_m/cell_m)));outer=2*rpx
    yy,xx=np.ogrid[-outer:outer+1,-outer:outer+1];d2=xx*xx+yy*yy
    disk=d2<=rpx*rpx
    ann=(d2>(1.25*rpx)**2)&(d2<=(2*rpx)**2)
    return rpx,outer,disk,ann,(disk|ann)

def support_fraction(valid,tr,x,y,radius_m,cell_m):
    row,col=rasterio.transform.rowcol(tr,float(x),float(y))
    _,outer,_,_,need=patch_geometry(radius_m,cell_m)
    r0,r1=row-outer,row+outer+1;c0,c1=col-outer,col+outer+1
    if r0<0 or c0<0 or r1>valid.shape[0] or c1>valid.shape[1]:return 0.0,int(row),int(col)
    v=valid[r0:r1,c0:c1]
    return float(np.mean(v[need])) if np.sum(need) else 0.0,int(row),int(col)

def back_metrics(region,x,y,radius_m):
    arr=region["data"];valid=region["valid"];tr=region["transform"]
    row,col=rasterio.transform.rowcol(tr,float(x),float(y))
    _,outer,disk,ann,need=patch_geometry(radius_m,BACK_CELL)
    r0,r1=row-outer,row+outer+1;c0,c1=col-outer,col+outer+1
    if r0<0 or c0<0 or r1>arr.shape[0] or c1>arr.shape[1]:return None,0.0,int(row),int(col)
    a=arr[r0:r1,c0:c1];v=valid[r0:r1,c0:c1]
    mf=float(np.mean(v[need])) if np.sum(need) else 0.0
    if mf<BACK_MASK_MIN:return None,mf,int(row),int(col)
    dv=a[disk&v].astype(np.float64);av=a[ann&v].astype(np.float64)
    if len(dv)<10 or len(av)<10:return None,mf,int(row),int(col)
    dstd=float(np.std(dv))
    contrast=float(abs(np.mean(dv)-np.mean(av))/(np.std(av)+1e-6))
    drange=float(np.quantile(dv,.95)-np.quantile(dv,.05))
    core=a.astype(np.float64)
    gx=np.abs(np.diff(core,axis=1));gy=np.abs(np.diff(core,axis=0))
    vx=v[:,:-1]&v[:,1:];vy=v[:-1,:]&v[1:,:]
    mx=disk[:,:-1]&disk[:,1:]&vx;my=disk[:-1,:]&disk[1:,:]&vy
    edges=np.concatenate([gx[mx],gy[my]])
    edge=float(np.median(edges)) if len(edges) else float("nan")
    vals=[dstd,contrast,drange,edge]
    if not all(np.isfinite(vals)):return None,mf,int(row),int(col)
    return {
      "disk_standard_deviation":dstd,
      "absolute_center_vs_annulus_contrast_over_annulus_sd":contrast,
      "disk_q95_minus_q05_range":drange,
      "median_first_difference_edge_energy":edge
    },mf,int(row),int(col)

candj,artifact_sha,candidate_sha=get_candidate_json()
candidates=list(candj["frozen_candidates"])
bathy,bathy_meta=load_bathymetry_regions()
back,back_meta=load_backscatter()

frozen_xy=np.array([[float(c["x_m"]),float(c["y_m"])] for c in candidates],dtype=float)
def far_from_candidates(x,y):
    return bool(np.min((frozen_xy[:,0]-x)**2+(frozen_xy[:,1]-y)**2)>=MIN_DIST**2)

candidate_base=[]
by_radius=defaultdict(list)
for c in candidates:
    rid=c["region_id"];x=float(c["x_m"]);y=float(c["y_m"]);radius=int(c["radius_m"])
    bsup,br,bc=support_fraction(bathy[rid]["valid"],bathy[rid]["transform"],x,y,radius,BATH_CELL)
    metrics=None;backsup=0.0;rr=cc=None
    if bsup>=BATHY_SUPPORT_MIN:
        metrics,backsup,rr,cc=back_metrics(back[rid],x,y,radius)
    rec={
      "candidate_id":c["candidate_id"],"blind_rank":c["blind_rank"],"region_id":rid,
      "x_m":x,"y_m":y,"radius_m":radius,
      "bathymetry_support_fraction":bsup,"backscatter_native_support_fraction":backsup,
      "metrics":metrics
    }
    candidate_base.append(rec);by_radius[radius].append(rec)

# Generate one frozen control set shared by v2/v3.
control_rng=random.Random(CONTROL_SEED)
controls_by_radius={}
control_generation={}
regions=list(bathy.keys())
area_weights={}
for rid in regions:
    core=bathy[rid]["core_bounds"]
    area_weights[rid]=max(1.0,(core[2]-core[0])*(core[3]-core[1]))
total_area=sum(area_weights.values())
cum=[];acc=0.0
for rid in regions:
    acc+=area_weights[rid]/total_area;cum.append((acc,rid))

def choose_region(rng):
    u=rng.random()
    for p,rid in cum:
        if u<=p:return rid
    return cum[-1][1]

for radius in RADII:
    scored=[q for q in by_radius[radius] if q["metrics"] is not None and q["bathymetry_support_fraction"]>=BATHY_SUPPORT_MIN and q["backscatter_native_support_fraction"]>=BACK_MASK_MIN]
    target=max(MIN_CONTROLS,len(scored)*CONTROLS_PER)
    controls=[];tries=0
    while len(controls)<target and tries<target*2500:
        tries+=1;rid=choose_region(control_rng)
        core=bathy[rid]["core_bounds"]
        x=control_rng.uniform(core[0],core[2]);y=control_rng.uniform(core[1],core[3])
        if not far_from_candidates(x,y):continue
        bsup,_,_=support_fraction(bathy[rid]["valid"],bathy[rid]["transform"],x,y,radius,BATH_CELL)
        if bsup<BATHY_SUPPORT_MIN:continue
        m,backsup,_,_=back_metrics(back[rid],x,y,radius)
        if m is None or backsup<BACK_MASK_MIN:continue
        controls.append({"region_id":rid,"x_m":x,"y_m":y,"bathymetry_support_fraction":bsup,
                         "backscatter_native_support_fraction":backsup,"metrics":m})
    controls_by_radius[radius]=controls
    control_generation[str(radius)]={"scored_candidate_count":len(scored),"required_controls":target,
                                     "control_count":len(controls),"tries":tries}

perm_state=control_rng.getstate()

def score_method(method_id,relative_guard):
    strata={};blocked=[];all_ks=[];all_cs=[];method_candidates=[]
    for base in candidate_base:
        method_candidates.append({k:v for k,v in base.items()})
    cmap={x["candidate_id"]:x for x in method_candidates}

    for radius in RADII:
        frozen=by_radius[radius]
        scored=[q for q in frozen if q["metrics"] is not None and q["bathymetry_support_fraction"]>=BATHY_SUPPORT_MIN and q["backscatter_native_support_fraction"]>=BACK_MASK_MIN]
        controls=controls_by_radius[radius]
        target=max(MIN_CONTROLS,len(scored)*CONTROLS_PER)
        key=str(radius)
        if not scored:
            strata[key]={"radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":0,
                         "status":"COVERAGE_BLOCKED","active_metrics":[],"candidate_support_count":0}
            blocked.append({"stratum":key,"reason":"COVERAGE_BLOCKED"});continue
        if len(controls)<target:
            strata[key]={"radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":len(scored),
                         "control_count":len(controls),"required_control_count":target,"status":"CONTROL_SUPPORT_BLOCKED",
                         "active_metrics":[],"candidate_support_count":0}
            blocked.append({"stratum":key,"reason":"CONTROL_SUPPORT_BLOCKED"});continue

        cm=np.array([[q["metrics"][m] for m in METRICS] for q in controls],dtype=float)
        med=np.median(cm,axis=0);mad=np.median(np.abs(cm-med),axis=0)
        relative=np.array([float(mad[i]/max(abs(med[i]),1e-12)) for i in range(len(METRICS))])
        if relative_guard:
            active_idx=[i for i in range(len(METRICS)) if float(mad[i])>0.0 and float(relative[i])>=RELATIVE_MAD_FLOOR]
        else:
            active_idx=[i for i in range(len(METRICS)) if float(mad[i])>0.0]
        inactive_idx=[i for i in range(len(METRICS)) if i not in active_idx]
        active=[METRICS[i] for i in active_idx];inactive=[METRICS[i] for i in inactive_idx]
        if len(active)<MIN_ACTIVE:
            strata[key]={"radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":len(scored),
                         "control_count":len(controls),"status":"DEGENERATE_METRICS_BLOCKED","active_metrics":active,
                         "excluded_metrics":inactive,"relative_mad":{METRICS[i]:float(relative[i]) for i in range(len(METRICS))},
                         "candidate_support_count":0}
            blocked.append({"stratum":key,"reason":"DEGENERATE_METRICS_BLOCKED"});continue

        meda=med[active_idx];mada=mad[active_idx];scale=1.4826*mada
        ca=cm[:,active_idx];cz=np.maximum(0.0,(ca-meda)/scale)
        cs=0.5*np.max(cz,axis=1)+0.5*np.mean(cz,axis=1)
        q95=float(np.quantile(cs,.95))

        km=np.array([[q["metrics"][m] for m in METRICS] for q in scored],dtype=float)[:,active_idx]
        kz=np.maximum(0.0,(km-meda)/scale)
        ks=0.5*np.max(kz,axis=1)+0.5*np.mean(kz,axis=1)
        support=0
        for q,zv,score in zip(scored,kz,ks):
            o=cmap[q["candidate_id"]]
            o["active_metrics"]=active;o["excluded_metrics"]=inactive
            o["metric_positive_robust_z"]={m:float(zv[i]) for i,m in enumerate(active)}
            o["aggregate_score"]=float(score);o["same_stratum_control_q95"]=q95
            o["cross_channel_support"]=bool(score>=q95);support+=int(score>=q95)
        all_ks.extend(float(x) for x in ks);all_cs.extend(float(x) for x in cs)
        strata[key]={
          "radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":len(scored),
          "unscored_candidate_count":len(frozen)-len(scored),"control_count":len(controls),"status":"SCORABLE",
          "control_metric_median":{m:float(med[i]) for i,m in enumerate(METRICS)},
          "control_metric_mad":{m:float(mad[i]) for i,m in enumerate(METRICS)},
          "control_metric_relative_mad":{m:float(relative[i]) for i,m in enumerate(METRICS)},
          "active_metrics":active,"excluded_metrics":inactive,"control_score_q95":q95,
          "candidate_support_count":support,"_control_scores":[float(x) for x in cs]
        }

    scorable=[k for k,v in strata.items() if v["status"]=="SCORABLE"]
    observed=float(np.mean(all_ks)) if all_ks else float("nan")
    control_mean=float(np.mean(all_cs)) if all_cs else float("nan")
    perm_p=None
    if not blocked and scorable:
        prng=random.Random();prng.setstate(perm_state)
        ge=0
        for _ in range(ITERS):
            sim=[]
            for key in scorable:
                st=strata[key];n=int(st["scored_candidate_count"]);pool=st["_control_scores"]
                if n>len(pool):raise RuntimeError("permutation n > controls")
                sim.extend(prng.sample(pool,n))
            if float(np.mean(sim))>=observed:ge+=1
        perm_p=(ge+1)/(ITERS+1)
    if all_ks and all_cs:
        mw=mannwhitneyu(all_ks,all_cs,alternative="greater",method="asymptotic")
        mw_u=float(mw.statistic);mw_p=float(mw.pvalue)
    else:mw_u=mw_p=None
    support_count=sum(int(v.get("candidate_support_count",0)) for v in strata.values())
    scored_count=sum(int(v.get("scored_candidate_count",0)) for v in strata.values() if v["status"]=="SCORABLE")
    support_fraction=support_count/scored_count if scored_count else None
    enrichment=support_fraction/.05 if support_fraction is not None else None
    internal_pass=bool(not blocked and perm_p is not None and perm_p<INTERNAL_ALPHA and
                       np.isfinite(observed) and np.isfinite(control_mean) and observed>control_mean and
                       enrichment is not None and enrichment>1.0)
    study_promotion=bool(internal_pass and perm_p<PROMOTION_ALPHA)
    for v in strata.values():v.pop("_control_scores",None)
    return {
      "method_id":method_id,"relative_mad_guard":relative_guard,"strata":strata,
      "candidate_results":method_candidates,
      "population":{
        "candidate_mean_aggregate_score":observed if np.isfinite(observed) else None,
        "control_mean_aggregate_score":control_mean if np.isfinite(control_mean) else None,
        "candidate_support_count":support_count,"candidate_support_fraction":support_fraction,
        "nominal_control_tail_fraction":0.05,"support_enrichment_ratio":enrichment,
        "stratified_permutation_iterations":ITERS,"stratified_permutation_p":perm_p,
        "mann_whitney_u":mw_u,"mann_whitney_one_sided_p":mw_p,
        "internal_primary_pass_at_0_05":internal_pass,
        "study_bonferroni_p_threshold":PROMOTION_ALPHA,
        "study_scientific_promotion":study_promotion
      },
      "blocked_strata":blocked
    }

v2=score_method("PAIRED_SUPPORT_DOMAIN_V2_UNCHANGED",False)
v3=score_method("PAIRED_SUPPORT_DOMAIN_V3_RELATIVE_MAD_GUARD",True)

out={
 "artifact_id":"JANUS-KUSTO-ELIZABETH-MIDDLETON-EM2-DUAL-V2-V3-UNSEEN-BACKSCATTER-RUN-2026-09-24-v1.0",
 "dual_method_freeze":DUAL["artifact_id"],
 "candidate_authority":{
   "receipt":CAND_AUTH["artifact_id"],"artifact_id":ARTIFACT_ID,
   "artifact_zip_sha256_verified":artifact_sha,"candidate_json_sha256_verified":candidate_sha,
   "frozen_candidate_count":len(candidates)
 },
 "inputs":{"bathymetry":bathy_meta,"backscatter":back_meta},
 "shared_control_generation":control_generation,
 "methods":{"v2":v2,"v3":v3},
 "truth_firewall":{"auv_imagery_read":False,"sediment_labels_read":False,"bruv_imagery_read":False,
                    "candidate_coordinates_changed":False,"candidate_radii_changed":False},
 "study_interpretation":{
   "v2_promotion":v2["population"]["study_scientific_promotion"],
   "v3_promotion":v3["population"]["study_scientific_promotion"],
   "familywise_alpha":0.05,"bonferroni_per_method":PROMOTION_ALPHA
 },
 "claim_ceiling":"DUAL_PREDECLARED_UNSEEN_CROSSCHANNEL_COMPARISON_ON_INDEPENDENT_EM2040C_SENSOR_LINEAGE__NO_VISUAL_GROUNDTRUTH_YET"
}
p=OUT/"JANUS-KUSTO-ELIZABETH-MIDDLETON-EM2-DUAL-V2-V3-UNSEEN-BACKSCATTER-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"candidate_hash_verified":candidate_sha,
 "controls":control_generation,
 "v2":v2["population"],"v2_blocked":v2["blocked_strata"],
 "v3":v3["population"],"v3_blocked":v3["blocked_strata"],
 "output_json_sha256":hashlib.sha256(raw.encode()).hexdigest()
},indent=2))
