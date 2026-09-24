#!/usr/bin/env python3
from __future__ import annotations

import hashlib, io, json, math, os, random, tempfile, zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import requests
import rasterio
from pyproj import Transformer, CRS
from scipy.stats import mannwhitneyu

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
IMP=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1H-CB13_03-UNCHANGED-V2-IMPLEMENTATION-FREEZE-2026-09-25-v1.0.json").read_text())
BR=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1G-CB13_03-WGS84-BACKSCATTER-NATIVE-METRIC-BRIDGE-2026-09-25-v1.0.json").read_text())
FR={"selected":{}}
for _k,_src in [("bathymetry",IMP["inputs"]["bathymetry"]),("heldout_backscatter",IMP["inputs"]["heldout_backscatter"])]:
    FR["selected"][_k]={"archive_url":_src["archive_url"],"member":_src["member"],"crc32":_src["member_crc32"],"uncompressed_size":_src["member_bytes"]}
GEO={"bathymetry":IMP["inputs"]["bathymetry"],"heldout_backscatter":IMP["inputs"]["heldout_backscatter"]}
AUTH=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1F-CB13_03-BLIND-BATHYMETRY-RECEIPT-2026-09-25-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-INFOMAR-I1H-CB13_03-v2/1.0","Accept-Encoding":"identity"}
ARTIFACT_ID=int(AUTH["run"]["artifact_id"])
EXPECTED_ARTIFACT_SHA=AUTH["run"]["artifact_zip_sha256"]
EXPECTED_JSON_SHA=AUTH["run"]["candidate_json_sha256"]
EXPECTED_CANDIDATES=int(AUTH["blind_execution"]["frozen_candidate_count"])

CELL=float(IMP["inputs"]["bathymetry"]["native_cell_m"])
DX=float(BR["solution"]["reference_pixel_axis_metres"][0])
DY=float(BR["solution"]["reference_pixel_axis_metres"][1])
BATH_SUPPORT_MIN=.80
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
to_utm=Transformer.from_crs("EPSG:4326","EPSG:32629",always_xy=True)

def get(url,headers=None,timeout=600,stream=False):
    h=dict(UA)
    if headers: h.update(headers)
    r=requests.get(url,headers=h,timeout=timeout,allow_redirects=True,stream=stream)
    r.raise_for_status()
    return r

def download(url,path,headers=None):
    h=hashlib.sha256(); n=0
    with get(url,headers=headers,stream=True) as r:
        with open(path,"wb") as f:
            for b in r.iter_content(1024*1024):
                if not b: continue
                f.write(b); h.update(b); n+=len(b)
    return h.hexdigest(),n

def candidate_json():
    token=os.environ.get("GITHUB_TOKEN")
    repo=os.environ.get("GITHUB_REPOSITORY","hawkar-usls/Janus-Cosmos")
    if not token: raise RuntimeError("GITHUB_TOKEN required")
    url=f"https://api.github.com/repos/{repo}/actions/artifacts/{ARTIFACT_ID}/zip"
    raw=get(url,headers={
      "Authorization":f"Bearer {token}",
      "Accept":"application/vnd.github+json",
      "X-GitHub-Api-Version":"2022-11-28"
    }).content
    zsha=hashlib.sha256(raw).hexdigest()
    if zsha!=EXPECTED_ARTIFACT_SHA: raise RuntimeError(f"candidate artifact SHA mismatch {zsha}")
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        js=[n for n in zf.namelist() if n.lower().endswith(".json") and not n.startswith("__MACOSX/")]
        if len(js)!=1: raise RuntimeError(f"expected one candidate JSON, found {js}")
        jraw=zf.read(js[0])
    jsha=hashlib.sha256(jraw).hexdigest()
    if jsha!=EXPECTED_JSON_SHA: raise RuntimeError(f"candidate JSON SHA mismatch {jsha}")
    j=json.loads(jraw)
    if int(j["frozen_candidate_count"])!=EXPECTED_CANDIDATES: raise RuntimeError("candidate count mismatch")
    return j,zsha,jsha

def extract_verified(label,spec,expected,tmp):
    import shutil, subprocess
    zpath=tmp/f"{label}.zip"; outdir=tmp/f"{label}_extract"; outdir.mkdir()
    asha,abytes=download(spec["archive_url"],zpath)
    if asha.lower()!=expected["archive_sha256"].lower(): raise RuntimeError(f"{label}: archive SHA drift")
    if abytes!=int(expected["archive_bytes"]): raise RuntimeError(f"{label}: archive size drift")
    seven=shutil.which("7z")
    if not seven: raise RuntimeError("7z required for Deflate64 source archives")
    subprocess.run([seven,"e","-y",f"-o{outdir}",str(zpath),spec["member"]],check=True,stdout=subprocess.DEVNULL)
    tpath=outdir/Path(spec["member"]).name
    msha=hashlib.sha256(); mbytes=0
    with open(tpath,"rb") as src:
        while True:
            b=src.read(8*1024*1024)
            if not b: break
            msha.update(b); mbytes+=len(b)
    if mbytes!=int(spec["uncompressed_size"]): raise RuntimeError(f"{label}: member size drift")
    if msha.hexdigest().lower()!=expected["member_sha256"].lower(): raise RuntimeError(f"{label}: member SHA drift")
    return tpath,{"archive_sha256":asha,"archive_bytes":abytes,"member":spec["member"],"member_sha256":msha.hexdigest(),"member_bytes":mbytes}

def bathy_geometry(radius):
    rpx=max(1,int(round(float(radius)/CELL))); outer=2*rpx
    yy,xx=np.ogrid[-outer:outer+1,-outer:outer+1]
    d2=xx*xx+yy*yy
    disk=d2<=rpx*rpx
    ann=(d2>(1.25*rpx)**2)&(d2<=(2*rpx)**2)
    return outer,disk,ann,(disk|ann)

def back_geometry(radius):
    ox=int(math.ceil(2.0*float(radius)/DX))
    oy=int(math.ceil(2.0*float(radius)/DY))
    rr=np.arange(-oy,oy+1,dtype=np.float64)[:,None]
    cc=np.arange(-ox,ox+1,dtype=np.float64)[None,:]
    d=np.sqrt((cc*DX)**2+(rr*DY)**2)
    disk=d<=float(radius)
    ann=(d>1.25*float(radius))&(d<=2.0*float(radius))
    return oy,ox,disk,ann,(disk|ann)

def horizontal_epsg(crs):
    c=CRS.from_user_input(crs)
    e=c.to_epsg()
    if e: return e
    subs=c.sub_crs_list
    return subs[0].to_epsg() if subs else None

candj,artifact_sha,candidate_sha=candidate_json()
candidates=list(candj["frozen_candidates"])
frozen_xy=np.array([[float(c["x_m"]),float(c["y_m"])] for c in candidates],dtype=np.float64)

with tempfile.TemporaryDirectory(prefix="kusto_infomar_i1h_cb1303_") as td:
    tmp=Path(td)
    bath_path,bmeta=extract_verified("bathymetry",FR["selected"]["bathymetry"],GEO["bathymetry"],tmp)
    back_path,kmeta=extract_verified("backscatter",FR["selected"]["heldout_backscatter"],GEO["heldout_backscatter"],tmp)

    with rasterio.open(bath_path) as ds:
        if str(ds.crs)!="EPSG:32629": raise RuntimeError(f"unexpected bathy CRS {ds.crs}")
        ma=ds.read(1,masked=True)
        bath=np.asarray(ma.data,dtype=np.float32)
        bvalid=(~np.ma.getmaskarray(ma))
        bvalid &= np.isfinite(bath)
        btransform=ds.transform
        bmeta.update({"shape":[int(ds.height),int(ds.width)],"crs":str(ds.crs),"transform":list(ds.transform)[:6],"valid_pixels":int(np.sum(bvalid))})
        del ma

    with rasterio.open(back_path) as ds:
        if horizontal_epsg(ds.crs)!=4326: raise RuntimeError(f"unexpected backscatter horizontal CRS {ds.crs}")
        ma=ds.read(1,masked=True)
        back=np.asarray(ma.data,dtype=np.float32)
        kvalid=(~np.ma.getmaskarray(ma))
        kvalid &= np.isfinite(back)
        ktransform=ds.transform
        kmeta.update({"shape":[int(ds.height),int(ds.width)],"crs":str(ds.crs),"transform":list(ds.transform)[:6],
                      "valid_pixels":int(np.sum(kvalid)),"explicit_nodata":None if ds.nodata is None else float(ds.nodata)})
        del ma

    geom_b={float(r):bathy_geometry(float(r)) for r in sorted({float(c["radius_m"]) for c in candidates})}
    geom_k={float(r):back_geometry(float(r)) for r in sorted({float(c["radius_m"]) for c in candidates})}

    def bathy_support(x,y,radius):
        row,col=rasterio.transform.rowcol(btransform,float(x),float(y))
        outer,disk,ann,need=geom_b[float(radius)]
        r0,r1=row-outer,row+outer+1; c0,c1=col-outer,col+outer+1
        if r0<0 or c0<0 or r1>bvalid.shape[0] or c1>bvalid.shape[1]: return 0.0
        v=bvalid[r0:r1,c0:c1]
        return float(np.mean(v[need])) if np.any(need) else 0.0

    def back_metrics(row,col,radius):
        oy,ox,disk,ann,need=geom_k[float(radius)]
        r0,r1=row-oy,row+oy+1; c0,c1=col-ox,col+ox+1
        if r0<0 or c0<0 or r1>back.shape[0] or c1>back.shape[1]: return None,0.0
        a=back[r0:r1,c0:c1]; v=kvalid[r0:r1,c0:c1]
        frac=float(np.mean(v[need])) if np.any(need) else 0.0
        if frac<BACK_SUPPORT_MIN: return None,frac
        dv=a[disk&v].astype(np.float64); av=a[ann&v].astype(np.float64)
        if len(dv)<10 or len(av)<10: return None,frac
        dstd=float(np.std(dv))
        contrast=float(abs(np.mean(dv)-np.mean(av))/(np.std(av)+1e-6))
        drange=float(np.quantile(dv,.95)-np.quantile(dv,.05))
        gx=np.abs(np.diff(a.astype(np.float64),axis=1))
        gy=np.abs(np.diff(a.astype(np.float64),axis=0))
        vx=v[:,:-1]&v[:,1:]; vy=v[:-1,:]&v[1:,:]
        maskx=disk[:,:-1]&disk[:,1:]&vx
        masky=disk[:-1,:]&disk[1:,:]&vy
        edges=np.concatenate([gx[maskx],gy[masky]])
        edge=float(np.median(edges)) if len(edges) else float("nan")
        vals=[dstd,contrast,drange,edge]
        if not all(np.isfinite(vals)): return None,frac
        return {
          "disk_standard_deviation":dstd,
          "absolute_center_vs_annulus_contrast_over_annulus_sd":contrast,
          "disk_q95_minus_q05_range":drange,
          "median_first_difference_edge_energy":edge
        },frac

    def far(x,y):
        d2=(frozen_xy[:,0]-x)**2+(frozen_xy[:,1]-y)**2
        return bool(np.min(d2)>=MIN_DIST**2)

    candidate_results=[]; by_radius=defaultdict(list)
    for c in candidates:
        radius=float(c["radius_m"])
        x,y=float(c["x_m"]),float(c["y_m"])
        bsup=bathy_support(x,y,radius)
        row,col=rasterio.transform.rowcol(ktransform,float(c["lon"]),float(c["lat"]))
        metrics=None; ksup=0.0
        if bsup>=BATH_SUPPORT_MIN:
            metrics,ksup=back_metrics(int(row),int(col),radius)
        rec={
          "candidate_id":c["candidate_id"],"blind_rank":int(c["blind_rank"]),
          "x_m":x,"y_m":y,"lon":float(c["lon"]),"lat":float(c["lat"]),
          "radius_m":radius,"radius_cells":int(c["radius_cells"]),
          "bathymetry_support_fraction":bsup,
          "backscatter_native_support_fraction":ksup,
          "backscatter_row":int(row),"backscatter_col":int(col),
          "metrics":metrics
        }
        candidate_results.append(rec); by_radius[radius].append(rec)

    strata={}; blocked=[]; all_candidate_scores=[]; all_control_scores=[]
    h,w=back.shape
    for radius in sorted(by_radius):
        key=f"{radius:.12g}"
        frozen=by_radius[radius]
        scored=[q for q in frozen if q["metrics"] is not None and q["bathymetry_support_fraction"]>=BATH_SUPPORT_MIN and q["backscatter_native_support_fraction"]>=BACK_SUPPORT_MIN]
        if not scored:
            strata[key]={"radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":0,"status":"COVERAGE_BLOCKED","candidate_support_count":0,"active_metrics":[],"excluded_zero_mad_metrics":[]}
            blocked.append({"stratum":key,"reason":"COVERAGE_BLOCKED"}); continue

        target=max(MIN_CONTROLS,len(scored)*CONTROLS_PER)
        oy,ox,_,_,_=geom_k[radius]
        controls=[]; tries=0
        while len(controls)<target and tries<target*3000:
            tries+=1
            row=RNG.randint(oy,h-oy-1); col=RNG.randint(ox,w-ox-1)
            lon,lat=rasterio.transform.xy(ktransform,row,col,offset="center")
            x,y=to_utm.transform(float(lon),float(lat)); x=float(x); y=float(y)
            if not far(x,y): continue
            bsup=bathy_support(x,y,radius)
            if bsup<BATH_SUPPORT_MIN: continue
            metrics,ksup=back_metrics(row,col,radius)
            if metrics is None or ksup<BACK_SUPPORT_MIN: continue
            controls.append({"x_m":x,"y_m":y,"row":int(row),"col":int(col),"bathymetry_support_fraction":bsup,"backscatter_native_support_fraction":ksup,"metrics":metrics})

        if len(controls)<target:
            strata[key]={"radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":len(scored),
                         "control_count":len(controls),"required_control_count":target,"control_attempts":tries,
                         "status":"CONTROL_SUPPORT_BLOCKED","candidate_support_count":0,"active_metrics":[],"excluded_zero_mad_metrics":[]}
            blocked.append({"stratum":key,"reason":"CONTROL_SUPPORT_BLOCKED"}); continue

        cm=np.array([[q["metrics"][m] for m in METRICS] for q in controls],dtype=np.float64)
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
            blocked.append({"stratum":key,"reason":"DEGENERATE_METRICS_BLOCKED"}); continue

        meda=med[active_idx]; mada=mad[active_idx]; scale=1.4826*mada
        ca=cm[:,active_idx]; cz=np.maximum(0.0,(ca-meda)/scale)
        cscores=.5*np.max(cz,axis=1)+.5*np.mean(cz,axis=1)
        q95=float(np.quantile(cscores,.95))
        km=np.array([[q["metrics"][m] for m in METRICS] for q in scored],dtype=np.float64)[:,active_idx]
        kz=np.maximum(0.0,(km-meda)/scale)
        kscores=.5*np.max(kz,axis=1)+.5*np.mean(kz,axis=1)

        support=0
        for q,zv,score in zip(scored,kz,kscores):
            q["active_metrics"]=active; q["excluded_zero_mad_metrics"]=excluded
            q["metric_positive_robust_z"]={m:float(zv[i]) for i,m in enumerate(active)}
            q["aggregate_score"]=float(score); q["same_stratum_control_q95"]=q95
            q["cross_channel_support"]=bool(score>=q95); support+=int(score>=q95)

        all_candidate_scores.extend(float(x) for x in kscores)
        all_control_scores.extend(float(x) for x in cscores)
        strata[key]={
          "radius_m":radius,"frozen_candidate_count":len(frozen),"scored_candidate_count":len(scored),
          "unscored_candidate_count":len(frozen)-len(scored),"control_count":len(controls),"control_attempts":tries,
          "status":"SCORABLE",
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
                s=strata[key]; n=int(s["scored_candidate_count"]); pool=s["_control_scores"]
                sim.extend(RNG.sample(pool,n))
            if float(np.mean(sim))>=observed: ge+=1
        perm_p=(ge+1)/(ITERS+1)

    if all_candidate_scores and all_control_scores:
        mw=mannwhitneyu(all_candidate_scores,all_control_scores,alternative="greater",method="asymptotic")
        mw_u=float(mw.statistic); mw_p=float(mw.pvalue)
    else:
        mw_u=mw_p=None

    support_count=sum(int(v.get("candidate_support_count",0)) for v in strata.values())
    scored_count=sum(int(v.get("scored_candidate_count",0)) for v in strata.values() if v["status"]=="SCORABLE")
    support_fraction=support_count/scored_count if scored_count else None
    enrichment=support_fraction/.05 if support_fraction is not None else None
    primary_pass=bool(
      not blocked and perm_p is not None and perm_p<ALPHA and
      np.isfinite(observed) and np.isfinite(control_mean) and observed>control_mean and
      enrichment is not None and enrichment>1.0
    )
    for v in strata.values(): v.pop("_control_scores",None)

    status=("PASS_PREREGISTERED_UNSEEN_INFOMAR_CB13_03_PAIRED_SUPPORT_DOMAIN_V2_VALIDATION"
            if primary_pass else
            "PRIMARY_VALIDATION_BLOCKED_BY_SUPPORT_COVERAGE_OR_DEGENERATE_STRATUM"
            if blocked else
            "NO_PREREGISTERED_UNSEEN_INFOMAR_CB13_03_PAIRED_SUPPORT_DOMAIN_V2_VALIDATION")

    out={
      "artifact_id":"JANUS-KUSTO-INFOMAR-I1H-CB13_03-UNSEEN-BACKSCATTER-V2-RUN-2026-09-25-v1.0",
      "implementation_freeze":IMP["artifact_id"],"representation_bridge":BR["artifact_id"],
      "candidate_authority":{
        "receipt":AUTH["artifact_id"],"artifact_id":ARTIFACT_ID,
        "artifact_zip_sha256_verified":artifact_sha,"candidate_json_sha256_verified":candidate_sha,
        "frozen_candidate_count":len(candidates)
      },
      "paired_inputs":{"bathymetry":bmeta,"backscatter":kmeta},
      "support_domain_v2":{
        "bathymetry_support_min_fraction":BATH_SUPPORT_MIN,
        "backscatter_native_mask_min_fraction":BACK_SUPPORT_MIN,
        "candidate_and_control_rule_identical":True,
        "backscatter_resampled":False,
        "backscatter_value_based_support_inference_used":False,
        "native_wgs84_physical_patch_axis_metres":[DX,DY]
      },
      "strata":strata,"candidate_results":candidate_results,
      "coverage":{"frozen_candidate_count":len(candidates),"scorable_candidate_count":scored_count,"blocked_strata":blocked},
      "population":{
        "candidate_mean_aggregate_score":observed if np.isfinite(observed) else None,
        "control_mean_aggregate_score":control_mean if np.isfinite(control_mean) else None,
        "candidate_support_count":support_count,"candidate_support_fraction":support_fraction,
        "nominal_control_tail_fraction":.05,"support_enrichment_ratio":enrichment,
        "stratified_permutation_iterations":ITERS,"stratified_permutation_p":perm_p,
        "mann_whitney_u":mw_u,"mann_whitney_one_sided_p":mw_p,
        "primary_pass":primary_pass,"scientific_promotion":primary_pass
      },
      "truth_firewall":{
        "candidate_coordinates_changed_after_freeze":False,
        "candidate_radii_changed_after_freeze":False,
        "groundtruth_or_object_labels_read":False,
        "backscatter_values_first_read_only_after_implementation_freeze":True
      },
      "status":status,
      "does_not_establish":[
        "OBJECT_IDENTITY","ARTIFICIALITY","ARCHAEOLOGICAL_OR_WRECK_CLASSIFICATION",
        "GLOBAL_FALSE_DISCOVERY_RATE","INDEPENDENT_SURVEY_REPLICATION","DIRECT_VISUAL_CONFIRMATION"
      ],
      "claim_ceiling":"PREREGISTERED_UNSEEN_INFOMAR_CB13_03_SAME_SURVEY_CROSSCHANNEL_V2_VALIDATION"
    }
    p=OUT/"JANUS-KUSTO-INFOMAR-I0L-CB13_02-UNSEEN-BACKSCATTER-V2-RUN-2026-09-24-v1.0.json"
    raw=json.dumps(out,indent=2); p.write_text(raw)
    print(json.dumps({
      "artifact_id":out["artifact_id"],"status":status,
      "frozen_candidates":len(candidates),"scorable_candidates":scored_count,
      "blocked_strata":blocked,"candidate_support_count":support_count,
      "candidate_support_fraction":support_fraction,"support_enrichment_ratio":enrichment,
      "candidate_mean_aggregate_score":out["population"]["candidate_mean_aggregate_score"],
      "control_mean_aggregate_score":out["population"]["control_mean_aggregate_score"],
      "stratified_permutation_p":perm_p,"mann_whitney_one_sided_p":mw_p,
      "primary_pass":primary_pass,
      "output_json_sha256":hashlib.sha256(raw.encode()).hexdigest()
    },indent=2))
