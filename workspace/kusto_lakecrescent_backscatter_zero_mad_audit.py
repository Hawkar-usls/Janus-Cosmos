#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import io
import json
import math
import os
import random
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import requests
import rasterio
from affine import Affine
from scipy.stats import mannwhitneyu

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"
OUT.mkdir(parents=True,exist_ok=True)

PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-LAKE-CRESCENT-3M-BACKSCATTER-REPLAY-PREREG-2026-09-23-v1.0.json").read_text())
VALIDITY=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-LAKE-CRESCENT-3M-BACKSCATTER-VALIDITY-ADDENDUM-2026-09-23-v1.0.json").read_text())
ROBUSTNESS=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-LAKE-CRESCENT-BACKSCATTER-ZERO-MAD-ROBUSTNESS-AUDIT-2026-09-23-v1.0.json").read_text())
GEOREF_REPAIR=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-LAKE-CRESCENT-BACKSCATTER-GEOREF-REPAIR-ADDENDUM-2026-09-23-v1.0.json").read_text())

ARTIFACT_ID=int(PRE["candidate_authority"]["artifact_id"])
EXPECTED_ARTIFACT_SHA=PRE["candidate_authority"]["artifact_zip_sha256"]
EXPECTED_JSON_SHA=PRE["candidate_authority"]["candidate_json_sha256"]
BACK_URL="https://www.sciencebase.gov/catalog/file/get/5eea74d882ce3bd58d8572dd?name=LakeCrescent_backscatter_3m_UTM10_NAD83.zip"
UA={"User-Agent":"JANUS-KUSTO-LakeCrescent-backscatter-zeroMAD-audit/1.0"}
RNG=random.Random(int(PRE["control_design"]["rng_seed"]))

def get(url,headers=None,timeout=180):
    h=dict(UA)
    if headers:h.update(headers)
    r=requests.get(url,headers=h,timeout=timeout,allow_redirects=True)
    r.raise_for_status()
    return r

def get_frozen_candidate_json():
    token=os.environ.get("GITHUB_TOKEN")
    repo=os.environ.get("GITHUB_REPOSITORY","hawkar-usls/Janus-Cosmos")
    if not token:
        raise RuntimeError("GITHUB_TOKEN required")
    api=f"https://api.github.com/repos/{repo}/actions/artifacts/{ARTIFACT_ID}/zip"
    r=get(api,headers={
        "Authorization":f"Bearer {token}",
        "Accept":"application/vnd.github+json",
        "X-GitHub-Api-Version":"2022-11-28"
    })
    blob=r.content
    zsha=hashlib.sha256(blob).hexdigest()
    if zsha!=EXPECTED_ARTIFACT_SHA:
        raise RuntimeError(f"Frozen artifact ZIP sha mismatch {zsha} != {EXPECTED_ARTIFACT_SHA}")
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        js=[n for n in zf.namelist() if n.lower().endswith(".json") and not n.startswith("__MACOSX/")]
        if len(js)!=1:
            raise RuntimeError(f"Expected one frozen candidate JSON, found {js}")
        raw=zf.read(js[0])
    jsha=hashlib.sha256(raw).hexdigest()
    if jsha!=EXPECTED_JSON_SHA:
        raise RuntimeError(f"Frozen candidate JSON sha mismatch {jsha} != {EXPECTED_JSON_SHA}")
    j=json.loads(raw)
    if int(j["frozen_candidate_count"])!=int(PRE["candidate_authority"]["frozen_candidate_count"]):
        raise RuntimeError("Frozen candidate count mismatch")
    return j,zsha,jsha

def get_backscatter():
    blob=get(BACK_URL).content
    sha=hashlib.sha256(blob).hexdigest()
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        names=zf.namelist()
        tifs=[n for n in names if n.lower().endswith((".tif",".tiff")) and not n.startswith("__MACOSX/") and "/._" not in n]
        tfws=[n for n in names if n.lower().endswith((".tfw",".tifw",".wld")) and not n.startswith("__MACOSX/") and "/._" not in n]
        xmls=[n for n in names if n.lower().endswith(".xml") and not n.startswith("__MACOSX/") and "/._" not in n]
        if len(tifs)!=1:
            raise RuntimeError(f"Expected one real TIFF, found {tifs}")
        if len(tfws)!=1:
            raise RuntimeError(f"Expected one externally documented world file, found {tfws}")
        tif_bytes=zf.read(tifs[0])
        tfw_text=zf.read(tfws[0]).decode("utf-8",errors="replace")
    vals=[float(x.strip()) for x in tfw_text.splitlines() if x.strip()]
    if len(vals)!=6:
        raise RuntimeError(f"World file must contain six numeric lines, found {vals}")
    A,D,B,E,C,F=vals
    # World-file C/F are the center of the upper-left pixel; Rasterio Affine uses its outer corner.
    world_transform=Affine(A,B,C-0.5*A-0.5*B,D,E,F-0.5*D-0.5*E)
    tif_path=OUT/"_lakecrescent_backscatter_tmp.tif"
    tif_path.write_bytes(tif_bytes)
    return blob,sha,tifs[0],tif_path,names,tfws[0],vals,world_transform,xmls

def patch_metrics(arr,valid,row,col,rpx):
    outer=2*rpx
    r0=row-outer;r1=row+outer+1;c0=col-outer;c1=col+outer+1
    if r0<0 or c0<0 or r1>arr.shape[0] or c1>arr.shape[1]:
        return None
    a=arr[r0:r1,c0:c1]
    v=valid[r0:r1,c0:c1]
    yy,xx=np.ogrid[-outer:outer+1,-outer:outer+1]
    d2=xx*xx+yy*yy
    disk=d2<=rpx*rpx
    ann=(d2>(1.25*rpx)**2)&(d2<=(2*rpx)**2)
    need=disk|ann
    if float(np.mean(v[need]))<float(PRE["candidate_rules"]["minimum_patch_valid_fraction"]):
        return None
    dv=a[disk&v].astype(np.float64)
    av=a[ann&v].astype(np.float64)
    if len(dv)<10 or len(av)<10:
        return None
    dstd=float(np.std(dv))
    contrast=float(abs(np.mean(dv)-np.mean(av))/(np.std(av)+1e-6))
    drange=float(np.quantile(dv,0.95)-np.quantile(dv,0.05))

    # First-difference energy inside the disk; only valid adjacent pairs contribute.
    core=a.astype(np.float64)
    gx=np.abs(np.diff(core,axis=1)); gy=np.abs(np.diff(core,axis=0))
    vx=v[:,:-1]&v[:,1:]
    vy=v[:-1,:]&v[1:,:]
    maskx=disk[:,:-1]&disk[:,1:]&vx
    masky=disk[:-1,:]&disk[1:,:]&vy
    edges=np.concatenate([gx[maskx],gy[masky]])
    edge=float(np.median(edges)) if len(edges) else float("nan")
    return {
        "disk_standard_deviation":dstd,
        "absolute_center_vs_annulus_contrast_over_annulus_sd":contrast,
        "disk_q95_minus_q05_range":drange,
        "median_first_difference_edge_energy":edge
    }

candj,artifact_sha,candidate_sha=get_frozen_candidate_json()
candidates=candj["frozen_candidates"]

back_blob,back_sha,tif_name,tif_path,zip_members,tfw_name,tfw_values,world_transform,xml_names=get_backscatter()
with rasterio.open(tif_path) as src:
    raw=src.read(1,masked=True)
    data=np.asarray(raw.filled(np.nan),dtype=np.float32)
    native_valid=(~np.ma.getmaskarray(raw))
    valid=native_valid & np.isfinite(data) & (data>0)
    embedded_transform=src.transform
    transform=embedded_transform
    embedded_res=(abs(float(embedded_transform.a)),abs(float(embedded_transform.e)))
    use_world_file=(src.crs is None) or max(abs(embedded_res[0]-3.0),abs(embedded_res[1]-3.0))>1e-6
    if use_world_file:
        transform=world_transform
    crs=str(src.crs) if src.crs is not None else "EPSG:26910"
    res=(abs(float(transform.a)),abs(float(transform.e)))
    shape=[src.height,src.width]
    if max(abs(res[0]-3.0),abs(res[1]-3.0))>1e-6:
        raise RuntimeError(f"Unexpected backscatter cell size after documented world-file repair {res}")

    cand_records=[]
    by_radius_candidates=defaultdict(list)
    for c in candidates:
        row,col=rasterio.transform.rowcol(transform,float(c["x_utm10_nad83_m"]),float(c["y_utm10_nad83_m"]))
        r=int(round(float(c["radius_m"])/res[0]))
        m=patch_metrics(data,valid,row,col,r)
        rec={
            "candidate_id":c["candidate_id"],
            "x_utm10_nad83_m":c["x_utm10_nad83_m"],
            "y_utm10_nad83_m":c["y_utm10_nad83_m"],
            "radius_m":c["radius_m"],
            "row":int(row),"col":int(col),
            "metrics":m
        }
        cand_records.append(rec)
        if m is not None:
            by_radius_candidates[int(c["radius_m"])].append(rec)

    frozen_xy=np.array([[float(c["x_utm10_nad83_m"]),float(c["y_utm10_nad83_m"])] for c in candidates],dtype=np.float64)

    def far_from_candidates(x,y):
        d2=(frozen_xy[:,0]-x)**2+(frozen_xy[:,1]-y)**2
        return bool(np.min(d2)>=float(PRE["control_design"]["minimum_distance_from_any_frozen_candidate_m"])**2)

    controls=[]
    by_radius_controls=defaultdict(list)
    height,width=data.shape
    for radius in sorted(by_radius_candidates):
        target=max(
            int(PRE["control_design"]["minimum_controls_per_radius"]),
            len(by_radius_candidates[radius])*int(PRE["control_design"]["controls_per_candidate"])
        )
        rpx=int(round(radius/res[0]))
        outer=2*rpx
        tries=0
        while len(by_radius_controls[radius])<target and tries<target*500:
            tries+=1
            row=RNG.randint(outer,height-outer-1)
            col=RNG.randint(outer,width-outer-1)
            if not valid[row,col]:
                continue
            x,y=rasterio.transform.xy(transform,row,col,offset="center")
            if not far_from_candidates(float(x),float(y)):
                continue
            m=patch_metrics(data,valid,row,col,rpx)
            if m is None:
                continue
            rec={
                "control_id":f"R{radius}_CTRL_{len(by_radius_controls[radius])+1:04d}",
                "radius_m":radius,
                "x_utm10_nad83_m":float(x),
                "y_utm10_nad83_m":float(y),
                "row":int(row),"col":int(col),
                "metrics":m
            }
            by_radius_controls[radius].append(rec)
            controls.append(rec)
        if len(by_radius_controls[radius])<target:
            raise RuntimeError(f"Could not generate enough controls at radius {radius}: {len(by_radius_controls[radius])}/{target}")

metric_names=list(PRE["backscatter_metrics"])
radius_models={}
all_candidate_scores=[]
all_control_scores=[]
support_count=0
scored_candidate_count=0

for radius in sorted(by_radius_candidates):
    cs=by_radius_candidates[radius]
    ks=by_radius_controls[radius]
    control_mat=np.array([[q["metrics"][m] for m in metric_names] for q in ks],dtype=np.float64)
    med=np.median(control_mat,axis=0)
    mad=np.median(np.abs(control_mat-med),axis=0)
    active_idx=[i for i,mad_i in enumerate(mad) if float(mad_i)>0.0]
    excluded_idx=[i for i,mad_i in enumerate(mad) if float(mad_i)==0.0]
    if not active_idx:
        raise RuntimeError(f"No nondegenerate backscatter metrics remain at radius {radius}")
    active_names=[metric_names[i] for i in active_idx]
    excluded_names=[metric_names[i] for i in excluded_idx]
    med_active=med[active_idx]
    mad_active=mad[active_idx]
    scale=1.4826*mad_active
    control_active=control_mat[:,active_idx]
    cz=np.maximum(0.0,(control_active-med_active)/scale)
    cscore=0.5*np.max(cz,axis=1)+0.5*np.mean(cz,axis=1)
    q95=float(np.quantile(cscore,0.95))
    all_control_scores.extend(float(x) for x in cscore)

    cand_mat=np.array([[q["metrics"][m] for m in metric_names] for q in cs],dtype=np.float64)
    cand_active=cand_mat[:,active_idx]
    kz=np.maximum(0.0,(cand_active-med_active)/scale)
    kscore=0.5*np.max(kz,axis=1)+0.5*np.mean(kz,axis=1)

    for q,zv,sc in zip(cs,kz,kscore):
        q["metric_positive_robust_z"]={m:float(zv[i]) for i,m in enumerate(active_names)}
        q["excluded_zero_mad_metrics"]=excluded_names
        q["aggregate_score"]=float(sc)
        q["same_radius_control_q95"]=q95
        q["cross_channel_support"]=bool(sc>=q95)
        scored_candidate_count+=1
        support_count+=int(sc>=q95)
        all_candidate_scores.append(float(sc))

    radius_models[str(radius)]={
        "candidate_count":len(cs),
        "control_count":len(ks),
        "control_metric_median":{m:float(med[i]) for i,m in enumerate(metric_names)},
        "control_metric_mad":{m:float(mad[i]) for i,m in enumerate(metric_names)},
        "active_metrics":active_names,
        "excluded_zero_mad_metrics":excluded_names,
        "control_score_q95":q95,
        "candidate_support_count":sum(1 for q in cs if q.get("cross_channel_support"))
    }

# Radius-matched permutation on candidate mean score.
observed=float(np.mean(all_candidate_scores)) if all_candidate_scores else float("nan")
iters=int(PRE["population_test"]["iterations"])
ge=0
for _ in range(iters):
    sim=[]
    for radius in sorted(by_radius_candidates):
        n=len(by_radius_candidates[radius])
        pool=[q for q in by_radius_controls[radius]]
        # scores are reconstructed from radius model to avoid storing a second large table in output
        rm=radius_models[str(radius)]
        active=rm["active_metrics"]
        med=np.array([rm["control_metric_median"][m] for m in active])
        mad=np.array([rm["control_metric_mad"][m] for m in active])
        scale=1.4826*mad
        chosen=RNG.sample(pool,n)
        mat=np.array([[q["metrics"][m] for m in active] for q in chosen],dtype=np.float64)
        z=np.maximum(0.0,(mat-med)/scale)
        sc=0.5*np.max(z,axis=1)+0.5*np.mean(z,axis=1)
        sim.extend(float(x) for x in sc)
    if float(np.mean(sim))>=observed:
        ge+=1
perm_p=(ge+1)/(iters+1)

if all_candidate_scores and all_control_scores:
    mw=mannwhitneyu(all_candidate_scores,all_control_scores,alternative="greater",method="asymptotic")
    mw_u=float(mw.statistic); mw_p=float(mw.pvalue)
else:
    mw_u=float("nan");mw_p=float("nan")

support_fraction=(support_count/scored_candidate_count) if scored_candidate_count else float("nan")
control_tail=0.05
enrichment=(support_fraction/control_tail) if scored_candidate_count else float("nan")
robustness_signal=bool(scored_candidate_count>0 and perm_p<float(PRE["population_test"]["alpha"]) and observed>float(np.mean(all_control_scores)))

out={
    "artifact_id":"JANUS-KUSTO-LAKE-CRESCENT-BACKSCATTER-ZERO-MAD-ROBUSTNESS-RUN-2026-09-23-v1.0",
    "prereg":PRE["artifact_id"],
    "validity_addendum":VALIDITY["artifact_id"],
    "robustness_audit":ROBUSTNESS["artifact_id"],
    "georef_repair_addendum":GEOREF_REPAIR["artifact_id"],
    "candidate_authority":{
        "artifact_id":ARTIFACT_ID,
        "artifact_zip_sha256_verified":artifact_sha,
        "candidate_json_sha256_verified":candidate_sha,
        "candidate_count":len(candidates)
    },
    "backscatter_input":{
        "url":BACK_URL,
        "zip_bytes":len(back_blob),
        "zip_sha256":back_sha,
        "tif_name":tif_name,
        "zip_members":zip_members,
        "shape":shape,
        "crs":crs,
        "resolution_m":list(res),
        "world_file_name":tfw_name,
        "world_file_values":tfw_values,
        "world_file_applied":bool(use_world_file),
        "xml_members":xml_names,
        "embedded_transform":list(embedded_transform)[:6],
        "effective_transform":list(transform)[:6],
        "valid_pixels":int(np.sum(valid))
    },
    "coverage":{
        "frozen_candidate_count":len(candidates),
        "scored_candidate_count":scored_candidate_count,
        "unscored_candidate_count":len(candidates)-scored_candidate_count
    },
    "radius_models":radius_models,
    "candidate_results":cand_records,
    "population":{
        "candidate_mean_aggregate_score":observed,
        "control_mean_aggregate_score":float(np.mean(all_control_scores)) if all_control_scores else float("nan"),
        "candidate_support_count":support_count,
        "candidate_support_fraction":support_fraction,
        "nominal_control_tail_fraction":control_tail,
        "support_enrichment_ratio":enrichment,
        "radius_matched_permutation_iterations":iters,
        "radius_matched_permutation_p":perm_p,
        "mann_whitney_u":mw_u,
        "mann_whitney_one_sided_p":mw_p,
        "robustness_signal":robustness_signal,
        "scientific_promotion":False
    },
    "truth_firewall":{
        "geologic_feature_labels_read":False,
        "landslide_or_fault_coordinates_read":False,
        "candidate_coordinates_changed":False,
        "candidate_radii_changed":False
    },
    "interpretation":(
        "PASS_SAME_SURVEY_CROSS_CHANNEL_ENRICHMENT"
        if pass_primary else
        "NO_DEMONSTRATED_SAME_SURVEY_CROSS_CHANNEL_ENRICHMENT"
    ),
    "claim_ceiling":"SAME_SURVEY_CROSS_CHANNEL_ANOMALY_ENRICHMENT_ONLY__NOT_INDEPENDENT_SURVEY_VALIDATION"
}
p=OUT/"JANUS-KUSTO-LAKE-CRESCENT-BACKSCATTER-ZERO-MAD-ROBUSTNESS-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
    "artifact_id":out["artifact_id"],
    "candidate_hash_verified":candidate_sha,
    "backscatter_zip_sha256":back_sha,
    "shape":shape,
    "valid_pixels":out["backscatter_input"]["valid_pixels"],
    "scored_candidates":scored_candidate_count,
    "support_count":support_count,
    "support_fraction":support_fraction,
    "support_enrichment_ratio":enrichment,
    "candidate_mean_score":out["population"]["candidate_mean_aggregate_score"],
    "control_mean_score":out["population"]["control_mean_aggregate_score"],
    "permutation_p":perm_p,
    "mann_whitney_p":mw_p,
    "robustness_signal":robustness_signal,
    "scientific_promotion":False,
    "interpretation":out["interpretation"]
},indent=2))
