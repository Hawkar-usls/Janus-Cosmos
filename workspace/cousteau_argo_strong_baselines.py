#!/usr/bin/env python3
"""Falsify Argo Ridge against seasonal, latitude, and geographic baselines.

Comparators fixed in the source-pinned preregistration before first execution.
Uses original Argo90 without new data, and makes no automatic model promotion.
"""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).parent))
import cousteau_argo_external_holdout as prior

FIELDS=("temp_c","psal_psu")
DEPTHS=(10,100,300,600,1000)

def longitude(p):
    return ((float(p["longitude_deg"])+180)%360)-180

def latband(p):
    return min(3,max(0,int((-float(p["latitude_deg"])-50)//10)))

def month(p):
    return int(p["index_date_utc"][4:6])

def haversine_deg(p,q):
    a,b=map(math.radians,(float(p["latitude_deg"]),float(q["latitude_deg"])))
    dlat=b-a
    dlon=math.radians(longitude(q)-longitude(p))
    v=math.sin(dlat/2)**2+math.cos(a)*math.cos(b)*math.sin(dlon/2)**2
    return 2*math.asin(math.sqrt(max(0,min(v,1))))

def same_data_rows(floats, predicate, field):
    train=[(p,z) for p in floats if not predicate(p) for z in p["sample_levels"] if field in z]
    test=[(p,z) for p in floats if predicate(p) for z in p["sample_levels"] if field in z]
    return train,test

def metric(pred,truth):
    delta=np.asarray(pred,dtype=float)-np.asarray(truth,dtype=float)
    return {"rmse":round(float(np.sqrt(np.mean(delta**2))),6),
            "mae":round(float(np.mean(np.abs(delta))),6)}

def eval_comparators(floats,predicate,field):
    train,test=same_data_rows(floats,predicate,field)
    if len({p["wmo"] for p,z in train})<20 or len({p["wmo"] for p,z in test})<8:
        return {"status":"INSUFFICIENT_SAMPLE"}
    lookup_depth=defaultdict(list);lookup_month=defaultdict(list);lookup_lat=defaultdict(list)
    rows_by_depth=defaultdict(list)
    for p,z in train:
        d=int(z["pressure_dbar"])
        val=float(z[field])
        lookup_depth[d].append(val)
        lookup_month[d,month(p)].append(val)
        lookup_lat[d,latband(p)].append(val)
        rows_by_depth[d].append((p,val))
    methods={"depth_only":[],"month_depth":[],"latitude_depth":[],"nearest5_geospatial":[]}
    truth=[]
    fallback_counts={"month_depth":0,"latitude_depth":0}
    for p,z in test:
        d=int(z["pressure_dbar"])
        if d not in lookup_depth: return {"status":"MISSING_TRAINING_DEPTH"}
        truth.append(float(z[field]))
        depth_mean=float(np.mean(lookup_depth[d]))
        methods["depth_only"].append(depth_mean)
        arr=lookup_month.get((d,month(p)))
        if not arr:fallback_counts["month_depth"]+=1
        methods["month_depth"].append(float(np.mean(arr)) if arr else depth_mean)
        arr=lookup_lat.get((d,latband(p)))
        if not arr:fallback_counts["latitude_depth"]+=1
        methods["latitude_depth"].append(float(np.mean(arr)) if arr else depth_mean)
        nearest=sorted(rows_by_depth[d],key=lambda q:(haversine_deg(p,q[0]),q[0]["wmo"]))
        methods["nearest5_geospatial"].append(float(np.mean([y for _,y in nearest[:5]])))
    return {"status":"EVALUATED",
            "n_train_levels":len(train),"n_test_levels":len(test),
            "metrics":{k:metric(v,truth) for k,v in methods.items()},
            "fallback_counts":fallback_counts}

def scenarios():
    opts={
      "future_year_2024_vs_2022_2023":lambda f:f["index_date_utc"][:4]=="2024",
      "winter_JJA_vs_other_seasons":lambda f:month(f) in (6,7,8),
      "lat_60S_or_south_vs_north_of_60S":lambda f:float(f["latitude_deg"])<=-60,
      "lat_50_to_60S_vs_south_of_60S":lambda f:float(f["latitude_deg"])>-60
    }
    for s in range(4):
        lo=-180+s*90;hi=lo+90
        opts[f"longitude_{lo}_{hi}E_vs_others"]=lambda f,lo=lo,hi=hi:lo<=longitude(f)<hi
    return opts

def evaluate(obs,receipt,raw):
    floats=prior.validate(obs,receipt,raw)
    results=[]
    weaknesses=[]
    for name,split in scenarios().items():
        previous=prior.assess(floats,split,name)
        r={"split":name,"n_train_floats":previous["train_float_count"],
           "n_test_floats":previous["test_float_count"],"model_fields":{}}
        for field in FIELDS:
            old=previous["models"].get(field,{})
            comparison=eval_comparators(floats,split,field)
            if old.get("status")!="EVALUATED" or comparison.get("status")!="EVALUATED":
                r["model_fields"][field]={"status":"NOT_COMPARABLE"}
                continue
            baseline=old["depth_only_training_mean"]["rmse"]
            newdepth=comparison["metrics"]["depth_only"]["rmse"]
            if abs(baseline-newdepth)>1e-6:
                raise ValueError("Comparison mismatch against frozen prior depth-only test")
            ridge=old["ridge"]["rmse"]
            beats={key:bool(ridge<m["rmse"]) for key,m in comparison["metrics"].items()}
            if not all(beats.values()):
                weaknesses.append({"split":name,"field":field,"ridge_rmse":ridge,
                                   "unbeaten_comparators":[k for k,b in beats.items() if not b]})
            r["model_fields"][field]={
                "status":"EVALUATED","ridge_rmse":ridge,
                "baselines":comparison["metrics"],"ridge_beats":beats,
                "fallbacks":comparison["fallback_counts"]}
        results.append(r)
    return {
      "schema":"janus.cousteau.argo90_strong_baselines.v1",
      "source_training_run":receipt["run_id"],
      "source_sha256":prior.sha(raw),
      "status":"EXPLORATORY_POSTHOC_STRONG_BASELINE_EVALUATION",
      "all_source_hashes_verified":True,
      "comparator_protocol":"COUSTEAU-ARGO-90-STRONG-CONTROLS-PREREG-2026-11-09.json",
      "scenarios":results,
      "negative_findings":weaknesses,
      "promote_model":False,
      "scientific_ceiling":"Posthoc research on 90-profile convenience sample; cannot validate prospective prediction, observed heat transport, under-ice location or glacier melting."
    }

def selftest():
    prof=[]
    for i in range(40):
        prof.append({"wmo":str(400000+i),"data_mode":"D","position_qc":"1",
                     "latitude_deg":-51-(i%4)*10,
                     "longitude_deg":-140+(i%8)*40,
                     "index_date_utc":str(2022+i%3)+"%02d01000000"%(1+i%12),
                     "sample_levels":[{"pressure_dbar":d,
                                       "temp_c":-1+0.002*d+i*.04,
                                       "psal_psu":34+.0002*d+i*.002} for d in DEPTHS]})
    o={"profiles":prof,"admissible_profiles":40,"index_sha256":"fake"}
    raw=json.dumps(o).encode()
    p={"status":"REAL_ARGO_DATA_TRAINED_AND_EVALUATED",
       "observations_sha256":prior.sha(raw),"source_index_sha256":"fake","run_id":"synthetic"}
    r=evaluate(o,p,raw)
    assert len(r["scenarios"])==8 and r["promote_model"] is False
    assert r["scenarios"][0]["model_fields"]["temp_c"]["status"]=="EVALUATED"
    print("ARGO_STRONG_CONTROLS_SELFTEST_SYNTHETIC_PASS")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--observations")
    ap.add_argument("--receipt")
    ap.add_argument("--out")
    a=ap.parse_args()
    if a.self_test:return selftest()
    if not (a.observations and a.receipt and a.out):
        ap.error("source and output arguments required")
    raw=Path(a.observations).read_bytes()
    obs=json.loads(raw)
    receipt=json.loads(Path(a.receipt).read_bytes())
    report=evaluate(obs,receipt,raw)
    path=Path(a.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print(json.dumps({"status":report["status"],"negative_findings":report["negative_findings"],
                      "scenarios":len(report["scenarios"])},sort_keys=True))

if __name__=="__main__":main()
