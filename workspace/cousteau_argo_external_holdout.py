#!/usr/bin/env python3
"""Cross-year / latitude / longitude stress tests for frozen Argo model design.

Never tunes hyperparameters using holdout, never promotes model based on a
single convenience sample, and preserves missing-coverage/failed gates.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

FIELDS=("temp_c","psal_psu")
LEVELS=(10,100,300,600,1000)

def sha(b): return hashlib.sha256(b).hexdigest()

def features(row, pressure):
    month=int(row["index_date_utc"][4:6])
    lat=float(row["latitude_deg"])
    lon=math.radians(float(row["longitude_deg"]))
    p=math.log1p(pressure)/7.0
    x=(lat+60.0)/20.0
    return [p,x,math.sin(lon),math.cos(lon),
            math.sin(month*math.tau/12),math.cos(month*math.tau/12),
            p*x,p*math.sin(lon),p*math.cos(lon)]

def validate(obs,receipt,raw):
    assert receipt["status"]=="REAL_ARGO_DATA_TRAINED_AND_EVALUATED", "No real-data source receipt"
    assert sha(raw)==receipt["observations_sha256"], "Observation source hash mismatch"
    assert obs["index_sha256"]==receipt["source_index_sha256"], "Index SHA mismatch"
    floats=obs["profiles"]
    if len(floats)!=obs["admissible_profiles"] or len(floats)<24: raise ValueError("Incomplete sample")
    if len(set(f["wmo"] for f in floats))!=len(floats): raise ValueError("Duplicated float WMO")
    for f in floats:
        if f["data_mode"]!="D" or f["position_qc"]!="1":
            raise ValueError("Unqualified data mode/position")
        if not -90<=float(f["latitude_deg"])<=-50: raise ValueError("latitude outside Southern Ocean")
        if len(set(int(p["pressure_dbar"]) for p in f["sample_levels"]))!=len(f["sample_levels"]):
            raise ValueError("duplicate pressure level")
        if any(int(p["pressure_dbar"]) not in LEVELS for p in f["sample_levels"]):
            raise ValueError("unexpected pressure level")
    return floats

def assess(rows, is_test, split_name):
    training=[f for f in rows if not is_test(f)]
    test=[f for f in rows if is_test(f)]
    base={"name":split_name,"train_float_count":len(training),"test_float_count":len(test)}
    if len(training)<20 or len(test)<8:
        return dict(base,status="INSUFFICIENT_GROUP_POWER",models={})
    per_field={}
    for field in FIELDS:
        train=[(f,p) for f in training for p in f["sample_levels"] if field in p]
        held=[(f,p) for f in test for p in f["sample_levels"] if field in p]
        if len(train)<30 or len(held)<12:
            per_field[field]={"status":"INSUFFICIENT_VALID_LEVELS"};continue
        y=np.array([p[field] for f,p in train],dtype=float)
        truth=np.array([p[field] for f,p in held],dtype=float)
        tr=np.array([features(f,int(p["pressure_dbar"])) for f,p in train],dtype=float)
        te=np.array([features(f,int(p["pressure_dbar"])) for f,p in held],dtype=float)
        if not np.all(np.isfinite(y)) or not np.all(np.isfinite(tr)) or not np.all(np.isfinite(te)):
            raise ValueError("Nonfinite training or evaluation data")
        mu=tr.mean(axis=0); std=tr.std(axis=0); std[std<1e-9]=1.0
        z=(tr-mu)/std
        alpha=10.0
        weights=np.linalg.solve(z.T@z+alpha*np.eye(z.shape[1]),z.T@(y-y.mean()))
        pred=float(y.mean())+((te-mu)/std)@weights
        by_depth={}
        for depth in LEVELS:
            values=[p[field] for f,p in train if int(p["pressure_dbar"])==depth]
            if values: by_depth[depth]=float(np.mean(values))
        if any(int(p["pressure_dbar"]) not in by_depth for f,p in held):
            per_field[field]={"status":"BASELINE_MISSING_TRAIN_DEPTH"};continue
        baseline=np.array([by_depth[int(p["pressure_dbar"])] for f,p in held])
        def error(est):
            d=est-truth
            return {"rmse":round(float(np.sqrt(np.mean(d*d))),6),
                    "mae":round(float(np.mean(np.abs(d))),6)}
        ridge=error(pred); base_metric=error(baseline)
        per_field[field]={"status":"EVALUATED",
             "n_train_levels":len(train),"n_test_levels":len(held),
             "ridge":ridge,"depth_only_training_mean":base_metric,
             "ridge_beats_baseline":bool(ridge["rmse"]<base_metric["rmse"]),
             "rmse_delta_baseline_minus_ridge":round(base_metric["rmse"]-ridge["rmse"],6),
             "fixed_alpha":10.0}
    return dict(base,status="EVALUATED",models=per_field,
                claims="RETROSPECTIVE_GEO_TEMPORAL_STRESS_TEST_ONLY")

def run(obs,receipt,raw):
    floats=validate(obs,receipt,raw)
    def longitude(f):
        return ((float(f["longitude_deg"])+180)%360)-180
    scenarios={
       "future_year_2024_vs_2022_2023":lambda f:f["index_date_utc"][:4]=="2024",
       "winter_JJA_vs_other_seasons":lambda f:int(f["index_date_utc"][4:6]) in (6,7,8),
       "lat_60S_or_south_vs_north_of_60S":lambda f:float(f["latitude_deg"])<=-60.0,
       "lat_50_to_60S_vs_south_of_60S":lambda f:float(f["latitude_deg"])>-60.0,
    }
    for sector in range(4):
        lo=-180+sector*90;hi=lo+90
        scenarios[f"longitude_{lo}_{hi}E_vs_others"]=(
            lambda f,lo=lo,hi=hi: lo<=longitude(f)<hi)
    out=[assess(floats,fn,name) for name,fn in scenarios.items()]
    not_evaluable=[v["name"] for v in out if v["status"]!="EVALUATED"]
    negative=[{"split":v["name"],"field":k,"baseline_rmse":m["depth_only_training_mean"]["rmse"],
               "ridge_rmse":m["ridge"]["rmse"]}
              for v in out if v["status"]=="EVALUATED" for k,m in v["models"].items()
              if m["status"]=="EVALUATED" and not m["ridge_beats_baseline"]]
    return {
      "schema":"janus.cousteau.argo_stress_holdout.v1",
      "status":"REPRODUCED_STRATIFIED_RESTROSPECTIVE_QC_STRESS_TEST",
      "source_run":receipt["run_id"],
      "source_obs_sha256":sha(raw),
      "validation_scopes":"No tuning on any holdout, distinct WMO floats across train/test",
      "sample_float_count":len(floats),
      "model_spec":"unchanged fixed alpha=10 Ridge, 9 predictors; training-only standardization; heldout group",
      "scenarios":out,
      "unevaluable_splits":not_evaluable,
      "model_losses_against_baseline":negative,
      "model_promotion_authorized":False,
      "limitations":[
        "Convenience sample of 90 floats, not representative of all Southern Ocean.",
        "Year-heldout test is a retrospective period split, not immutable future observations.",
        "Geographic heldout tests measure extrapolation across sparse sampled regions, not glacier melt.",
        "No source data from IBCSO/MEOP/BGC are used in this benchmark.",
        "Seasonal group holdout is a deliberate distribution shift and must not be described as one-month future forecast.",
        "Individual depth levels interpolated between observed pressure levels per predeclared QC rule."
      ]
    }

def selftest():
    profiles=[]
    for n in range(40):
        profiles.append({"wmo":str(70000+n),"data_mode":"D","position_qc":"1",
            "latitude_deg":-51-(n%4)*8,"longitude_deg":-150+(n%8)*40,
            "index_date_utc":str(2022+(n%3))+"%02d01000000"%(1+n%12),
            "sample_levels":[{"pressure_dbar":d,
                 "temp_c":round(1+n*.02-.001*d,5),
                 "psal_psu":round(34+n*.001+.0001*d,5)} for d in LEVELS]})
    obs={"profiles":profiles,"admissible_profiles":len(profiles),"index_sha256":"synthetic"}
    raw=json.dumps(obs).encode()
    rec={"status":"REAL_ARGO_DATA_TRAINED_AND_EVALUATED",
         "observations_sha256":sha(raw),"source_index_sha256":"synthetic","run_id":"selftest_only"}
    output=run(obs,rec,raw)
    assert len(output["scenarios"])==8 and output["model_promotion_authorized"] is False
    assert output["sample_float_count"]==40
    try: validate(obs,rec,raw+b"tamper")
    except AssertionError: pass
    else: raise RuntimeError("SHA validation bypassed")
    print("STRESS_HOLDOUT_SELFTEST_PASS_SYNTHETIC_ONLY")

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--self-test",action="store_true")
    parser.add_argument("--observations")
    parser.add_argument("--receipt")
    parser.add_argument("--out")
    args=parser.parse_args()
    if args.self_test:return selftest()
    if not (args.observations and args.receipt and args.out):
        parser.error("Need --observations --receipt --out")
    raw=Path(args.observations).read_bytes()
    result=run(json.loads(raw),json.loads(Path(args.receipt).read_text()),raw)
    path=Path(args.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"sample":result["sample_float_count"],
        "negative_scenarios":result["model_losses_against_baseline"],
        "unevaluable_scenarios":result["unevaluable_splits"]},sort_keys=True))

if __name__=="__main__":main()
