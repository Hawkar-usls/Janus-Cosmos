#!/usr/bin/env python3
"""Source-year gated Argo T/S selector, never chosen on test-year scores.

This does not provide an independent prospective test: 2024 was already
studied by previous exploratory experiments and remains explicitly labelled.
"""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import cousteau_argo_external_holdout as ridge
import cousteau_argo_strong_baselines as strong

def year(float_row):
    return str(float_row["index_date_utc"])[:4]

def policy(observations, receipt, source_bytes):
    floats=ridge.validate(observations,receipt,source_bytes)
    train_22=[f for f in floats if year(f)=="2022"]
    calibrate_23=[f for f in floats if year(f)=="2023"]
    test_24=[f for f in floats if year(f)=="2024"]
    cohorts=[("2022",train_22),("2023",calibrate_23),("2024",test_24)]
    if len(train_22)<15 or len(calibrate_23)<10 or len(test_24)<30:
        return {
          "schema":"janus.cousteau.argo_year_selector_receipt.v1",
          "status":"REFUSE_INSUFFICIENT_COHORT_POWER",
          "cohorts":{k:len(v) for k,v in cohorts},
          "model_promotion_authorized":False}
    if len({f["wmo"] for f in floats})!=len(floats):
        raise ValueError("Duplicate float crossed year partitions")
    base_floats=train_22+calibrate_23
    cal_test=lambda f:year(f)=="2023"
    future_test=lambda f:year(f)=="2024"
    cal_ridge=ridge.assess(base_floats,cal_test,"calibration_2023")
    test_ridge=ridge.assess(floats,future_test,"holdout_2024")
    per_field={}
    for field in ("temp_c","psal_psu"):
        c=cal_ridge["models"][field]
        t=test_ridge["models"][field]
        near_cal=strong.eval_comparators(base_floats,cal_test,field)
        near_test=strong.eval_comparators(floats,future_test,field)
        if any(d.get("status")!="EVALUATED" for d in (c,t,near_cal,near_test)):
            per_field[field]={"status":"REFUSE_INCOMPLETE_EVAL"}
            continue
        c_r=c["ridge"]["rmse"]
        c_k=near_cal["metrics"]["nearest5_geospatial"]["rmse"]
        selected="ridge" if c_r<c_k else "nearest5_geospatial"
        test_rmse={"ridge":t["ridge"]["rmse"],
                   "nearest5_geospatial":near_test["metrics"]["nearest5_geospatial"]["rmse"],
                   "depth_only_baseline":t["depth_only_training_mean"]["rmse"]}
        selected_error=test_rmse[selected]
        per_field[field]={
           "status":"EVALUATED_POSTHOC_ONLY",
           "calibration_2023_rmse":{"ridge":c_r,"nearest5_geospatial":c_k},
           "selected_without_using_2024":selected,
           "test_2024_rmse":test_rmse,
           "selected_2024_rmse":selected_error,
           "selected_vs_depth_baseline":selected_error<test_rmse["depth_only_baseline"],
           "selected_vs_best_candidate":selected_error<=min(test_rmse["ridge"],test_rmse["nearest5_geospatial"]),
           "n_calibration_pressure_levels":c["n_test_levels"],
           "n_holdout_pressure_levels":t["n_test_levels"]
        }
    return {
      "schema":"janus.cousteau.argo_year_selector_receipt.v1",
      "status":"EXPLORATORY_YEAR_GATED_TWO_PHASE_SELECTOR_EVALUATED",
      "frozen_source_run":receipt["run_id"],
      "frozen_source_sha256":ridge.sha(source_bytes),
      "cohorts":{k:len(v) for k,v in cohorts},
      "separate_float_ids_by_year":True,
      "model_selector_rule":"Tune choice on 2023 only, train 2022, refit 2022+2023, evaluate year 2024",
      "fixed_candidate_models":["ridge_alpha10","nearest5_haversine_at_same_pressure"],
      "variables":per_field,
      "model_promotion_authorized":False,
      "limitations":[
         "Only 90 selected D-mode Argo floats south of 50S, not the complete Argo Southern Ocean.",
         "Posthoc: year 2024 dataset already studied in earlier experiments, despite 2024 being excluded from this selector's choice.",
         "Historical backtest does not prove future generalization, uncertainty calibration, density, ocean heat flux, under-ice navigation, or glacier melt.",
         "Candidate choices and baseline are conditional on severe spatiotemporal observation gaps."
      ]
    }

def selftest():
    samples=[]
    for i in range(100):
        yr=("2022","2023","2024")[i%3]
        samples.append({"wmo":str(50000+i),"data_mode":"D","position_qc":"1",
                        "latitude_deg":-51-(i%4)*10,
                        "longitude_deg":-160+(i%8)*40,
                        "index_date_utc":yr+"%02d12000000"%(1+i%12),
                        "sample_levels":[{"pressure_dbar":d,"temp_c":round(2-0.002*d+i*.013,6),
                                          "psal_psu":round(34+0.0001*d+i*.003,6)}
                                         for d in (10,100,300,600,1000)]})
    obs={"profiles":samples,"admissible_profiles":len(samples),"index_sha256":"synthetic"}
    raw=json.dumps(obs).encode()
    receipt={"status":"REAL_ARGO_DATA_TRAINED_AND_EVALUATED",
             "observations_sha256":ridge.sha(raw),
             "source_index_sha256":"synthetic","run_id":"selftest_only"}
    res=policy(obs,receipt,raw)
    assert res["cohorts"]=={"2022":34,"2023":33,"2024":33}
    assert res["status"]=="EXPLORATORY_YEAR_GATED_TWO_PHASE_SELECTOR_EVALUATED"
    assert res["model_promotion_authorized"] is False
    assert all(v["selected_without_using_2024"] in ("ridge","nearest5_geospatial") for v in res["variables"].values())
    print("ARGO_YEAR_SELECTOR_SYNTHETIC_SELFTEST_PASS")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--observations")
    ap.add_argument("--receipt")
    ap.add_argument("--out")
    args=ap.parse_args()
    if args.self_test:return selftest()
    if not (args.observations and args.receipt and args.out):
        ap.error("Missing source/receipt/out")
    raw=Path(args.observations).read_bytes()
    obs=json.loads(raw)
    receipt=json.loads(Path(args.receipt).read_bytes())
    answer=policy(obs,receipt,raw)
    dst=Path(args.out)
    dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(answer,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":answer["status"],"cohorts":answer["cohorts"],
        "variables":answer.get("variables",{})},sort_keys=True))

if __name__=="__main__":main()
