#!/usr/bin/env python3
"""Recover and compare two full Argo GDAC acquisition receipts.

This is an *observability* comparison, not a same-input model A/B test:
GDAC source index hashes differ. No assertion that improved completeness
caused improved prediction performance is allowed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import cousteau_southern_observability_audit as observer

def sha(data):
    return hashlib.sha256(data).hexdigest()

def validate(obs_bytes,rec_bytes):
    obs=json.loads(obs_bytes);rec=json.loads(rec_bytes)
    if rec.get("status")!="REAL_ARGO_DATA_TRAINED_AND_EVALUATED":
        raise ValueError("No valid successful real-data training receipt")
    if sha(obs_bytes)!=rec["observations_sha256"]:
        raise ValueError("Source observation hash mismatch")
    if obs["index_sha256"]!=rec["source_index_sha256"]:
        raise ValueError("GDAC index SHA mismatch")
    if len(obs["profiles"])!=obs["admissible_profiles"]:
        raise ValueError("Incomplete profile list")
    if len({p["wmo"] for p in obs["profiles"]})!=len(obs["profiles"]):
        raise ValueError("WMO duplicated")
    return obs,rec

def summarize(obs,rec,bytes_):
    old=observer.build(obs,rec,bytes_,source_run=rec["run_id"])
    sample=old["sample"]
    return {
        "run_id":rec["run_id"],
        "observation_sha256":sha(bytes_),
        "GDAC_index_sha256":obs["index_sha256"],
        "selection_rule":obs.get("index_selection_stats",{}).get("selection_rule",
                            "Legacy newest-per-WMO (before lat-month balancing)"),
        "accepted_float_profiles":sample["unique_floats"],
        "attempts":sample["download_attempts"],
        "rejects":sample["rejections"],
        "by_year":sample["by_year"],
        "by_month":sample["by_month"],
        "by_season":sample["by_season"],
        "by_latitude_band":sample["by_latitude_band"],
        "nonempty_lat_month_cells":sum(c["count"]>0 for c in sample["latitude_month_cells"]),
        "empty_lat_month_cells":sum(c["count"]==0 for c in sample["latitude_month_cells"]),
        "pressure_depth_coverage":sample["pressure_coverage_dbar"],
        "ridge_rmse_temp_c":rec["result"]["model_benchmark"]["results"]["temp_c"]["ridge"]["rmse"],
        "ridge_rmse_psal":rec["result"]["model_benchmark"]["results"]["psal_psu"]["ridge"]["rmse"]
    }

def compare(old_obs_bytes,old_rec_bytes,new_obs_bytes,new_rec_bytes):
    old_o,old_r=validate(old_obs_bytes,old_rec_bytes)
    new_o,new_r=validate(new_obs_bytes,new_rec_bytes)
    o=summarize(old_o,old_r,old_obs_bytes)
    n=summarize(new_o,new_r,new_obs_bytes)
    shifted_indices=o["GDAC_index_sha256"]!=n["GDAC_index_sha256"]
    return {
        "schema":"janus.cousteau.argo_source_snapshots_coverage_comparison.v1",
        "status":"TWO_SOURCE_SNAPSHOTS_VALIDATED__MODEL_SUPERIORITY_NOT_INFERABLE",
        "source_runs":[o["run_id"],n["run_id"]],
        "old":o,"new":n,
        "coverage_difference":{
            "accepted_profile_count_change":n["accepted_float_profiles"]-o["accepted_float_profiles"],
            "nonempty_lat_month_cells_change":n["nonempty_lat_month_cells"]-o["nonempty_lat_month_cells"],
            "high_latitude_70_80S_change":n["by_latitude_band"]["70-80S"]-o["by_latitude_band"]["70-80S"],
            "december_profile_count_change":n["by_month"].get("12",0)-o["by_month"].get("12",0),
            "winter_JJA_profile_count_change":n["by_season"].get("austral_winter_JJA",0)-o["by_season"].get("austral_winter_JJA",0),
            "source_GDAC_index_changed":shifted_indices,
            "apples_to_apples_same_input_model_comparison":not shifted_indices and False,
            "source_snapshots_immutable_sha_verified":True,
            "training_sample_bias_not_eliminated":True,
            "model_automatic_promotion":False
        },
        "interpretation":"Observed differences are jointly affected by sample-selection policy, source-index snapshot and QC rejection rates; cannot causally assign all changes to stratification.",
        "model_performance_comparison_prohibited":True,
        "falsifiability_note":"A future frozen identical GDAC index with both sampling policies and independent high-latitude holdout is required to attribute coverage and model gains to selection.",
        "underice_coordinates_unknown":True,
        "raw_sonar_or_ice_shelf_melt_inference_allowed":False
    }

def selftest():
    def make(prefix,month,latitude):
        o={"profiles":[],"admissible_profiles":2,
            "profiles_download_attempted":2,
            "index_sha256":prefix,
            "rejections":{},
            "index_selection_stats":{}}
        for i in range(2):
            o["profiles"].append({"wmo":str(900000+i),"data_mode":"D","position_qc":"1",
               "latitude_deg":latitude,"longitude_deg":20,
               "index_date_utc":"2023%02d02000000"%month,
               "sample_levels":[{"pressure_dbar":100,"temp_c":0.5,"psal_psu":34.5}]})
        raw=json.dumps(o).encode()
        r={"run_id":prefix,"status":"REAL_ARGO_DATA_TRAINED_AND_EVALUATED",
           "observations_sha256":sha(raw),"source_index_sha256":prefix,
           "result":{"model_benchmark":{"results":{
             "temp_c":{"ridge":{"rmse":1.1}},
             "psal_psu":{"ridge":{"rmse":.2}}
           }}}}
        return raw,json.dumps(r).encode()
    result=compare(*make("A",12,-55),*make("B",6,-75))
    assert result["coverage_difference"]["source_GDAC_index_changed"] is True
    assert result["coverage_difference"]["nonempty_lat_month_cells_change"]==0
    assert result["model_automatic_promotion"] is False if "model_automatic_promotion" in result else result["coverage_difference"]["model_automatic_promotion"] is False
    try:
        o,r=make("A",12,-55)
        validate(o+b"tamper",r)
    except ValueError: pass
    else: raise AssertionError("bad SHA was not refused")
    print("ARGO_SNAPSHOT_COMPARE_SELFTEST_SYNTHETIC_PASS")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--self-test",action="store_true")
    for a in ("old-observations","old-receipt","new-observations","new-receipt","out"):
        p.add_argument("--"+a)
    x=p.parse_args()
    if x.self_test:return selftest()
    if not all(getattr(x,k.replace("-","_")) for k in ("old-observations","old-receipt","new-observations","new-receipt","out")):
        p.error("all four source files and output required")
    paths=[getattr(x,k.replace("-","_")) for k in ("old-observations","old-receipt","new-observations","new-receipt")]
    result=compare(*(Path(f).read_bytes() for f in paths))
    output=Path(x.out);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":result["status"],"coverage":result["coverage_difference"]},sort_keys=True))

if __name__=="__main__":main()
