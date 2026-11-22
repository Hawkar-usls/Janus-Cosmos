"""JANUS SIT2C — fixed RFI-aware temporal score, X-development → L-band test.

SOURCES: true public GBT 3I/ATLAS X and L receiver HDF5 *power arrays*.
TARGET TRUTH: only synthetic post-extraction weak/strong line and RFI overlays.
No external pointing witness, no sky truth, no retrocausal communication.

Frozen source and comparator:
  - existing SIT1 real HDF5 decoder and its 16x128 measured panel extraction
  - original COSMOS R42A8A blind ridge path/score
  - SIT2B exact 672-path direct matched filters and synthetic generator
  - equal-information direct RFI-aware search and pinned R5F mirrored replay

Test L-band scores and labels never change the development-frozen X thresholds.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import janus_sit2b_temporal_recovery as previous
from tools.tachyon_r5f_stock_mirror_audit import mirror_case, mirrored_replay
from tools.tachyon_r5f_causal_memory_audit import ABSTAIN

ROOT=Path(__file__).resolve().parents[2]
Q=json.loads((ROOT/"research/seti/janus_sit2c_protocol.json").read_text())
SIT1=previous.SIT1
OUT=ROOT/"runtime/sit2c"
SEED=int(Q["dataset"]["seed"])
PATHS=previous.PATHS
T=np.arange(16,dtype=np.int32)[None,:]
LABELS=tuple(Q["dataset"]["classes"])
ARMS=("A_NATIVE_COSMOS_RIDGE",
      "B_STRONG_DIRECT_SIDEBAND",
      "C_STRONG_DIRECT_OFF_AWARE",
      "D_RFI_AWARE_DIRECT",
      "E_TACHYON_REVERSE_RFI_AWARE")


def off_reference(rows):
    """Unlabelled actual OFF power, before any synthetic injections; same for D/E/C."""
    source=[r["panel"] for r in rows if r["role"]=="OFF"]
    if not source:raise ValueError("NO_NOMINAL_OFF_REFERENCE")
    references={}
    for idx in range(4):
        zz=[]
        for r in rows:
            if r["role"]=="OFF" and r["panel_idx"]==idx:
                z,_=previous.robust_normalize(r["panel"])
                zz.append(z)
        if len(zz) not in (1,3):
            raise ValueError("OFF_REFERENCE_COUNT_NOT_ONE_OR_THREE")
        refs=np.median(np.stack(zz),axis=(0,1))
        assert refs.shape==(128,) and np.isfinite(refs).all()
        references[idx]=np.array(refs,dtype=np.float64)
    return references


def score_bank(panel,off_ref,method,backward=False):
    """Same 672 candidate paths and 16 time rows for all matched objective arms.

    One receiver's original, unmodified nominal-OFF reference is shared
    with all OFF-aware arms, including the strong direct competitor C.
    """
    z,_=previous.robust_normalize(panel)
    target=z[T,PATHS]
    side=target-.5*(z[T,PATHS-3]+z[T,PATHS+3])
    if backward:
        side=side[:,::-1]
    if method=="sideband":
        return side.mean(axis=1)
    off=np.clip(off_ref[PATHS]-1.5,0,10).mean(axis=1)
    if method=="offaware":
        return side.mean(axis=1)-.35*off
    if method=="rfi":
        persistence=(side<.2).mean(axis=1)
        spread=np.maximum(side.std(axis=1)-2.0,0)
        return (np.median(side,axis=1)
                +.25*np.clip(side,0,5).mean(axis=1)
                -.30*persistence-.15*spread-.35*off)
    raise ValueError("UNKNOWN_PREREGISTERED_SCORER")


def choose(bank):
    # Deterministic tie resolution; never depends on route scan order.
    return previous.best_path(bank)


def score_all(panel,off_ref,model):
    t=time.perf_counter()
    native_path,unused_shadow,native_ridge=SIT1.score_one(panel,model)
    native_time=time.perf_counter()-t
    out={ARMS[0]:(native_path.copy(),native_ridge)}
    costs={"native_seconds":native_time}
    for label,method in ((ARMS[1],"sideband"),(ARMS[2],"offaware"),(ARMS[3],"rfi")):
        t=time.perf_counter()
        path,score,_=choose(score_bank(panel,off_ref,method,False))
        costs[label+"_seconds"]=time.perf_counter()-t
        out[label]=(path,score)
    t=time.perf_counter()
    path,score,_=choose(score_bank(panel,off_ref,"rfi",True))
    costs[ARMS[4]+"_backward_seconds"]=time.perf_counter()-t
    direct_path,direct_score=out[ARMS[3]]
    if not np.array_equal(path,direct_path) or abs(score-direct_score)>1e-9:
        raise AssertionError("REVERSE_RFI_OBJECTIVE_IS_NOT_FORWARD_EQUIVALENT")
    out[ARMS[4]]=(path,score)
    return out,costs


def trace_actual_source(row,path,score,source_sha,model_lock):
    """Pinned causal graph traces *computation provenance*, not sky events."""
    if row["source_file_sha256"]!=source_sha:
        raise AssertionError("SOURCE_SHA_WAS_SUBSTITUTED")
    record=SIT1.manifest(row,path,float(score),model_lock)
    graph=mirror_case(record["graph"])
    r=mirrored_replay(graph,budget=500)
    if r["prediction"]==ABSTAIN:
        raise AssertionError("TACHYON_REVERSE_PROVENANCE_ABSTAINED")
    if record["panel_sha256"]!=SIT1.hash_array(row["panel"],dtype="<f4"):
        raise AssertionError("CASUAL_TEST_PANEL_DIGEST_INCORRECT")
    if record["candidate_path_sha256"]!=SIT1.hash_array(path,dtype="<i4"):
        raise AssertionError("CASUAL_PATH_DIGEST_INCORRECT")
    return int(r["operations"])


def download_group(rows,group,prior_x,folder):
    receipts=[];panels=[];errors=[]
    for scan in rows:
        info=SIT1.retrieve({"number":scan["scan"],"role":scan["role"],"path":scan["filename"]},folder)
        info["group"]=group
        receipts.append(info)
        if not info.get("ok"):
            errors.append({"group":group,"scan":scan["scan"],"kind":"DOWNLOAD","attempts":info.get("retry_errors")})
            continue
        if group=="X_DEV" and info["first_seen_sha256"]!=prior_x.get(scan["scan"]):
            errors.append({"group":group,"scan":scan["scan"],"kind":"SOURCE_SHA_MISMATCH"})
            (folder/scan["filename"]).unlink(missing_ok=True)
            continue
        try:
            decoded=SIT1.reduce_real_h5(folder/scan["filename"],
                                      {"scan":scan["scan"],"role":scan["role"],
                                       "first_seen_sha256":info["first_seen_sha256"]})
            panels.extend(decoded["panels"])
            info["header"]=decoded["header"]
            info["shape"]=decoded["dataset_shape"]
            info["panel_count"]=len(decoded["panels"])
        except Exception as e:
            errors.append({"group":group,"scan":scan["scan"],
                           "kind":"HDF5_FORMAT_OR_PHYSICAL_AXIS","message":str(e)[:300]})
        finally:
            (folder/scan["filename"]).unlink(missing_ok=True)
    expected=8 if group=="X_DEV" else 24
    if len(panels)!=expected:
        errors.append({"group":group,"kind":"MISSING_REAL_PANELS","got":len(panels),"expected":expected})
    if len(panels)==expected:
        base=panels[0]
        if any(abs(x["native_freq_start_mhz"]-base["native_freq_start_mhz"])>1e-5
               for x in panels if x["panel_idx"]==base["panel_idx"]):
            errors.append({"group":group,"kind":"MISMATCHED_PHYSICAL_FREQ_GRID"})
    return receipts,panels,errors


def make_cases(panels,rng):
    cases=[]
    for base in panels:
        for replicate in range(4):
            for kind_idx,kind in enumerate(LABELS):
                # Reuse the frozen SIT2B injection generator, fixed ranges,
                # separate source path labels inaccessible to scorers.
                value,truth=previous.generate_overlay(
                    base["panel"],previous.CATEGORIES[kind_idx],rng)
                row=previous.prev.matched_row(value,base)
                cases.append({"row":row,"kind":kind,
                              "truth_for_eval_only":truth,
                              "base_source_sha":base["source_file_sha256"],
                              "replicate":replicate})
    return cases


def evaluate_rows(rows):
    dev=[x for x in rows if x["group"]=="X_DEV"]
    eva=[x for x in rows if x["group"]=="L_EVAL"]
    if len(dev)!=128 or len(eva)!=384:
        raise AssertionError("SIT2C_SPLIT_LEAK_OR_COUNT_ERROR")
    outcome={}
    internal={}
    for arm in ARMS:
        dev_neg=[z["scores"][arm] for z in dev if z["kind"] in LABELS[:2]]
        if len(dev_neg)!=64:
            raise AssertionError("WRONG_DEVELOPMENT_CONTROL_COUNT")
        threshold=float(np.quantile(np.array(dev_neg),.95,method="higher"))
        weak=[z for z in eva if z["kind"]==LABELS[2]]
        strong=[z for z in eva if z["kind"]==LABELS[3]]
        controls=[z for z in eva if z["kind"] in LABELS[:2]]
        weak_ok={z["id"] for z in weak if z["scores"][arm]>threshold and z["localized"][arm]}
        strong_ok={z["id"] for z in strong if z["scores"][arm]>threshold and z["localized"][arm]}
        control_fail={z["id"] for z in controls if z["scores"][arm]>threshold}
        good=len(strong_ok)>=87 and len(control_fail)<=9
        outcome[arm]={
            "development_only_95th_percentile_threshold":threshold,
            "heldout_weak_count":len(weak),
            "heldout_strong_count":len(strong),
            "heldout_control_count":len(controls),
            "weak_path_localized":sum(z["localized"][arm] for z in weak),
            "strong_path_localized":sum(z["localized"][arm] for z in strong),
            "weak_detected":len(weak_ok),
            "strong_detected":len(strong_ok),
            "controls_above_development_threshold":len(control_fail),
            "strong_and_false_control_gate_pass":good,
            "test_set_operating_point_is_cross_receiver_without_readaptation":True,
        }
        internal[arm]=weak_ok
    direct=outcome[ARMS[3]]
    reverse=outcome[ARMS[4]]
    if direct!=reverse or internal[ARMS[3]]!=internal[ARMS[4]]:
        raise AssertionError("D_DIRECT_AND_E_REVERSE_NOT_EXACTLY_EQUAL")
    baseline_weak=max(outcome[ARMS[1]]["weak_detected"],outcome[ARMS[2]]["weak_detected"])
    baseline_mincontrols=min(outcome[ARMS[1]]["controls_above_development_threshold"],
                             outcome[ARMS[2]]["controls_above_development_threshold"])
    gain={
      "rfi_aware_extra_weak_vs_best_direct_matched":direct["weak_detected"]-baseline_weak,
      "rfi_aware_extra_control_exceedances_vs_least_direct_matched":
          direct["controls_above_development_threshold"]-baseline_mincontrols,
      "tachyon_reverse_extra_weak_vs_same_objective_direct":0,
      "reverse_same_objective_exactly_equal_in_all_384_cases":True,
      "rfi_aware_candidate_threshold_qualification":(
          direct["strong_and_false_control_gate_pass"] and
          direct["weak_detected"]-baseline_weak>=10 and
          direct["controls_above_development_threshold"]<=baseline_mincontrols
      ),
      "external_scientific_sensitivity_gain_confirmed":False,
    }
    return outcome,gain


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    folder=OUT/"source_hdf5";folder.mkdir(exist_ok=True)
    original=json.loads((ROOT/"research/seti/janus_sit1_open_result.json").read_text())
    prior={x["scan"]:x["first_seen_sha256"] for x in original["runtime_inputs"]["hdf5_files"]}
    dev_a,dev_b,dev_errors=download_group(Q["inputs"]["dev_files"],"X_DEV",prior,folder)
    eval_a,eval_b,eval_errors=download_group(Q["inputs"]["heldout_files"],"L_EVAL",prior,folder)
    errors=dev_errors+eval_errors
    result={
      "schema":"janus.cosmos.sit2c.real_radio_crossband_rfi_score_result.v1",
      "status":"REAL_SOURCE_DATA_OR_FORMAT_HOLD" if errors else "REAL_GBT_CROSSBAND_RFI_AWARE_SCORE_REPRODUCED__SCIENCE_HOLD",
      "source_scans":[{key:x.get(key) for key in ("scan","role","group","filename","bytes","first_seen_sha256","header","shape","error","retry_errors")}
                      for x in dev_a+eval_a],
      "error_receipts":errors,
      "true_X_development_panels":len(dev_b),
      "true_L_evaluation_panels":len(eval_b),
      "development_archive":"X receiver GBT 3I/ATLAS, previously seen, 2 scans",
      "evaluation_archive":"L receiver GBT 3I/ATLAS, 6 NEW and different receiver/cadence scans, same date/target/telescope",
      "original_L_source_sha256_published_by_provider":False,
      "external_independent_observation_dates":False,
      "independent_antenna_pointing":False,
      "native_model_weights_modified":False,
      "production_detector_or_threshold_modified":False,
      "physical_time_reversal":False,
      "real_ETI_or_natural_radio_ground_truth":False,
      "scientific_gate":"SCIENCE_HOLD_ONLY_REAL_BG_PLUS_SYNTHETIC_LABELS__ONE_TELESCOPE_ONE_DAY__NO_INDEPENDENT_POINTING"
    }
    if errors:return result
    dev_off=off_reference(dev_b)
    eva_off=off_reference(eval_b)
    model=SIT1.cosmos.frozen_model()
    lock=SIT1.cosmos.Q["frozen_model"]["fit_lock_sha256"]
    rng=np.random.default_rng(SEED)
    records=[]
    cost={name:0.0 for name in ARMS}
    cost["R5F_forward_trace_seconds"]=0.0
    count_ops=0
    for group,rows,off in (("X_DEV",dev_b,dev_off),("L_EVAL",eval_b,eva_off)):
        for item in make_cases(rows,rng):
            row=item["row"]
            all_scores,elapsed=score_all(row["panel"],off[row["panel_idx"]],model)
            for key,t in elapsed.items():
                prefix=key.replace("_seconds","").replace("_backward","")
                arm=key if key in cost else prefix
                if prefix=="native":arm=ARMS[0]
                if prefix in ARMS:cost[prefix]+=t
                if arm in cost and arm!=prefix:cost[arm]+=t
            start=time.perf_counter()
            path,score=all_scores[ARMS[4]]
            count_ops+=trace_actual_source(row,path,score,item["base_source_sha"],lock)
            cost["R5F_forward_trace_seconds"]+=time.perf_counter()-start
            truth=item["truth_for_eval_only"]
            records.append({
               "id":len(records),"group":group,"scan":row["scan"],
               "kind":item["kind"],
               "scores":{k:float(v[1]) for k,v in all_scores.items()},
               "localized":{k:previous.localizes(v[0],truth) for k,v in all_scores.items()},
               "real_parent_scan_sha":item["base_source_sha"],
               "post_injection_test_panel_sha":row["measured_panel_sha256"]
            })
    assert len(records)==512
    outcome,gain=evaluate_rows(records)
    result["status"]="REAL_GBT_X_TO_L_CROSSBAND_EXAM_COMPUTED__NO_EXTERNAL_SCIENTIFIC_SENSITIVITY_VALIDATION"
    result["dataset"]={
       "total_cases":512,"development_cases":128,"heldout_cases":384,
       "heldout_weak_injections":96,"heldout_strong_injections":96,
       "heldout_control_cases":192,"one_L_cadence_not_independent_records":True,
       "no_threshold_fit_on_L_evaluation":True,"artificial_labels_only":True,
       "actual_L_receiver_file_count":len(eval_a)
    }
    result["outcomes"]=outcome
    result["paired_gains"]=gain
    result["wall_seconds_exploratory"]={k:round(v,5) for k,v in cost.items()}
    result["R5F_extra_graph_operations"]=count_ops
    result["all_exact_forward_reverse_equal"]=True
    result["source_manifest_hash"]=hashlib.sha256(json.dumps(
        result["source_scans"],sort_keys=True).encode()).hexdigest()
    result["score_manifest_hash"]=hashlib.sha256(json.dumps(
        [{"scores":x["scores"],"kind":x["kind"],"scan":x["scan"]} for x in records],
        sort_keys=True).encode()).hexdigest()
    return result


def main():
    pa=argparse.ArgumentParser()
    pa.add_argument("--audit",action="store_true")
    args=pa.parse_args()
    if not args.audit:pa.error("--audit required")
    result=run()
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"sit2c_result.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("SIT2C_METRICS",json.dumps({
        "status":result["status"],"errors":result["error_receipts"],
        "outcomes":result.get("outcomes"),"paired_gains":result.get("paired_gains")
    },sort_keys=True),flush=True)
    if result.get("outcomes"):
        assert result["dataset"]["heldout_cases"]==384
        assert result["all_exact_forward_reverse_equal"]
        assert result["paired_gains"]["tachyon_reverse_extra_weak_vs_same_objective_direct"]==0
        print("JANUS_SIT2C_REAL_L_RECEIVER_CROSSBAND_SCORE_TEST_PASS__SCIENCE_HOLD")
    else:
        print("JANUS_SIT2C_NEW_SOURCE_FORMAT_OR_ACCESS_HOLD")
if __name__=="__main__":
    main()
