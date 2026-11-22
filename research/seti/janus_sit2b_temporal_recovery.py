"""JANUS SIT2B — temporal coherent-stacking and reversed path tracing.

Data: all six SHA256-verified actual GBT 3I/ATLAS HDF5 files used in SIT1
as an observational *background*. The labelled line and RFI perturbations
are SIMULATED and added only after acquisition. Results do not prove real
astronomical sensitivity or any physical reversal of time.

Five predeclared scoring arms:
- unchanged native COSMOS ridge and frozen shadow score;
- exact 672-path forward matched filter;
- exact 672-path forward symmetric sideband matched filter;
- same 672-path sideband score traversed BACKWARDS, then independently
  forward-replayed and source-checked using pinned DemiHead R5F mirror.
All thresholds frozen on the first ON/OFF development scans only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
import janus_sit2_weak_rfi as prev
from tools.tachyon_r5f_stock_mirror_audit import mirror_case, mirrored_replay
from tools.tachyon_r5f_causal_memory_audit import ABSTAIN

ROOT=Path(__file__).resolve().parents[2]
Q=json.loads((ROOT/"research/seti/janus_sit2b_protocol.json").read_text())
OUT=ROOT/"runtime/sit2b"
SIT1=prev.sit1
SEED=int(Q["preregistered_dataset"]["seed"])
STARTS=np.arange(16,112,dtype=np.int32)
SLOPES=(-1.5,-1.,-.5,0.,.5,1.,1.5)
TIMES=np.arange(16,dtype=np.int32)
PATHS=np.stack([
    np.rint(st+slope*(TIMES-7.5)).astype(np.int32)
    for slope in SLOPES for st in STARTS
],axis=0)
assert PATHS.shape==(672,16) and PATHS.min()>=3 and PATHS.max()<=124
T=TIMES[None,:]
CATEGORIES=tuple(Q["preregistered_dataset"]["conditions"])
DEVELOPMENT=set(Q["preregistered_dataset"]["development_scans"])
EVALUATION=set(Q["preregistered_dataset"]["evaluation_scans"])


def robust_normalize(panel):
    x=np.asarray(panel,dtype=np.float64)
    if x.shape!=(16,128) or not np.isfinite(x).all():
        raise ValueError("SIT2B_INVALID_MEASURED_OR_INJECTED_PANEL")
    median=np.median(x,axis=1)
    mad=np.maximum(1.4826*np.median(np.abs(x-median[:,None]),axis=1),.02)
    return (x-median[:,None])/mad[:,None],mad


def generate_overlay(source,kind,rng):
    v=np.array(source,dtype=np.float32,copy=True)
    _,mad=robust_normalize(v)
    truth=None
    if kind==CATEGORIES[0]:
        return v,truth
    if kind==CATEGORIES[1]:
        # Clearly artificial RFI-like perturbation, not independently known RFI.
        gap=int(rng.integers(6,13))
        offset=int(rng.integers(2,gap+1))
        strength=float(rng.uniform(1,2.5))
        v[:,offset::gap]+=np.asarray(strength*mad[:,None],dtype=np.float32)
        for ti in rng.choice(16,3,replace=False):
            v[ti,:]+=float(rng.uniform(1,2))*mad[ti]
        return v,truth
    if kind not in (CATEGORIES[2],CATEGORIES[3]):
        raise ValueError("UNREGISTERED_CATEGORY")
    first=int(rng.integers(34,97))
    slope=float((-1,-.5,0,.5,1)[int(rng.integers(5))])
    truth=np.rint(first+slope*(TIMES-7.5)).astype(np.int32)
    bounds=Q["preregistered_dataset"]["weak_amplitude_range_robust_row_mad"] if kind==CATEGORIES[2] else Q["preregistered_dataset"]["strong_amplitude_range_robust_row_mad"]
    amp=rng.uniform(*bounds,size=16)*mad
    for t,k in enumerate(truth):
        v[t,k]+=amp[t]
        v[t,k-1]+=.2*amp[t]
        v[t,k+1]+=.2*amp[t]
    return v,truth


def all_cases(authentic_panels):
    rng=np.random.default_rng(SEED)
    items=[]
    for source in authentic_panels:
        base=np.asarray(source["panel"],dtype=np.float32)
        for replicate in range(4):
            for category in CATEGORIES:
                v,truth=generate_overlay(base,category,rng)
                row=prev.matched_row(v,source)
                items.append({"row":row,"kind":category,
                              "ground_truth_for_evaluation_only":truth,
                              "source_panel_sha256":source["measured_panel_sha256"],
                              "replicate":replicate})
    assert len(items)==384
    return items


def scores_for_paths(panel,sideband=False,backwards=False):
    norm,_=robust_normalize(panel)
    center=norm[T,PATHS]
    sample=center
    if sideband:
        sample=center-.5*(norm[T,PATHS-3]+norm[T,PATHS+3])
    if backwards:
        # Candidate hypotheses indexed by their terminal frequency in
        # reverse graph; every hypothesis is still scored exactly once.
        order=np.lexsort((np.arange(len(PATHS)),PATHS[:,-1]))
        scores=np.zeros(len(PATHS),dtype=np.float64)
        for index in order:
            # Accumulate from observed endpoint C (t=15) to earliest t=0.
            scores[index]=float(sample[index,::-1].sum()/16.)
    else:
        scores=sample.mean(axis=1)
    return scores


def best_path(scores):
    # Stable canonical tie rule independent of search order and end-frequency.
    # Rounding only resolves floating-point associativity at <=1e-10.
    idx=int(np.argmax(np.round(np.asarray(scores,dtype=np.float64),10)))
    return PATHS[idx].copy(),float(scores[idx]),idx


def direct_matched(panel,sideband=False):
    return best_path(scores_for_paths(panel,sideband=sideband,backwards=False))


def reverse_sideband_verified(row,model_lock):
    start=time.perf_counter()
    inverse_path,inverse_score,index=best_path(scores_for_paths(row["panel"],True,True))
    time_reverse=time.perf_counter()-start
    start=time.perf_counter()
    # Independent forward-time replay against identical astronomical pixels.
    fpath,fscore,fidx=direct_matched(row["panel"],True)
    if fidx!=index or not np.array_equal(fpath,inverse_path) or not np.isclose(fscore,inverse_score,rtol=0,atol=1e-9):
        raise AssertionError("SIT2B_BACKWARD_FORWARD_SIDE_BAND_DISAGREEMENT")
    m=SIT1.manifest(row,inverse_path,inverse_score,model_lock)
    mirror=mirror_case(m["graph"])
    verified=mirrored_replay(mirror,budget=500)
    if verified["prediction"]==ABSTAIN:
        raise AssertionError("SIT2B_TACHYON_PROVENANCE_GRAPH_FAILED")
    if (m["source_file_sha256"]!=row["source_file_sha256"]
        or m["panel_sha256"]!=SIT1.hash_array(row["panel"],dtype="<f4")
        or m["candidate_path_sha256"]!=SIT1.hash_array(inverse_path,dtype="<i4")):
        raise AssertionError("SIT2B_INPUT_PROVENANCE_DIGEST_FAILED")
    replay_time=time.perf_counter()-start
    return inverse_path,inverse_score,{"reverse_seconds":time_reverse,
                                      "forward_and_R5F_replay_seconds":replay_time,
                                      "graph_operations":verified["operations"]}


def decide_threshold(developer_negatives):
    values=np.asarray(developer_negatives,dtype=np.float64)
    if values.size!=64:raise ValueError("SIT2B_DEVELOPER_NEGATIVE_COUNT_MISMATCH")
    return float(np.quantile(values,.95,method="higher"))


def localizes(path,truth):
    if truth is None:return False
    return bool((np.abs(path-truth)<=2).sum()>=12)


def evaluate(rows):
    """Frozen thresholds from dev negative only; no heldout adaptation."""
    labels=(
      "A_NATIVE_COSMOS_RIDGE","B_FROZEN_COSMOS_SHADOW",
      "C_STRONG_DIRECT_MATCHED_FILTER","D_STRONG_DIRECT_SIDEBAND",
      "E_TACHYON_REVERSE_SIDEBAND"
    )
    dev=[x for x in rows if x["split"]=="development"]
    eva=[x for x in rows if x["split"]=="evaluation"]
    assert len(dev)==128 and len(eva)==256
    out={}
    for arm in labels:
        n=[z["scores"][arm] for z in dev if z["kind"] in CATEGORIES[:2]]
        threshold=decide_threshold(n)
        weak=[z for z in eva if z["kind"]==CATEGORIES[2]]
        strong=[z for z in eva if z["kind"]==CATEGORIES[3]]
        ctrl=[z for z in eva if z["kind"] in CATEGORIES[:2]]
        detected_weak={z["case_id"] for z in weak if z["scores"][arm]>threshold and z["located"][arm]}
        detected_strong={z["case_id"] for z in strong if z["scores"][arm]>threshold and z["located"][arm]}
        raised_ctrl={z["case_id"] for z in ctrl if z["scores"][arm]>threshold}
        path_weak=sum(z["located"][arm] for z in weak)
        path_strong=sum(z["located"][arm] for z in strong)
        qualified=len(detected_strong)>=56 and len(raised_ctrl)<=8
        out[arm]={
            "frozen_development_negative_p95_threshold":threshold,
            "heldout_weak_count":len(weak),"heldout_strong_count":len(strong),
            "heldout_controls_count":len(ctrl),
            "heldout_weak_path_localized":int(path_weak),
            "heldout_strong_path_localized":int(path_strong),
            "heldout_weak_recovered":len(detected_weak),
            "heldout_strong_recovered":len(detected_strong),
            "heldout_synthetic_controls_above_threshold":len(raised_ctrl),
            "positive_guard_qualified":qualified,
            "weak_positive_id_digest":hashlib.sha256(json.dumps(sorted(detected_weak)).encode()).hexdigest(),
        }
        out[arm]["_internal_weak_detected_ids"]=detected_weak
    # Internally compare paired exact heldout identifiers but do not export IDs.
    a=out["A_NATIVE_COSMOS_RIDGE"]["_internal_weak_detected_ids"]
    b=out["B_FROZEN_COSMOS_SHADOW"]["_internal_weak_detected_ids"]
    c=out["C_STRONG_DIRECT_MATCHED_FILTER"]["_internal_weak_detected_ids"]
    d=out["D_STRONG_DIRECT_SIDEBAND"]["_internal_weak_detected_ids"]
    e=out["E_TACHYON_REVERSE_SIDEBAND"]["_internal_weak_detected_ids"]
    gains={
        "reverse_extra_weak_vs_native_ridge":len(e-a),
        "reverse_extra_weak_vs_frozen_shadow":len(e-b),
        "reverse_extra_weak_vs_exact_matched_filter":len(e-c),
        "reverse_extra_weak_vs_same_objective_strong_direct_sideband":len(e-d),
        "reverse_missing_weak_vs_strong_direct_sideband":len(d-e),
        "reverse_vs_strong_direct_same_objective_identical_weak_ids":d==e,
        "all_arms_positive_guards_pass":all(z["positive_guard_qualified"] for z in out.values()),
        "reverse_scientific_claim_allowed":False,
    }
    for z in out.values():
        z.pop("_internal_weak_detected_ids")
    # A true TACHYON-specific advantage cannot be announced without the
    # preregistered >=5 new weak and equal controls / independent evidence.
    r=out["E_TACHYON_REVERSE_SIDEBAND"]
    strong=out["D_STRONG_DIRECT_SIDEBAND"]
    gains["tachyon_specific_benchmark_gain_gate"]=bool(
        gains["reverse_extra_weak_vs_same_objective_strong_direct_sideband"]>=5 and
        r["heldout_synthetic_controls_above_threshold"]<=strong["heldout_synthetic_controls_above_threshold"] and
        r["positive_guard_qualified"] and strong["positive_guard_qualified"])
    return out,gains


def run():
    receipts,authentic,errors,download_seconds=prev.cases_from_authentic_sources()
    result={
      "schema":"janus.cosmos.sit2b.real_gbt_temporal_accumulation_result.v1",
      "status":"DATA_SOURCE_HOLD" if errors else "REAL_GBT_BACKGROUNDS_WITH_SYNTHETIC_LINES_SIT2B_RESEARCH_COMPUTED",
      "actual_radio_observatory":"Green Bank Telescope",
      "obs_date":"2025-12-18",
      "real_gbt_sources_hash_verified":not bool(errors),
      "downloads":[{"scan":x["scan"],"role":x["role"],"sha256":x.get("first_seen_sha256"),
                    "bytes":x.get("bytes"),"error":x.get("error")} for x in receipts],
      "errors":errors,"real_background_count":len(authentic),
      "archive_download_seconds":download_seconds,
      "original_measurements_modified":False,
      "production_model_or_threshold_modified":False,
      "model_training":False,
      "independent_sky_signal_truth":False,
      "independent_observing_cadences":1,
      "real_false_alarm_reduction":None,
      "independent_actual_antenna_pointing":False,
      "physical_time_reversal":False,
      "scientific_gate":"HOLD_ONLY_ONE_CADENCE_SYNTHETIC_LABELS_OOD_NATIVE_LOW_RES_AND_NO_INDEPENDENT_SKY_TRUTH",
    }
    if errors:
        return result
    native_model=SIT1.cosmos.frozen_model()
    model_lock=SIT1.cosmos.Q["frozen_model"]["fit_lock_sha256"]
    raw=all_cases(authentic)
    t0=time.perf_counter()
    rows=[]
    summed_times={"native_cosmos_seconds":0.,"strong_matched_seconds":0.,
                  "strong_direct_sideband_seconds":0.,
                  "tachyon_reverse_sideband_seconds":0.,
                  "tachyon_forward_provenance_replay_seconds":0.}
    graph_ops=0
    for case_id,item in enumerate(raw):
        source=item["row"];v=source["panel"]
        now=time.perf_counter()
        npth,ns,rg=SIT1.score_one(v,native_model)
        summed_times["native_cosmos_seconds"]+=time.perf_counter()-now
        now=time.perf_counter()
        mp,ms,m_idx=direct_matched(v,False)
        summed_times["strong_matched_seconds"]+=time.perf_counter()-now
        now=time.perf_counter()
        dp,ds,d_idx=direct_matched(v,True)
        summed_times["strong_direct_sideband_seconds"]+=time.perf_counter()-now
        ep,es,verification=reverse_sideband_verified(source,model_lock)
        summed_times["tachyon_reverse_sideband_seconds"]+=verification["reverse_seconds"]
        summed_times["tachyon_forward_provenance_replay_seconds"]+=verification["forward_and_R5F_replay_seconds"]
        graph_ops+=verification["graph_operations"]
        # All computations operate on exactly one identical blinded panel.
        if not np.array_equal(dp,ep) or not np.isclose(ds,es,rtol=0,atol=1e-9):
            raise AssertionError("REVERSE_AND_SAME_OBJECTIVE_DIRECT_MATCHED_FILTER_DISAGREE")
        if not np.array_equal(npth,mp) or not np.isclose(rg,ms,rtol=0,atol=1e-9):
            raise AssertionError("NATIVE_RIDGE_AND_IDENTICAL_MF_HYPOTHESIS_GRID_DISAGREE")
        kind=item["kind"]
        truth=item["ground_truth_for_evaluation_only"]
        scores={
          "A_NATIVE_COSMOS_RIDGE":rg,
          "B_FROZEN_COSMOS_SHADOW":ns,
          "C_STRONG_DIRECT_MATCHED_FILTER":ms,
          "D_STRONG_DIRECT_SIDEBAND":ds,
          "E_TACHYON_REVERSE_SIDEBAND":es,
        }
        paths={
          "A_NATIVE_COSMOS_RIDGE":npth,
          "B_FROZEN_COSMOS_SHADOW":npth,
          "C_STRONG_DIRECT_MATCHED_FILTER":mp,
          "D_STRONG_DIRECT_SIDEBAND":dp,
          "E_TACHYON_REVERSE_SIDEBAND":ep,
        }
        # Endpoints & labels not consumed by any preceding score operation.
        rows.append({
          "case_id":case_id,
          "scan":source["scan"],
          "split":"development" if source["scan"] in DEVELOPMENT else "evaluation",
          "kind":kind,"scores":scores,
          "located":{k:localizes(path,truth) for k,path in paths.items()},
        })
    s,g=evaluate(rows)
    result["dataset"]={
      "real_backgrounds":len(authentic),"cases":len(rows),
      "development_cases":sum(x["split"]=="development" for x in rows),
      "heldout_cases":sum(x["split"]=="evaluation" for x in rows),
      "weak_heldout":64,"strong_heldout":64,"controls_heldout":128,
      "injected_signal_truth_real_sky":False,
      "per_panel_repetitions":4,"same_cadence_dependent_samples":True,
    }
    result["outcomes"]=s
    result["paired_gains"]=g
    result["cost"]={
      "same_hypotheses_per_match_filter":672,
      "same_spectral_samples_per_path_direct_and_reverse":16,
      "direct_and_reverse_sideband_ops_equal_excluding_reverse_provenance":True,
      "extra_R5F_mirrored_graph_ops_total":graph_ops,
      "wall_seconds_exploratory":{k:round(v,4) for k,v in summed_times.items()},
      "complete_model_and_matched_filter_benchmark_time_excludes_file_download":round(time.perf_counter()-t0,4),
    }
    result["score_predictions_independently_verified_all_cases"]=True
    result["original_native_candidate_and_frozen_model_unmodified"]=True
    result["status"]=(
        "REAL_GBT_INJECTION_EXAM_POSITIVE_CONTROLS_PASS__NO_TACHYON_SPECIFIC_GAIN__SCIENCE_HOLD"
        if g["all_arms_positive_guards_pass"] and not g["tachyon_specific_benchmark_gain_gate"]
        else "REAL_GBT_INJECTION_EXAM_COMPUTED__POSITIVE_OR_FALSE_CONTROL_QUALIFICATION_HOLD__NO_TACHYON_SPECIFIC_GAIN"
    )
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--audit",action="store_true")
    a=ap.parse_args()
    if not a.audit:ap.error("--audit required")
    OUT.mkdir(parents=True,exist_ok=True)
    result=run()
    (OUT/"sit2b_result.json").write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("SIT2B_RESULT",json.dumps({
       "status":result["status"],"errors":result.get("errors",[]),
       "outcomes":result.get("outcomes",{}),"paired_gains":result.get("paired_gains",{})},
       sort_keys=True,ensure_ascii=False),flush=True)
    if result.get("outcomes"):
        assert result["dataset"]["cases"]==384
        assert result["dataset"]["heldout_cases"]==256
        assert result["score_predictions_independently_verified_all_cases"]
        assert result["paired_gains"]["reverse_vs_strong_direct_same_objective_identical_weak_ids"]
        assert result["paired_gains"]["reverse_extra_weak_vs_same_objective_strong_direct_sideband"]==0
        print("SIT2B_REAL_GBT_TEMPORAL_INJECTION_EXAM_REPRODUCIBLE__SCIENTIFIC_HOLD",flush=True)
    else:
        print("SIT2B_REAL_SOURCE_FAIL_CLOSED_HOLD",flush=True)


if __name__=="__main__":
    main()
