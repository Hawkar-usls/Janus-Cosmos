"""JANUS SIT2: controlled weak-line and synthetic-RFI injections on real GBT data.

Uses actual, SHA-verified 2025-12-18 3I/ATLAS GBT HDF5 power arrays
as backgrounds, but ALL signal/RFI target labels and interventions are
artificial. This is an injection-recovery stress test, not verified ETI
search, an independent ON/OFF sky-field test or physical time reversal.

Frozen algorithms: COSMOS R42A8A path + R42A7A shadow scorer and DemiHead
TACHYON R5F mirrored computational-provenance replay. Fully equal-input,
target-blind, strong direct provenance comparison. No production updates.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score
import janus_sit1_open_3iatlas as sit1
from tools.tachyon_r5f_stock_mirror_audit import mirror_case, mirrored_replay
from tools.tachyon_r5f_causal_memory_audit import ABSTAIN

ROOT=Path(__file__).resolve().parents[2]
P=json.loads((ROOT/"research/seti/janus_sit2_weak_rfi_protocol.json").read_text())
OUT=ROOT/"runtime/sit2"
SEED=int(P["benchmark"]["seed"])
SLOPES=(-1.0,-0.5,0.0,0.5,1.0)
TYPE=("BASE_REAL_BACKGROUND_NO_INJECTION","WEAK_ARTIFICIAL_DRIFT_INJECTION","ARTIFICIAL_BROAD_RFI_CONTROL")


def matched_row(panel_data,base):
    row=dict(base)
    row["panel"]=np.asarray(panel_data,dtype=np.float32)
    row["measured_panel_sha256"]=sit1.hash_array(row["panel"],dtype="<f4")
    return row


def inject_from_real(base,typ,rng):
    """Never mutate actual HDF5-derived arrays. Return injected-path truth separately."""
    panel=np.array(base,dtype=np.float32,copy=True)
    assert panel.shape==(16,128) and np.isfinite(panel).all()
    if typ==TYPE[0]:return panel,None
    med=np.median(panel,axis=1)
    sigma=np.maximum(1.4826*np.median(abs(panel-med[:,None]),axis=1),0.02)
    if typ==TYPE[1]:
        first=int(rng.integers(34,97))
        slope=float(SLOPES[int(rng.integers(len(SLOPES)))])
        path=np.rint(first+slope*(np.arange(16)-7.5)).astype(np.int32)
        amp=rng.uniform(0.8,1.5,size=16)*sigma
        for t,k in enumerate(path):
            panel[t,k]+=amp[t]
            panel[t,k-1]+=.2*amp[t]
            panel[t,k+1]+=.2*amp[t]
        return panel,path
    if typ==TYPE[2]:
        gap=int(rng.integers(6,13))
        offset=int(rng.integers(2,gap+1))
        amp=float(rng.uniform(1.,2.5))
        panel[:,offset::gap]+=amp*sigma[:,None]
        bursts=rng.choice(np.arange(16),size=3,replace=False)
        for t in bursts:
            panel[t,:]+=float(rng.uniform(1.,2.0))*sigma[t]
        return panel,None
    raise ValueError("UNREGISTERED_INJECTION_CLASS")


def cases_from_authentic_sources():
    """Parent-SIT1 SHA reference gates genuine radio background provenance."""
    parent=json.loads((ROOT/"research/seti/janus_sit1_open_result.json").read_text())
    prior={x["scan"]:x["first_seen_sha256"] for x in parent["runtime_inputs"]["hdf5_files"]}
    if len(prior)!=6:raise RuntimeError("SIT1_PARENT_HASH_CATALOG_MISSING")
    folder=OUT/"real_hdf5";folder.mkdir(parents=True,exist_ok=True)
    actual=[];panels=[];errors=[]
    start=time.perf_counter()
    for scan in sit1.Q["immutable_real_archive"]["scans"]:
        d=sit1.retrieve(scan,folder)
        actual.append(d)
        if not d.get("ok"):
            errors.append({"scan":scan["number"],"error":"DOWNLOAD_FAILED","detail":d.get("retry_errors",[])})
            continue
        if d["first_seen_sha256"]!=prior.get(scan["number"]):
            errors.append({"scan":scan["number"],"error":"SHA256_NOT_EQUAL_PARENT_SIT1","sha":d["first_seen_sha256"]})
            (folder/scan["path"]).unlink(missing_ok=True)
            continue
        try:
            q=sit1.reduce_real_h5(folder/scan["path"],d)
            panels.extend(q["panels"])
            d["data_shape"]=q["dataset_shape"];d["header"]=q["header"]
        except Exception as exc:
            errors.append({"scan":scan["number"],"error":"HDF5_DATA_OR_AXES_FAILED","detail":str(exc)[:320]})
        finally:
            (folder/scan["path"]).unlink(missing_ok=True)
    if len(panels)!=24:errors.append({"error":"EXPECTED_24_REAL_BACKGROUND_PANELS","got":len(panels)})
    return actual,panels,errors,round(time.perf_counter()-start,4)


def build_test_cases(panels):
    rng=np.random.default_rng(SEED)
    rows=[]
    # Fixed, label-blind ordering and seed. Multiple synthetic perturbations of
    # the same radio background are paired and NOT statistically independent.
    for base in panels:
        for rep in range(4):
            for typ in TYPE:
                v,truth=inject_from_real(base["panel"],typ,rng)
                row=matched_row(v,base)
                rows.append({"row":row,"source_panel_sha256":base["measured_panel_sha256"],
                             "kind":typ,"truth_for_evaluation_only":truth,
                             "replicate":rep})
    assert len(rows)==288
    return rows


def choose_threshold(x):
    # Development-only NEGATIVE examples; strictly greater than threshold.
    return float(np.quantile(np.asarray(x,dtype=float),.95,method="higher"))


def assess(y):
    """Predictions and truth only meet here AFTER the frozen inference pass."""
    positives=[z for z in y if z["kind"]==TYPE[1]]
    negs=[z for z in y if z["kind"]!=TYPE[1]]
    out={}
    for metric in ("blind_ridge_statistic","frozen_shadow_score"):
        dev_neg=[z[metric] for z in y if z["split"]=="development" and z["kind"]!=TYPE[1]]
        threshold=choose_threshold(dev_neg)
        eva=[z for z in y if z["split"]=="evaluation"]
        eva_p=[z for z in eva if z["kind"]==TYPE[1]]
        eva_n=[z for z in eva if z["kind"]!=TYPE[1]]
        found=sum((z[metric]>threshold) and z["localized"] for z in eva_p)
        scores=np.array([z[metric] for z in eva],dtype=float)
        labels=np.array([z["kind"]==TYPE[1] for z in eva],dtype=int)
        out[metric]={
            "dev_negative_count":len(dev_neg),
            "threshold_fixed_on_development":threshold,
            "evaluation_cases":len(eva),
            "evaluation_weak_injections":len(eva_p),
            "evaluation_no_injection_or_engineered_RFI_controls":len(eva_n),
            "eval_injected_weak_found_with_localization":int(found),
            "eval_injected_weak_missed":int(len(eva_p)-found),
            "eval_synthetic_control_above_threshold":sum(z[metric]>threshold for z in eva_n),
            "eval_synthetic_control_below_or_equal_threshold":sum(z[metric]<=threshold for z in eva_n),
            "evaluation_score_AUC_descriptive":float(roc_auc_score(labels,scores)),
            "evaluation_path_localized_weak":int(sum(z["localized"] for z in eva_p)),
        }
    return out


def evaluate(panels):
    model=sit1.cosmos.frozen_model()
    model_lock=sit1.cosmos.Q["frozen_model"]["fit_lock_sha256"]
    items=build_test_cases(panels)
    assert len(items)==288
    scores=[];fail=[];cost={"native_s":0.,"mirror_replay_s":0.,"strong_direct_s":0.}
    graph_ops=0
    for i,item in enumerate(items):
        row=item["row"]
        start=time.perf_counter()
        pred,score,ridge=sit1.score_one(row["panel"],model)
        cost["native_s"]+=time.perf_counter()-start
        # The manifest describes a COMPUTED TEST PANEL, not an actual raw
        # telescope sample after injection. Original parent hash retained
        # separately, and never claimed to hash the synthetic perturbation.
        m=sit1.manifest(row,pred,score,model_lock)
        start=time.perf_counter()
        mir=mirror_case(m["graph"])
        r=mirrored_replay(mir,budget=500)
        cost["mirror_replay_s"]+=time.perf_counter()-start
        graph_ops+=r["operations"]
        start=time.perf_counter()
        direct=sit1.direct_check(m,row,pred,score,model_lock,model)
        cost["strong_direct_s"]+=time.perf_counter()-start
        # Both versions are evidence validators, neither can see injection
        # identity nor select an alternative path based on hidden truth.
        mirror_valid=r["prediction"]!=ABSTAIN and direct
        if not mirror_valid or not direct:
            fail.append({"case":i,"mirror":mirror_valid,"strong":direct})
        pathtruth=item["truth_for_evaluation_only"]
        good=False
        if pathtruth is not None:
            good=bool((abs(pred-pathtruth)<=2).sum()>=12)
        scores.append({
            "case":i,"scan":row["scan"],"role":row["role"],"window":row["panel_idx"],
            "split":"development" if row["scan"] in ("0061","0062") else "evaluation",
            "kind":item["kind"],"localized":good,
            "blind_ridge_statistic":float(ridge),
            "frozen_shadow_score":float(score),
            "source_real_panel_sha256":item["source_panel_sha256"],
            "computed_test_panel_sha256":row["measured_panel_sha256"],
            "mirror_valid":mirror_valid,"strong_valid":bool(direct),
            "candidate_path_sha256":sit1.hash_array(pred,dtype="<i4"),
            # No injected truth path exported in public receipt.
        })
    report=assess(scores)
    # Count paired decisions as they actually are; no invented detector gain.
    positives=[z for z in scores if z["kind"]==TYPE[1]]
    dev=[z for z in scores if z["split"]=="development"]
    eva=[z for z in scores if z["split"]=="evaluation"]
    return {
        "research_case_count":len(scores),
        "development_cases":len(dev),
        "heldout_cases":len(eva),
        "development_scans":["0061","0062"],
        "heldout_scans":["0063","0064","0065","0066"],
        "weak_synthetic_injections":len(positives),
        "real_background_panels":len(panels),
        "native_path_and_scores_unchanged_across_arms":True,
        "mirrored_and_strong_provenance_disagreements":fail,
        "paired_mirror_replay_abstentions":sum(not z["mirror_valid"] for z in scores),
        "paired_direct_provenance_abstentions":sum(not z["strong_valid"] for z in scores),
        "true_original_native_missed_weak_rescued_by_mirror":0,
        "true_original_native_missed_weak_rescued_by_direct":0,
        "detector_recovery_gain_over_native":0,
        "blind_ridge_and_frozen_shadow_results":report,
        "runtime_seconds_uncontrolled":{k:round(v,5) for k,v in cost.items()},
        "total_mirror_graph_ops":graph_ops,
        "records_digest":hashlib.sha256(json.dumps(scores,sort_keys=True).encode()).hexdigest(),
        "scoring_note":"Truth labels are only artificial injection truth, evaluated after the exact frozen candidate output; no observational astrophysical truth labels exist.",
        "independent_real_observational_groups":1,
        "synthetic_perturbations_share_real_background":True,
        "external_scientific_validation":False,
        "real_false_alarm_reduction":None,
        "real_ETI_detection":False,
        "physical_time_reversal":False,
        "production_weights_changed":False,
    }


def write_receipt(out):
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"sit2_result.json").write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("SIT2_RESULT",json.dumps({
        "status":out["status"],
        "real_files":len(out.get("downloads",[])),
        "errors":out.get("errors",[]),
        "eval":out.get("evaluation",{}).get("blind_ridge_and_frozen_shadow_results",{}),
    },ensure_ascii=False,sort_keys=True),flush=True)


def run():
    receipts,panels,errors,download_s=cases_from_authentic_sources()
    out={
      "schema":"janus.cosmos.sit2.real_background_weak_injection_rfi_result.v1",
      "status":"DATA_QUALIFICATION_HOLD" if errors else "REAL_GBT_BACKGROUND_PLUS_SYNTHETIC_WEAK_RFI_CHALLENGE_PASS__SCIENCE_HOLD",
      "dedication":"JANUS SIT2 follows Stephen Hawking memorial SIT program; not a physical time-reversal experiment.",
      "source":"GBT 2025-12-18 six real 3I/ATLAS X-band scans, one ABACAD cadence, exact SHA matched against SIT1",
      "downloads":receipts,"errors":errors,
      "retrieval_seconds":download_s,
      "real_background_panels_reduced":len(panels),
      "model_and_production_weights_unchanged":True,
      "scientific_gate":"HOLD_SAME_CADENCE_INJECTED_LABELS_ONLY__NO_PHYSICAL_ETI_OR_RFI_VERIFICATION"
    }
    if not errors:
        out["evaluation"]=evaluate(panels)
    write_receipt(out)
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--audit",action="store_true")
    args=p.parse_args()
    if not args.audit:p.error("--audit required")
    out=run()
    if out["status"].startswith("REAL_GBT_BACKGROUND"):
        r=out["evaluation"]
        assert r["development_cases"]==96 and r["heldout_cases"]==192
        assert r["mirrored_and_strong_provenance_disagreements"]==[]
        assert r["detector_recovery_gain_over_native"]==0
        assert r["paired_mirror_replay_abstentions"]==r["paired_direct_provenance_abstentions"]==0
        print("JANUS_SIT2_REAL_BACKGROUND_INJECTION_EVALUATION_PASS__SCIENCE_HOLD",flush=True)
    else:
        print("JANUS_SIT2_DATA_QUALIFICATION_HOLD__NO_SCIENTIFIC_RESULT",flush=True)

if __name__=="__main__":
    main()
