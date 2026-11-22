"""JANUS COSMOS R42A8A × TACHYON R5F STOCK MIRROR — research sidecar.

Real reusable code:
- frozen COSMOS R42A7A shadow scoring and blind R42A8A line search
- pinned DemiHead R5F mirror_case / mirrored_exact / mirrored_replay

Test inputs are 1024 internally SIMULATED time-frequency arrays, not telescope
observations. A→B→C is data/computation lineage: simulated spectrum, blind
candidate trajectory, frozen score. C→B→A is backward GRAPH TRACING, not an
astronomical or physical backward-in-time signal.

The COSMOS model, candidate trajectories and decisions are NEVER modified.
Metadata faults are injected after the frozen scores exist; these faults are
only test fixtures and are not proof of real telescope RFI or data tampering.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import random
import struct
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(HERE / "_tachyon"))
sys.path.insert(0,str(HERE / "research" / "seti"))
import r42a8_ood_stress as cosmos
import r42a7_shadow_rfi as shadow
from tools.tachyon_r5f_stock_mirror_audit import (
    mirror_case, mirrored_exact, mirrored_replay, transpose,
)
from tools.tachyon_r5f_causal_memory_audit import ABSTAIN

PROTOCOL_PATH=Path("research/seti/cosmos_tachyon_r5f_bridge_protocol.json")
PROTOCOL=json.loads(PROTOCOL_PATH.read_text())
CATEGORIES=tuple(PROTOCOL["adversarial_overlay"]["types"])
ASSERTED_MODEL_SHA=PROTOCOL["source_pins"]["frozen_model_source"]
BUDGET=500
OUTPUT=Path("runtime/cosmos_tachyon_r5f")
OUTPUT.mkdir(parents=True,exist_ok=True)


def sha(data:bytes)->str:
    return hashlib.sha256(data).hexdigest()


def spectrum_sha(spectrum:np.ndarray)->str:
    return sha(np.ascontiguousarray(spectrum,dtype=np.float32).tobytes())


def path_sha(path:np.ndarray)->str:
    return sha(np.ascontiguousarray(path,dtype=np.int32).tobytes())


def score_sha(score:float)->str:
    return sha(struct.pack(">d",float(score)))


def original_manifest(i:int,spectrum:np.ndarray,path:np.ndarray,score:float,lock:str)->dict:
    """No simulated injection/control labels or expected path in inference input."""
    ident=f"{i:04d}"
    raw_node="raw:"+ident
    path_node="path:"+ident
    score_node="score:"+ident
    source_graph={
        "schema":"janus.tachyon.r5f.synthetic_reverse_road_graph.v1",
        "case_ref":ident,
        "terminal":score_node,
        # Observed here means an OBSERVED COMPUTED OUTPUT, never a sky event.
        "terminal_kind":"OBSERVED",
        "forward_edges":{
            "A_B":{raw_node:[path_node]},
            "B_C":{path_node:[score_node]},
        },
        "reverse_proposals":{
            "C_B":{score_node:[path_node]},
            "B_A":{path_node:[raw_node]},
            "complete_certified":True,
        },
        "witness_reports":[],
        "witness_window_seconds":120,
        "synthetic_only":True,
    }
    return {
        "schema":"janus.cosmos.tachyon_r5f.synthetic_lineage_manifest.v1",
        "case_ref":ident,
        "raw_sha256":spectrum_sha(spectrum),
        "candidate_path_sha256":path_sha(path),
        "frozen_score_sha256":score_sha(score),
        "frozen_model_lock_sha256":lock,
        "graph":source_graph,
        "synthetic_only":True,
    }


def tamper_metadata(original:dict,category:str)->dict:
    rec=json.loads(json.dumps(original))
    if category=="INTACT":
        pass
    elif category=="BAD_RAW_DIGEST":
        rec["raw_sha256"]="0"*64
    elif category=="BAD_PATH_DIGEST":
        rec["candidate_path_sha256"]="0"*64
    elif category=="BAD_SCORE_DIGEST":
        rec["frozen_score_sha256"]="0"*64
    elif category=="BAD_MODEL_LOCK":
        rec["frozen_model_lock_sha256"]="0"*64
    elif category=="MISSING_RAW_DIGEST":
        rec.pop("raw_sha256")
    elif category=="INCOMPLETE_PREDECESSOR_INDEX":
        rec["graph"]["reverse_proposals"]["complete_certified"]=False
    elif category=="GHOST_REVERSE_PREDECESSOR":
        node=next(iter(rec["graph"]["reverse_proposals"]["B_A"]))
        rec["graph"]["reverse_proposals"]["B_A"][node].append("raw:unverified-extraneous")
    else:
        raise ValueError("unexpected tampering fixture category")
    # The category is not in rec. Actual labels and true injection paths
    # are also not present anywhere in rec.
    return rec


def overlay_assignments(n:int=1024)->list[str]:
    if n%len(CATEGORIES):
        raise ValueError("balanced metadata overlay requires divisible rows")
    names=[x for x in CATEGORIES for _ in range(n//len(CATEGORIES))]
    random.Random(int(PROTOCOL["adversarial_overlay"]["seed"])).shuffle(names)
    return names


def _graph_consistency(graph:dict)->bool:
    if graph.get("synthetic_only") is not True:
        return False
    ab=graph.get("forward_edges",{}).get("A_B",{})
    bc=graph.get("forward_edges",{}).get("B_C",{})
    proposed=graph.get("reverse_proposals",{})
    if not proposed.get("complete_certified"):
        return False
    return (proposed.get("C_B")==transpose(bc) and
            proposed.get("B_A")==transpose(ab))


def _metadata_check(record:dict,spectrum:np.ndarray,path:np.ndarray,score:float,lock:str)->bool:
    if record.get("synthetic_only") is not True:
        return False
    if record.get("raw_sha256")!=spectrum_sha(spectrum):
        return False
    if record.get("candidate_path_sha256")!=path_sha(path):
        return False
    if record.get("frozen_score_sha256")!=score_sha(score):
        return False
    if record.get("frozen_model_lock_sha256")!=lock:
        return False
    return _graph_consistency(record.get("graph",{}))


def cosmos_only(rec,arr,path,score,lock)->tuple[bool,int]:
    """Baseline: evaluate spectrum and score; no provenance requirements."""
    return True,0


def mirror_graph_only(rec,arr,path,score,lock)->tuple[bool,int]:
    """Weak ablation: graph isomorphism without source/score authenticity."""
    mir=mirror_case(rec["graph"])
    result=mirrored_exact(mir,budget=BUDGET)
    return result["prediction"]!=ABSTAIN,result["operations"]


def mirror_replay(rec,arr,path,score,lock)->tuple[bool,int]:
    mir=mirror_case(rec["graph"])
    outcome=mirrored_replay(mir,budget=BUDGET)
    if outcome["prediction"]==ABSTAIN:
        return False,outcome["operations"]
    return _metadata_check(rec,arr,path,score,lock),outcome["operations"]


def strong_direct(rec,arr,path,score,lock)->tuple[bool,int]:
    """Direct full audit, same accessible evidence, no TACHYON helper."""
    graph=rec["graph"]
    if not _metadata_check(rec,arr,path,score,lock):
        return False,0
    src=next(iter(graph["forward_edges"]["A_B"]))
    ab=graph["forward_edges"]["A_B"].get(src,[])
    terminal=graph["terminal"]
    verified=0
    for mid in ab:
        for end in graph["forward_edges"]["B_C"].get(mid,[]):
            if end==terminal:
                verified+=1
    return verified==1,len(ab)+sum(len(graph["forward_edges"]["B_C"].get(mid,[])) for mid in ab)


def always_abstain(rec,arr,path,score,lock)->tuple[bool,int]:
    return False,0


ARMS={
    "A_COSMOS_ONLY_NO_PROVENANCE":cosmos_only,
    "B_TACHYON_MIRROR_GRAPH_ONLY":mirror_graph_only,
    "C_TACHYON_MIRROR_FORWARD_REPLAY":mirror_replay,
    "D_STRONG_DIRECT_PROVENANCE":strong_direct,
    "E_ALWAYS_ABSTAIN":always_abstain,
}


def run():
    if PROTOCOL["source_pins"]["tachyon_operation_budget_per_case"]!=BUDGET:
        raise AssertionError("TACHYON base resource budget drifted")
    # The archive is pinned by SHA inside cosmos.frozen_model(), not our code.
    start_native=time.perf_counter()
    mean,scale,coef,intercept=cosmos.frozen_model()
    import r42a8_ood_stress as frozen_source
    frozen_source.Q["panel"]["seed"] == PROTOCOL["dataset"]["from_original_r42a8_seed"] or (_ for _ in ()).throw(ValueError("seed drift"))
    rng=np.random.default_rng(int(cosmos.Q["panel"]["seed"]))
    data=cosmos.batch(rng)
    v=data["V"]
    candidate_path,ridge=cosmos.blind_line_candidate(v)
    features=shadow.features(v,candidate_path)
    frozen_score=((features-mean)/scale)@coef+intercept
    native_duration=time.perf_counter()-start_native
    assert v.shape==(1024,16,128) and candidate_path.shape==(1024,16)
    assert frozen_score.shape==(1024,) and np.isfinite(frozen_score).all()
    # Eval labels are not used until AFTER all algorithms' accept/reject calls.
    del features
    conditions=overlay_assignments(len(v))
    stats={a:Counter() for a in ARMS}
    elapsed={a:0. for a in ARMS}
    mismatches={a:[] for a in ARMS}
    source_lock=cosmos.Q["frozen_model"]["fit_lock_sha256"]
    valid_scores_pre=sha(np.ascontiguousarray(frozen_score,dtype=np.float64).tobytes())
    valid_paths_pre=sha(np.ascontiguousarray(candidate_path,dtype=np.int32).tobytes())
    raw_input_stream=hashlib.sha256()
    for i,kind in enumerate(conditions):
        raw=original_manifest(i,v[i],candidate_path[i],frozen_score[i],source_lock)
        rec=tamper_metadata(raw,kind)
        assert "label" not in rec and "truth" not in rec and "attack_category" not in rec
        raw_input_stream.update(json.dumps(rec,sort_keys=True,separators=(",",":")).encode())
        accepted_gold=(kind=="INTACT")
        for name,fn in ARMS.items():
            t0=time.perf_counter()
            ok,ops=fn(rec,v[i],candidate_path[i],frozen_score[i],source_lock)
            elapsed[name]+=time.perf_counter()-t0
            s=stats[name]
            s["total"]+=1
            s["accepted"]+=int(ok)
            s["intact_accepted"]+=int(ok and accepted_gold)
            s["intact_rejected"]+=int(not ok and accepted_gold)
            s["corrupt_accepted"]+=int(ok and not accepted_gold)
            s["corrupt_rejected"]+=int(not ok and not accepted_gold)
            s["graph_operations"]+=ops
            if ok!=accepted_gold and len(mismatches[name])<6:
                mismatches[name].append(raw["case_ref"])
    assert sha(np.ascontiguousarray(frozen_score,dtype=np.float64).tobytes())==valid_scores_pre
    assert sha(np.ascontiguousarray(candidate_path,dtype=np.int32).tobytes())==valid_paths_pre
    # Only now consult simulated category labels for descriptive native AUC.
    aucs=[]
    for k in (0,1):
        ix=data["group"]==k
        aucs.append({
            "regime":cosmos.Q["panel"]["regimes"][k],
            "rows":int(ix.sum()),
            "synthetic_positives":int(data["label"][ix].sum()),
            "synthetic_controls":int((~data["label"][ix]).sum()),
            "frozen_R42A7A_shadow_auc":cosmos.calc_auc(data["label"][ix],frozen_score[ix]),
            "ridge_baseline_auc":cosmos.calc_auc(data["label"][ix],ridge[ix]),
        })
    paired={}
    for name,s in stats.items():
        paired[name]={
            "cases":s["total"],
            "intact_accepted":s["intact_accepted"],
            "intact_rejected":s["intact_rejected"],
            "corrupt_accepted":s["corrupt_accepted"],
            "corrupt_rejected":s["corrupt_rejected"],
            "metadata_acceptance_coverage":round(s["accepted"]/s["total"],6),
            "mean_graph_operations_per_case":round(s["graph_operations"]/s["total"],4),
            "postprocess_seconds":round(elapsed[name],4),
            "error_case_refs":mismatches[name],
        }
    c=paired["C_TACHYON_MIRROR_FORWARD_REPLAY"]
    d=paired["D_STRONG_DIRECT_PROVENANCE"]
    return {
        "schema":"janus.cosmos.tachyon_r5f_mirror_sidecar_result.v1",
        "status":"SYNTHETIC_COSMOS_TACHYON_PROVENANCE_ENGINEERING_ONLY__SCIENCE_HOLD",
        "source_model":"Actual R42A7A frozen shadow model through verified R42A8A artifact lock; 1024 native simulated spectrogram candidates",
        "tachyon_code":"DemiHead R5F STOCK MIRROR pinned at 91ef855a6f7119a733990c8f3adeb478e79beb05",
        "source_r42a8_seed":cosmos.Q["panel"]["seed"],
        "test_cases":len(v),
        "original_synthetic_positive_count":int(data["label"].sum()),
        "original_synthetic_control_count":int((~data["label"]).sum()),
        "metadata_overlay_types":{q:conditions.count(q) for q in CATEGORIES},
        "native_regimes":aucs,
        "native_model_and_candidate_path_arrays_unchanged":True,
        "native_model_pipeline_seconds_including_authenticated_model_artifact_fetch":round(native_duration,4),
        "source_lock_sha256":source_lock,
        "frozen_score_array_sha256":valid_scores_pre,
        "candidate_path_array_sha256":valid_paths_pre,
        "manifest_stream_sha256":raw_input_stream.hexdigest(),
        "mirror_only_first_party_source_functions":["mirror_case","mirrored_exact","mirrored_replay"],
        "arms":paired,
        "mirror_plus_replay_detects_all_injected_bad_metadata":c["corrupt_accepted"]==0,
        "mirror_plus_replay_retains_all_intact":c["intact_accepted"]==128,
        "mirror_plus_replay_equivalent_to_strong_direct":all(c[k]==d[k] for k in ("intact_accepted","intact_rejected","corrupt_accepted","corrupt_rejected")),
        "algorithmic_gain_over_strong_direct":False,
        "production_decisions_changed":0,
        "historical_318_localized_69_weak_tested":False,
        "actual_false_alarms_reduced":0,
        "actual_protected_signals_lost":0,
        "safety_note":"The 0 loss/reduction counters mean no production decisions were edited, not that 318/69 historical detections were freshly examined.",
        "real_telescope_evidence":False,
        "real_ONOFF_controls":False,
        "external_custodian":False,
        "trained_model_modified":False,
        "physical_time_reversal":False,
        "scientific_gate":"HOLD_UNTIL_REAL_INDEPENDENT_PROVENANCE_AND_END_TO_END_COSMOS_TEST",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--audit",action="store_true")
    args=p.parse_args()
    if not args.audit:p.error("--audit required")
    result=run()
    print("COSMOS_TACHYON_R5F_RESULT",json.dumps(result,sort_keys=True,allow_nan=False),flush=True)
    OUTPUT.joinpath("cosmos_tachyon_r5f_result.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    a=result["arms"]["C_TACHYON_MIRROR_FORWARD_REPLAY"]
    b=result["arms"]["D_STRONG_DIRECT_PROVENANCE"]
    assert result["test_cases"]==1024
    assert a["intact_accepted"]==128 and a["corrupt_accepted"]==0
    assert a["intact_rejected"]==0 and a["corrupt_rejected"]==896
    assert result["mirror_plus_replay_equivalent_to_strong_direct"]
    assert b["intact_accepted"]==128 and b["corrupt_accepted"]==0
    assert result["native_model_and_candidate_path_arrays_unchanged"]
    print("COSMOS_X_TACHYON_R5F_INTEGRATION_SYNTHETIC_PROVENANCE_PASS__SCIENCE_HOLD",flush=True)


if __name__=="__main__":
    main()
