"""JANUS SIT0 — Stephen Hawking dedication.

Real-observation-derived natural methanol radio profiles W3OH and G42.03+0.19,
from the SHA256-pinned R37F2D3/R37F3A JANUS COSMOS archive.
TACHYON R5F stock mirror traces an *analysis graph* C->B->A, then checks
actual native R37A rebin and R37F3A spectral metrics forward A->B->C.

The archived inputs are **one-dimensional calibrated-derived averaged
spectra**, NOT native time-frequency radio candidates, not authenticated
spatial ON/OFF controls, and cannot be fed to R42A8A's synthetic 16x128
classifier without an unvalidated adapter. This sidecar does not do that.

All adversarial overlays are synthetic METADATA ONLY. No astrophysical source,
scientific line power, product decision, model weight, or physical timeline
is changed.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

BASE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(BASE / "research" / "seti"))
sys.path.insert(0,str(BASE / "_tachyon"))
import r37f3a_natural_preservation as native
from tools.tachyon_r5f_stock_mirror_audit import mirror_case, mirrored_exact, mirrored_replay, transpose
from tools.tachyon_r5f_causal_memory_audit import ABSTAIN

P=BASE/"research/seti/janus_sit0_hawking_protocol.json"
Q=json.loads(P.read_text())
OVERLAYS=tuple(Q["adversarial"]["overlays"])
OUT=BASE/"runtime/janus_sit0"
BUDGET=500

def digest(b:bytes)->str:
    return hashlib.sha256(b).hexdigest()

def arrays_sha(*arrays:np.ndarray)->str:
    h=hashlib.sha256()
    for a in arrays:
        a=np.ascontiguousarray(a,dtype="<f8")
        h.update(len(a.shape).to_bytes(1,"big"))
        for n in a.shape:h.update(n.to_bytes(8,"big"))
        h.update(a.tobytes())
    return h.hexdigest()

def jsonsha(data:dict)->str:
    return digest(json.dumps(data,sort_keys=True,allow_nan=False,separators=(",",":")).encode())

def edge_rebin(native_edges, native_profile, transform, native_rebin):
    delta=float(np.median(np.diff(native_edges)))*transform["factor"]
    start=native_edges[0]+delta*transform["phase"]
    inner=np.arange(start,native_edges[-1],delta)
    inner=inner[(inner>native_edges[0]+1e-9)&(inner<native_edges[-1]-1e-9)]
    edges=np.r_[native_edges[0],inner,native_edges[-1]]
    rebinned,err=native_rebin(native_profile[None,:],np.array([0.,1.]),native_edges,np.array([0.,1.]),edges)
    assert rebinned.shape==(1,len(edges)-1)
    assert np.isfinite(rebinned).all()
    assert err<=native.L["conservation_rel"]
    return edges,rebinned[0].copy(),float(err)

def manifest(ctx:dict)->dict:
    tag=ctx["tag"]
    A="profile:"+tag
    B="conservative_rebin:"+tag
    C="natural_line_metrics:"+tag
    return {
      "schema":"janus.sit0.hawking.real_profile_lineage_manifest.v1",
      "tag":tag,
      "source_id":ctx["source"],
      "mjd_day":ctx["mjd_day"],
      "profile_name":ctx["profile_name"],
      "source_archive_artifact_id":Q["source"]["parent_artifact_id"],
      "original_profile_npz_sha256":ctx["npz_sha"],
      "original_numeric_spectrum_sha256":arrays_sha(ctx["v"],ctx["y"]),
      "native_rebin_arrays_sha256":arrays_sha(ctx["f"],ctx["z"]),
      "native_preservation_result_sha256":jsonsha(ctx["compare"]),
      "pinned_rebin_source_blob_sha1":Q["source"]["native_rebin_source_expected_sha"],
      "transform_settings":dict(ctx["transform"]),
      "forward_lineage_graph":{
        "schema":"janus.tachyon.r5f.synthetic_reverse_road_graph.v1",
        "case_ref":tag,
        "terminal":C,
        "terminal_kind":"OBSERVED",
        "forward_edges":{"A_B":{A:[B]},"B_C":{B:[C]}},
        "reverse_proposals":{"C_B":{C:[B]},"B_A":{B:[A]},"complete_certified":True},
        "witness_reports":[],  # NO independent physical witness roots invented.
        "witness_window_seconds":120,
        "synthetic_only":True, # graph construction, NOT radio measurement
      },
    }

def corrupt(m:dict,kind:str)->dict:
    x=copy.deepcopy(m)
    if kind=="INTACT":return x
    if kind=="BAD_RAW_SHA":x["original_numeric_spectrum_sha256"]="0"*64
    elif kind=="BAD_REBIN_SHA":x["native_rebin_arrays_sha256"]="0"*64
    elif kind=="BAD_METRIC_SHA":x["native_preservation_result_sha256"]="0"*64
    elif kind=="WRONG_SOURCE_ID":x["source_id"]="WRONG_SKY_SOURCE"
    elif kind=="MISSING_SOURCE_ID":x.pop("source_id")
    elif kind=="INCOMPLETE_REVERSE_INDEX":x["forward_lineage_graph"]["reverse_proposals"]["complete_certified"]=False
    elif kind=="GHOST_REVERSE_EDGE":
        props=x["forward_lineage_graph"]["reverse_proposals"]["B_A"]
        b=next(iter(props))
        props[b].append("profile:UNVERIFIED_DECOY")
    else:raise ValueError("Unknown preregistered corruption")
    return x

def trusted_provenance_valid(m:dict,ctx:dict,rebin):
    """Authoritative direct provenance & TRUE forward computation, no labels."""
    if m.get("source_archive_artifact_id")!=Q["source"]["parent_artifact_id"]:
        return False
    for field,key in [("source_id","source"),("mjd_day","mjd_day"),("profile_name","profile_name"),
                      ("original_profile_npz_sha256","npz_sha")]:
        if m.get(field)!=ctx[key]:
            return False
    if m.get("pinned_rebin_source_blob_sha1")!=Q["source"]["native_rebin_source_expected_sha"]:
        return False
    if m.get("transform_settings")!=ctx["transform"]:
        return False
    if m.get("original_numeric_spectrum_sha256")!=arrays_sha(ctx["v"],ctx["y"]):
        return False
    # Independently *reexecute* native COSMOS rebin and preservation comparison.
    f,z,err=edge_rebin(ctx["e"],ctx["y"],ctx["transform"],rebin)
    if m.get("native_rebin_arrays_sha256")!=arrays_sha(f,z):
        return False
    bg=ctx["bg"]
    comp=native.compare(ctx["e"],ctx["y"],f,z,bg,ctx["line"])
    if m.get("native_preservation_result_sha256")!=jsonsha(comp):
        return False
    if comp!=ctx["compare"] or not np.array_equal(f,ctx["f"]) or not np.array_equal(z,ctx["z"]):
        return False
    g=m.get("forward_lineage_graph")
    if not isinstance(g,dict) or g.get("terminal_kind")!="OBSERVED":return False
    if not g.get("reverse_proposals",{}).get("complete_certified"):
        return False
    fw=g.get("forward_edges",{})
    props=g.get("reverse_proposals",{})
    if props.get("C_B")!=transpose(fw.get("B_C",{})) or props.get("B_A")!=transpose(fw.get("A_B",{})):
        return False
    anchor=manifest(ctx)["forward_lineage_graph"]
    # Original physical computation lineage must not be replaced or expanded.
    return fw==anchor["forward_edges"] and g.get("terminal")==anchor["terminal"]

def cosmos_native_unmodified(m,ctx,rebin):
    # Existing native R37F3A calculations are executed before audit; this
    # ablation does not verify *new overlay metadata* and accepts it.
    return True,0

def tachyon_mirror_only(m,ctx,rebin):
    graph=mirror_case(m["forward_lineage_graph"])
    d=mirrored_exact(graph,budget=BUDGET)
    return d["prediction"]!=ABSTAIN,d["operations"]

def tachyon_mirror_replay(m,ctx,rebin):
    graph=mirror_case(m["forward_lineage_graph"])
    d=mirrored_replay(graph,budget=BUDGET)
    if d["prediction"]==ABSTAIN:return False,d["operations"]
    return bool(trusted_provenance_valid(m,ctx,rebin)),d["operations"]

def strong_direct_replay(m,ctx,rebin):
    return bool(trusted_provenance_valid(m,ctx,rebin)),2

def always_abstain(m,ctx,rebin):
    return False,0

ARMS={
 "A_COSMOS_NATIVE_NO_OVERLAY_AUDIT":cosmos_native_unmodified,
 "B_TACHYON_GRAPH_ONLY":tachyon_mirror_only,
 "C_TACHYON_MIRROR_PLUS_FORWARD_REPLAY":tachyon_mirror_replay,
 "D_STRONG_DIRECT_FULL_REPLAY":strong_direct_replay,
 "E_ALWAYS_ABSTAIN":always_abstain,
}

def load_real_chains():
    rebin,sha_native=native.pinned_rebin()
    assert sha_native==Q["source"]["native_rebin_source_expected_sha"]
    controls=Q["observational_sample"]["physical_sources"]
    chains=[];profiles=[]
    for index,(source,receipt,v,y,e,step,npz_sha) in enumerate(native.fetch()):
        control=controls[index]
        assert source["source"]==control["source"] and source["day"]==control["mjd_day"]
        assert np.isfinite(v).all() and np.isfinite(y).all() and y.min()>=0
        assert len(v)==len(y) and len(v)>100
        assert receipt["status"]=="CONDITIONAL_NATURAL_LINE_SUPPORT"
        assert receipt["profile"]==control["native_profile"]
        assert npz_sha and len(npz_sha)==64
        original_profile_sha=arrays_sha(v,y)
        v.setflags(write=False);y.setflags(write=False);e.setflags(write=False)
        line=list(source["line"])
        lo,hi=line
        bg=float(np.median(y[(v<lo-10)|(v>hi+10)]))
        assert math.isfinite(bg) and bg>0
        source_summary={"source":source["source"],"mjd_day":source["day"],
                        "profile_npz_sha256":npz_sha,
                        "numeric_arrays_sha256":original_profile_sha,
                        "native_profile_bins":int(len(v)),
                        "input_resolution_kind":"R37F2D3 averaged natural line profile (not raw filterbank)",
                        "physical_spatial_OFF_verified":False}
        profiles.append(source_summary)
        for tix,transform in enumerate(native.S["transforms"]):
            f,z,err=edge_rebin(e,y,transform,rebin)
            comp=native.compare(e,y,f,z,bg,line)
            c={
                "tag":f"{source['source']}:{source['day']}:T{tix}",
                "source":source["source"],"mjd_day":source["day"],
                "profile_name":source["profile"],"npz_sha":npz_sha,
                "v":v,"y":y,"e":e,"f":f,"z":z,
                "bg":bg,"line":line,"transform":dict(transform),
                "compare":comp,"whole_support_conservation_rel":err,
            }
            chains.append(c)
    if len(profiles)!=2 or len(chains)!=16:
        raise AssertionError("Expected exactly two genuine archived real profiles ×8 transforms")
    return rebin,profiles,chains

def run():
    start=time.perf_counter()
    rebin,profiles,chains=load_real_chains()
    load_s=time.perf_counter()-start
    scores={name:Counter() for name in ARMS}
    elapsed={name:0. for name in ARMS}
    observed_canary=Counter()
    native_chain_summaries=[]
    merged_hash=hashlib.sha256()
    originals={x["source"]:x["npz_sha"] for x in chains}
    for ctx in chains:
        native_chain_summaries.append({
            "source":ctx["source"],"mjd_day":ctx["mjd_day"],
            "transform":ctx["transform"],
            "preservation_status":ctx["compare"]["status"],
            "whole_support_conservation_rel":ctx["whole_support_conservation_rel"],
            "spectrum_numeric_sha256":arrays_sha(ctx["v"],ctx["y"]),
            "rebin_arrays_sha256":arrays_sha(ctx["f"],ctx["z"]),
            "review_checks":ctx["compare"].get("failed",[]),
        })
        for kind in OVERLAYS:
            obs=corrupt(manifest(ctx),kind)
            # Category withheld from every audit function. It is only used for
            # evaluation *after* all algorithms have issued their verdicts.
            assert "attack_type" not in obs and "kind" not in obs
            merged_hash.update(json.dumps(obs,sort_keys=True,separators=(",",":")).encode())
            decisions={}
            for name,func in ARMS.items():
                start=time.perf_counter()
                passed,operations=func(obs,ctx,rebin)
                elapsed[name]+=time.perf_counter()-start
                decisions[name]=(passed,operations)
            is_intact=kind=="INTACT"
            observed_canary[kind]+=1
            for name,(passed,operations) in decisions.items():
                m=scores[name]
                m["total"]+=1
                m["intact_accepted"]+=int(passed and is_intact)
                m["intact_rejected"]+=int(not passed and is_intact)
                m["altered_accepted"]+=int(passed and not is_intact)
                m["altered_rejected"]+=int(not passed and not is_intact)
                m["operations"]+=operations
    assert {x["source"]:x["npz_sha"] for x in chains}==originals
    arms={}
    for name,m in scores.items():
        arms[name]={
            "total":m["total"],"intact_accepted":m["intact_accepted"],
            "intact_rejected":m["intact_rejected"],
            "altered_accepted":m["altered_accepted"],
            "altered_rejected":m["altered_rejected"],
            "postprocessing_seconds_exploratory":round(elapsed[name],6),
            "mean_graph_operations_per_trial":round(m["operations"]/m["total"],4),
        }
    c=arms["C_TACHYON_MIRROR_PLUS_FORWARD_REPLAY"]
    d=arms["D_STRONG_DIRECT_FULL_REPLAY"]
    return {
      "schema":"janus.cosmos.sit0.hawking.real_radio_profile_mirror_audit_result.v1",
      "experiment_name":"JANUS SIT0 — Hawking Dedication",
      "status":"REAL_DERIVED_RADIO_PROFILES_REPLAYED__METADATA_INTEGRITY_PILOT__SCIENCE_HOLD",
      "physical_sources_count":len(profiles),
      "real_spectrum_profiles":profiles,
      "real_native_spectral_rebin_chains":len(chains),
      "metadata_overlay_trials":len(chains)*len(OVERLAYS),
      "overlay_counts":dict(observed_canary),
      "original_derived_science_chain_summaries":native_chain_summaries,
      "source_artifact_verified_sha256":Q["source"]["parent_artifact_archive_sha256"],
      "source_artifact_id":Q["source"]["parent_artifact_id"],
      "source_model_sha1":Q["source"]["native_rebin_source_expected_sha"],
      "tachyon_stock_mirror_pin":Q["source"]["tachyon_pin"],
      "manifest_stream_sha256":merged_hash.hexdigest(),
      "native_archive_and_profile_load_seconds":round(load_s,4),
      "arms":arms,
      "C_mirror_plus_replay_passed_all_16_intact":c["intact_accepted"]==16 and c["intact_rejected"]==0,
      "C_mirror_plus_replay_rejected_all_112_metadata_mutations":c["altered_accepted"]==0 and c["altered_rejected"]==112,
      "equal_to_strong_direct_exact_check":all(c[k]==d[k] for k in ("intact_accepted","intact_rejected","altered_accepted","altered_rejected")),
      "proved_gain_over_strong_exact_baseline":False,
      "original_telemetry_arrays_modified":False,
      "original_native_spectral_decisions_modified":False,
      "full_COSMOS_R28_or_R42_detector_executed":False,
      "original_native_2p794Hz_filterbank_reloaded":False,
      "real_2D_spectrogram_tested":False,
      "spatial_OFF_measured":False,
      "real_technosignature_detected":False,
      "physical_retrocausality_or_backward_time_proven":False,
      "independent_external_holdout":False,
      "production_model_weights_changed":False,
      "dedication":"To Stephen Hawking; research tribute, not an endorsement.",
      "scientific_gate":"HOLD_REAL_DERIVED_PROFILE_NOT_FULL_RADIO_DETECTOR_OR_INDEPENDENT_VALIDATION",
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--audit",action="store_true")
    args=parser.parse_args()
    if not args.audit:parser.error("--audit required")
    OUT.mkdir(parents=True,exist_ok=True)
    result=run()
    (OUT/"result.json").write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("JANUS_SIT0_RESULT",json.dumps({
        "status":result["status"],
        "sources":result["physical_sources_count"],
        "native_chains":result["real_native_spectral_rebin_chains"],
        "metadata_trials":result["metadata_overlay_trials"],
        "arms":result["arms"],
        "strong_baseline_equal":result["equal_to_strong_direct_exact_check"],
        "original_native_spectral_decisions_modified":result["original_native_spectral_decisions_modified"]
    },sort_keys=True),flush=True)
    assert result["C_mirror_plus_replay_passed_all_16_intact"]
    assert result["C_mirror_plus_replay_rejected_all_112_metadata_mutations"]
    assert result["equal_to_strong_direct_exact_check"]
    assert not result["original_native_spectral_decisions_modified"]
    print("JANUS_SIT0_REAL_DERIVED_PROFILE_AUDIT_PASS__PHYSICS_AND_DETECTOR_SCIENCE_HOLD",flush=True)

if __name__=="__main__":
    main()
