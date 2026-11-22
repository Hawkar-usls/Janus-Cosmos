"""JANUS SIT2D: exact archived R41S1H1 physics morphology donor on real GBT.

No data substitution: retrieve 6 X-development and 6 new blc22-evaluation
public rawspec HDF5 files; verify X SHA256 against frozen SIT1, compute
and record evaluation SHA256 from actual bytes, preserve the original data.

Artificial injections are applied only to temporary derived panels, into
both ON scans of a 4-role test unit. All labels used after inference only.
Compare the exact byte-for-byte R41S1H1 physical profile and its original
fixed-Gaussian control to strong same-input direct 4-role narrowband contrast.
Pinned R5F TACHYON replays computational provenance, not lost observations.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path
import numpy as np

import janus_sit1_open_3iatlas as sit1
import janus_sit2b_temporal_recovery as sit2b
import r41s1h1_physical_profile as donor
from tools.tachyon_r5f_stock_mirror_audit import mirror_case, mirrored_replay
from tools.tachyon_r5f_causal_memory_audit import ABSTAIN

ROOT=Path(__file__).resolve().parents[2]
Q=json.loads((ROOT/"research/seti/janus_sit2d_archive_donors_protocol.json").read_text())
OUT=ROOT/"runtime/sit2d"
T=np.arange(16,dtype=np.int32)
VALID=np.all((sit2b.PATHS>=24)&(sit2b.PATHS<=103),axis=1)
BANK=sit2b.PATHS[VALID]
KIND=tuple(Q["study"]["classes"])
SCORE_NAMES=("A_STRONG_DIRECT_FOUR_ROLE","B_R41_H1_FIXED_ORIGINAL",
             "C_R41_H1_PHYSICAL_ORIGINAL","D_TACHYON_R41_H1_PHYSICAL")
SEED=int(Q["study"]["seed"])
assert BANK.shape[1]==16 and len(BANK)>0

def sha(x):
    return hashlib.sha256(x).hexdigest()

def acquire(files,group,old_sha,folder):
    evidence=[];panels=[];errors=[]
    for spec in files:
        scan={"number":spec["scan"],"role":spec["role"],"path":spec["filename"]}
        item=sit1.retrieve(scan,folder)
        item["split"]=group
        evidence.append(item)
        if not item.get("ok"):
            errors.append({"scan":scan["number"],"error":"DOWNLOAD","attempts":item.get("retry_errors")})
            continue
        if group=="X_DEV" and item["first_seen_sha256"]!=old_sha.get(scan["number"]):
            errors.append({"scan":scan["number"],"error":"FROZEN_ORIGINAL_X_SHA256_MISMATCH"})
            (folder/scan["path"]).unlink(missing_ok=True)
            continue
        try:
            q=sit1.reduce_real_h5(folder/scan["path"],item)
            item["header"]=q["header"];item["shape"]=q["dataset_shape"]
            panels.extend(q["panels"])
        except Exception as exc:
            errors.append({"scan":scan["number"],"error":"REAL_HDF5_AXES_OR_FORMAT",
                           "message":str(exc)[:280]})
        finally:
            (folder/scan["path"]).unlink(missing_ok=True)
    if len(panels)!=24:
        errors.append({"group":group,"error":"EXPECTED_6_SCANS_X_4_REAL_PANELS","got":len(panels)})
    if not errors:
        # Mandatory same physical frequency bin for all four roles.
        for window in range(4):
            rows=[x for x in panels if x["panel_idx"]==window]
            if len(rows)!=6 or not all(
                abs(x["native_freq_start_mhz"]-rows[0]["native_freq_start_mhz"])<1e-5 and
                abs(x["native_freq_end_mhz"]-rows[0]["native_freq_end_mhz"])<1e-5
                for x in rows):
                errors.append({"group":group,"window":window,"error":"FREQUENCY_GRID_MISMATCH"})
    return evidence,panels,errors

def role_groups(panels):
    byscan={}
    for x in panels:
        byscan[(x["scan"],x["panel_idx"])]=x
    order=list(dict.fromkeys(x["scan"] for x in panels))
    if len(order)!=6:raise ValueError("SIX_SCANS_REQUIRED")
    out=[]
    # Three overlapping pairings from the one six-scan ON/OFF cadence.
    for i,j in itertools.combinations((0,1,2),2):
        for win in range(4):
            spec=[byscan[(order[i*2],win)],byscan[(order[i*2+1],win)],
                  byscan[(order[j*2],win)],byscan[(order[j*2+1],win)]]
            if [x["role"] for x in spec]!=["ON","OFF","ON","OFF"]:
                raise ValueError("WRONG_4_ROLE_ORDER")
            out.append({"scans":[x["scan"] for x in spec],"window":win,
                        "sources":spec,
                        "four_original_sha256":[x["source_file_sha256"] for x in spec],
                        "matrix":[np.array(x["panel"],dtype=np.float32,copy=True) for x in spec]})
    assert len(out)==12
    return out

def overlay(mats,kind,rng):
    v=[np.array(z,dtype=np.float32,copy=True) for z in mats]
    truepath=None
    if kind==KIND[0]:return v,truepath
    if kind==KIND[1]:
        gap=int(rng.integers(6,13))
        off=int(rng.integers(2,gap+1))
        selected=rng.choice(T,size=3,replace=False)
        for i in range(4):
            _,mad=sit2b.robust_normalize(v[i])
            v[i][:,off::gap]+=(float(rng.uniform(1,2.5))*mad[:,None]).astype(np.float32)
            for t in selected:
                v[i][t,:]+=(float(rng.uniform(1,2))*mad[t])
        return v,truepath
    if kind not in KIND[2:]:
        raise ValueError("NOT_PREREGISTERED")
    first=int(rng.integers(34,97))
    slope=float((-1,-.5,0,.5,1)[int(rng.integers(5))])
    truepath=np.rint(first+slope*(T-7.5)).astype(np.int32)
    amps=(1.2,2.) if kind==KIND[2] else (8.,11.)
    # Both independent ON input records have the same synthetic emitter path.
    # OFFs are untouched. The injection is NOT a real telescope emission.
    for ix in (0,2):
        _,mad=sit2b.robust_normalize(v[ix])
        a=rng.uniform(*amps,size=16)*mad
        for t,k in enumerate(truepath):
            v[ix][t,k]+=float(a[t])
            v[ix][t,k-1]+=.2*float(a[t])
            v[ix][t,k+1]+=.2*float(a[t])
    return v,truepath

def select_blind_path(mats):
    # Same selected path for all methods; known synthetic path not visible.
    mean_on=(mats[0].astype(np.float64)+mats[2].astype(np.float64))/2
    norm,_=sit2b.robust_normalize(mean_on)
    z=norm[T[None,:],BANK]
    side=z-.5*(norm[T[None,:],BANK-3]+norm[T[None,:],BANK+3])
    i=int(np.argmax(np.round(side.mean(axis=1),10)))
    return BANK[i].copy(),float(side.mean(axis=1)[i])

def donor_cutout(mats,path):
    if np.any(path<24) or np.any(path>103):raise ValueError("PATH_OUTSIDE_REAL_H1_SUPPORT")
    return np.stack([
        np.stack([v[t,path[t]-24:path[t]+25] for t in T],axis=0)
        for v in mats
    ],axis=0)

def direct_four_role(mats,path):
    def q(v):
        norm,_=sit2b.robust_normalize(v)
        return float((norm[T,path]-.5*(norm[T,path-3]+norm[T,path+3])).mean())
    on1,off1,on2,off2=[q(z) for z in mats]
    return min(on1,on2)-max(0.,off1,off2)

def donor_evidence(mats,path,template,df):
    cube=donor_cutout(mats,path)
    return donor.score(cube,*template,df)

def tachyon_evidence(source_sha,selected_path,score,original_quad_sha):
    # Root fingerprint includes all four source observations. Graph is data
    # provenance, NOT a physical cause for the celestial emission.
    combo=sha(json.dumps(original_quad_sha,sort_keys=True).encode())
    root="four_real_source_scans:"+combo[:18]
    intermediate="blind_candidates:"+sit1.hash_array(selected_path,dtype="<i4")[:18]
    result="four_role_physical_score:"+sha(np.float64(score).tobytes())[:18]
    graph={
      "schema":"janus.tachyon.r5f.synthetic_reverse_road_graph.v1",
      "case_ref":combo[:16],"terminal":result,"terminal_kind":"OBSERVED",
      "forward_edges":{"A_B":{root:[intermediate]},"B_C":{intermediate:[result]}},
      "reverse_proposals":{"C_B":{result:[intermediate]},
                           "B_A":{intermediate:[root]},"complete_certified":True},
      "witness_reports":[],"witness_window_seconds":120,"synthetic_only":True,
    }
    mirrored=mirrored_replay(mirror_case(graph),budget=500)
    if mirrored["prediction"]==ABSTAIN:raise AssertionError("R5F_CANNOT_REPLAY_CORRECT_FOUR_ROLE_GRAPH")
    return mirrored["operations"]

def create_cases(groups,rng):
    cases=[]
    for group in groups:
        for rep in range(3):
            for label in KIND:
                v,truth=overlay(group["matrix"],label,rng)
                cases.append({"group":group,"rep":rep,"role_mats":v,
                              "kind_for_evaluation_only":label,"truth_for_evaluation_only":truth})
    assert len(cases)==144
    return cases

def score_cases(cases,group,template,df):
    out=[];costs={name:0. for name in SCORE_NAMES};graphops=0
    for item in cases:
        v=item["role_mats"]
        path,selector_score=select_blind_path(v)
        start=time.perf_counter()
        direct_score=direct_four_role(v,path)
        costs[SCORE_NAMES[0]]+=time.perf_counter()-start
        start=time.perf_counter()
        d=donor_evidence(v,path,template,df)
        costs[SCORE_NAMES[1]]+=time.perf_counter()-start
        costs[SCORE_NAMES[2]]+=0.  # donor computes both physical+fixed together
        start=time.perf_counter()
        physical2=donor_evidence(v,path,template,df)["physical"]
        graphops+=tachyon_evidence("",path,physical2,item["group"]["four_original_sha256"])
        costs[SCORE_NAMES[3]]+=time.perf_counter()-start
        if abs(physical2-d["physical"])>1e-10:
            raise AssertionError("SOURCE_DONOR_FORWARD_REPLAY_NOT_DETERMINISTIC")
        truth=item["truth_for_evaluation_only"]
        located=truth is not None and int(np.sum(np.abs(path-truth)<=2))>=12
        out.append({
            "split":group,
            "kind":item["kind_for_evaluation_only"],
            "localized":bool(located),
            "scores":{
                SCORE_NAMES[0]:direct_score,
                SCORE_NAMES[1]:float(d["fixed"]),
                SCORE_NAMES[2]:float(d["physical"]),
                SCORE_NAMES[3]:float(physical2),
            },
            "source_sha256_combination":sha(json.dumps(item["group"]["four_original_sha256"]).encode()),
        })
    return out,costs,graphops

def frozen_eval(dev,exam):
    out={};ids={}
    for name in SCORE_NAMES:
        neg=[x["scores"][name] for x in dev if x["kind"] in KIND[:2]]
        if len(neg)!=72:raise AssertionError("NOT_72_DEVELOPMENT_CONTROLS")
        cutoff=float(np.quantile(np.asarray(neg),.95,method="higher"))
        weak=[x for x in exam if x["kind"]==KIND[2]]
        strong=[x for x in exam if x["kind"]==KIND[3]]
        null=[x for x in exam if x["kind"] in KIND[:2]]
        recovered={i for i,x in enumerate(exam) if x["kind"]==KIND[2] and x["localized"] and x["scores"][name]>cutoff}
        strong_count=sum(x["localized"] and x["scores"][name]>cutoff for x in strong)
        false=sum(x["scores"][name]>cutoff for x in null)
        out[name]={
            "development_negative_p95":cutoff,
            "heldout_weak_positive":len(weak),"heldout_strong_positive":len(strong),
            "heldout_artificial_controls":len(null),
            "heldout_weak_path_localized":sum(x["localized"] for x in weak),
            "heldout_weak_recovered":len(recovered),
            "heldout_strong_recovered":int(strong_count),
            "heldout_control_exceedances":int(false),
            "strong_and_control_qualified":strong_count>=32 and false<=3,
        }
        ids[name]=recovered
    d=out[SCORE_NAMES[2]]
    direct=out[SCORE_NAMES[0]]
    fixed=out[SCORE_NAMES[1]]
    reverse=out[SCORE_NAMES[3]]
    if d!=reverse or ids[SCORE_NAMES[2]]!=ids[SCORE_NAMES[3]]:
        raise AssertionError("R5F_REPLAY_GAIN_FABRICATED")
    gain={
      "donor_extra_weak_over_strong_direct_4_role":len(ids[SCORE_NAMES[2]]-ids[SCORE_NAMES[0]]),
      "donor_extra_weak_over_original_fixed_Gaussian":len(ids[SCORE_NAMES[2]]-ids[SCORE_NAMES[1]]),
      "tachyon_extra_weak_over_identical_donor_direct":0,
      "donor_strong_control_gate":d["strong_and_control_qualified"],
      "strict_donor_weak_gain_gate":bool(
          d["strong_and_control_qualified"] and
          len(ids[SCORE_NAMES[2]]-ids[SCORE_NAMES[0]])>=4 and
          len(ids[SCORE_NAMES[2]]-ids[SCORE_NAMES[1]])>=4 and
          d["heldout_control_exceedances"]<=min(
             direct["heldout_control_exceedances"],fixed["heldout_control_exceedances"])),
      "astrophysical_independent_gain_proved":False
    }
    return out,gain

def run():
    OUT.mkdir(parents=True,exist_ok=True)
    folder=OUT/"source_hdf5";folder.mkdir(exist_ok=True)
    sit1report=json.loads((ROOT/"research/seti/janus_sit1_open_result.json").read_text())
    prior={r["scan"]:r["first_seen_sha256"] for r in sit1report["runtime_inputs"]["hdf5_files"]}
    x,xb,xerr=acquire(Q["real_radio"]["dev_files"],"X_DEV",prior,folder)
    e,eb,eerr=acquire(Q["real_radio"]["evaluation_files"],"NEW_BLC22",prior,folder)
    errors=xerr+eerr
    report={
       "schema":"janus.cosmos.sit2d.archived_h1_morphology_real_gbt.v1",
       "status":"DATA_SOURCE_OR_FORMAT_HOLD" if errors else "REAL_GBT_BL_C22_DONOR_EVALUATION_ENGINEERING_COMPLETE__SCIENCE_HOLD",
       "source_name":"GBT 3I/ATLAS public rawspec.0002 HDF5; X dev vs new blc22 evaluation one telescope/day",
       "source_receipts":[{k:z.get(k) for k in ("split","scan","role","filename","first_seen_sha256","bytes","header","shape","error","retry_errors")} for z in x+e],
       "errors":errors,
       "x_development_real_panels":len(xb),
       "new_blc22_real_panels":len(eb),
       "donor_exact_source_blob_sha1":Q["archive_audit"]["exact_donor_blob_sha1"],
       "real_sky_truth":False,
       "synthetic_labels_added_to_real_backgrounds":True,
       "external_independent_pointing":False,
       "source_date_independent":False,
       "production_modified":False,
       "astrophysical_sensitivity_gain":None,
       "physical_retrocausality":False,
       "scientific_gate":"HOLD_EXTERNAL_NEW_DAY_HIGH_RES_AND_AUTHENTICATED_RFI_SOURCE_CLASSES"
    }
    if errors:return report
    xg,eg=role_groups(xb),role_groups(eb)
    df_x=abs(float(x[0]["header"]["foff"]))*1e6
    df_e=abs(float(e[0]["header"]["foff"]))*1e6
    # Approximate real temporal spacing after 279 observed time bins -> 16 means.
    t_x=float(x[0]["header"]["tsamp"])*int(x[0]["shape"][0])/16
    t_e=float(e[0]["header"]["tsamp"])*int(e[0]["shape"][0])/16
    t0=time.perf_counter()
    templ_x=donor.templates(df_x,t_x)
    templ_e=donor.templates(df_e,t_e)
    rng=np.random.default_rng(SEED)
    dr,dc,dops=score_cases(create_cases(xg,rng),"X_DEV",templ_x,df_x)
    er,ec,eops=score_cases(create_cases(eg,rng),"NEW_BLC22",templ_e,df_e)
    metrics,gain=frozen_eval(dr,er)
    report.update({
        "dataset":{"dev_cases":len(dr),"eval_cases":len(er),
                   "eval_weak":36,"eval_strong":36,"eval_controls":72,
                   "real_four_role_base_groups_per_receiver":len(xg),
                   "same_date_and_target_across_receivers":True,
                   "synthetic_injection_truth_not_known_physical_sky":True,
                   "all_cases_share_dependent_original_scans":True},
        "receiver_physical_header":{"x_fch1":x[0]["header"]["fch1"],
                                    "blc22_fch1":e[0]["header"]["fch1"],
                                    "x_foff":x[0]["header"]["foff"],
                                    "blc22_foff":e[0]["header"]["foff"],
                                    "x_df_hz":df_x,"blc22_df_hz":df_e,
                                    "x_aggregated_dt_seconds":t_x,
                                    "blc22_aggregated_dt_seconds":t_e},
        "donor_exact_number_of_physical_templates":len(templ_x[0]),
        "matched_four_role_methods":metrics,
        "comparisons":gain,
        "experimental_wall_seconds_excluding_download":round(time.perf_counter()-t0,4),
        "algorithm_seconds_development":{k:round(v,5) for k,v in dc.items()},
        "algorithm_seconds_evaluation":{k:round(v,5) for k,v in ec.items()},
        "tachyon_mirror_graph_operations_total":dops+eops,
        "real_gbt_sha256_dev_matches_prior_sit1":True,
        "new_evaluation_sha256_provider_published_reference":False
    })
    report["status"]=(
        "DONOR_H1_REAL_GBT_ENGINEERING_PASS_POSITIVE_GUARD_PASS__BUT_EXTERNAL_SCIENCE_HOLD"
        if gain["strict_donor_weak_gain_gate"]
        else "DONOR_H1_REAL_GBT_ENGINEERING_PASS__NO_QUALIFIED_WEAK_GAIN__SCIENCE_HOLD"
    )
    return report

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--audit",action="store_true")
    args=ap.parse_args()
    if not args.audit:ap.error("--audit required")
    result=run()
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"result.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("SIT2D_DONOR_RESULT",json.dumps({
        "status":result["status"],"errors":result["errors"],
        "outcomes":result.get("matched_four_role_methods"),
        "gains":result.get("comparisons")},sort_keys=True),flush=True)
    if result.get("matched_four_role_methods"):
        assert result["dataset"]["eval_cases"]==144
        assert result["comparisons"]["tachyon_extra_weak_over_identical_donor_direct"]==0
        print("SIT2D_ARCHIVE_EXACT_R41_H1_REAL_GBT_SOURCE_TEST_DONE__SCIENTIFIC_HOLD")
    else:
        print("SIT2D_ARCHIVE_OR_INPUT_HOLD_NO_SCIENTIFIC_RESULT")

if __name__=="__main__":
    main()
