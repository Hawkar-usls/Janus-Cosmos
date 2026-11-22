"""JANUS SIT2E-E. Preregistered research-only Source/Reference soft contrast.

Inputs: four archived, *full SHA256 verified*, July and September Parkes files.
No October bytes, no production changes, no free parameter search.
All science labels are synthetic injections on real radio backgrounds.
September p95 is retrospective in-sample only; July was previously inspected.
A real physical detection, a verified FPR or improved JANUS sensitivity is NOT claimed.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import requests

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"research"/"seti"))
import janus_sit2eb_power_calibration as b
import janus_sit2ec_shift_diagnostics as c
import janus_sit2eb_full_integrity as integrity
import janus_sit2ed_polarization_806_witness as d

OUT=ROOT/"runtime/sit2ee"
DATES=("2025-07-31","2025-09-12")
METHODS=("JANUS_NATIVE_RIDGE","DIRECT_672_MATCHED","DIRECT_672_SIDEBAND",
         "DIRECT_ON_OFF_CONTRAST")
CATEGORIES=("NO_ADDED_SIGNAL","SYNTHETIC_ON_ONLY_RFI","SYNTHETIC_SHARED_ON_OFF_RFI",
            "WEAK_ON_ONLY","STRONG_ON_ONLY","WEAK_ON_REF_COLLISION")
NEGATIVE=set(CATEGORIES[:3])
SEED=22061123
REPLICATES=16
OFF_COEFFICIENT=1.0  # fixed unit-coefficient ON minus OFF morphology; not fitted
LOCK=json.loads((ROOT/"research/seti/janus_sit2eb_july_september_calibration_receipt.json").read_text())
assert set(LOCK["numeric_thresholds"])==set(b.METHODS)
T=b.TIX
PATHS=b.PATHS

def reference_profile(panel):
    z,_=b.robust_normalize(panel)
    # Fixed OFF frequency profile only; OFF is not fitted to ON pixels or labels.
    return np.median(z,axis=0)

def sideband_bank(z):
    return (z[T,PATHS]-.5*(z[T,PATHS-3]+z[T,PATHS+3])).mean(axis=1)

def off_contrast(on_panel,off_panel):
    z,_=b.robust_normalize(on_panel)
    r=reference_profile(off_panel)
    on=sideband_bank(z)
    off=(r[PATHS]-.5*(r[PATHS-3]+r[PATHS+3])).mean(axis=1)
    bank=on-OFF_COEFFICIENT*off
    index=int(np.argmax(np.round(bank,10)))
    return PATHS[index].copy(),float(bank[index]),{
        "max_reference_matched_occupancy":float(off[index]),
        "raw_selected_on_score":float(on[index])}

def all_scores(on,ref):
    frozen=b.score(on)
    path,score,audit=off_contrast(on,ref)
    return {**frozen,"DIRECT_ON_OFF_CONTRAST":(path,score)},audit

def overlay_path(base,path,rng,strength):
    out=np.asarray(base,dtype=np.float32).copy()
    _,mad=b.robust_normalize(out)
    amps=rng.uniform(*strength,size=16)*mad
    for t,f in enumerate(path):
        ix=int(f)
        out[t,ix]+=amps[t]
        out[t,ix-1]+=.2*amps[t]
        out[t,ix+1]+=.2*amps[t]
    return out

def case_overlay(base_on,base_ref,kind,rng):
    on=np.asarray(base_on,dtype=np.float32).copy()
    ref=np.asarray(base_ref,dtype=np.float32).copy()
    truth=None
    if kind=="NO_ADDED_SIGNAL":return on,ref,truth
    if kind=="SYNTHETIC_ON_ONLY_RFI":
        on,_=b.create_overlay(on,"SYNTHETIC_RFI",rng)
    elif kind=="SYNTHETIC_SHARED_ON_OFF_RFI":
        # Same narrowband frequency path with independent source/ref power units.
        path=PATHS[int(rng.integers(0,len(PATHS)))]
        a=float(rng.uniform(5,8))
        on=overlay_path(on,path,rng,(a,a+.01))
        ref=overlay_path(ref,path,rng,(a,a+.01))
    elif kind in ("WEAK_ON_ONLY","STRONG_ON_ONLY"):
        amp="WEAK_LINE" if kind=="WEAK_ON_ONLY" else "STRONG_LINE"
        on,truth=b.create_overlay(on,amp,rng)
    elif kind=="WEAK_ON_REF_COLLISION":
        # Evaluation-only truth: actual strongest OFF sideband candidate.
        # Fixed criterion, independent of changed ON and its injected outcome.
        truth=b.direct_matched(ref,sideband=True)[0]
        on=overlay_path(on,truth,rng,(1.2,2.0))
    else:raise ValueError("UNREGISTERED_TEST_KIND")
    return on,ref,truth

def calc_cases(pairs):
    rng=np.random.default_rng(SEED)
    out=[]
    for date in DATES:
        for band in range(4):
            on0=pairs[(date,"source")][band]["panel"]
            ref0=pairs[(date,"reference")][band]["panel"]
            assert pairs[(date,"source")][band]["start_bin"]==pairs[(date,"reference")][band]["start_bin"]
            for rep in range(REPLICATES):
                for kind in CATEGORIES:
                    on,off,truth=case_overlay(on0,ref0,kind,rng)
                    scores,audit=all_scores(on,off)
                    out.append({
                        "id":f"{date}:{band}:{rep}:{kind}",
                        "date":date,"band":band,"rep":rep,"kind":kind,
                        "scores":{m:float(v[1]) for m,v in scores.items()},
                        "localized":{m:b.localizes(v[0],truth) for m,v in scores.items()},
                        "ref_selected_occupancy":audit["max_reference_matched_occupancy"]
                    })
    assert len(out)==2*4*REPLICATES*len(CATEGORIES)
    return out

def result_metrics(rows,threshold,method,date):
    group=[r for r in rows if r["date"]==date]
    neg=[r for r in group if r["kind"] in NEGATIVE]
    categories={}
    for kind in CATEGORIES:
        sample=[r for r in group if r["kind"]==kind]
        ispos=kind not in NEGATIVE
        hits=[r for r in sample if r["scores"][method]>threshold
              and (r["localized"][method] if ispos else True)]
        categories[kind]={
            "cases":len(sample),"above_threshold_and_localized_if_positive":len(hits),
            "per_parent_band_hits":[sum(r["id"] in {x["id"] for x in hits}
                                        for r in sample if r["band"]==band) for band in range(4)]
        }
    positives=[r for r in group if r["kind"] not in NEGATIVE]
    return {"controls":len(neg),"control_exceedances":sum(r["scores"][method]>threshold for r in neg),
            "parent_band_groups":4,"categories":categories,
            "strong_detected":sum(r["scores"][method]>threshold and r["localized"][method]
                                  for r in positives if r["kind"]=="STRONG_ON_ONLY"),
            "weak_detected":sum(r["scores"][method]>threshold and r["localized"][method]
                                for r in positives if r["kind"]=="WEAK_ON_ONLY"),
            "collision_weak_detected":sum(r["scores"][method]>threshold and r["localized"][method]
                                for r in positives if r["kind"]=="WEAK_ON_REF_COLLISION")}

def evaluate(rows):
    outcome={}
    for method in METHODS:
        sept_neg=[r["scores"][method] for r in rows
                  if r["date"]==DATES[1] and r["kind"] in NEGATIVE]
        assert len(sept_neg)==4*REPLICATES*3
        diagnostic_p95=float(np.quantile(sept_neg,.95,method="higher"))
        out={"retrospective_sept_p95_NOT_FROZEN":diagnostic_p95,
             "threshold_basis":"same September controls, different synthetic scenario vs SIT2E-B",
             "july":result_metrics(rows,diagnostic_p95,method,DATES[0]),
             "september":result_metrics(rows,diagnostic_p95,method,DATES[1])}
        # Preserve an old-threshold stress test separately, without adopting
        # it as a calibrated operating point for the changed OFF-aware scorer.
        if method in b.METHODS:
            t=LOCK["numeric_thresholds"][method]
            out["original_SIT2EB_threshold"]=t
            out["july_at_original_SIT2EB_threshold"]=result_metrics(rows,t,method,DATES[0])
            out["sept_at_original_SIT2EB_threshold"]=result_metrics(rows,t,method,DATES[1])
        outcome[method]=out
    direct=outcome["DIRECT_672_SIDEBAND"]
    off=outcome["DIRECT_ON_OFF_CONTRAST"]
    jdir=direct["july"];joff=off["july"]
    # Different empirical in-sample thresholds: use as exploratory comparisons,
    # not scientific gain. Include collisions and strong guard to catch losses.
    guard={
       "july_direct_sideband_shared_rfi":jdir["categories"]["SYNTHETIC_SHARED_ON_OFF_RFI"]["above_threshold_and_localized_if_positive"],
       "july_contrast_shared_rfi":joff["categories"]["SYNTHETIC_SHARED_ON_OFF_RFI"]["above_threshold_and_localized_if_positive"],
       "july_direct_sideband_strong":jdir["strong_detected"],
       "july_contrast_strong":joff["strong_detected"],
       "july_direct_sideband_weak":jdir["weak_detected"],
       "july_contrast_weak":joff["weak_detected"],
       "july_direct_sideband_collision_weak":jdir["collision_weak_detected"],
       "july_contrast_collision_weak":joff["collision_weak_detected"],
       "strong_guard_no_losses":joff["strong_detected"]>=jdir["strong_detected"],
       "weak_guard_no_losses":joff["weak_detected"]>=jdir["weak_detected"],
       "collision_guard_no_losses":joff["collision_weak_detected"]>=jdir["collision_weak_detected"],
       "july_5pct_controls_guard":joff["control_exceedances"]<=int(.05*joff["controls"]),
       "independent_astronomy_validation":False,
       "never_auto_adopt":True
    }
    return outcome,guard

def exact_806150_channel_audit(cubes):
    """Probe actual 806.150 MHz in all four archived files, IF0 and IF1.

    Static 128-bin window 806.000–806.254 MHz extracted from genuine 300-row
    spectra, no signal-dependent frequency selection. Fixed zero-drift path at
    channel 75 = 806.150 MHz. No threshold fit, hard veto, or cosmic claim.
    """
    rows=[]
    fixed_path=np.full(16,75,dtype=np.int32)
    choices=np.where(np.all(PATHS==fixed_path[None,:],axis=1))[0]
    if len(choices)!=1:raise AssertionError("FIXED_806150_PATH_NOT_IN_ORIGINAL_672_BANK")
    idx=int(choices[0])
    for date in DATES:
        for iff in (0,1):
            cube_src=cubes[(date,"source")]
            cube_ref=cubes[(date,"reference")]
            assert cube_src.shape==cube_ref.shape==(300,4,376)
            # d.CSTART=51000 => 806.000 MHz and 806.150 is channel 75.
            s=np.stack([x.mean(axis=0) for x in np.array_split(cube_src[:,iff,:128],16)]).astype(np.float32)
            ref=np.stack([x.mean(axis=0) for x in np.array_split(cube_ref[:,iff,:128],16)]).astype(np.float32)
            zs,_=b.robust_normalize(s)
            zr=reference_profile(ref)
            on=float(sideband_bank(zs)[idx])
            off=float((zr[PATHS]-.5*(zr[PATHS-3]+zr[PATHS+3])).mean(axis=1)[idx])
            contrast=on-off
            best_path,best_val,aud=off_contrast(s,ref)
            rng=np.random.default_rng(SEED+iff+(0 if date==DATES[0] else 100))
            injection_diagnostics={}
            for kind,strength in (("WEAK_ON_COLLISION",(1.2,2.0)),
                                  ("STRONG_ON_COLLISION",(8.0,11.0))):
                deltas=[]
                selected=[]
                for rep in range(32):
                    augmented=overlay_path(s,fixed_path,rng,strength)
                    zaug,_=b.robust_normalize(augmented)
                    newer=float(sideband_bank(zaug)[idx])-off
                    deltas.append(newer-contrast)
                    pa,_,_=off_contrast(augmented,ref)
                    selected.append(bool(b.localizes(pa,fixed_path)))
                injection_diagnostics[kind]={
                   "cases":32,
                   "fixed_path_signal_increment_median":float(np.median(deltas)),
                   "fixed_path_signal_increment_min":float(min(deltas)),
                   "positive_increment_count":int(sum(x>0 for x in deltas)),
                   "chosen_path_localized_count":int(sum(selected)),
                   "ground_truth":"SYNTHETIC_ADDED_ON_ONLY",
                   "pre_existing_common_line_not_an_on_target_detection":True
                }
            rows.append({
                "date":date,"IF_index_not_certified_polarization":iff,
                "frequency_mhz":806.150,
                "window_mhz":[806.0,806.254],
                "source_reference_fixed_path_on_sideband_score":on,
                "reference_fixed_path_sideband_score":off,
                "source_minus_reference_fixed_path_contrast_score":contrast,
                "fixed_path_signed_magnitude_suppression":float(abs(on)-abs(contrast)),
                "best_contrast_score_after_full_672_bank_search":float(best_val),
                "best_contrast_selected_fixed_806150_path":bool(np.array_equal(best_path,fixed_path)),
                "injections":injection_diagnostics,
                "same_300s_source_and_reference_time_cadence":False,
                "original_radio_power_verified_by_full_file_sha256":True
            })
    return rows

def run():
    OUT.mkdir(parents=True,exist_ok=True)
    archive=OUT/"archive";archive.mkdir(exist_ok=True)
    selected=integrity.selected_files(json.loads(integrity.MANIFEST.read_text()))
    if len(selected)!=4 or tuple(dict.fromkeys(s[0] for s in selected))!=DATES:
        raise ValueError("OCTOBER_OR_UNREGISTERED_DATE_SELECTED")
    record={"schema":"janus.cosmos.sit2ee.off_aware_signal_protection.v1",
       "status":"SOURCE_DATA_OR_SCIENCE_HOLD","method":"no fitting; fixed unit ON-OFF spectral sideband contrast",
       "source_dates":list(DATES),"science_gate":"HOLD",
       "previous_sit2eb_lock_sha256":LOCK["computation_sha256"],
       "previous_thresholds_unchanged":LOCK["numeric_thresholds"],
       "new_thresholds_approved":False,"production_changed":False,
       "october_files_requested":False,"october_power_read":False,
       "real_astronomical_ground_truth":False,"TCS_encoder_pointing_verified":False,
       "physical_polarization_IF_mapping_verified":False,
       "methods":list(METHODS),"test_classes":list(CATEGORIES),
       "seed":SEED,"replicates_per_parent":REPLICATES,
       "parent_groups_per_date":4,"synthetic_repeats_non_independent":True,
       "files":[],"panel_hashes":[],"exact_806150_test":[],"methods_outcomes":{},
       "protection_gates":{},"errors":[]}
    pairs={}
    narrow_cubes={}
    try:
        with requests.Session() as sess:
            for date,role,spec in selected:
                p=None
                try:
                    p,receipt=c.download_verified(sess,date,role,spec,archive)
                    panels=b.read_panels(p,date,spec)
                    cube,header=d.power_cube(p,date)
                    narrow_cubes[(date,role)]=cube
                    pairs[(date,role)]=panels
                    record["files"].append(receipt)
                    record["panel_hashes"].extend([{
                        "date":date,"role":role,"band":i,
                        "sha256":x["panel_sha256"],"channel_start":x["start_bin"]
                    } for i,x in enumerate(panels)])
                finally:
                    if p is not None:p.unlink(missing_ok=True)
        record["exact_806150_test"]=exact_806150_channel_audit(narrow_cubes)
        rows=calc_cases(pairs)
        record["dataset_cases"]=len(rows)
        record["methods_outcomes"],record["protection_gates"]=evaluate(rows)
        record["status"]="OFF_CONTRAST_SAME_INFO_PAIRED_EXPERIMENT_COMPUTED__PROSPECTIVE_VALIDATION_HOLD"
    except Exception as exc:
        record["errors"].append({"type":type(exc).__name__,"message":str(exc)[:450]})
    finally:
        for f in archive.glob("*"):f.unlink(missing_ok=True)
        output=OUT/"sit2ee_result.json"
        output.write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+"\n")
        digest=hashlib.sha256(output.read_bytes()).hexdigest()
        print("SIT2EE_RESULT",json.dumps({"status":record["status"],"errors":record["errors"],
            "files_verified":len(record["files"]),"cases":record.get("dataset_cases"),
            "guard":record["protection_gates"],"exact_806150_test":record["exact_806150_test"],"report_sha256":digest,
            "methods_outcomes":record["methods_outcomes"]},sort_keys=True),flush=True)
    if record["errors"]:raise RuntimeError(record["errors"][0]["message"])
    return record

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--audit",required=True,action="store_true")
    parser.parse_args()
    run()
