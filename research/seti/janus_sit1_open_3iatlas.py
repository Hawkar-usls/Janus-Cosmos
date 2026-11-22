"""JANUS SIT1-OPEN: genuine GBT 3I/ATLAS X-band 2025-12-18 ON/OFF pilot.

Research-only:
- download exact public HDF5 rawspec.0002 subband files without replacements
- hash real bytes, record original time/frequency header metadata
- form 16x128 *exploratory OOD* real time-frequency panels from measured arrays
- run actual frozen COSMOS R42A8A blind-line search and R42A7A shadow score
- compare frozen COSMOS vs pinned TACHYON R5F stock mirror+forward replay
  vs strong exact direct provenance; do not change native scores
- do not claim independent actual antenna pointing, SETI sensitivity, ETI
  detection, real false alarm reduction, or physical retrocausal information.

The on/off labels come from the archive, not independent Antenna FITS.
Every source SHA256 is first computed from downloaded bytes; no published
provider checksum is claimed unless separately verified.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import os
import pathlib
import struct
import sys
import time
import traceback
from collections import Counter

import h5py
import numpy as np
import requests

ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"research"/"seti"))
sys.path.insert(0,str(ROOT/"_tachyon"))
import r42a8_ood_stress as cosmos
import r42a7_shadow_rfi as shadow
from tools.tachyon_r5f_stock_mirror_audit import mirror_case, mirrored_replay, transpose
from tools.tachyon_r5f_causal_memory_audit import ABSTAIN

Q=json.loads((ROOT/"research/seti/janus_sit1_open_protocol.json").read_text())
OUT=ROOT/"runtime/sit1_open"
BASE_URL="https://bldata.berkeley.edu/ATLAS/GB_ATLAS/"
MAX_BYTES=90_000_000
GROUPS=(0.20,0.40,0.60,0.80)
FAIL={"schema":"janus.cosmos.sit1_open.real_gbt_result.v1",
      "status":"DATA_ACQUISITION_OR_COMPATIBILITY_HOLD",
      "stage":"DATA_QUALIFICATION",
      "source_name":"GBT 3I/ATLAS X band 2025-12-18 blc20 lowest-resolution HDF5",
      "replayed_real_files":0,
      "no_synthetic_substitution":True,
      "production_models_modified":False,
      "scientific_gain":None,
      "false_alarm_reduction":None,
      "actual_pointing_independently_authenticated":False,
      "independent_test_groups":0,
      "downloads":[],
      "errors":[]}


def hash_bytes(b:bytes)->str:
    return hashlib.sha256(b).hexdigest()


def hash_array(a:np.ndarray, dtype="<f8")->str:
    v=np.ascontiguousarray(a,dtype=dtype)
    h=hashlib.sha256()
    h.update(json.dumps(v.shape).encode())
    h.update(v.tobytes())
    return h.hexdigest()


def retrieve(scan:dict,folder:pathlib.Path)->dict:
    path=folder/scan["path"]
    result={"scan":scan["number"],"role":scan["role"],"filename":scan["path"],
            "url":BASE_URL+scan["path"],"provider_sha256_known":False}
    h=hashlib.sha256()
    size=0
    try:
        with requests.get(result["url"],timeout=(18,100),stream=True,headers={"User-Agent":"JANUS-SIT1-Open-Academic-Provenance/1.0"}) as resp:
            resp.raise_for_status()
            if resp.headers.get("content-length") and int(resp.headers["content-length"])>MAX_BYTES:
                raise ValueError("DOWNLOAD_EXCEEDS_PRESET_MAX_BYTES")
            with path.open("wb") as output:
                for block in resp.iter_content(chunk_size=1024*1024):
                    if not block:continue
                    size+=len(block)
                    if size>MAX_BYTES:
                        raise ValueError("DOWNLOAD_EXCEEDS_PRESET_MAX_BYTES")
                    h.update(block)
                    output.write(block)
        if size<4096:
            raise ValueError("DOWNLOAD_TOO_SMALL_FOR_HDF5")
        result.update({"ok":True,"bytes":size,"first_seen_sha256":h.hexdigest()})
        return result
    except Exception as exc:
        path.unlink(missing_ok=True)
        result.update({"ok":False,"bytes_received_before_error":size,
                       "error_type":type(exc).__name__,"error":str(exc)[:300]})
        return result


def python_scalar(x):
    if isinstance(x,bytes):return x.decode("utf-8",errors="replace")
    if isinstance(x,np.ndarray):
        if x.size==1:return python_scalar(x.ravel()[0])
        return [python_scalar(y) for y in x.ravel()[:8]]
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,(str,int,float,bool)) or x is None:return x
    return str(x)


def source_header(f:h5py.File)->dict:
    a={str(k).lower():python_scalar(v) for k,v in f.attrs.items()}
    if "data" not in f:
        raise ValueError("DATASET_DATA_NOT_FOUND")
    for k,v in f["data"].attrs.items():
        a.setdefault(str(k).lower(),python_scalar(v))
    needed=Q["numerical_adapter"]["require_source_header_fields"]
    missing=[k for k in needed if k not in a]
    if missing:raise ValueError("MISSING_HDF5_PHYSICAL_AXES:"+",".join(missing))
    for k in needed:
        a[k]=float(a[k])
        if not np.isfinite(a[k]):raise ValueError("NONFINITE_HDF5_AXIS:"+k)
    if a["tsamp"]<=0 or a["foff"]==0:
        raise ValueError("INVALID_HDF5_TIME_FREQUENCY_AXIS")
    # Pointing fields are reported, but NEVER equated with measured antenna direction.
    wanted=("source_name","src_raj","src_dej","tstart","tsamp","fch1","foff",
            "nchans","nbits","nifs","telescope_id","machine_id","rawdatafile")
    return {k:a[k] for k in wanted if k in a}


def reduce_real_h5(path:pathlib.Path,record:dict):
    with h5py.File(path,"r") as f:
        d=f["data"]
        if d.ndim!=3 or d.shape[1]<1 or d.shape[2]<128 or d.shape[0]<16:
            raise ValueError("INCOMPATIBLE_HDF5_SHAPE:"+str(d.shape))
        head=source_header(f)
        ntime,nbeam,nfreq=tuple(map(int,d.shape))
        panels=[]
        for fraction in GROUPS:
            start=int((nfreq-128)*fraction)
            # Real unchanged 128-channel contiguous frequency slice.
            values=np.asarray(d[:,0,start:start+128],dtype=np.float64)
            if values.shape!=(ntime,128) or not np.isfinite(values).all():
                raise ValueError("NONFINITE_OR_MISSING_MEASURED_PANEL:"+str(fraction))
            # Pure temporal mean aggregation of real rows; never synthesizes
            # extra time resolution, inpaints NaNs or shifts frequencies.
            parts=np.array_split(values,16,axis=0)
            panel=np.stack([np.mean(p,axis=0) for p in parts]).astype(np.float32)
            if not np.isfinite(panel).all():
                raise ValueError("NONFINITE_TIME_AGGREGATION")
            f_start=float(head["fch1"]+head["foff"]*start)
            f_stop=float(head["fch1"]+head["foff"]*(start+127))
            panels.append({
                "scan":record["scan"],"role":record["role"],"panel_idx":len(panels),
                "frac":fraction,"start_bin":start,"native_bin_count":128,
                "native_freq_start_mhz":f_start,
                "native_freq_end_mhz":f_stop,
                "scan_first_mjd":float(head["tstart"]),
                "source_file_sha256":record["first_seen_sha256"],
                "measured_panel_sha256":hash_array(panel,dtype="<f4"),
                "panel":panel,
            })
        return {"record":record,"header":head,
                "dataset_shape":[ntime,nbeam,nfreq],
                "time_rows_aggregated":ntime,
                "selected_frequency_segments":len(panels),
                "panels":panels}


def score_one(panel:np.ndarray,model):
    mean,scale,coef,intercept=model
    v=panel[None,:,:]
    path,ranking=cosmos.blind_line_candidate(v)
    feat=shadow.features(v,path)
    val=float(((feat-mean)/scale)@coef+intercept)
    if not np.isfinite(val):raise ValueError("NONFINITE_FROZEN_COSMOS_SCORE")
    return path[0],val,float(ranking[0])


def manifest(row,path,score,model_lock):
    A="measured_scan:"+row["source_file_sha256"][:16]
    B="real_16x128_panel:"+row["measured_panel_sha256"][:16]
    C="cosmos_frozen_score:"+hash_bytes(struct.pack(">d",score))[:16]
    return {
        "source_file_sha256":row["source_file_sha256"],
        "panel_sha256":row["measured_panel_sha256"],
        "candidate_path_sha256":hash_array(path,dtype="<i4"),
        "frozen_score_sha256":hash_bytes(struct.pack(">d",score)),
        "frozen_model_lock_sha256":model_lock,
        "window_id":str(row["scan"])+":"+str(row["panel_idx"]),
        "graph":{
            "schema":"janus.tachyon.r5f.synthetic_reverse_road_graph.v1",
            "case_ref":str(row["scan"])+":"+str(row["panel_idx"]),
            "terminal":C,
            "terminal_kind":"OBSERVED",  # computed output, not astronomical truth
            "forward_edges":{"A_B":{A:[B]},"B_C":{B:[C]}},
            "reverse_proposals":{"C_B":{C:[B]},"B_A":{B:[A]},"complete_certified":True},
            "witness_reports":[],
            "witness_window_seconds":120,
            "synthetic_only":True,
        },
    }


def mutate(m,kind):
    new=copy.deepcopy(m)
    if kind=="INTACT":pass
    elif kind=="RAW_SHA":new["source_file_sha256"]="0"*64
    elif kind=="PANEL_SHA":new["panel_sha256"]="0"*64
    elif kind=="SCORE_SHA":new["frozen_score_sha256"]="0"*64
    elif kind=="MODEL_LOCK":new["frozen_model_lock_sha256"]="0"*64
    elif kind=="GHOST_EDGE":
        b=next(iter(new["graph"]["reverse_proposals"]["B_A"]))
        new["graph"]["reverse_proposals"]["B_A"][b].append("unverified_unknown")
    elif kind=="INCOMPLETE_INDEX":
        new["graph"]["reverse_proposals"]["complete_certified"]=False
    else:raise ValueError("UNDECLARED_OVERLAY")
    return new


def direct_check(m,row,path,score,model_lock,model):
    # Independent forward replay of real measurements with the *same*
    # locked native search and scorer; avoids comparing to a mutable guess.
    ret_path,ret_score,_=score_one(row["panel"],model)
    if not np.array_equal(path,ret_path) or score!=ret_score:return False
    expected=manifest(row,path,score,model_lock)
    if any(m.get(k)!=expected.get(k) for k in
           ("source_file_sha256","panel_sha256","candidate_path_sha256",
            "frozen_score_sha256","frozen_model_lock_sha256","window_id")):
        return False
    g=m.get("graph",{})
    if (g.get("forward_edges")!=expected["graph"]["forward_edges"] or
        g.get("terminal")!=expected["graph"]["terminal"] or
        not g.get("reverse_proposals",{}).get("complete_certified")):
        return False
    p=g["reverse_proposals"]
    return p.get("C_B")==transpose(g["forward_edges"]["B_C"]) and p.get("B_A")==transpose(g["forward_edges"]["A_B"])


def paired_audit(rows,model,lock):
    methods=("COSMOS_NATIVE_FROZEN","COSMOS_PLUS_TACHYON_MIRROR_REPLAY",
             "COSMOS_PLUS_STRONG_DIRECT_PROVENANCE")
    stats={name:Counter() for name in methods}
    scans=[]
    overlay_cases=("INTACT","RAW_SHA","PANEL_SHA","SCORE_SHA","MODEL_LOCK",
                   "GHOST_EDGE","INCOMPLETE_INDEX")
    for row in rows:
        path,score,ranking=score_one(row["panel"],model)
        scans.append({"scan":row["scan"],"role":row["role"],"panel_idx":row["panel_idx"],
                      "native_line_path":path.tolist(),"frozen_shadow_score":score,
                      "blind_candidate_ranking":ranking,
                      "frequency_range_mhz":[row["native_freq_start_mhz"],
                                              row["native_freq_end_mhz"]],
                      "measured_panel_sha256":row["measured_panel_sha256"]})
        original=manifest(row,path,score,lock)
        for kind in overlay_cases:
            manifest_candidate=mutate(original,kind)
            label=kind=="INTACT"
            # Each algorithm sees exactly the same real panel and manifest.
            v_direct=direct_check(manifest_candidate,row,path,score,lock,model)
            native_accept=True  # no extra metadata guard in stock scorer
            mirror=mirror_case(manifest_candidate["graph"])
            mirror_out=mirrored_replay(mirror,budget=500)
            mirror_accept=mirror_out["prediction"]!=ABSTAIN and v_direct
            accepts=(native_accept,mirror_accept,v_direct)
            for name,accept in zip(methods,accepts):
                s=stats[name]
                s["cases"]+=1;s["intact_accepted"]+=int(label and accept)
                s["intact_rejected"]+=int(label and not accept)
                s["corrupted_accepted"]+=int(not label and accept)
                s["corrupted_rejected"]+=int(not label and not accept)
    return {
        "overlay_scenarios":list(overlay_cases),
        "paired_arms":{n:dict(s) for n,s in stats.items()},
        "real_candidate_panels":scans,
        "same_native_spectrum_path_score_retained":True,
        "strong_direct_vs_mirror_metadata_verdict_equal":
           all(stats["COSMOS_PLUS_TACHYON_MIRROR_REPLAY"][k]==stats["COSMOS_PLUS_STRONG_DIRECT_PROVENANCE"][k]
               for k in ("intact_accepted","corrupted_rejected",
                         "corrupted_accepted","intact_rejected")),
    }


def record_result(r):
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"sit1_open_receipt.json").write_text(json.dumps(r,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("SIT1_OPEN_RECEIPT",json.dumps({"status":r["status"],"real_files":r.get("replayed_real_files"),
          "errors":r.get("errors",[]),"source_files":len(r.get("downloads",[]))},
          ensure_ascii=False),flush=True)


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    folder=OUT/"real_hdf5"
    folder.mkdir(exist_ok=True)
    r=copy.deepcopy(FAIL)
    r["source_catalog"]=Q["immutable_real_archive"]["public_index"]
    r["expected_scans"]=6
    r["sha256_provider_status"]="FIRST_SEEN_SHA256_COMPUTED_NOT_PUBLISHED_REFERENCE"
    all_scans=[]
    for scan in Q["immutable_real_archive"]["scans"]:
        log=retrieve(scan,folder)
        r["downloads"].append(log)
        if not log.get("ok"):
            r["errors"].append({"scan":scan["number"],"error":log.get("error","HTTP_FAILURE")})
            continue
        try:
            q=reduce_real_h5(folder/scan["path"],log)
            all_scans.append(q)
            log["data_shape"]=q["dataset_shape"]
            log["header"]=q["header"]
            log["panels"]=q["selected_frequency_segments"]
        except Exception as exc:
            r["errors"].append({"scan":scan["number"],"error_type":type(exc).__name__,
                                "error":str(exc)[:320]})
        finally:
            (folder/scan["path"]).unlink(missing_ok=True)
    r["replayed_real_files"]=len(all_scans)
    if len(all_scans)!=6:
        r["status"]="PUBLIC_3IATLAS_DATA_PARTIAL_OR_INCOMPATIBLE__SIT1_OPEN_HOLD"
        record_result(r)
        return r
    # Archive roles are provided, and not independently measured pointing.
    r["archive_on_count"]=sum(x["record"]["role"]=="ON" for x in all_scans)
    r["archive_off_count"]=sum(x["record"]["role"]=="OFF" for x in all_scans)
    if (r["archive_on_count"],r["archive_off_count"])!=(3,3):
        r["errors"].append({"error":"WRONG_ABACAD_ON_OFF_COUNT"})
        record_result(r)
        return r
    # The same physical freq bins must actually overlap for paired ON/OFF.
    ref=all_scans[0]["panels"]
    aligned=all(all(
        abs(x["native_freq_start_mhz"]-y["native_freq_start_mhz"])<1e-5 and
        abs(x["native_freq_end_mhz"]-y["native_freq_end_mhz"])<1e-5
        for x,y in zip(ref,z["panels"])) for z in all_scans[1:])
    r["same_native_physical_freq_grid_across_six"]=aligned
    if not aligned:
        r["errors"].append({"error":"FREQUENCY_GRIDS_NOT_MATCHED_TO_FIXED_1E-5_MHZ_TOLERANCE"})
        record_result(r)
        return r
    rows=[p for scan in all_scans for p in scan["panels"]]
    r["real_panels_count"]=len(rows)
    if len(rows)!=24:raise AssertionError("wrong preregistered panel count")
    model=cosmos.frozen_model()
    lock=cosmos.Q["frozen_model"]["fit_lock_sha256"]
    check=paired_audit(rows,model,lock)
    r.update(check)
    r["status"]="REAL_GBT_3IATLAS_TIME_FREQUENCY_ONOFF_OOD_PILOT_COMPUTED__SCIENCE_HOLD"
    r["native_R42A8A_16x128_algorithm_invoked"]=True
    r["synthetic_truth_labels_used"]=False
    r["production_decisions_changed"]=0
    r["independent_actual_pointing_proved"]=False
    r["external_custodian"]=False
    r["real_ETI_or_false_alarm_impact_proven"]=False
    r["full_band_sensitivity_tested"]=False
    r["scan_groups_independent"]=1
    record_result(r)
    return r


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--audit",action="store_true")
    args=p.parse_args()
    if not args.audit:p.error("--audit required")
    outcome=run()
    if outcome["status"].startswith("REAL_GBT_3IATLAS"):
        x=outcome["paired_arms"]["COSMOS_PLUS_TACHYON_MIRROR_REPLAY"]
        y=outcome["paired_arms"]["COSMOS_PLUS_STRONG_DIRECT_PROVENANCE"]
        assert x==y, "MIRROR_DIRECT_DISAGREEMENT"
        assert x["intact_accepted"]==24 and x["corrupted_rejected"]==144
        assert outcome["production_decisions_changed"]==0
        print("JANUS_SIT1_OPEN_REAL_GBT_ENGINEERING_PILOT_PASS__SCIENTIFIC_HOLD",flush=True)
    else:
        print("JANUS_SIT1_OPEN_REAL_SOURCE_OR_FORMAT_HOLD__NO_SCIENTIFIC_RESULT",flush=True)


if __name__=="__main__":
    main()
