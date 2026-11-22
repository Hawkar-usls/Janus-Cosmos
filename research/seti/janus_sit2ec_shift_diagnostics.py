"""SIT2E-C forensic domain-shift audit; read July/Sep ONLY, October is sealed.

No threshold changes, no search over alternate operating points. Instrument IF
channels are not called Stokes/polarization without verified header semantics.
Reports physical sample summaries, original fixed threshold replay, stratified
controls and effective sample counts. This is observational-background analysis
with artificial evaluation overlays, NOT real-source astrophysical ground truth.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import pathlib
import sys
from collections import defaultdict

import h5py
import hdf5plugin  # Required for archived BL compression filters.
import numpy as np
import requests

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "research/seti"))
import janus_sit2eb_full_integrity as integrity
import janus_sit2eb_power_calibration as baseline
import janus_sit2e_parkes_range_headers as sigproc

DATES = ("2025-07-31", "2025-09-12")
OUT = ROOT / "runtime/sit2ec"
RECEIPT = ROOT / "research/seti/janus_sit2eb_july_september_calibration_receipt.json"
Q = json.loads(RECEIPT.read_text())
LOCKED = Q["numeric_thresholds"]
METHODS = baseline.METHODS
assert set(LOCKED) == set(METHODS)
FRACTIONS = baseline.FRACTIONS
SEED = baseline.SEED
EXPECTED_CONTROLS = {
    "JANUS_NATIVE_RIDGE": {"2025-07-31": 74, "2025-09-12": 9},
    "DIRECT_672_MATCHED": {"2025-07-31": 74, "2025-09-12": 9},
    "DIRECT_672_SIDEBAND": {"2025-07-31": 63, "2025-09-12": 9},
}

def nums(a):
    a=np.asarray(a,dtype=np.float64)
    if not a.size or not np.isfinite(a).all():
        raise ValueError("NONFINITE_OR_EMPTY_DIAGNOSTIC_VALUES")
    quantiles=np.quantile(a,[0,.05,.25,.5,.75,.95,.99,1])
    return {"n":int(a.size),"min":float(quantiles[0]),
            "q05":float(quantiles[1]),"q25":float(quantiles[2]),
            "median":float(quantiles[3]),"q75":float(quantiles[4]),
            "q95":float(quantiles[5]),"q99":float(quantiles[6]),
            "max":float(quantiles[7]),"mean":float(a.mean()),
            "std":float(a.std())}

def attr(v):
    if isinstance(v,bytes):return v.decode("utf-8","replace")
    if isinstance(v,np.ndarray):
        if v.size==1:return attr(v.reshape(-1)[0])
        return [attr(x) for x in v.reshape(-1)[:4]]
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,(bool,int,float,str)) or v is None:return v
    return repr(v)[:100]

def extract_all_if(path,date):
    panels={}
    if date == DATES[0]:
        h=sigproc.parse_sigproc_header(path.open("rb").read(8192))
        if h["nifs"]!=4 or h["nbits"]!=32 or h["nchans"]!=64000:
            raise ValueError("JULY_IF_LAYOUT_CHANGED")
        if path.stat().st_size != h["header_bytes"]+300*4*64000*4:
            raise ValueError("JULY_BYTES_AND_IF_OFFSET_DISAGREE")
        data=np.memmap(path,dtype="<f4",mode="r",offset=h["header_bytes"],shape=(300,4,64000))
        head={k:attr(h.get(k)) for k in ("source_name","fch1","foff","tsamp","tstart","nifs","nbits","src_raj","src_dej")}
        head["container"]="SIGPROC_FIL"; head["layout"]="300x4x64000"; head["native_dtype"]="<f4"
        def get_block(i,a,b):
            return np.asarray(data[:,i,a:b],dtype=np.float64)
    else:
        f=h5py.File(path,"r")
        data=f["data"]
        if tuple(data.shape)!=(300,4,64000):raise ValueError("SEPT_H5_ARRAY_LAYOUT_CHANGED")
        attrs={str(k).lower():attr(v) for k,v in f.attrs.items()}
        attrs.update({str(k).lower():attr(v) for k,v in data.attrs.items()})
        head={k:attrs.get(k) for k in ("source_name","fch1","foff","tsamp","tstart","nifs","nbits","src_raj","src_dej")}
        head["container"]="HDF5";head["layout"]=list(data.shape)
        head["native_dtype"]=str(data.dtype)
        head["compression"]=str(data.compression)
        def get_block(i,a,b):
            return np.asarray(data[:,i,a:b],dtype=np.float64)
    try:
        if abs(float(head["fch1"])-704)>1e-6 or abs(float(head["foff"])-.002)>1e-9:
            raise ValueError("PHYSICAL_FREQUENCY_GRID_CHANGED")
        if abs(float(head["tsamp"])-1.)>1e-9:raise ValueError("NATIVE_TIME_GRID_CHANGED")
        for i in range(4):
            panels[i]=[]
            for frac in FRACTIONS:
                start=int((64000-128)*frac)
                block=get_block(i,start,start+128)
                if block.shape!=(300,128) or not np.isfinite(block).all():
                    raise ValueError("BAD_NATIVE_SPECTRAL_BLOCK")
                # Identical time averaging as locked SIT2E-B; no channel rebin.
                panel=np.stack([x.mean(axis=0) for x in np.array_split(block,16)]).astype(np.float32)
                z,mad=baseline.robust_normalize(panel)
                med=np.median(panel,axis=1)
                pol={"if_index":i,"band_index":len(panels[i]),
                     "channel_start":start,"freq_mhz":round(704+start*.002,6),
                     "panel_sha256":hashlib.sha256(panel.tobytes()).hexdigest(),
                     "raw_power":nums(panel.ravel()),
                     "native_raw_nonpositive_fraction":float(np.mean(block<=0)),
                     "native_raw_negative_fraction":float(np.mean(block<0)),
                     "native_raw_zero_fraction":float(np.mean(block==0)),
                     "normalization":{
                         "row_median":nums(med),
                         "row_robust_sigma_MAD":nums(mad),
                         "sigma_floor_fraction":float(np.mean(mad<=.020000001)),
                         "normalized_z":nums(z.ravel()),
                         "fraction_normalized_z_gt5":float(np.mean(z>5)),
                         "spectral_persistent_channel_occupancy_gt3":
                             float(np.mean((z>3).sum(axis=0)>=8)),
                         "row_median_span":float(np.ptp(med))
                     },
                     "frozen_uninjected_max_path_scores":{
                         k:float(t[1]) for k,t in baseline.score(panel).items()
                     }}
                panels[i].append({"panel":panel,"audit":pol})
        # Binary reproduction against SIT2E-B IF0 decoder, without changing it.
        original=baseline.read_panels(path,date,None)
        for i in range(4):
            if not np.array_equal(original[i]["panel"],panels[0][i]["panel"]):
                raise AssertionError("IF0_NOT_BYTE_EXACT_WITH_SIT2EB")
        return head,panels
    finally:
        if date==DATES[0]:del data
        else:f.close()

def download_verified(session,date,role,spec,folder):
    dst=folder/spec["filename"]
    n=0
    try:
        with session.get(spec["url"],stream=True,timeout=(20,100),headers={
                "Accept-Encoding":"identity"},allow_redirects=False) as res:
            res.raise_for_status()
            if res.status_code!=200 or res.headers.get("Content-Encoding","identity").lower() not in ("","identity"):
                raise ValueError("UNEXPECTED_HTTP_RESPONSE")
            if res.headers.get("Content-Length") and int(res.headers["Content-Length"])!=spec["bytes"]:
                raise ValueError("CHANGED_CONTENT_LENGTH")
            with dst.open("wb") as f:
                for c in res.iter_content(chunk_size=4*1024*1024):
                    if not c:continue
                    n+=len(c)
                    if n>spec["bytes"]:raise ValueError("OVERLONG_ARCHIVE_SOURCE")
                    f.write(c)
        digest=baseline.assert_file(dst,spec)
        return dst,{"date":date,"role":role,"filename":spec["filename"],"bytes":n,"sha256":digest}
    except Exception:
        dst.unlink(missing_ok=True)
        raise

def summarize_synthetic_controls(all_panels):
    # Preserve original SIT2E-B RNG ordering, data and unchanged Sept thresholds.
    # No alternative thresholds, tuning or October access.
    rng=np.random.default_rng(SEED)
    case_rows=[]
    for date in DATES:
        for band in range(4):
            for role in ("source","reference"):
                base=all_panels[(date,role)][0][band]["panel"]
                for rep in range(12):
                    for category in baseline.ALL_CLASSES:
                        new,truth=baseline.create_overlay(base,category,rng)
                        r=baseline.score(new)
                        case_rows.append({"date":date,"band":band,"role":role,
                                          "rep":rep,"kind":category,
                                          "scores":{k:float(v[1]) for k,v in r.items()},
                                          "located":{k:baseline.localizes(v[0],truth) for k,v in r.items()}})
    if len(case_rows)!=768:raise AssertionError("CASE_COUNT_DRIFT")
    breakdown={};date_metrics={}
    for method in METHODS:
        breakdown[method]={};date_metrics[method]={}
        for date in DATES:
            rows=[x for x in case_rows if x["date"]==date]
            neg=[x for x in rows if x["kind"] in baseline.CONTROL_CLASSES]
            pos=[x for x in rows if x["kind"] in baseline.POS_CLASSES]
            threshold=float(LOCKED[method])
            exceeds=sum(x["scores"][method]>threshold for x in neg)
            if exceeds!=EXPECTED_CONTROLS[method][date]:
                raise AssertionError("PRIOR_FROZEN_SIT2EB_EXCEEDANCES_NOT_REPRODUCED")
            by_group=[]
            for band in range(4):
                for role in ("source","reference"):
                    for category in baseline.CONTROL_CLASSES:
                        samples=[x["scores"][method] for x in neg if x["band"]==band and x["role"]==role and x["kind"]==category]
                        by_group.append({"band":band,"role":role,"class":category,
                            "cases":len(samples),"exceedances":sum(x>threshold for x in samples),
                            "raw_score":nums(samples)})
            date_metrics[method][date]={
                "calibration_threshold_unchanged":threshold,
                "negative_counts":len(neg),
                "negative_exceedances":exceeds,
                "negative_fraction":exceeds/len(neg),
                "neg_scores":nums([x["scores"][method] for x in neg]),
                "no_injection_scores":nums([x["scores"][method] for x in neg if x["kind"]=="NO_INJECTION"]),
                "engineered_RFI_scores":nums([x["scores"][method] for x in neg if x["kind"]=="SYNTHETIC_RFI"]),
                "strong_detected":sum(x["scores"][method]>threshold and x["located"][method] for x in pos if x["kind"]=="STRONG_LINE"),
                "weak_detected":sum(x["scores"][method]>threshold and x["located"][method] for x in pos if x["kind"]=="WEAK_LINE"),
                "diagnostic_p95_NOT_AN_ADOPTABLE_THRESHOLD":float(np.quantile(
                   [x["scores"][method] for x in neg],.95,method="higher")),
                "by_band_role_class":by_group,
                "independent_base_panels":8,
                "twelve_duplicate_no_injection_scores_per_base_panel":True
            }
        breakdown[method]["July_to_September_control_median_shift"] = (
            date_metrics[method][DATES[0]]["neg_scores"]["median"] -
            date_metrics[method][DATES[1]]["neg_scores"]["median"])
    return date_metrics,breakdown

def run():
    OUT.mkdir(parents=True,exist_ok=True)
    folder=OUT/"download"
    folder.mkdir(exist_ok=True)
    specs=integrity.selected_files(json.loads(integrity.MANIFEST.read_text()))
    if len(specs)!=4 or tuple(dict.fromkeys(x[0] for x in specs))!=DATES:
        raise ValueError("SPLIT_INVALID_OR_OCTOBER_SELECTED")
    report={"schema":"janus.cosmos.sit2ec.parkes_shift_forensic_v1",
            "status":"SIT2EC_DIAGNOSTIC_INCOMPLETE",
            "october_power_read":False,"october_source_requested":False,
            "thresholds_changed":False,"production_changed":False,
            "prior_sit2eb_calibration_lock":Q["computation_sha256"],
            "frozen_thresholds":LOCKED,
            "physical_IF_semantics":"UNVERIFIED: 4 native IFs may contain derived polarizations; do not assert Stokes labels",
            "frequency_coverage":"4 fixed windows of 128 channels = 1.024 MHz of 128 MHz band, per IF",
            "files":[],"headers":[],"panel_metrics":[],"IF_coherence":[],
            "paired_source_reference":[],"control_replay":{},
            "errors":[],"scientific_gate":"HOLD"}
    all_panels={}
    try:
        with requests.Session() as s:
            for date,role,spec in specs:
                path=None
                try:
                    path,verified=download_verified(s,date,role,spec,folder)
                    head,per_if=extract_all_if(path,date)
                    report["files"].append(verified)
                    report["headers"].append({"date":date,"role":role,**head})
                    all_panels[(date,role)]=per_if
                    for i in range(4):
                        report["panel_metrics"].extend([{
                            "date":date,"role":role,**p["audit"]} for p in per_if[i]])
                    for band in range(4):
                        z0,_=baseline.robust_normalize(per_if[0][band]["panel"])
                        for i in range(1,4):
                            zz,_=baseline.robust_normalize(per_if[i][band]["panel"])
                            x=z0.ravel();y=zz.ravel()
                            corr=(float(np.corrcoef(x,y)[0,1]) if x.std()>1e-9 and y.std()>1e-9 else None)
                            report["IF_coherence"].append({"date":date,"role":role,
                                "band":band,"if_index":i,"z_correlation_vs_IF0":corr,
                                "byte_exact_equal_to_IF0":bool(np.array_equal(per_if[0][band]["panel"],per_if[i][band]["panel"]))})
                finally:
                    if path is not None:path.unlink(missing_ok=True)
        for date in DATES:
            for i in range(4):
                for band in range(4):
                    a=all_panels[(date,"source")][i][band]["panel"]
                    b=all_panels[(date,"reference")][i][band]["panel"]
                    az,_=baseline.robust_normalize(a);bz,_=baseline.robust_normalize(b)
                    report["paired_source_reference"].append({
                        "date":date,"if_index":i,"band":band,
                        "raw_median_source":float(np.median(a)),
                        "raw_median_reference":float(np.median(b)),
                        "raw_median_ratio":float(np.median(a)/np.median(b)) if np.median(b)>0 else None,
                        "source_reference_normalized_z_RMS_difference":float(np.sqrt(np.mean((az-bz)**2)))})
        metrics,shift=summarize_synthetic_controls(all_panels)
        report["control_replay"]=metrics
        report["score_shift_summary"]=shift
        report["status"]="SIT2EC_REAL_POWER_IF_AND_CONTROL_DIAGNOSTICS_COMPLETE__NO_THRESHOLD_CHANGE__SCIENCE_HOLD"
        report["total_panel_audits"]=len(report["panel_metrics"])
        report["files_verified"]=len(report["files"])
        report["october_power_read"]=False
    except Exception as exc:
        report["errors"].append({"type":type(exc).__name__,"message":str(exc)[:500]})
    finally:
        (OUT/"sit2ec_shift_audit.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
        for f in folder.glob("*"):f.unlink(missing_ok=True)
    digest=hashlib.sha256((OUT/"sit2ec_shift_audit.json").read_bytes()).hexdigest()
    print("SIT2EC_RESULT",json.dumps({"status":report["status"],"files":len(report["files"]),
        "panel_audits":len(report["panel_metrics"]), "errors":report["errors"],
        "october_power_read":report["october_power_read"],"report_sha256":digest,
        "score_shift_summary":report.get("score_shift_summary",{})},sort_keys=True),flush=True)
    if report["errors"]:raise RuntimeError(report["errors"][0]["message"])

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--audit",action="store_true",required=True)
    p.parse_args()
    run()
