"""JANUS SIT2E-D: full-byte verified July/Sep polarization and 806 MHz witness.

FAIL CLOSED. These observations cannot determine IF/Stokes ordering by themselves.
Never download October or fit, replace, or relax a detector threshold.
Source/Reference are separated pointings, not simultaneous sky labels.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import sys
import numpy as np
import h5py
import hdf5plugin
import requests

ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"research"/"seti"))
import janus_sit2ec_shift_diagnostics as c
import janus_sit2eb_power_calibration as b
import janus_sit2e_parkes_range_headers as hdr
import janus_sit2eb_full_integrity as integrity

OUT=ROOT/"runtime"/"sit2ed"
DATES=("2025-07-31","2025-09-12")
F0=704.
DF=.002
CSTART=51000  # 806.000 MHz; separate from July's previously flagged 806.194 MHz
CSTOP=51376  # end exclusive, 806.750 MHz
REGION=[F0+CSTART*DF,F0+(CSTOP-1)*DF]
STUDY_START=51097 # 806.194 MHz
STUDY_STOP=51225 # end exclusive, 806.448 MHz
FRACTIONS=b.FRACTIONS
FROZEN_RECEIPT=ROOT/"research"/"seti"/"janus_sit2eb_july_september_calibration_receipt.json"
LOCK=json.loads(FROZEN_RECEIPT.read_text())
POLARIZATION_DOCS={
 "3i_atlas":"https://seti.berkeley.edu/atlas/",
 "rawspec":"https://github.com/UCBerkeleySETI/rawspec",
 "note":"Provider: July full polarization; later dates total intensity. rawspec -p 1 total, -p 4 cross-correlations, -p -4 full Stokes; NIFS=4 cannot identify mode."
}

def scalar(v):
    if isinstance(v,bytes):return v.decode("utf-8","replace")
    if isinstance(v,np.ndarray):
        if v.size==1:return scalar(v.reshape(-1)[0])
        return [scalar(z) for z in v.reshape(-1)[:8]]
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,(str,float,int,bool)) or v is None:return v
    return str(v)[:120]

def quant(a):
    x=np.asarray(a,dtype=np.float64)
    if not np.isfinite(x).all() or not x.size:raise ValueError("NONFINITE_OR_MISSING_POWER")
    return {"median":float(np.median(x)),
            "p05":float(np.quantile(x,.05)),"p95":float(np.quantile(x,.95)),
            "minimum":float(np.min(x)),"maximum":float(np.max(x)),
            "negative_fraction":float(np.mean(x<0)),
            "zero_fraction":float(np.mean(x==0))}

def safe_corr(a,bb):
    x=np.asarray(a,dtype=np.float64).ravel();y=np.asarray(bb,dtype=np.float64).ravel()
    if x.size!=y.size or x.std()<1e-12 or y.std()<1e-12:return None
    return float(np.corrcoef(x,y)[0,1])

def robust_rows(block):
    med=np.median(block,axis=1)
    sig=np.maximum(1.4826*np.median(np.abs(block-med[:,None]),axis=1),.02)
    return (block-med[:,None])/sig[:,None]

def power_cube(path,date):
    """Return one 300x4x376 *actual* frequency grid after full SHA check."""
    if date==DATES[0]:
        with path.open("rb") as stream:h=hdr.parse_sigproc_header(stream.read(8192))
        if (h["nifs"],h["nchans"],h["nbits"])!=(4,64000,32):
            raise ValueError("JULY_POLARIZATION_AXES_CHANGED")
        if path.stat().st_size!=h["header_bytes"]+300*4*64000*4:
            raise ValueError("BAD_JULY_SIGPROC_PAYLOAD")
        arr=np.memmap(path,dtype="<f4",mode="r",offset=h["header_bytes"],shape=(300,4,64000))
        cube=np.asarray(arr[:,:,CSTART:CSTOP],dtype=np.float64).copy()
        keys={k:scalar(v) for k,v in h.items() if k in ("nifs","nbits","telescope_id","machine_id","source_name","rawdatafile","fch1","foff","tsamp","tstart","src_raj","src_dej")}
        keys["format"]="SIGPROC_FIL";keys["pol_type_header_absent"]="pol_type" not in h
        del arr
    else:
        with h5py.File(path,"r") as handle:
            d=handle["data"]
            if tuple(d.shape)!=(300,4,64000):raise ValueError("HDF5_POLARIZATION_AXES_CHANGED")
            attrs={str(k).lower():scalar(v) for k,v in handle.attrs.items()}
            attrs.update({str(k).lower():scalar(v) for k,v in d.attrs.items()})
            keys={k:attrs.get(k) for k in ("nifs","nbits","source_name","pol_type","polarization","stokes","npol","foff","fch1","tsamp","tstart","src_raj","src_dej","telescope_id","machine_id")}
            keys["all_attribute_keys"]=sorted(attrs.keys())
            keys["format"]="HDF5"
            cube=np.asarray(d[:,:,CSTART:CSTOP],dtype=np.float64)
    if cube.shape!=(300,4,CSTOP-CSTART) or not np.isfinite(cube).all():
        raise ValueError("INVALID_806_POWER_ARRAY")
    return cube,keys

def candidate_profiles(cube, role):
    """Pure physical-frequency line and IF diagnostics. Nothing is fitted."""
    out=[]
    for iff in range(4):
        a=cube[:,iff,:]
        z=robust_rows(a)
        spec=np.median(z,axis=0)
        peaks=np.argsort(spec)[-8:][::-1]
        strong=spec>5
        out.append({
          "if_index":iff,"role":role,
          "raw_distribution":quant(a),
          "normalized_z_distribution":quant(z),
          "median_frequency_profile_top_8":[
              {"frequency_mhz":round(F0+(CSTART+int(idx))*DF,6),
               "frequency_channel":CSTART+int(idx),
               "median_time_robust_z":float(spec[idx]),
               "time_fraction_z_gt_5":float(np.mean(z[:,idx]>5))}
              for idx in peaks],
          "channels_with_median_time_robust_z_gt5":int(strong.sum()),
          "region806_194_to806_448":{
            "median_robust_z_peak":float(np.max(spec[STUDY_START-CSTART:STUDY_STOP-CSTART])),
            "above5_bins":int(np.sum(strong[STUDY_START-CSTART:STUDY_STOP-CSTART]))}
        })
    return out

def common_line(cubes):
    # Frequency is the same, but the pointings are ~342/347 s apart.
    rows=[]
    for date in DATES:
        source=cubes[(date,"source")]
        off=cubes[(date,"reference")]
        for iff in range(4):
            za=robust_rows(source[:,iff,:])
            zb=robust_rows(off[:,iff,:])
            a=np.median(za,axis=0)
            d=np.median(zb,axis=0)
            shared=np.where((a>5)&(d>5))[0]
            union=np.where((a>5)|(d>5))[0]
            best_both=int(np.argmax(np.minimum(a,d)))
            rows.append({"date":date,"if_index":iff,
                "shared_z_gt5_frequency_bins":int(shared.size),
                "union_z_gt5_frequency_bins":int(union.size),
                "jaccard":float(shared.size/union.size) if union.size else None,
                "median_spectral_shape_correlation":safe_corr(a,d),
                "max_min_two_role_line_z":float(np.minimum(a,d)[best_both]),
                "max_min_two_role_frequency_mhz":round(F0+(CSTART+best_both)*DF,6),
                "two_roles_separated_not_simultaneous":True})
    return rows

def polarization_hypothesis_diagnostics(cubes):
    # Both are *probes*, not claims about backend polarization mode or Stokes.
    results=[]
    for date in DATES:
        for role in ("source","reference"):
            a=cubes[(date,role)]
            # Hypothesis H1 = [I,Q,U,V] index order, physical I >= sqrt(Q^2+U^2+V^2).
            # Noise and uncalibrated channels can violate; it is NOT a conclusive test.
            valid_stokes=np.mean(a[:,0,:]+1e-6>=np.sqrt(np.sum(a[:,1:,:]**2,axis=1)))
            # Hypothesis H2 = [XX,YY,ReXY,ImXY] (one possible product order).
            # Cauchy-Schwarz: XX*YY >= ReXY^2 + ImXY^2.
            xx=a[:,0,:];yy=a[:,1,:]
            valid_cross=np.mean((xx>=0)&(yy>=0)&(xx*yy+1e-6>=a[:,2,:]**2+a[:,3,:]**2))
            results.append({
                "date":date,"role":role,
                "stokes_IQUV_inequality_pass_fraction_if_assumed_order_0_1_2_3":float(valid_stokes),
                "cross_XX_YY_ReIm_inequality_pass_fraction_if_assumed_order_0_1_2_3":float(valid_cross),
                "sampled_cells":int(a.shape[0]*a.shape[2]),
                "conditional_not_proof_of_polarization_order":True,
                "per_IF_806_distributions":[quant(a[:,j,:]) for j in range(4)],
                "per_IF_806_robust_correlation_vs_IF0":[safe_corr(a[:,0,:],a[:,j,:]) for j in range(1,4)]
                })
    return results

def comparable_processing_probe(per_if):
    """Actually compute frozen native/matched scores under explicit basis hypotheses.

    No new thresholds or pretense of known physical polarization mapping.
    """
    records=[]
    for date in DATES:
        for band in range(4):
            for mode in (("IF0_only",) if date==DATES[1]
                         else ("IF0_only","IF0_plus_IF1_conditional")):
                for role in ("source","reference"):
                    x=per_if[(date,role)][0][band]["panel"].copy()
                    if mode=="IF0_plus_IF1_conditional":
                        x=x+per_if[(date,role)][1][band]["panel"]
                    scores=b.score(x)
                    records.append({
                        "date":date,"role":role,"band":band,
                        "frequency_start_mhz":round(704+int((64000-128)*FRACTIONS[band])*.002,6),
                        "hypothesis":mode,
                        "frozen_original_objective_scores":{
                            k:float(v[1]) for k,v in scores.items()},
                        "panel_sha256":hashlib.sha256(x.tobytes()).hexdigest(),
                        "polarization_basis_verified":False,
                        "threshold_fitted_or_changed":False,
                        "eligible_for_science_or_october":False
                    })
    assert len(records)==24
    return records

def run():
    OUT.mkdir(parents=True,exist_ok=True)
    folder=OUT/"verified_orig";folder.mkdir(exist_ok=True)
    manifest=json.loads(integrity.MANIFEST.read_text())
    selected=integrity.selected_files(manifest)
    if len(selected)!=4 or tuple(dict.fromkeys(x[0] for x in selected))!=DATES:
        raise ValueError("ONLY_JULY_SEPT_PERMITTED")
    result={"schema":"janus.cosmos.sit2ed.polarization_frequency_witness.v1",
       "status":"INCOMPLETE_HOLD","dates":list(DATES),
       "original_sit2eb_thresholds":LOCK["numeric_thresholds"],
       "calibration_lock_sha256":LOCK["computation_sha256"],
       "october_files_read":False,"october_files_requested":False,
       "thresholds_changed":False,"production_models_changed":False,
       "polarization_docs":POLARIZATION_DOCS,
       "candidate_region_mhz":REGION,"target_806_mhz_window_mhz":[806.194,806.448],
       "hash_receipts":[],"per_scan_headers":[],"per_role_if_profiles":[],
       "if_order_hypotheses":[],"paired_frequency_structure":[],
       "conditional_total_intensity_mapping":[],"errors":[],
       "proven_stokes_channel_mapping":False,
       "independent_tcs_encoder_pointing":False,
       "real_transmitter_identified":False,
       "science_gate":"HOLD",
       "allow_october_exam":False}
    cubes={}
    all_panels={}
    try:
        with requests.Session() as sess:
            for date,role,spec in selected:
                file=None
                try:
                    file,receipt=c.download_verified(sess,date,role,spec,folder)
                    cube,h=power_cube(file,date)
                    _,per_if=c.extract_all_if(file,date)
                    all_panels[(date,role)]=per_if
                    cubes[(date,role)]=cube
                    result["hash_receipts"].append(receipt)
                    result["per_scan_headers"].append({"date":date,"role":role,"header":h})
                    result["per_role_if_profiles"].append({
                       "date":date,"role":role,"profiles":candidate_profiles(cube,role)})
                finally:
                    if file is not None:file.unlink(missing_ok=True)
        result["paired_frequency_structure"]=common_line(cubes)
        result["if_order_hypotheses"]=polarization_hypothesis_diagnostics(cubes)
        result["conditional_total_intensity_mapping"]=comparable_processing_probe(all_panels)
        result["status"]="SIT2ED_PHYSICAL_806_FREQ_AND_IF_HYPOTHESES_PASS__MAPPING_AND_POINTING_HOLD"
    except Exception as exc:
        result["errors"].append({"type":type(exc).__name__,"message":str(exc)[:400]})
    finally:
        for x in folder.glob("*"):x.unlink(missing_ok=True)
        output=OUT/"sit2ed_forensics.json"
        output.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+"\n")
        result["report_sha256"]=hashlib.sha256(output.read_bytes()).hexdigest()
        print("SIT2ED_RESULT",json.dumps({
            "status":result["status"],"errors":result["errors"],
            "files":len(result["hash_receipts"]),
            "report_sha256":result["report_sha256"],
            "october_files_read":result["october_files_read"],
            "role_pair_frequency_structure":result["paired_frequency_structure"]},
            sort_keys=True),flush=True)
    if result["errors"]:raise RuntimeError(result["errors"][0]["message"])
    return result

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--audit",action="store_true",required=True)
    parser.parse_args()
    run()
