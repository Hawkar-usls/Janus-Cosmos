#!/usr/bin/env python3
"""Bounded on-archive *feasibility* probe for MOCSS open-shutter Spec2/x1d.

IMPORTANT: A virtual single MSA slit is not the published 4-dither,
open-shutter, full-row-background MOCSS reduction. Any returned x1d spectrum
is a *diagnostic pipeline product*, not a valid flux limit or source spectrum.
Actual science extraction requires MOCSS-specific 2D stacked calibration,
3-shutter throughput, variance and artificial-source recovery.

This probe preserves original RATE and MSA files and hashes the archive
sources. It records a receipt even if the official Spec2 pipeline rejects
the unusual test geometry.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import traceback

from fetch_jwst_mocss_mast_subset import read_and_store, BASE, ROOT_URI
from prepare_mocss_wcs_probe import create_probe

REPO=Path(__file__).resolve().parents[1]
RECEIPT=REPO/"domains/space/MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json"
MSA={
"A":"jw03290001001_05_msa.fits",
"B":"jw03290001001_05_msa.fits",
"C":"jw03290002001_01_msa.fits",
"D":"jw03290002001_01_msa.fits"
}
MSA_SHA={
"jw03290001001_05_msa.fits":"fb6239b287d25c478d2e22737ff956850a9bc1e0f365500cfc5661abbfe17a44",
"jw03290002001_01_msa.fits":"25a061a226521a8bbf2f30d3f6922a22fb13437800047f12c8511e354820c91e"
}
ORIGINAL_CONTEXT="jwst_1413.pmap"
MODALITY="VIRTUAL_SINGLE_BACKGROUND_SHUTTER_DIAGNOSTIC_NOT_SCIENCE"
def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
    return h.hexdigest()

def fetch(name, output_dir, expected_sha):
    from urllib.parse import quote
    p=output_dir/name
    read_and_store(BASE+quote(ROOT_URI+name,safe=""),p)
    actual=sha(p)
    if actual!=expected_sha: raise ValueError(f"Archive SHA mismatch {name}: {actual}")
    return p

def attempt(entry,base,model_name="TEST_PILOT"):
    from jwst.pipeline import Spec2Pipeline
    from astropy.io import fits
    import jwst
    os.environ["CRDS_CONTEXT"]=ORIGINAL_CONTEXT
    assert jwst.__version__=="1.19.1",jwst.__version__
    group=entry["group"]
    det=entry["expected_science_detector_from_WCS"]
    d=entry["detectors"][det]
    input_dir=base/"original_inputs";input_dir.mkdir(exist_ok=True,parents=True)
    rate=fetch(d["rate_file"],input_dir,d["sha256"])
    msa_name=MSA[group]; msa=fetch(msa_name,input_dir,MSA_SHA[msa_name])
    slit=d["selected_virtual_shutter"]
    assert slit is not None
    tmp=base/f"diag_{group}_{entry['exposure_index']}"
    tmp.mkdir(exist_ok=True,parents=True)
    probe=create_probe(rate,msa,tmp,grid=1,target={
        "quadrant":slit["quadrant"],"row":slit["row"],"column":slit["column"]},
        half_width=0)
    assert probe["num_virtual_slits"]==1
    out={
        "group":group,
        "exposure_index":entry["exposure_index"],
        "detector":det,
        "original_rate_file":rate.name,
        "original_rate_sha256":sha(rate),
        "original_msa_file":msa.name,
        "original_msa_sha256":sha(msa),
        "original_msametid":probe["original_msa_metadata_id"],
        "virtual_msametid_signed16":probe["normalized_signed16_id"],
        "virtual_shutter":slit,
        "mode":MODALITY,
        "pipeline_version":jwst.__version__,
        "crds_context":ORIGINAL_CONTEXT,
        "test_rate_sha256":sha(Path(probe["probe_rate"])),
        "test_msa_sha256":sha(Path(probe["probe_msa"])),
        "stage2_result":"NOT_ATTEMPTED",
        "calibrated_line_flux":None,
        "stage2_x1d_scientifically_valid":False,
        "required_next_mocss_reduction":"Group-aligned four-dither stack, custom full-row sky subtraction, MOCSS 3-adjacent-shutter photometry and injection."
    }
    try:
        # A deliberately bounded test, not a re-extraction: critical MOCSS-
        # specific reduction steps are intentionally not completed here.
        pipeline=Spec2Pipeline()
        pipeline.output_dir=str(tmp)
        pipeline.save_results=True
        out["available_stage2_steps"]=sorted(pipeline.step_defs)
        out["skipped_stage2_steps"]=[]
        out["missing_requested_skip_steps"]=[]
        for name in ("master_background_mos","bkg_subtract","pathloss","barshadow","resample_spec"):
            if hasattr(pipeline,name):
                getattr(pipeline,name).skip=True
                out["skipped_stage2_steps"].append(name)
            else:
                out["missing_requested_skip_steps"].append(name)
        # extract_1d is left enabled to test whether the pipeline actually
        # creates any formally readable x1d from a virtual background shutter.
        result=pipeline.run(probe["probe_rate"])
        out["stage2_result"]="CALL_RETURNED"
        out["returned_model_type"]=type(result).__name__
        models=result if isinstance(result,(list,tuple)) else ([] if result is None else [result])
        out["returned_model_count"]=len(models)
        out["returned_model_types"]=[type(m).__name__ for m in models]
        out["calibration_step_states_by_model"]=[]
        for model in models:
            meta=getattr(model,"meta",None)
            cal_step=getattr(meta,"cal_step",None)
            states={}
            for step in ("assign_wcs","extract_2d","srctype","wavecorr","flat_field","photom","extract_1d"):
                v=getattr(cal_step,step,None)
                if v is not None:states[step]=str(v)
            out["calibration_step_states_by_model"].append(states)
            if hasattr(model,"close"):model.close()
    except Exception as exc:
        out["stage2_result"]="PIPELINE_OR_CONFIGURATION_FAILURE_NOT_NEGATIVE_SPECTRAL_RESULT"
        out["exception_class"]=type(exc).__name__
        out["exception_short"]=str(exc)[:1200]
        out["traceback_tail"]=traceback.format_exc()[-4000:]
    files=[]
    for p in tmp.rglob("*.fits"):
        if p.name.startswith("TEST_"):continue
        try:
            with fits.open(p,memmap=True) as f:
                spec=[h for h in f if getattr(h,"name","").upper()=="EXTRACT1D"]
            files.append({"name":p.name,"bytes":p.stat().st_size,"sha256":sha(p),"has_extract1d":bool(spec)})
        except Exception as exc:
            files.append({"name":p.name,"inspection_error":str(exc)[:100]})
    out["stage2_output_fits"]=files
    out["x1d_table_count"]=sum(1 for f in files if f.get("has_extract1d"))
    out["astrophysical_calibrated_x1d"]=False
    if out["stage2_result"]=="CALL_RETURNED" and not out["x1d_table_count"]:
        out["stage2_result"]="NO_DIAGNOSTIC_X1D_REPRODUCED"
    out["source_redshift"]=None
    out["directional_trail_detected"]=False
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--data-dir",default="/tmp/janus_mocss_stage2_pilot")
    ap.add_argument("--output",default="/tmp/janus_mocss_stage2_pilot_receipt.json")
    ap.add_argument("--exposures",default="8,13")
    args=ap.parse_args()
    import jwst
    base=Path(args.data_dir);base.mkdir(parents=True,exist_ok=True)
    receipt=json.loads(RECEIPT.read_text())
    assert receipt["jwst_pipeline"]=="1.19.1" and receipt["crds"]==ORIGINAL_CONTEXT
    idx=[int(z) for z in args.exposures.split(",")]
    assert 1<=len(idx)<=2 and len(set(idx))==len(idx)
    mapped={e["exposure_index"]:e for e in receipt["exposures"]}
    assert all(i in mapped for i in idx)
    data={"schema":"janus-mocss-spec2-virtual-shutter-feasibility-v1",
          "method":MODALITY,
          "created_context":"2026-10-18 GO-3290",
          "pipeline_version":jwst.__version__,
          "CRDS_CONTEXT":ORIGINAL_CONTEXT,
          "original_source_receipt":str(RECEIPT.relative_to(REPO)),
          "results":[],
          "scientific_authority":"NONE_NOT_A_REDUCTION_OF_PUBLISHED_MOCSS_STACKS",
          "calibrated_flux_by_group":None,
          "oiii_companion_detected":None,"pa_gamma_detected":None,
          "redshift":None,"source_class":None,"gas_trail":None,
          "next_step":"MOCSS-aware full-stack custom sky, 2D calibration, PSF-throughput and injection-tested x1d. Reject single virtual background slit as astrophysical measurement."}
    output=Path(args.output);output.parent.mkdir(parents=True,exist_ok=True)
    for i in idx:
        try:
            r=attempt(mapped[i],base)
        except Exception as exc:
            r={"exposure_index":i,"stage2_result":"FAILED_BEFORE_PIPELINE",
               "exception_class":type(exc).__name__,"exception_short":str(exc)[:1000],
               "traceback_tail":traceback.format_exc()[-4000:]}
        data["results"].append(r)
        output.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
        print("STAGE2_PILOT",i,r["stage2_result"],r.get("exception_class",""),
              "OUTPUT",str(output),flush=True)
    print("SCIENTIFIC_VERDICT=NO_CALIBRATED_MOCSS_LINE_FLUX_CLAIM",flush=True)
    if any(r.get("stage2_result") not in ("CALL_RETURNED","NO_DIAGNOSTIC_X1D_REPRODUCED") for r in data["results"]):
        raise SystemExit("CI_FAIL: Stage2 pipeline exception or malformed data; saved receipt is NOT a successful calibration")
    if any(r.get("stage2_result") == "NO_DIAGNOSTIC_X1D_REPRODUCED" for r in data["results"]):
        print("STAGE2_PILOT_NOT_SCIENTIFIC_SUCCESS: No virtual-background-shutter x1d; need MOCSS source-aware extraction.",flush=True)

if __name__=="__main__":main()
