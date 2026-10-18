#!/usr/bin/env python3
"""Diagnostic SOURCE_INFO-aware NIRSpec MOCSS three-shutter pilot.

This uses the original GO-3290 NIRSpec RATE and MSA data and a synthetic,
non-authoritative *source-marked* test MSA generated in a temporary directory.
The inserted source does NOT exist in the original MPT catalog, and centering
at 0.5 is a geometry trial, not a verified exact PSF position. Never use the
outputs for line significance/flux without full 4-dither MOCSS custom sky,
throughput, slit/PSF corrections and realistic injection/recovery.
"""
import argparse
import json
import os
import traceback
from pathlib import Path

import numpy as np
from astropy.io import fits

from attempt_mocss_stage2_open_shutter_pilot import (
    REPO, RECEIPT, ORIGINAL_CONTEXT, MSA, MSA_SHA, sha, fetch
)
from prepare_mocss_wcs_probe import create_probe

SOURCE_ID=16004
RA,DEC=53.13127,-27.78908

def prepare_source_aware(entry,base):
    group=entry["group"]
    detector=entry["expected_science_detector_from_WCS"]
    det=entry["detectors"][detector]
    archive=base/"raw";archive.mkdir(exist_ok=True,parents=True)
    original_rate=fetch(det["rate_file"],archive,det["sha256"])
    msa_name=MSA[group]
    original_msa=fetch(msa_name,archive,MSA_SHA[msa_name])
    where=det["selected_virtual_shutter"]
    assert where is not None
    tmp=base/f"source_test_{group}_{entry['exposure_index']}"
    tmp.mkdir(parents=True,exist_ok=True)
    probe=create_probe(original_rate,original_msa,tmp,grid=1,
        target={"quadrant":where["quadrant"],"row":where["row"],
                "column":where["column"]},half_width=1)
    target_r,target_c=where["row"],where["column"]
    source_msa=Path(probe["probe_msa"])
    sky_trace=det["sky_registered_trace"]["16004A"]
    sky_position_x=0.5+float(sky_trace["slit_x"])
    sky_position_y=0.5+float(sky_trace["slit_y"])
    if not (0<sky_position_x<1 and 0<sky_position_y<1):
        raise ValueError("Sky WCS outside nominal microshutter fractional bounds")
    with fits.open(source_msa,memmap=False,mode="update") as f:
        shutter=f["SHUTTER_INFO"]
        rows=shutter.data
        selected=rows[(rows["shutter_quadrant"]==where["quadrant"])&
                      (rows["shutter_row"]==target_r)&
                      (np.abs(rows["shutter_column"].astype(int)-target_c)<=1)]
        assert len(selected)==3,[(v["shutter_row"],v["shutter_column"]) for v in selected]
        assert set(map(int,selected["shutter_column"]))==set(range(target_c-1,target_c+2))
        new=np.array(selected,dtype=rows.dtype,copy=True)
        # Three adjacent columns cross the dispersion/slit axis in MSA
        # metadata: only the central shutter contains the alleged source.
        new["slitlet_id"]=1
        new["source_id"]=SOURCE_ID
        new["estimated_source_in_shutter_x"]=np.nan
        new["estimated_source_in_shutter_y"]=np.nan
        for i in range(len(new)):
            central=int(new[i]["shutter_column"])==target_c
            new[i]["background"]="N" if central else "Y"
            new[i]["primary_source"]="Y" if central else "N"
            if central:
                new[i]["estimated_source_in_shutter_x"]=sky_position_x
                new[i]["estimated_source_in_shutter_y"]=sky_position_y
        f["SHUTTER_INFO"]=fits.BinTableHDU(
            data=new,header=shutter.header,name="SHUTTER_INFO")
        info=f["SOURCE_INFO"]
        source=np.zeros(1,dtype=info.data.dtype)
        source[0]["program"]=3290
        source[0]["source_id"]=SOURCE_ID
        source[0]["source_name"]="JANUS16004"
        source[0]["alias"]="TEST_ONLY"
        source[0]["ra"]=RA
        source[0]["dec"]=DEC
        source[0]["preimage_id"]=""
        source[0]["stellarity"]=1.0
        f["SOURCE_INFO"]=fits.BinTableHDU(
            data=source,header=info.header,name="SOURCE_INFO")
        f.flush(output_verify="silentfix")
    with fits.open(source_msa,memmap=False) as f:
        sh=f["SHUTTER_INFO"].data
        sources=f["SOURCE_INFO"].data
        assert len(sh)==3 and len(sources)==1
        assert sum(sh["primary_source"]=="Y")==1
        assert sum(sh["background"]=="N")==1
        center=sh[sh["background"]=="N"][0]
        assert abs(float(center["estimated_source_in_shutter_x"])-sky_position_x)<1e-6
        assert abs(float(center["estimated_source_in_shutter_y"])-sky_position_y)<1e-6
        assert int(sources[0]["source_id"])==SOURCE_ID
        assert abs(float(sources[0]["ra"])-RA)<1e-8
        assert len(set(sh["slitlet_id"]))==1
    return {
        "group":group,"detector":detector,
        "exposure_index":entry["exposure_index"],
        "original_rate":str(original_rate),"original_rate_sha256":sha(original_rate),
        "original_msa":str(original_msa),"original_msa_sha256":sha(original_msa),
        "test_rate":probe["probe_rate"],"test_msa":str(source_msa),
        "test_rate_sha256":sha(Path(probe["probe_rate"])),
        "test_msa_sha256":sha(source_msa),
        "source_id_diagnostic_only":SOURCE_ID,
        "source_ra_dec_deg":[RA,DEC],
        "synthetic_source_position_WCS_fractional":[sky_position_x,sky_position_y],
        "source_position_method":"0.5 + sky WCS slit_x/slit_y from pinned original jwst_1413.pmap projection at 16004 Angstrom; this is not measured PSF centering",
        "three_shutter_cols":list(range(target_c-1,target_c+2)),
        "shutter_row":target_r,"shutter_quadrant":where["quadrant"],
        "source_metadata_inserted_in_test_copy_only":True,
        "group_wcs_16004_original_slit_x":det["sky_registered_trace"]["16004A"]["slit_x"],
        "group_wcs_16004_original_slit_y":det["sky_registered_trace"]["16004A"]["slit_y"],
        "astrophysical_calibration":False,
        "output_directory":str(tmp)}

def test_pipeline(row):
    import jwst
    from jwst.pipeline import Spec2Pipeline
    from jwst.extract_1d import Extract1dStep
    assert jwst.__version__=="1.19.1",jwst.__version__
    os.environ["CRDS_CONTEXT"]=ORIGINAL_CONTEXT
    tmp=Path(row["output_directory"])
    data={"stage2_status":"NOT_ATTEMPTED","returned_types":[],
          "diagnostic_x1d_count":0,"source_flux_valid":False,
          "source_identifier_verified_by_real_MSA":False,
          "science_x1d":None,"calibrated_flux":None,
          "independent_line_detection":None,"redshift":None,"gas_trail":None}
    models=[]
    try:
        pipeline=Spec2Pipeline()
        pipeline.output_dir=str(tmp)
        pipeline.save_results=False
        for name in ("master_background_mos","bkg_subtract",
                     "pathloss","barshadow","resample_spec"):
            if hasattr(pipeline,name):getattr(pipeline,name).skip=True
        raw=pipeline.run(row["test_rate"])
        models=raw if isinstance(raw,(list,tuple)) else ([] if raw is None else [raw])
        data["stage2_status"]="CALL_RETURNED"
        data["returned_types"]=[type(m).__name__ for m in models]
        for i,model in enumerate(models):
            data.setdefault("stage2_results",[]).append({
                "model_type":type(model).__name__,
                "slit_count":len(getattr(model,"slits",[])),
                "meta_calibration":{
                    name:str(getattr(getattr(getattr(model,"meta",None),"cal_step",None),name,None))
                    for name in ("assign_wcs","extract_2d","srctype","wavecorr",
                                 "flat_field","photom","extract_1d")
                }
            })
            try:
                if getattr(model,"slits",None):
                    # Source metadata must be in slit for any diagnostic
                    # extraction; no physical source is thereby proven.
                    slits=list(model.slits)
                    data["stage2_results"][-1]["slit_source_ids"]=[
                        int(getattr(s,"source_id",-999)) for s in slits]
                    if SOURCE_ID not in data["stage2_results"][-1]["slit_source_ids"]:
                        data["stage2_results"][-1]["source_id_mismatch"]=True
                        continue
                # Trigger extract_1d directly even when pipeline skipped it.
                spec=Extract1dStep.call(model,save_results=False)
                has=len(getattr(spec,"spec",[])) if spec is not None else 0
                data["diagnostic_x1d_count"]+=has
                data.setdefault("manual_extract_1d_steps",[]).append({
                    "source_calibrated_false":True,
                    "returned_type":type(spec).__name__,
                    "number_of_spectra":has
                })
                if has:
                    output=tmp/f"JANUS_DIAGNOSTIC_{row['group']}_{i}_x1d.fits"
                    spec.save(output)
                    data.setdefault("diagnostic_fits",[]).append({
                        "file":output.name,"bytes":output.stat().st_size,
                        "sha256":sha(output),"source_claim":False})
                if hasattr(spec,"close"):spec.close()
            except Exception as exc:
                data.setdefault("manual_extraction_failures",[]).append({
                    "exception_class":type(exc).__name__,
                    "message":str(exc)[:900]})
        if not data["diagnostic_x1d_count"]:
            data["stage2_status"]="NO_DIAGNOSTIC_SOURCE_X1D"
    except Exception as exc:
        data["stage2_status"]="PIPELINE_ERROR"
        data["pipeline_error_class"]=type(exc).__name__
        data["pipeline_error"]=str(exc)[:1500]
        data["traceback_tail"]=traceback.format_exc()[-5000:]
    finally:
        for m in models:
            if hasattr(m,"close"):m.close()
    return data

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",default="/tmp/janus_mocss_source_aware_pilot")
    ap.add_argument("--output",default="/tmp/janus_mocss_source_aware_pilot.json")
    ap.add_argument("--indices",default="8,13")
    ap.add_argument("--geometry-only",action="store_true")
    a=ap.parse_args()
    root=Path(a.root);root.mkdir(parents=True,exist_ok=True)
    frozen=json.loads(RECEIPT.read_text())
    mapped={e["exposure_index"]:e for e in frozen["exposures"]}
    ids=[int(t) for t in a.indices.split(",")]
    assert ids and len(ids)<=2 and all(x in mapped for x in ids)
    receipt={
        "schema":"janus-mocss-source-aware-diagnostic-only-v1",
        "source_ra_dec_deg":[RA,DEC],
        "JWST_version":"1.19.1","CRDS_CONTEXT":ORIGINAL_CONTEXT,
        "method":"SYNTHETIC_SOURCE_INFO_COPIED_MSA_THREE_ADJACENT_SHUTTERS_NOT_REAL_MPT_SOURCE",
        "original_rate_or_msa_mutated":False,
        "group_A_to_D_independent_spectra":None,
        "scientifically_valid_x1d":False,
        "source_identity":None,"redshift":None,"physical_gas_trail":None,
        "results":[],
        "still_needed":[
            "All 16 original exposures in four registered stacks",
            "Original open-MSA full-row sky subtraction before stack",
            "Explicit wavelength-varying bar/transmission correction and physical PSF",
            "Independent PSF-shaped injected-source flux recovery with variance",
            "Cross-angle same-sky null controls and off-line trials"]}
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    for i in ids:
        r={"exposure_index":i}
        try:
            prep=prepare_source_aware(mapped[i],root)
            r["geometry"]=prep
            r["calibration"]=({"stage2_status":"SKIPPED_GEOMETRY_ONLY"}
                              if a.geometry_only else test_pipeline(prep))
        except Exception as exc:
            r["fatal_error"]=type(exc).__name__+": "+str(exc)
            r["traceback_tail"]=traceback.format_exc()[-5000:]
        receipt["results"].append(r)
        out.write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
        print("SOURCE_AWARE_TRIAL",i,
              r.get("calibration",{}).get("stage2_status","FAIL_BEFORE_CALIBRATION"),
              "spectra",r.get("calibration",{}).get("diagnostic_x1d_count",0),flush=True)
    failures=[r for r in receipt["results"] if "fatal_error" in r or
              r.get("calibration",{}).get("stage2_status")=="PIPELINE_ERROR"]
    if failures:
        raise SystemExit("FAILED source-aware test MSA geometry or Stage2: see JSON for full evidence")
    print("SCIENCE_GATE=NO_ASTROPHYSICALLY_CALIBRATED_X1D_FROM_ARTIFICIAL_SOURCE_METADATA",flush=True)

if __name__=="__main__":main()
