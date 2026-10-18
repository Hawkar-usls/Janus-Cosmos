#!/usr/bin/env python3
"""Non-authoritative spectroscopy audit of *actual* Stage2 diagnostic x1d FITS.

Reads FITS generated from original verified JWST RATE + artificially source-
labelled three-shutter MSA copies. Tests spectral coverage of three frozen
wavelength windows. Optional Gaussian/continuum coefficients use x1d-provided
per-bin errors but MUST NOT be called detections, physical fluxes or limits.
The source is absent from the original MPT; custom full-row sky subtraction,
four-dither stacks, bar/PSF throughput and correlated uncertainty are absent.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits

from fit_mocss_calibrated_three_lines import (
    CENTER_ANGSTROM, line_fit, read_x1d_fits, sha256_file
)

def audit(src, root):
    assert src["scientifically_valid_x1d"] is False
    assert src["group_A_to_D_independent_spectra"] is None
    assert src["original_rate_or_msa_mutated"] is False
    results=[]
    for e in src["results"]:
        group=e["geometry"]["group"]
        exp=e["exposure_index"]
        cal=e["calibration"]
        row={"group":group,"exposure_index":exp,
             "stage2_status":cal["stage2_status"],
             "diagnostic_x1d_count":cal["diagnostic_x1d_count"],
             "source_flux_valid":False,
             "lines":{}}
        for d in cal.get("diagnostic_fits",[]):
            file=Path(e["geometry"]["output_directory"])/d["file"]
            assert file.is_file(),str(file)
            assert sha256_file(file)==d["sha256"], "FITS SHA mismatch"
            detail={"basename":file.name,"sha256":d["sha256"],"nbytes":file.stat().st_size}
            with fits.open(file,memmap=False) as hdul:
                detail["table_extensions"]=[
                    {"name":h.name,"columns":list(h.columns.names),
                     "units":{c.name:str(c.unit) for c in h.columns}}
                    for h in hdul if isinstance(h,fits.BinTableHDU)]
            try:
                wave,flux,err=read_x1d_fits(file)
                usable=np.isfinite(wave)&np.isfinite(flux)&np.isfinite(err)&(err>0)
                detail["numeric_quality"]={
                    "n_bins":len(wave),
                    "valid_positive_error_bins":int(usable.sum()),
                    "min_lambda_angstrom":float(wave.min()) if len(wave) else None,
                    "max_lambda_angstrom":float(wave.max()) if len(wave) else None,
                    "median_lambda_step_angstrom":float(np.median(np.diff(wave))) if len(wave)>1 else None}
                for name,center in CENTER_ANGSTROM.items():
                    around=np.abs(wave-center)<=75
                    n=int(np.count_nonzero(usable&around))
                    entry={"wavelength_angstrom":center,
                           "valid_bins_within_plusminus_75A":n,
                           "science_valid":False,
                           "line_detected":None,
                           "calibrated_flux_or_limit":None}
                    fit=line_fit(wave,flux,err,center)
                    if fit["status"].startswith("X1D_LINEAR"):
                        entry["conditional_diagnostic_continuum_line_fit"]=fit
                    else:
                        entry["unfitted_reason"]=fit["status"]
                    row["lines"][name]=entry
            except (ValueError,KeyError,TypeError,IndexError,ZeroDivisionError) as ex:
                detail["strict_x1d_reader_error"]=type(ex).__name__+": "+str(ex)[:300]
                for name,center in CENTER_ANGSTROM.items():
                    row["lines"][name]={
                        "wavelength_angstrom":center,
                        "status":"UNVERIFIED_SCHEMA_OR_UNUSABLE_SPECTRUM",
                        "science_valid":False,"line_detected":None}
            row.setdefault("diagnostic_products",[]).append(detail)
        results.append(row)
    return {
        "schema":"janus-mocss-source-aware-x1d-three-window-diagnostic-v1",
        "science_status":"REAL_RATE_STAGE2_DIAGNOSTIC_SPECTRA_NOT_ASTROPHYSICAL_FLUX",
        "data_path":"Original MAST GO3290 NIRSpec RATE+MSA -> copied test MSA with invented SOURCE_INFO -> calibrated 2D and manual Extract1dStep",
        "astrophysical_source_associated_with_original_MPT":False,
        "photometric_slit_bar_throughput_validated":False,
        "custom_mocss_open_shutter_sky_background_validated":False,
        "all_16_rate_frames_extracted":False,
        "source_detected":None,"source_class":None,"redshift":None,"directional_gas_trail":None,
        "physical_line_fluxes_or_limits":None,
        "analysis_windows_angstrom":CENTER_ANGSTROM,
        "results":results,
        "guardrail":"Formal x1d Gaussian coefficients lack valid background, throughput and correlated-noise model; NEVER use their apparent flux/SNR as science significance."
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-receipt",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    rec=json.loads(args.source_receipt.read_text())
    result=audit(rec,args.source_receipt.parent)
    args.output.parent.mkdir(exist_ok=True,parents=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("X1D_FITS_SHA_AND_THREE_LINES_DIAGNOSTIC_AUDITED",flush=True)
    for g in result["results"]:
        print(g["group"],g["stage2_status"],g["diagnostic_x1d_count"],
              {k:v.get("valid_bins_within_plusminus_75A","SCHEMA_UNREADABLE")
               for k,v in g["lines"].items()},flush=True)
    print("NO_PHYSICAL_LINE_FLUX_OR_REDSHIFT_CLAIM",flush=True)

if __name__=="__main__":
    main()
