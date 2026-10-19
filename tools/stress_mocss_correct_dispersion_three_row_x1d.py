#!/usr/bin/env python3
"""Falsification-first diagnostic X1D off-line controls and spectral injections.

INPUTS: SHA-pinned real JWST/NIRSpec Stage2 x1d outputs (A/B/C/D) from
2026-10-18 GitHub Actions MOCSS source-aware trial, itself using temporary
SOURCE_INFO metadata not present in MPT. No valid MOCSS custom sky subtraction,
4-dither stack, PSF/bar throughput or independent noise covariance exists here.

OUTPUT: empirical nearby off-line responses plus mismatched 1D spectral
injection-recovery. This cannot measure telescope throughput nor calibrated
line flux/significance; no astrophysical claims are allowed.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np
from fit_mocss_calibrated_three_lines import CENTER_ANGSTROM, line_fit, read_x1d_fits

EXPECTED_SHA={
  "A":"79eb00de3556ee94ecb67901cf621479371a26467603fc1bfeb1251071a6bcbb",
  "B":"47430bc196e9562a21ab9617b7e79d01c47ba9b9fa0a801994251335b56db6d7",
  "C":"e18a6afa000e006985dcc8c6f582500d4ebff262badff0de17149fff382839a8",
  "D":"dd780aacbff407dfe0a03a209be6a464a51816904caca3cf487a3c88f4d93fec",
}
ORIGINAL_RUN_ID=39223778736
C_A_PER_S=2.99792458e18
LINE_FLUX_CHALLENGE={
    "16004A":2.5e-18,
    "OIII_4960":2.5e-18/3.,
    "Pa_gamma":2.5e-18/1.6,
}

def digest(path):
    sha=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):sha.update(block)
    return sha.hexdigest()

def safe_float(x):
    x=float(x)
    return x if math.isfinite(x) else None

def robust_mad(values):
    a=np.asarray(values,dtype=float)
    center=float(np.median(a))
    sigma=float(1.4826*np.median(np.abs(a-center)))
    return center,sigma

def matched_offline_controls(wave,fnu,error,center,all_centers):
    # Pre-register local controls at ±125-425 A every 25 A. Adjacent +/-75 A
    # fitting windows overlap; these are not independent statistical trials.
    offsets=np.arange(-425,426,25)
    excluded=[]; measured=[]
    for dx in offsets:
        target=center+float(dx)
        if abs(dx)<125:continue
        if any(abs(target-c)<125 for c in all_centers):
            excluded.append(target);continue
        q=line_fit(wave,fnu,error,target)
        if q["status"]!="X1D_LINEAR_TEMPLATE_MEASURED_ERROR_UNVALIDATED":
            continue
        measured.append({"center_A":target,"flux_like_coefficient_cgs":q["fitted_flux_cgs"],
                         "formal_fit_ratio_NOT_SIGNIFICANCE":q["formal_template_amplitude_over_error_NOT_DETECTION_SIGNIFICANCE"],
                         "reduced_chi2":q["chi2"]/q["dof"]})
    if len(measured)<12:
        return {"status":"INSUFFICIENT_OFFLINE_SKY_CONTROLS","n":len(measured)}
    med,mad=robust_mad([q["flux_like_coefficient_cgs"] for q in measured])
    formal=[abs(q["formal_fit_ratio_NOT_SIGNIFICANCE"]) for q in measured]
    chi2=[q["reduced_chi2"] for q in measured]
    return {
        "status":"EMPIRICAL_DIAGNOSTIC_CONTROLS_NOT_SCIENCE_NULL",
        "n":len(measured),
        "nearby_control_centers_A":[q["center_A"] for q in measured],
        "controls":measured,
        "control_flux_median_cgs":med,
        "control_flux_robust_mad_cgs":mad,
        "maximum_absolute_formal_template_ratio":max(formal),
        "fraction_control_abs_formal_ratio_ge_3":sum(v>=3 for v in formal)/len(formal),
        "reduced_chi2_median":float(np.median(chi2)),
        "overlapping_controls_not_independent":True,
        "same_sky_spectral_real_emission_or_artifacts_possible":True
    }

def inject_spectral_line(wave,fnu,error,nominal_center,flux_cgs,actual_fwhm_A,offset_A):
    # Spectral 1D injection is not a detector/PSF/slit-thruput injection.
    w=np.asarray(wave,dtype=float)
    sigma=actual_fwhm_A/(2*math.sqrt(2*math.log(2)))
    continuous=np.exp(-.5*((w-nominal_center-offset_A)/sigma)**2)/(math.sqrt(2*math.pi)*sigma)
    binwidth=np.gradient(w)
    # Discrete-area normalize to exactly the injected integrated cgs flux.
    area=float(np.sum(continuous*binwidth))
    if area<0.1:return {"status":"NO_PROFILE_COVERAGE"}
    flam=flux_cgs*continuous/area
    fnu_delta=flam*w*w/(1e-23*C_A_PER_S)
    before=line_fit(w,fnu,error,nominal_center)
    after=line_fit(w,fnu+fnu_delta,error,nominal_center)
    if any(q["status"]!="X1D_LINEAR_TEMPLATE_MEASURED_ERROR_UNVALIDATED" for q in [before,after]):
        return {"status":"FIT_INSUFFICIENT_COVERAGE"}
    recovered=float(after["fitted_flux_cgs"]-before["fitted_flux_cgs"])
    return {
      "status":"DETECTOR_THROUGHPUT_NOT_TESTED",
      "flux_injected_cgs":flux_cgs,
      "observed_fwhm_A":actual_fwhm_A,
      "centroid_offset_A":offset_A,
      "flux_recovered_difference_cgs":recovered,
      "response_fraction":recovered/flux_cgs,
      "underlying_real_x1d_unchanged":True,
      "matched_template_model_is_fwhm_13p2":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    files=list(args.root.rglob("JANUS_DIAGNOSTIC_*_0_x1d.fits"))
    groups={}
    for file in files:
        m=re.search(r"JANUS_DIAGNOSTIC_([ABCD])_0_x1d[.]fits$",file.name)
        if not m:continue
        g=m.group(1)
        if g in groups:raise RuntimeError(f"Duplicate group {g}")
        if digest(file)!=EXPECTED_SHA[g]:
            raise RuntimeError(f"WRONG_ORIGINAL_X1D_SHA {g} {file.name}")
        groups[g]=file
    if set(groups)!=set("ABCD"):
        raise RuntimeError(f"MISSING_REAL_X1D_GROUPS found={sorted(groups)}")
    d={
      "schema":"janus-mocss-correct-dispersion-three-row-x1d-local-null-mismatch-injection-v1",
      "status":"CORRECT_MSA_THREE_ROW_DIAGNOSTIC_X1D_CONTROLS_NOT_ASTROPHYSICAL_FLUX",
      "upstream_action_run_id":ORIGINAL_RUN_ID,
      "upstream_action":"https://github.com/Hawkar-usls/Janus-Cosmos/actions/runs/39223778736",
      "original_jwst_rate_MSA_data_tested":True,
      "test_source_info_invented":True,
      "all_16_exposures_stacked":False,
      "valid_mocss_row_background_subtraction":False,
      "valid_three_shutter_PSF_bar_throughput":False,
      "injection_stage":"DIAGNOSTIC_X1D_ONLY_NOT_2D_DETECTOR_OR_MSA",
      "line_flux_cgs_from_this_run":None,
      "line_detection_from_this_run":None,
      "source_redshift":None,"source_class":None,"gas_trail":None,
      "three_preregistered_windows_angstrom":CENTER_ANGSTROM,
      "spectral_injected_flux_challenge_cgs":LINE_FLUX_CHALLENGE,
      "controls_spacing_angstrom":25,
      "controls_search_offsets_angstrom":[-425,425],
      "controls_exclusion_distance_from_any_target_angstrom":125,
      "results":{}
    }
    centers=list(CENTER_ANGSTROM.values())
    for g in "ABCD":
        wave,flux,err=read_x1d_fits(groups[g])
        assert len(wave)>1000 and np.all(np.diff(wave)>0)
        row={"x1d_filename":groups[g].name,"x1d_sha256":EXPECTED_SHA[g],
             "n_bins":len(wave),"spectral_wavelength_minmax_A":[float(wave.min()),float(wave.max())],
             "line_windows":{}}
        for label,center in CENTER_ANGSTROM.items():
            fit=line_fit(wave,flux,err,center)
            assert fit["status"]=="X1D_LINEAR_TEMPLATE_MEASURED_ERROR_UNVALIDATED", (g,label,fit)
            reduced=fit["chi2"]/fit["dof"]
            control=matched_offline_controls(wave,flux,err,center,centers)
            row["line_windows"][label]={
              "lambda_A":center,
              "published_flux_model_NOT_THIS_DATA":LINE_FLUX_CHALLENGE[label],
              "diagnostic_unvalidated_coefficient_cgs":fit["fitted_flux_cgs"],
              "formal_error_IGNORES_NOISE_COVARIANCE_cgs":fit["formal_x1d_independent_bin_sigma_cgs"],
              "formal_fit_ratio_NOT_DETECTION_SIGNIFICANCE":fit["formal_template_amplitude_over_error_NOT_DETECTION_SIGNIFICANCE"],
              "reduced_chi2":reduced,
              "chi2_model_rejected_for_unit_variance_if_gt_3":bool(reduced>3),
              "off_line":control,
              "one_dimensional_injection": [
                  inject_spectral_line(wave,flux,err,center,LINE_FLUX_CHALLENGE[label],w,shift)
                  for w in (7.0,13.2,22.0) for shift in (-3.0,0.,3.0)
              ],
              "detected":None,"calibrated_flux_cgs":None,"calibrated_upper_limit_cgs":None
            }
            print(g,label,"chi2/dof",round(reduced,3),
                  "formal_R_NOT_SIGMA",round(fit["formal_template_amplitude_over_error_NOT_DETECTION_SIGNIFICANCE"],3),
                  "control_n",control.get("n"),"control_absRmax",round(control.get("maximum_absolute_formal_template_ratio",float("nan")),2),flush=True)
        d["results"][g]=row
    # Deliberately prevent apparent huge formal peak -> discovery.
    assert all(d["results"][g]["line_windows"][name]["detected"] is None
               for g in "ABCD" for name in CENTER_ANGSTROM)
    assert len(d["results"])==4
    args.output.parent.mkdir(exist_ok=True,parents=True)
    args.output.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    print("PASS: FOUR_REAL_STAGE2_X1D_FITS_SHA_VERIFIED; EMPIRICAL_NULL_AND_1D_INJECTION_STRESS_ONLY")
    print("ASTROPHYSICAL_CALIBRATION=NOT_OBTAINED; NEW_LINE_DETECTION=NOT_ESTABLISHED")

if __name__=="__main__":
    main()
