#!/usr/bin/env python3
"""Three pre-registered line fits to independently calibrated JWST x1d spectra.

Requires per-group *real* fully calibrated 1D FITS and explicit provenance.
No RATE pixel contrast or author-reported flux may be supplied instead.
Statistical errors are conditional on x1d ERROR and independent pixels;
covariance, shutter transmission and injection/recovery MUST be validated
before accepting a line detection or limit as science-grade.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

CENTER_ANGSTROM = {
    "16004A":16004.,
    "OIII_4960":16004.*4960.295/5008.239,
    "Pa_gamma":16004.*10941.1/12821.6,
}
C_A_PER_S = 2.99792458e18
JY_CGS = 1e-23

def sha256_file(path):
    h=hashlib.sha256()
    with open(path,"rb") as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()

def line_fit(wavelength_A, fnu_jy, err_jy, center_A, observed_fwhm_A=13.2,
             window_halfwidth_A=75.0):
    wave=np.asarray(wavelength_A,dtype=float)
    fnu=np.asarray(fnu_jy,dtype=float)
    efnu=np.asarray(err_jy,dtype=float)
    mask=(np.isfinite(wave)&np.isfinite(fnu)&np.isfinite(efnu)&(efnu>0)
          &(np.abs(wave-center_A)<=window_halfwidth_A))
    wave,fnu,efnu=wave[mask],fnu[mask],efnu[mask]
    if len(wave)<12:
        return {"status":"INSUFFICIENT_WAVELENGTH_COVERAGE","n_bins":len(wave)}
    # WAVELENGTH measured in Angstrom. fnu Jy converts exactly to
    # erg/(s cm2 Angstrom) as Jy*1e-23*c_A_per_s / lambda_A**2.
    conversion=JY_CGS*C_A_PER_S/wave**2
    flam=fnu*conversion
    eflam=efnu*conversion
    fwhm=observed_fwhm_A*(center_A/16004.)
    profile_sigma=fwhm/(2*math.sqrt(2*math.log(2)))
    dx=wave-center_A
    profile=np.exp(-.5*(dx/profile_sigma)**2)/(math.sqrt(2*math.pi)*profile_sigma)
    # Normalize amplitude to units of 10^-18 cgs line flux; orthogonalize
    # continuum position/slope to reduce numerical conditioning issues.
    response=flam/1e-20
    errors=eflam/1e-20
    X=np.column_stack([np.ones_like(dx),dx/75,100*profile])
    weight=1/errors
    coeff,_,rank,_=np.linalg.lstsq(X*weight[:,None],response*weight,rcond=None)
    if rank<3: return {"status":"SINGULAR_LINE_CONTINUUM_FIT"}
    xtwx=(X*weight[:,None]).T@(X*weight[:,None])
    cov=np.linalg.inv(xtwx)
    # coefficient[2] * (1e-20 / 100) = 1e-22 cgs? NO: since
    # 100*profile * coeff[2] * 1e-20 gives integrated flux
    # coeff[2]*100*1e-20 = coeff[2]*1e-18 cgs.
    flux=coeff[2]*1e-18
    err=math.sqrt(cov[2,2])*1e-18
    resid=response-X@coeff
    chi2=float(np.sum((resid/errors)**2))
    return {
        "status":"X1D_LINEAR_TEMPLATE_MEASURED_ERROR_UNVALIDATED",
        "n_bins":len(wave),
        "fitted_flux_cgs":float(flux),
        "formal_x1d_independent_bin_sigma_cgs":float(err),
        "formal_template_amplitude_over_error_NOT_DETECTION_SIGNIFICANCE":float(flux/err),
        "observed_profile_fwhm_A_assumed":float(fwhm),
        "chi2":chi2,
        "dof":len(wave)-3,
        "window_A":[float(center_A-window_halfwidth_A),float(center_A+window_halfwidth_A)]
    }

def read_x1d_fits(path):
    from astropy.io import fits
    with fits.open(path,memmap=False) as hdul:
        tables=[h for h in hdul if isinstance(h,fits.BinTableHDU) and
                {"WAVELENGTH","FLUX","FLUX_ERROR"}.issubset(set(h.columns.names))]
        if len(tables)!=1:
            raise ValueError("Require exactly one x1d table with WAVELENGTH, FLUX, FLUX_ERROR; no ambiguous slit selection")
        tab=tables[0]
        units={n:str(tab.columns[n].unit).replace(" ","").lower() for n in ("WAVELENGTH","FLUX","FLUX_ERROR")}
        if units["WAVELENGTH"] not in ("um","micron","microns") or units["FLUX"]!="jy" or units["FLUX_ERROR"]!="jy":
            raise ValueError(f"Explicit FITS column units must be um/Jy/Jy, received {units}")
        d=tab.data
        wave=np.asarray(d["WAVELENGTH"],dtype=float)*1e4
        flux=np.asarray(d["FLUX"],dtype=float)
        err=np.asarray(d["FLUX_ERROR"],dtype=float)
        if "DQ" in tab.columns.names:
            keep=(np.asarray(d["DQ"])==0)
            wave,flux,err=wave[keep],flux[keep],err[keep]
        if np.any(np.diff(wave)<=0):
            raise ValueError("Wavelength grid must be strictly increasing")
        return wave,flux,err

def selftest():
    wave=np.arange(13500,16100,3.0)
    continuum=1e-20
    known_flux=2.5e-18
    sigma=13.2/(2*math.sqrt(2*math.log(2)))
    flam=continuum+known_flux*np.exp(-.5*((wave-16004)/sigma)**2)/(sigma*math.sqrt(2*math.pi))
    conversion=JY_CGS*C_A_PER_S/wave**2
    fnu=flam/conversion
    e=1e-20*np.ones(len(wave))/conversion
    q=line_fit(wave,fnu,e,16004.)
    assert q["status"]=="X1D_LINEAR_TEMPLATE_MEASURED_ERROR_UNVALIDATED"
    assert abs(q["fitted_flux_cgs"]/known_flux-1)<1e-10,q
    assert set(CENTER_ANGSTROM)=={"16004A","OIII_4960","Pa_gamma"}
    print("PASS: synthetic unit-conversion and Gaussian area only; NO REAL JWST X1D INPUT")
    print("NIRSPEC_INDEPENDENT_STAGE2_LINE_FLUX=NOT_PERFORMED")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--manifest",type=Path,help="User-provided calibrated x1d A/B/C/D manifest (required outside --selftest)")
    ap.add_argument("--output",type=Path,default=Path("/tmp/janus_mocss_three_line_x1d.json"))
    ap.add_argument("--selftest",action="store_true")
    args=ap.parse_args()
    if args.selftest:
        selftest();return
    if args.manifest is None:
        ap.error("real calibrated FITS and manifest are REQUIRED; cannot infer from RATE contrasts")
    man=json.loads(args.manifest.read_text())
    assert man["source_ra_dec_deg"]==[53.13127,-27.78908]
    assert man["pipeline_version"]=="1.19.1" and man["crds_context"]=="jwst_1413.pmap"
    required=["custom_open_msa_full_row_background_verified",
              "detector_flat_sflat_photom_completed",
              "throughput_geometry_validated",
              "x1d_error_and_dq_arrays_checked"]
    if any(man.get(field) is not True for field in required):
        raise RuntimeError("Physical calibration prerequisites incomplete: "+repr(required))
    groups=man["groups"]
    assert set(groups)==set("ABCD")
    result={
       "status":"FITS_FLUX_FIT_CONDITIONAL_ERROR_NOT_INJECTION_CALIBRATED",
       "source_ra_dec_deg":man["source_ra_dec_deg"],
       "calibration_provenance":{k:man[k] for k in ["pipeline_version","crds_context"]+required},
       "groups":{},
       "line_wavelengths_angstrom":CENTER_ANGSTROM,
       "not_performed":["independent_pixel_noise_covariance","line_injection_recovery",
                        "spectrograph_slitloss_stress_test","source_redshift_identification"],
       "warning":"Do not call fitted formal x1d error a validated significance or companion 3sigma limit."
    }
    for group in "ABCD":
        item=groups[group]
        p=Path(item["x1d_fits"])
        if sha256_file(p)!=item["x1d_sha256"]:
            raise ValueError("x1d file digest mismatch: "+str(p))
        wave,fnu,err=read_x1d_fits(p)
        result["groups"][group]={"source_file":p.name,"source_sha256":item["x1d_sha256"],
            "lines":{name:line_fit(wave,fnu,err,cen) for name,cen in CENTER_ANGSTROM.items()}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("X1D_DIAGNOSTIC_MEASUREMENTS_NOT_SCIENCE_GRADE",str(args.output))

if __name__=="__main__":
    main()
