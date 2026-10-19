#!/usr/bin/env python3
"""Monochromatic NIRSpec three-adjacent-MSA-shutter optical PSF stress experiment.

STPSF supplies MONOCHROMATIC IMAGING-PLANE PSFs, NOT dispersed detector PSFs.
This compares 'Three adjacent MSA open shutters' with 'MSA all open' under
identical STPSF configurations at THREE FROZEN WAVELENGTHS and two source
position offsets. Ratios are *model optical flux fractions within common
finite FOV*, never a calibrated point-source slit-throughput for MOCSS.
Physical detector-space injection requires mapping through WCS/dispersion,
flat/photom, and independently measured source position and source flux.
"""
import argparse, hashlib, json, math, os
from pathlib import Path
import numpy as np

WAVES_A={"Pa_gamma":13656.7483309,"OIII_4960":15850.7932988,"primary":16004.0}
MASKS=["MSA all open","Three adjacent MSA open shutters"]
# These offsets are *sensitivity probes*, not confirmed position-in-shutter.
OFFSETS_ARCSEC={"center":(0.,0.),"offcenter_stress":(0.05,0.12)}

def digest_array(a):
    return hashlib.sha256(np.ascontiguousarray(a,dtype="<f8").tobytes()).hexdigest()

def measure_one(inst, wave_a, offset, oversample, fov):
    inst.options["source_offset_x"]=offset[0]
    inst.options["source_offset_y"]=offset[1]
    vals={}
    for mask in MASKS:
        inst.image_mask=mask
        psf=inst.calc_psf(monochromatic=wave_a*1e-10,
                          oversample=oversample, fov_arcsec=fov)
        try:
            oversampled=np.asarray(psf["OVERSAMP"].data,dtype=float)
            assert np.isfinite(oversampled).all()
            assert np.min(oversampled)>-1e-10 and np.sum(oversampled)>0
            vals[mask]={"finite_fov_image_sum":float(np.sum(oversampled)),
                         "array_sha256":digest_array(oversampled),
                         "image_shape":list(oversampled.shape),
                         "header_pixscale":str(psf["OVERSAMP"].header.get("PIXELSCL","")),
                         "oversamp_primary_header":str(psf[0].header.get("OVERSAMP",""))}
        finally:
            psf.close()
    ratio=vals[MASKS[1]]["finite_fov_image_sum"]/vals[MASKS[0]]["finite_fov_image_sum"]
    assert 0 < ratio < 1.2, ratio
    return {"offset_source_arcsec":offset,
            "per_mask":vals,
            "finite_fov_three_shutter_to_all_open_sum_ratio":ratio,
            "finite_fov_three_shutter_loss_fraction":1-ratio}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("/tmp/janus_stpsf_three_wavelength_shutter_model.json"))
    p.add_argument("--oversample",type=int,default=4)
    p.add_argument("--fov",type=float,default=1.8)
    opt=p.parse_args()
    assert opt.oversample>=3 and opt.fov>=1.5
    import stpsf
    m=stpsf.NIRSpec()
    m.filter="F100LP"
    m.pupil_mask="NIRSpec grating"
    m.options["jitter"]="gaussian"
    m.options["jitter_sigma"]=0.007 # TEST choice, not measured for these observations
    receipt={"schema":"janus-nirspec-three-wavelength-stpsf-optical-mask-comparison-v1",
      "method":"STPSF_MONOCHROMATIC_OPTICAL_IMAGE_PSF_RATIO_NOT_DISPERSED_NIRSPEC_MOCSS_THROUGHPUT",
      "date_utc":"2026-10-19","stpsf_version":getattr(stpsf,"__version__","unavailable"),
      "config":{"filter":"F100LP","pupil_mask":"NIRSpec grating",
                "oversample":opt.oversample,"fov_arcsec":opt.fov,
                "jitter_sigma_arcsec_test_assumed":.007},
      "reference_data_dir":os.environ.get("STPSF_PATH","automatic_default"),
      "masks":MASKS,"results":{},"scientific_limits":{
        "true_four_dither_MOCSS_throughput_calibrated":False,
        "detector_2D_dispersed_psf_injection_done":False,
        "actual_source_centroid_per_dither_measured":False,
        "absolute_line_flux_calibration":False,
        "fov_comparison_ratio_is_approximation_only":True,
        "source_identity_or_redshift":None}}
    opt.output.parent.mkdir(parents=True,exist_ok=True)
    for name,wave in WAVES_A.items():
        values={}
        for label,offset in OFFSETS_ARCSEC.items():
            print("STPSF",name,wave,label,flush=True)
            values[label]=measure_one(m,wave,offset,opt.oversample,opt.fov)
            receipt["results"][name]={"wavelength_angstrom":wave,"optical_mask_results":values}
            opt.output.write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    print("SUCCESS: STPSF 3-lambda monochromatic optical mask ratio, not calibrated spectrograph throughput",flush=True)

if __name__=="__main__":main()
