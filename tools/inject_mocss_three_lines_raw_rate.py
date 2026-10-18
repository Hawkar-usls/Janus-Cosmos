#!/usr/bin/env python3
"""Real JWST/NIRSpec RATE-pixel three-window artificial-line recovery.

The injected PSF is an intentionally *parametric detector-space Gaussian*, not
STPSF, flat/photom calibration, sky-line completeness or MOCSS throughput.
Outputs are injection-response in raw DN/s ONLY, not flux, SNR, redshift, or a
science-grade nondetection limit. Four widely separated exposures A/B/C/D
are used to check basic detector-angle dependence, not a four-stack spectrum.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits
from fetch_jwst_mocss_mast_subset import read_and_store, BASE, ROOT_URI

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"domains/space/MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json"
LABELS=("16004A","15850p8A","13656p75A_Pagamma")
INDICES=(1,8,10,13)
AMP=(0.01,0.04,0.12)
SX,SY=1.15,1.30

def file_sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):h.update(part)
    return h.hexdigest()

def get(path,sha):
    from urllib.parse import quote
    path.parent.mkdir(parents=True,exist_ok=True)
    read_and_store(BASE+quote(ROOT_URI+path.name,safe=""),path)
    actual=file_sha(path)
    if actual != sha:raise ValueError(f"SHA changed: {path.name}")
    return path

def score_patch(patch,good,x_axis,y_axis,xc,yc):
    rows,cols=patch.shape
    xx=x_axis[None,:]-xc
    yy=y_axis[:,None]-yc
    g=np.exp(-.5*(xx/SX)**2-.5*(yy/SY)**2)
    # Local row-wise off-line pedestal. Mask sources and invalid pixels.
    sky=(np.abs(xx)>=18)&(np.abs(xx)<=45)&good
    pedestal=np.full(rows,np.nan)
    for k in range(rows):
        z=patch[k][sky[k]]
        if len(z)>=16:pedestal[k]=np.median(z)
    target=good&np.isfinite(pedestal[:,None])&(g>=.0001)
    if np.sum(target)<15:raise ValueError("Insufficient valid target pixels")
    residual=patch-pedestal[:,None]
    amp=float(np.sum((g*residual)[target])/np.sum((g*g)[target]))
    return amp

def measure_one(img,dq,cx,cy):
    ix=int(round(cx));iy=int(round(cy))
    halfw,halfh=105,12
    if ix-halfw<0 or ix+halfw>=img.shape[1] or iy-halfh<0 or iy+halfh>=img.shape[0]:
        return {"status":"CUTOUT_OFF_DETECTOR"}
    patch=np.asarray(img[iy-halfh:iy+halfh+1,ix-halfw:ix+halfw+1],dtype=float)
    dqpatch=np.asarray(dq[iy-halfh:iy+halfh+1,ix-halfw:ix+halfw+1],dtype=np.uint32)
    good=np.isfinite(patch)&((dqpatch&3)==0)
    x=np.arange(ix-halfw,ix+halfw+1,dtype=float)
    y=np.arange(iy-halfh,iy+halfh+1,dtype=float)
    if np.mean(good)<0.75:
        return {"status":"TOO_MANY_INVALID_PIXELS","valid_fraction":float(np.mean(good))}
    pre=score_patch(patch,good,x,y,cx,cy)
    # Source-free controls: center moved, unchanged data and same template.
    control_shifts=list(range(-78,-21,3))+list(range(24,79,3))
    null=np.array([score_patch(patch,good,x,y,cx+v,cy) for v in control_shifts])
    median=float(np.median(null))
    mad=float(1.4826*np.median(np.abs(null-median)))
    # Injection added *in-memory* to original measured detector cutout.
    g=np.exp(-.5*((x[None,:]-cx)/SX)**2-.5*((y[:,None]-cy)/SY)**2)
    inject=[]
    for amp in AMP:
        after=score_patch(patch+amp*g,good,x,y,cx,cy)
        inject.append({"injected_peak_dn_per_s":amp,
                       "recovered_change_peak_dn_per_s":float(after-pre),
                       "recovery_ratio_not_astrophysical_throughput":float((after-pre)/amp)})
    return {"status":"PARAMETRIC_GAUSSIAN_RAW_RATE_INJECTION_ONLY",
            "source_center_px":[cx,cy],
            "template_sigma_detector_px":[SX,SY],
            "original_rate_template_score_dn_per_s":pre,
            "off_line_control_n":len(null),
            "off_line_control_median_dn_per_s":median,
            "off_line_control_robust_mad_dn_per_s":mad,
            "control_response_raw_dn_s":null.tolist(),
            "injections":inject,
            "valid_pixel_fraction":float(np.mean(good))}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--data-dir",default="/tmp/janus_mocss_three_line_injection")
    parser.add_argument("--output",default="/tmp/janus_mocss_three_line_rate_injection.json")
    args=parser.parse_args()
    x=json.loads(SOURCE.read_text())
    assert x["jwst_pipeline"]=="1.19.1" and x["crds"]=="jwst_1413.pmap"
    entries={v["exposure_index"]:v for v in x["exposures"]}
    records=[]
    for ix in INDICES:
        v=entries[ix]; detector=v["expected_science_detector_from_WCS"]
        f=v["detectors"][detector]
        filename=f["rate_file"]
        p=get(Path(args.data_dir)/filename,f["sha256"])
        with fits.open(p,memmap=True) as hdul:
            sci=hdul["SCI"].data;dq=hdul["DQ"].data
            lines={}
            for label in LABELS:
                coord=f["sky_registered_trace"][label]
                assert coord["on_detector"] and coord["inside_shutter"]
                lines[label]=measure_one(sci,dq,coord["x_px"],coord["y_px"])
        records.append({"exposure_index":ix,"group":v["group"],
                        "detector":detector,
                        "original_rate_file":filename,"sha256":file_sha(p),
                        "wcs_version":"JWST 1.19.1 / jwst_1413.pmap",
                        "science_claim":False,"lines":lines})
        print("INJECTION_REAL_RATE",ix,v["group"],[(k,round(z.get("original_rate_template_score_dn_per_s",0),6),
                                                [round(j["recovery_ratio_not_astrophysical_throughput"],3) for j in z.get("injections",[])])
                                               for k,z in lines.items()],flush=True)
    output={"schema":"janus-mocss-three-line-rate-injection-v1",
            "result":"REAL_RATE_PIXELS_TESTED_GAUSSIAN_METHOD_NOT_SCIENCE_CALIBRATION",
            "date_utc":"2026-10-18",
            "original_rate_file_count":len(records),
            "wcs_geometry_receipt":str(SOURCE.relative_to(ROOT)),
            "source_ra_dec_deg":[53.13127,-27.78908],
            "labels":LABELS,"injected_peak_dn_per_s":AMP,
            "data_original_immutable":True,
            "records":records,
            "controls":"off-line detector x positions plus 3 detector-peak injection levels at frozen WCS predicted line centers",
            "no_claims":{
                "official_stage2_completed":False,
                "spectrophotometric_line_flux_measured":False,
                "line_detection_significance":None,
                "three_shutter_psf_bar_throughput_calibrated":False,
                "line_identity":None,
                "source_redshift":None,
                "source_class":None,
                "gas_trail":None},
            "remaining":"Use CRDS-pinned MOCSS custom stage2, blank-sky covariance, STPSF source-centered bar throughput, frame-level sky, PSF-injected physical line-flux recovery."}
    dest=Path(args.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
    print("PASS: real detector-pixel parametric injection; ABSOLUTE FLUX UNCALIBRATED",flush=True)

if __name__=="__main__":main()
