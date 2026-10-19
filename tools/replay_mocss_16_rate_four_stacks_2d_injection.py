#!/usr/bin/env python3
"""16 ORIGINAL JWST MOCSS RATE frames, 4 WCS-aligned line-specific stacks,
robust full-detector-row sky and 2D detector-pixel artificial line trials.

Every selected RATE and both MSA metatables are checked against earlier SHA
receipts. No synthetic x1d metadata is used in the INPUT sky/detector stacks.

CRITICAL LIMITS: sky uses robust whole-detector-row *proxy* and audits the
per-MSA-row open fraction, but cannot transfer that fraction to exact pixels
without full physical shutter-to-detector footprint and contamination mask;
it does NOT reproduce Huberty's custom illumination-aware full-row sky.
WCS maps only three point windows, not the full 2D dispersed spectrum. The
2D in-memory Gaussian injection is NOT STPSF-dispersed/throughput-calibrated.
Units DN/s; neither physical flux nor independent S/N is available.
"""
import argparse
import concurrent.futures as cf
import hashlib
import json
import math
import os
from pathlib import Path
from urllib.parse import quote
import numpy as np
from astropy.io import fits
from scipy.ndimage import map_coordinates

from fetch_jwst_mocss_mast_subset import read_and_store,BASE,ROOT_URI

ROOT=Path(__file__).resolve().parents[1]
WCS_RECEIPT=ROOT/"domains/space/MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json"
MSA={
  "jw03290001001_05_msa.fits":"fb6239b287d25c478d2e22737ff956850a9bc1e0f365500cfc5661abbfe17a44",
  "jw03290002001_01_msa.fits":"25a061a226521a8bbf2f30d3f6922a22fb13437800047f12c8511e354820c91e"
}
LINE_NAMES={"primary":"16004A","OIII_4960":"15850p8A","Pa_gamma":"13656p75A_Pagamma"}
HALF_X,HALF_Y=70,15
X,Y=np.meshgrid(np.arange(-HALF_X,HALF_X+1,dtype=float),
                np.arange(-HALF_Y,HALF_Y+1,dtype=float))
CAX=np.arange(-50,-19,4).tolist()+np.arange(20,51,4).tolist()

def digest(p):
    h=hashlib.sha256()
    with open(p,"rb") as fh:
        for chunk in iter(lambda:fh.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()

def fetch_expected(directory,name,sha):
    dest=directory/name;directory.mkdir(parents=True,exist_ok=True)
    read_and_store(BASE+quote(ROOT_URI+name,safe=""),dest)
    if digest(dest)!=sha:raise RuntimeError("ARCHIVE SHA MISMATCH "+name)
    return dest

def audit_msa_open_masks(metafiles,source):
    # SHUTTER_INFO contains all open shutters of each original MSA_METADATA_ID.
    # Comparing with 171 physical columns maps open-fraction per MSA row, NOT
    # pixel illumination. Distinguish the two physical arrays explicitly.
    records={}
    for index,entry in enumerate(source["exposures"]):
        name=("jw03290001001_05_msa.fits" if entry["orientation"]==1
              else "jw03290002001_01_msa.fits")
        file=metafiles[name]
        idvalue=entry["source_msametid_signed16"]
        with fits.open(file,memmap=False) as f:
            rows=f["SHUTTER_INFO"].data
            names=rows.dtype.names
            # Rows in these metatables are defined per original shutter IDs.
            key="msa_metadata_id"
            filtered=rows[rows[key]==idvalue]
            if not len(filtered):
                raise RuntimeError(f"NO_ORIGINAL_MSA_METADATA_ID {idvalue}")
            quadrants={}
            for q in range(1,5):
                sh=filtered[filtered["shutter_quadrant"]==q]
                counts=np.bincount(np.asarray(sh["shutter_row"],dtype=int),minlength=366)[1:366]
                # Some MSA rows may have duplicate shutters in metatable:
                # unique coordinates prevent double-counted opening.
                coords=set((int(v["shutter_row"]),int(v["shutter_column"])) for v in sh)
                unique_counts=np.bincount([row for row,col in coords],minlength=366)[1:366]
                if np.any(unique_counts>171):raise RuntimeError("MSA ROW >171 unique columns")
                quadrants[str(q)]={
                    "record_rows":len(sh),
                    "unique_open_shutter_count":len(coords),
                    "per_msa_row_fraction_open_minmax":[float(unique_counts.min()/171),float(unique_counts.max()/171)],
                    "median_fraction_open":float(np.median(unique_counts/171)),
                    "q5_q95_fraction_open":np.quantile(unique_counts/171,[.05,.95]).tolist(),
                    "valid_original_msa_identifier":int(idvalue)
                }
            records[str(entry["exposure_index"])]=quadrants
    return records

def robust_row_background(sci,valid,exclude_x):
    # Target windows are a negligible fraction of ~2048 pixel rows, but
    # exclude them explicitly and use sigma clipping on the remaining rows.
    allow=valid.copy()
    for xc in exclude_x:
        x0=max(0,int(xc)-65);x1=min(sci.shape[1],int(xc)+66)
        allow[:,x0:x1]=False
    good=np.where(allow,sci,np.nan)
    with np.errstate(all="ignore"):
        med=np.nanmedian(good,axis=1)
        mad=1.4826*np.nanmedian(np.abs(good-med[:,None]),axis=1)
        good=np.where(allow&(np.abs(sci-med[:,None])<4*np.maximum(mad[:,None],1e-8)),sci,np.nan)
        sky=np.nanmedian(good,axis=1)
    if np.mean(np.isfinite(sky))<.95:raise RuntimeError("Too many sky rows undefined")
    sky=np.where(np.isfinite(sky),sky,med)
    return sky, {"method":"FULL_DETECTOR_ROW_SIGMA_CLIPPED_MEDIAN_4MAD_PROXY",
                 "finite_fraction":float(np.mean(np.isfinite(sky))),
                 "median_sky_dn_per_s":float(np.median(sky)),
                 "median_row_mad_dn_per_s":float(np.median(mad[np.isfinite(mad)])),
                 "absolute_sky_peak_to_peak_dn_per_s":float(np.ptp(sky)),
                 "n_sky_columns_max":sci.shape[1],
                 "physical_msa_row_projection_applied":False}

def sample(image,good,x0,y0):
    coords=[Y+y0,X+x0]
    numerator=map_coordinates(np.where(good,image,0),coords,order=1,mode="constant",cval=0)
    denom=map_coordinates(good.astype(float),coords,order=1,mode="constant",cval=0)
    return np.divide(numerator,denom,out=np.full_like(numerator,np.nan),where=denom>.95),denom

def per_exposure(rec,data_dir):
    index=rec["exposure_index"];g=rec["group"];det=rec["expected_science_detector_from_WCS"]
    f=rec["detectors"][det]
    p=data_dir/f["rate_file"]
    if digest(p)!=f["sha256"]:raise RuntimeError("Second SHA failure on "+p.name)
    with fits.open(p,memmap=False) as h:
        sci=np.asarray(h["SCI"].data,dtype=float)
        dq=np.asarray(h["DQ"].data,dtype=np.uint32)
        good=np.isfinite(sci)&((dq&3)==0)
        if "ERR" not in h:raise RuntimeError("No original RATE ERR extension")
        err=np.asarray(h["ERR"].data,dtype=float)
        centers={name:f["sky_registered_trace"][key] for name,key in LINE_NAMES.items()}
        assert all(v["on_detector"] for v in centers.values())
        row,sky=robust_row_background(sci,good,[v["x_px"] for v in centers.values()])
        sci_sub=sci-row[:,None]
        vgood=good&np.isfinite(err)&(err>0)
        stamps={}
        for name,center in centers.items():
            im,valid=sample(sci_sub,vgood,center["x_px"],center["y_px"])
            variance,_=sample(err*err,vgood,center["x_px"],center["y_px"])
            # Bilinear interpolation needs square-weight variance, so
            # interpolated ERR^2 is a CONSERVATIVE indicative proxy only.
            if not np.all(np.isfinite(im[np.isfinite(variance)])):raise RuntimeError("Bad resampling")
            stamps[name]={
              "image_dn_per_s":im,
              "variance_proxy":variance,
              "valid_fraction":float(np.isfinite(im).mean()),
              "source_detector_x_y":[center["x_px"],center["y_px"]],
              "source_inside_one_selected_shutter":center["inside_shutter"]
            }
    return {"group":g,"exposure_index":index,"rate_sha256":f["sha256"],
            "filename":p.name,"detector":det,"sky":sky,"stamps":stamps}

def score(z,dx=0,dy=0,gauss_sx=1.2,gauss_sy=1.5):
    template=np.exp(-.5*((X-dx)/gauss_sx)**2-.5*((Y-dy)/gauss_sy)**2)
    finite=np.isfinite(z)&(template>.0001)
    return float(np.sum(z[finite]*template[finite])/np.sum(template[finite]**2))

def robust_mad(nums):
    arr=np.asarray(nums,dtype=float)
    median=float(np.median(arr))
    return median,float(1.4826*np.median(np.abs(arr-median)))

def coadd(stamps):
    assert len(stamps)==4
    images=np.stack([s["image_dn_per_s"] for s in stamps])
    vars=np.stack([s["variance_proxy"] for s in stamps])
    good=np.isfinite(images)&np.isfinite(vars)&(vars>0)
    invvar=np.where(good,1/np.maximum(vars,1e-9),0)
    den=invvar.sum(axis=0)
    co=np.divide(np.sum(np.where(good,images,0)*invvar,axis=0),den,
                 out=np.full_like(images[0],np.nan),where=den>0)
    # Empirical between-dither sample variance, no covariance of detector
    # interpolation nor wavelength-dependent MSA throughput.
    scatter=np.sqrt(np.nanmean(np.where(good,(images-co)**2,np.nan),axis=0))
    controls=[score(co,dx) for dx in CAX]
    median,mad=robust_mad(controls)
    model=score(co)
    return co,{"method":"WCS_POINT_WINDOW_FOUR_DITHER_INV_RATE_ERR2_PROXY",
               "n_input":len(stamps),"coverage_at_target":int(good[:,HALF_Y,HALF_X].sum()),
               "center_template_response_dn_per_s":model,
               "off_line_n":len(controls),
               "off_line_median_dn_per_s":median,
               "off_line_MAD_dn_per_s":mad,
               "center_minus_control_median_in_MAD_units_NOT_FORMAL_SIGMA":(model-median)/mad if mad>0 else None,
               "dither_scatter_rms_dn_per_s":float(np.nanmedian(scatter)),
               "recovered_physical_flux":None,
               "covariance_and_shutter_throughput_valid":False,
               "off_line_centers_dependently_sampled":True}

def inject_coadd(stamps,amplitude,shift_px,profile_x,profile_y=1.6):
    # Add onto 2D detector WCS patches BEFORE common resampling/coaddition;
    # source is PSF-proxy Gaussian, not a model output from STPSF.
    profile=np.exp(-.5*((X-shift_px)/profile_x)**2-.5*(Y/profile_y)**2)
    injected=[dict(s,image_dn_per_s=s["image_dn_per_s"]+amplitude*profile) for s in stamps]
    return coadd(injected)[1]["center_template_response_dn_per_s"]

def process(rows):
    groups={k:[] for k in "ABCD"}
    for r in rows:groups[r["group"]].append(r)
    report={}
    images={}
    for g in "ABCD":
        assert len(groups[g])==4
        current={"indices":[r["exposure_index"] for r in groups[g]],
                 "detector":groups[g][0]["detector"],"lines":{}}
        for key in LINE_NAMES:
            ls=[r["stamps"][key] for r in groups[g]]
            co,info=coadd(ls)
            rec=[]
            for amp in (.01,.04,.12):
                for sx in (0.8,1.2,2.0):
                    for shift in (-.7,0,.7):
                        recovered=inject_coadd(ls,amp,shift,sx)
                        rec.append({"peak_injected_dn_per_s":amp,
                                    "dispersion_sigma_pixels_TEST":sx,
                                    "wcs_x_centroid_offset_pixels_TEST":shift,
                                    "recovered_delta_template_dn_per_s":recovered-info["center_template_response_dn_per_s"],
                                    "recovery_ratio_not_PSF_throughput":(recovered-info["center_template_response_dn_per_s"])/amp})
            current["lines"][key]=dict(info,
                injected_detector_patch_2D_real_rate_challenge=rec,
                injected_physical_source_flux_cgs=None)
            images[g+"_"+key]=co.astype(np.float32)
            print("STACK",g,key,"dithers",info["n_input"],"MAD",round(info["off_line_MAD_dn_per_s"],6),flush=True)
        report[g]=current
    return report,images

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--data-dir",type=Path,default=Path("/tmp/janus_16rate_row_sky"))
    p.add_argument("--output",type=Path,default=Path("/tmp/janus_16rate_sky_fourstacks.json"))
    p.add_argument("--images",type=Path,default=Path("/tmp/janus_16rate_sky_fourstacks.npz"))
    p.add_argument("--workers",type=int,default=4)
    a=p.parse_args()
    r=json.loads(WCS_RECEIPT.read_text())
    assert len(r["exposures"])==16 and r["jwst_pipeline"]=="1.19.1"
    a.data_dir.mkdir(parents=True,exist_ok=True)
    tasks=[]
    for e in r["exposures"]:
        d=e["detectors"][e["expected_science_detector_from_WCS"]]
        tasks.append((d["rate_file"],d["sha256"]))
    tasks.extend(MSA.items())
    with cf.ThreadPoolExecutor(max_workers=a.workers) as pool:
        list(pool.map(lambda e:fetch_expected(a.data_dir,e[0],e[1]),tasks))
    filenames={name:a.data_dir/name for name in MSA}
    opened=audit_msa_open_masks(filenames,r)
    report={"schema":"janus-cosmos-mocss-16-rate-four-dither-row-background-proxy-v1",
     "status":"ALL_16_REAL_RATE_FITS_FOUR_WCS_DITHER_STACKS_AND_2D_PARAMETRIC_INJECTION_NOT_PHYSICAL_SOURCE_FLUX",
     "input_date":"2023-08-23_to_2023-11-19","created_utc":"2026-10-19",
     "source_ra_dec_deg":[r["target_ra_deg"],r["target_dec_deg"]],
     "wavelengths_A":{"primary":16004,"OIII_4960":15850.7932988,"Pa_gamma":13656.7483309},
     "original_rate_and_msa_sha256_verified":True,
     "original_inputs":[{"exposure_index":e["exposure_index"],
                          "group":e["group"],
                          "file":e["detectors"][e["expected_science_detector_from_WCS"]]["rate_file"],
                          "sha256":e["detectors"][e["expected_science_detector_from_WCS"]]["sha256"]}
                         for e in r["exposures"]],
     "original_msa_sha256":MSA,
     "MSA_metadata_open_shutter_row_histograms_NOT_PIXEL_MAPPED":opened,
     "background_method":"clipped median of each full detector row, central spectral windows masked; not physical shutter illumination model",
     "physical_shutter_row_illumination_to_detector_mapping_executed":False,
     "wcs_3_small_windows_not_full_dispersed_spectrum":True,
     "rate_ERR_interpolated_as_variance_proxy_not_exact_covariance":True,
     "group_report":{},
     "scientific_limits":{"authors_full_row_MOCSS_sky_reproduced":False,
                          "calwebb_spec2_on_all_four_coadds":False,
                          "STPSF_in_detector_space":False,
                          "absolute_line_flux":None,
                          "realistic_dispersed_PSF_shutter_injection":False,
                          "independent_detection_significance":None,
                          "line_identification":None,"redshift":None,"gas_trail":None}}
    # Save early immutable source provenance even if a later stage fails.
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    frames=[]
    for e in r["exposures"]:
        f=per_exposure(e,a.data_dir)
        frames.append(f)
        report.setdefault("exposure_sky",[]).append({"index":f["exposure_index"],
                   "group":f["group"],"rate_sha256":f["rate_sha256"],
                   "sky":f["sky"]})
        a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
        print("FRAME",e["exposure_index"],e["group"],"sky",round(f["sky"]["median_sky_dn_per_s"],5),flush=True)
    gr,ims=process(frames)
    report["group_report"]=gr
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    a.images.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(a.images,**ims)
    report["output_2D_rate_stacks_npz_sha256"]=digest(a.images)
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("PASS: 16 RAW JWST RATE SHA, 4 empirical line window stacks with row-median SKY PROXY, 2D GAUSSIAN RATE INJECTIONS")
    print("ABSOLUTE_PHYSICAL_FLUX=NOT_MEASURED; STPSF_DETECTOR_INJECTION=NOT_RUN")

if __name__=="__main__":main()
