#!/usr/bin/env python3
"""Empirical *matched-sky aperture* covariance of October 2023 F150W mosaics.

Uses four original public JWST NIRCam i2d SCI FITS byte ranges and their
committed SHA-256 receipts. It does not use a PSF model. Cross-image
correlation at fixed sky positions includes astrophysical confusion and
shared-resampling structure: NOT a pure detector-noise covariance estimator.
Never convert the output to source class, redshift or calibrated line limits.
"""
import argparse
import hashlib
import json
import math
import statistics
from pathlib import Path

import numpy as np
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from recheck_nircam_time_domain import fetch_byte_range, parse_headers
from replay_jwst_mocss_mast_subset import RA, DEC, aperture_photometry

REPO = Path(__file__).resolve().parents[1]
RECEIPT = REPO / "domains/space/NIRCAM_F150W_TIME_DOMAIN_RECEIPT_2026-10-16.json"
NEIGHBORS = [(53.131165,-27.788930),(53.131584,-27.789105),(53.131301,-27.789367)]

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(1024*1024), b""):
            h.update(data)
    return h.hexdigest()

def dist(e1,n1,e2,n2):
    return math.hypot(e1-e2,n1-n2)

def to_en(ra,dec):
    return (ra-RA)*3600*math.cos(math.radians(DEC)),(dec-DEC)*3600

def materialize(row, root):
    filename = row["product"]
    prefix = filename.removesuffix("_i2d.fits")
    head = root / (prefix+"_header.bin")
    data = root / (prefix+"_sci_"+str(row["sci_subset_rows"])+"rows.bin")
    fetch_byte_range(filename,0,262143,head)
    lo,hi = row["sci_byte_range_inclusive"]
    fetch_byte_range(filename,lo,hi,data)
    assert digest(head) == row["header_sha256"], "MAST header changed from pinned receipt: "+filename
    assert digest(data) == row["sci_subset_sha256"], "MAST pixels changed from pinned receipt: "+filename
    primary,sci,offset = parse_headers(head.read_bytes())
    assert sci["BITPIX"] == -32 and sci["EXTNAME"]=="SCI" and sci.get("BUNIT","").replace(" ","").lower()=="mjy/sr"
    width = sci["NAXIS1"]
    assert (lo-offset)%(4*width)==0 and hi-lo+1==4*width*row["sci_subset_rows"]
    y0 = (lo-offset)//(4*width)
    image = np.memmap(data,dtype=">f4",mode="r",shape=(row["sci_subset_rows"],width))
    wcs = WCS(sci)
    target=aperture_photometry(image,None,wcs,RA,DEC,y_origin=y0)
    assert target is not None
    assert abs(target["aperture_njy"]-row["flux_aperture_njy"]) < 6, "Original row includes control-median subtraction"
    return {"name":filename,"image":image,"wcs":wcs,"y0":y0,"target":target,
            "hash":row["sci_subset_sha256"],"date":row["date_time_utc"][:10]}

def candidate_points(radius,step):
    extent = int(radius/step)
    points = []
    for i in range(-extent,extent+1):
        for j in range(-extent,extent+1):
            e,n=i*step,j*step
            if math.hypot(e,n) > radius or math.hypot(e,n) < 1.6: continue
            if any(dist(e,n,*to_en(ra,dec)) < .9 for ra,dec in NEIGHBORS): continue
            points.append((e,n))
    return points

def robust_clip_matrix(a,limit=3.5):
    # Drop a *position* across all four images if any epoch exceeds the robust
    # per-image envelope, conservatively masking bright common sky objects.
    center = np.median(a,axis=0)
    mad = 1.4826*np.median(np.abs(a-center),axis=0)
    if np.any(mad <= 0): raise RuntimeError("zero MAD in sky control set")
    good = np.all(np.abs(a-center) <= limit*mad,axis=1)
    return good,center,mad

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--data-dir",default="/tmp/janus_f150w_time_domain")
    ap.add_argument("--output",default="/tmp/janus_f150w_matched_covariance.json")
    ap.add_argument("--grid-radius",type=float,default=5.4)
    ap.add_argument("--grid-spacing",type=float,default=0.9)
    ap.add_argument("--min-controls",type=int,default=35)
    args=ap.parse_args()
    assert 3.0<=args.grid_radius<=6.0 and 0.8<=args.grid_spacing<=1.3
    root=Path(args.data_dir);root.mkdir(parents=True,exist_ok=True)
    saved=json.loads(RECEIPT.read_text())
    rows=[r for r in saved["results"] if r.get("coverage") and
          r.get("date_time_utc","").startswith("2023-10")]
    rows.sort(key=lambda r:r["date_time_utc"])
    assert [r["date_time_utc"][:10] for r in rows]==["2023-10-01","2023-10-02","2023-10-03","2023-10-09"]
    print("INPUTS",len(rows),"pinned original SCI byte ranges",flush=True)
    frames=[materialize(r,root) for r in rows]
    coords=candidate_points(args.grid_radius,args.grid_spacing)
    all_samples=[]
    use_coords=[]
    for east,north in coords:
        ra=RA+east/(3600*math.cos(math.radians(DEC)))
        dec=DEC+north/3600
        ps=[aperture_photometry(f["image"],None,f["wcs"],ra,dec,
                               y_origin=f["y0"]) for f in frames]
        if any(z is None for z in ps): continue
        vals=np.array([z["aperture_njy"] for z in ps],dtype=float)
        if np.all(np.isfinite(vals)):
            all_samples.append(vals)
            use_coords.append([round(east,5),round(north,5)])
    if len(all_samples)<args.min_controls:
        raise RuntimeError(f"Insufficient co-covered control points {len(all_samples)}")
    v=np.asarray(all_samples)
    keep,_,_=robust_clip_matrix(v)
    clean=v[keep]
    if len(clean)<args.min_controls:
        raise RuntimeError(f"Fewer than {args.min_controls} controls survive bright-source mask: {len(clean)}")
    # Center with clipped-position medians. Covariance is empirical over
    # matched sky controls and contains stationary astrophysical variance.
    control_medians=np.median(clean,axis=0)
    cov=np.cov(clean,rowvar=False,ddof=1)
    sig=np.sqrt(np.diag(cov))
    correlation=cov/np.outer(sig,sig)
    assert np.min(np.linalg.eigvalsh(cov))>=-1e-8 and np.all(sig>0)
    original=np.array([r["flux_aperture_njy"] for r in rows])
    original_sd=np.array([r["control_noise_1sigma_njy"] for r in rows])
    # Retain original fixed positive inverse-variance weights, which were
    # frozen before the covariance probe. Reject unknown correlations.
    weights=(1/original_sd**2);weights=weights/sum(weights)
    mu=float(weights@original)
    empirical_sigma=float(np.sqrt(weights@cov@weights))
    nominal_sigma=float(np.sqrt(np.sum((weights*original_sd)**2)))
    # GLS is diagnostic only, and its optimised weight may be negative in
    # highly correlated data; report fixed weights as primary.
    out={
        "schema":"janus-cosmos-october-nircam-matched-sky-covariance-v1",
        "status":"MEASURED_MATCHED_APERTURE_COVARIANCE_NOT_PURE_INSTRUMENT_NOISE",
        "source_ra_dec_deg":[RA,DEC],
        "datasets":[{"product":f["name"],"date":f["date"],"sci_subset_sha256":f["hash"]}
                    for f in frames],
        "method":{"aperture_arcsec":0.12,"annulus_arcsec":[0.24,0.38],
                  "grid_radius_arcsec":args.grid_radius,
                  "grid_spacing_arcsec":args.grid_spacing,
                  "native_aperture_background_method":"median annulus",
                  "target_exclusion_radius_arcsec":1.6,
                  "neighbor_exclusion_radius_arcsec":0.9,
                  "joint_clip_threshold_mad_units":3.5,
                  "candidate_controls":len(coords),
                  "all_coverage_controls":len(v),
                  "surviving_controls":len(clean),
                  "not_pure_noise_warning":"Same-sky static faint sources and correlated drizzle/interpolation may induce covariance; no source-free segmentation or independent calibration of systematic errors."},
        "control_medians_njy":control_medians.tolist(),
        "matched_aperture_empirical_sd_njy":sig.tolist(),
        "matched_aperture_covariance_njy_squared":cov.tolist(),
        "matched_aperture_correlation":correlation.tolist(),
        "fixed_original_weights":weights.tolist(),
        "target_four_epoch_old_receipt_njy":original.tolist(),
        "old_independent_mean_njy":mu,
        "old_independent_assumption_sigma_njy":nominal_sigma,
        "same_weights_empirical_matched_covariance_sigma_njy":empirical_sigma,
        "control_coordinate_grid_arcsec":[p for p,k in zip(use_coords,keep) if bool(k)],
        "scientific_limits":{
            "calibrated_psf_injection_performed":False,
            "actual_NIRSpec_spec2_calibration_performed":False,
            "cross_epoch_instrument_noise_covariance_isolated":False,
            "formal_p_value_or_sigma_from_this_receipt":None,
            "intrinsic_variability_confirmed":False,
            "line_identification":None,
            "redshift":None,
            "directed_gas_trail":None
        }
    }
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status":out["status"],"controls":len(clean),
                      "correlation":np.round(correlation,3).tolist(),
                      "fixed_weight_empirical_sigma_njy":round(empirical_sigma,4),
                      "nominal_sigma_njy":round(nominal_sigma,4),
                      "output":str(path)},indent=2),flush=True)

if __name__=="__main__":
    main()
