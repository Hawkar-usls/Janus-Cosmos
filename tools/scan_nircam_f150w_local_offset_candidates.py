#!/usr/bin/env python3
"""Blind-ish local 4-epoch F150W offset-aperture search, with field-max null.

Uses four independently SHA-verified original MAST i2d SCI strips; each actual
sky position is aperture measured after local median-sky subtraction. Repeat
identical offset-grid max search around the 95 pinned same-sky controls to avoid
promoting a peak selected after a spatial search into a single-pixel sigma.

This maps BROADBAND light, not the 16004A line. Static galaxies and correlated
sky control windows are included; the field-max distribution does not yield a
formal false-alarm probability or an astrophysical line-host identification.
"""
import argparse
import json
import math
import statistics
import warnings
from pathlib import Path

import numpy as np
from replay_jwst_mocss_mast_subset import RA,DEC,aperture_photometry
from measure_nircam_matched_covariance import (
    RECEIPT, NEIGHBORS, candidate_points, materialize, robust_clip_matrix,
)


def sky(center_ra,center_dec,east,north):
    return center_ra+east/(3600*math.cos(math.radians(center_dec))), center_dec+north/3600


def offset_grid(step=.075,max_radius=.45):
    offsets=[]
    count=math.ceil(max_radius/step)
    for a in range(-count,count+1):
        for b in range(-count,count+1):
            east,north=a*step,b*step
            if math.hypot(east,north)<max_radius+1e-9:
                offsets.append((east,north))
    assert (0,0) in offsets
    return offsets


def sample(frames,center,offsets,weights,medians):
    all_rows=[]
    for east,north in offsets:
        ra,dec=sky(center[0],center[1],east,north)
        measures=[aperture_photometry(f["image"],None,f["wcs"],ra,dec,y_origin=f["y0"])
                  for f in frames]
        if any(m is None or not math.isfinite(m["aperture_njy"]) for m in measures):
            continue
        flux=np.array([m["aperture_njy"] for m in measures])
        val=float(weights@(flux-medians))
        all_rows.append({"offset_east_arcsec":east,"offset_north_arcsec":north,
                         "weighted_aperture_njy":val,
                         "radius_arcsec":math.hypot(east,north)})
    return all_rows


def maximum(rows,radius):
    choices=[r for r in rows if r["radius_arcsec"]<=radius+1e-10]
    return max(choices,key=lambda x:x["weighted_aperture_njy"]) if choices else None


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--data-dir",type=Path,default=Path("/tmp/janus_f150w_time_domain"))
    ap.add_argument("--output",type=Path,default=Path("/tmp/janus_f150w_offset_search.json"))
    a=ap.parse_args()
    rows=sorted([r for r in json.loads(RECEIPT.read_text())["results"] if r.get("coverage") and
                 r.get("date_time_utc","").startswith("2023-10")],key=lambda r:r["date_time_utc"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        frames=[materialize(r,a.data_dir) for r in rows]
    errors=np.array([r["control_noise_1sigma_njy"] for r in rows])
    weights=errors**-2;weights/=weights.sum()
    coords=[]
    original=[]
    for e,n in candidate_points(5.4,.9):
        ra,dec=sky(RA,DEC,e,n)
        ms=[aperture_photometry(f["image"],None,f["wcs"],ra,dec,y_origin=f["y0"])
            for f in frames]
        if all(z is not None for z in ms):
            coords.append([ra,dec,e,n])
            original.append([z["aperture_njy"] for z in ms])
    original=np.array(original)
    good,_,_=robust_clip_matrix(original)
    coords=[x for x,k in zip(coords,good) if k]
    original=original[good]
    assert len(coords)==95
    medians=np.median(original,axis=0)
    offsets=offset_grid()
    target=sample(frames,(RA,DEC),offsets,weights,medians)
    assert len(target)==len(offsets)
    radii=[.15,.225,.30,.45]
    actual_max={str(rad):maximum(target,rad) for rad in radii}
    centers=[]
    for ra,dec,e,n in coords:
        control=sample(frames,(ra,dec),offsets,weights,medians)
        if len(control)==len(offsets):
            centers.append({"center_east_arcsec":e,"center_north_arcsec":n,
                            "maxima_njy":{str(rad):maximum(control,rad)["weighted_aperture_njy"] for rad in radii}})
    assert len(centers)>=80, "control scan grid requires >=80 complete coverage fields"
    controls={}
    for rad in radii:
        key=str(rad)
        mx=np.array([q["maxima_njy"][key] for q in centers])
        actual=actual_max[key]["weighted_aperture_njy"]
        controls[key]={
            "target_field_max_njy":actual,
            "same_radius_matched_control_fields":len(mx),
            "control_field_max_median_njy":float(np.median(mx)),
            "control_field_max_95pct_njy":float(np.percentile(mx,95)),
            "control_fields_with_max_ge_target_max":int(np.sum(mx>=actual)),
            "empirical_control_exceedance_fraction_NOT_a_pvalue":float(np.mean(mx>=actual)),
        }
    # The closest known catalog galaxy is ~0.635 arcsec away and can leak
    # into outer radius=.45 apertures; never classify a bright point as source.
    out={
        "schema":"janus-cosmos-actual-f150w-four-epoch-local-offset-aperture-scan-v1",
        "status":"BROADBAND_LOCAL_MAX_SCANNED_WITH_MATCHED_FIELD_MAX_CONTROLS_NO_LINE_HOST_IDENTIFICATION",
        "source_ra_dec_deg":[RA,DEC],
        "images":[{"date":f["date"],"mast_original_i2d_fits":f["name"],"sha256":f["hash"]} for f in frames],
        "method":{"aperture_radius_arcsec":.12,"annulus_arcsec":[.24,.38],
                  "east_north_offset_grid_arcsec":.075,"max_search_radius_arcsec":.45,
                  "num_scan_positions":len(target),"fixed_weights":weights.tolist(),
                  "shared_clipped_controls":95,"complete_control_field_grid_count":len(centers),
                  "post_selection_is_accounted_for_using_fieldwise_maxima":True,
                  "source_companions_ra_dec_deg":NEIGHBORS,
                  "warning":"This is observed F150W broad-band aperture light. Positional scan has multiple tested points, repeated sky backgrounds, potential galaxies and drizzled noise; no formal p-value or line flux inferred."},
        "replayed_control_frame_medians_njy":medians.tolist(),
        "target_map_all_offset_positions":target,
        "target_field_maxima_by_radius":actual_max,
        "matched_control_fieldmax_distribution":controls,
        "science_limits":{"line_identification":None,"source_class":None,"redshift":None,"sky_direction_of_emitting_gas":None,
                          "multi_aperture_extended_source_test_done":False,
                          "completely_independent_noise_realizations":False,
                          "true_science_target_psf_calibration":False,
                          "statistically_calibrated_false_alarm_p":None},
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status":out["status"],"source_aperture_centered_njy":next(r["weighted_aperture_njy"] for r in target if r["radius_arcsec"]==0),
                      "search_positions":len(target),"controls":len(centers),"fieldmax":controls,
                      "output":str(a.output)},indent=2))

if __name__=="__main__":
    main()
