#!/usr/bin/env python3
"""Reproducible, *approximate* F150W PSF injection into actual MAST SCI pixels.

Four SHA-256-pinned JW01180 October 2023 i2d mosaics are used without replacing
their observed pixels. Pixel-integrated synthetic PSFs are added to each local
measurement stamp, followed by the SAME 0.12 arcsec aperture / 0.24--0.38
arcsec median sky annulus as the existing archival analysis.

IMPORTANT: The two analytic PSFs are constrained by official STScI F150W i2d
EE50=.051 arcsec and EE80=.134 arcsec, NOT empirical PSFs extracted from these
four fields or STPSF optical simulations. This is sensitivity / recovery
bracketing, neither a calibrated source flux nor a formal variability test.

Replay: PYTHONPATH=<numpy+astropy+scipy env> python3
  tools/inject_nircam_f150w_empirical_psf_brackets.py
  --data-dir /tmp/janus_f150w_time_domain --output /tmp/f150w_injected.json

Original archive byte ranges and SHA256 live in
domains/space/NIRCAM_F150W_TIME_DOMAIN_RECEIPT_2026-10-16.json.
"""
import argparse
import json
import math
import warnings
from pathlib import Path

import numpy as np
from astropy.wcs.utils import proj_plane_pixel_area, proj_plane_pixel_scales
from scipy.optimize import root

from measure_nircam_matched_covariance import (
    RECEIPT, candidate_points, materialize, robust_clip_matrix,
)
from replay_jwst_mocss_mast_subset import RA, DEC, aperture_photometry

EE50_ARCSEC = 0.051
EE80_ARCSEC = 0.134
PSF_OVERSAMPLE = 7
AP_ARCSEC = 0.12
SKY_IN_ARCSEC = 0.24
SKY_OUT_ARCSEC = 0.38


def moffat_fit():
    # 2D Moffat EE(r)=1-[1+(r/alpha)^2]^(1-beta).
    def for_beta(beta):
        alpha = EE50_ARCSEC / math.sqrt(2 ** (1 / (beta - 1)) - 1)
        return alpha * math.sqrt(5 ** (1 / (beta - 1)) - 1), alpha

    lo, hi = 1.12, 12.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if for_beta(mid)[0] > EE80_ARCSEC:
            lo = mid
        else:
            hi = mid
    beta = (lo + hi) / 2
    return beta, for_beta(beta)[1]


def core_halo_fit(weight=0.72):
    # A structurally different profile satisfying the two same EE anchors.
    def ee(radius, core, halo):
        return 1 - weight * math.exp(-radius ** 2 / (2 * core ** 2)) - \
            (1 - weight) * math.exp(-radius ** 2 / (2 * halo ** 2))

    sol = root(lambda p: [ee(EE50_ARCSEC, *p) - 0.5,
                          ee(EE80_ARCSEC, *p) - 0.8], [.034, .16])
    if not sol.success or not (0 < sol.x[0] < sol.x[1]):
        raise RuntimeError("Two-Gaussian PSF-anchor solution failed")
    return float(sol.x[0]), float(sol.x[1]), float(weight)


def radial_density(r2, model, params):
    if model == "moffat":
        beta, alpha = params
        return (beta - 1) / (math.pi * alpha ** 2) * (1 + r2 / alpha ** 2) ** (-beta)
    core, halo, weight = params
    return (weight / (2 * math.pi * core ** 2) * np.exp(-r2 / (2 * core ** 2)) +
            (1 - weight) / (2 * math.pi * halo ** 2) * np.exp(-r2 / (2 * halo ** 2)))


def psf_ee(radius, model, params):
    if model == "moffat":
        beta, alpha = params
        return 1 - (1 + (radius / alpha) ** 2) ** (1 - beta)
    core, halo, weight = params
    return 1 - weight * math.exp(-radius ** 2 / (2 * core ** 2)) - \
        (1 - weight) * math.exp(-radius ** 2 / (2 * halo ** 2))


def measure_local(frame, ra, dec, injected_njy=0, model=None, params=None,
                  injection_offset_en_arcsec=(0,0)):
    """Run original estimator with PSF integrated into original local pixels.

    nJy/(MJy/sr pixel) is computed from the actual celestial WCS. No global
    image copies, no smoothed/bespoke mock noise, no normalization to a finite
    stamp, and no changed FITS source files.
    """
    image, wcs = frame["image"], frame["wcs"]
    cx, cy = wcs.all_world2pix([[ra, dec]], 0)[0]
    cy -= frame["y0"]
    scale = float(np.mean(proj_plane_pixel_scales(wcs.celestial)) * 3600)
    r_ap, r_in, r_out = (rad / scale for rad in
                          [AP_ARCSEC, SKY_IN_ARCSEC, SKY_OUT_ARCSEC])
    if cx < r_out + 2 or cy < r_out + 2 or cx > image.shape[1] - r_out - 2 or cy > image.shape[0] - r_out - 2:
        return None
    x0, x1 = math.floor(cx-r_out-1), math.ceil(cx+r_out+2)
    y0, y1 = math.floor(cy-r_out-1), math.ceil(cy+r_out+2)
    ys, xs = np.mgrid[y0:y1, x0:x1]
    rsq = (xs-cx)**2 + (ys-cy)**2
    patch = np.array(image[y0:y1, x0:x1], dtype=np.float64)
    good = np.isfinite(patch)
    ap = good & (rsq <= r_ap**2)
    sky = good & (rsq >= r_in**2) & (rsq <= r_out**2)
    if ap.sum() < .8*math.pi*r_ap**2 or sky.sum() < 35:
        return None

    pixel_area_sr = float(proj_plane_pixel_area(wcs.celestial)) * (math.pi/180)**2
    if injected_njy:
        assert model in {"moffat", "two_gaussian"} and params
        inj_e, inj_n = injection_offset_en_arcsec
        inj_ra = ra + inj_e/(3600*math.cos(math.radians(dec)))
        inj_dec = dec + inj_n/3600
        icx, icy = wcs.all_world2pix([[inj_ra,inj_dec]],0)[0]
        icy -= frame["y0"]
        # Uniform 7x7 midpoint quadrature integrates the radial surface
        # brightness across each pixel, including subpixel source centering.
        integrated = np.zeros_like(patch)
        n = PSF_OVERSAMPLE
        for iy in range(n):
            dy = (iy+.5)/n-.5
            for ix in range(n):
                dx = (ix+.5)/n-.5
                radius2 = ((xs+dx-icx)*scale)**2 + ((ys+dy-icy)*scale)**2
                integrated += radial_density(radius2, model, params)
        integrated *= scale**2/n**2
        # Units: [nJy per sky-pixel] / [1e15*nJy/(MJy/sr * sr)]
        patch += (injected_njy/(pixel_area_sr*1e15))*integrated

    background = float(np.median(patch[sky]))
    flux = float(np.sum(patch[ap]-background))*pixel_area_sr*1e15
    return {"flux_njy": flux, "background_mjysr":background,
            "x":float(cx),"y":float(cy),"aperture_pixels":int(ap.sum())}


def percentile(values):
    return {str(p):float(np.percentile(values,p)) for p in [0, 5, 16, 50, 84, 95, 100]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="/tmp/janus_f150w_time_domain")
    parser.add_argument("--output", default="/tmp/janus_f150w_psf_injection.json")
    args = parser.parse_args()
    saved = json.loads(RECEIPT.read_text())
    rows = sorted([r for r in saved["results"] if r.get("coverage") and
                  r.get("date_time_utc","").startswith("2023-10")],
                  key=lambda r:r["date_time_utc"])
    assert [r["date_time_utc"][:10] for r in rows] == [
        "2023-10-01", "2023-10-02", "2023-10-03", "2023-10-09"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        frames = [materialize(row, Path(args.data_dir)) for row in rows]

    all_coords = candidate_points(5.4, 0.9)
    coords, controls = [], []
    for east, north in all_coords:
        ra = RA + east/(3600*math.cos(math.radians(DEC)))
        dec = DEC + north/3600
        measured = [measure_local(f,ra,dec) for f in frames]
        if all(m is not None for m in measured):
            coords.append((ra,dec,east,north))
            controls.append([m["flux_njy"] for m in measured])
    good, _, _ = robust_clip_matrix(np.asarray(controls))
    coords = [c for c,k in zip(coords,good) if k]
    controls = np.asarray(controls)[good]
    assert len(coords) == 95, "Control grid/astrometry differs from covariance receipt"
    fixed_errors=np.array([r["control_noise_1sigma_njy"] for r in rows])
    weights=(fixed_errors**(-2)); weights/=weights.sum()
    median=np.median(controls,axis=0)
    cov=np.cov(controls,rowvar=False,ddof=1)
    cov_sigma=float(np.sqrt(weights@cov@weights))
    baseline=[measure_local(f,RA,DEC) for f in frames]
    for f,b in zip(frames,baseline):
        reference=aperture_photometry(f["image"],None,f["wcs"],RA,DEC,y_origin=f["y0"])
        # Original estimator sums float32 mosaic pixels; our injected stamps
        # are intentionally float64 to avoid quantizing faint injected flux.
        assert abs(b["flux_njy"]-reference["aperture_njy"]) < 1e-3
    baseline_array=np.array([b["flux_njy"] for b in baseline])
    observed_centered=float(weights@(baseline_array-median))
    assert abs(cov_sigma-.662980825637632) < 1e-5

    models={"moffat":moffat_fit(), "two_gaussian":core_halo_fit()}
    profiles={}
    for model,params in models.items():
        assert all(abs(psf_ee(rad,model,params)-fraction)<1e-10 for rad,fraction in
                   [(EE50_ARCSEC,.5),(EE80_ARCSEC,.8)])
        injections={}
        for total_flux in [1.1,2.5,5.5]:
            measured=[measure_local(f,RA,DEC,total_flux,model,params)["flux_njy"]
                      for f in frames]
            shifts=np.array(measured)-baseline_array
            effective=float(weights@shifts/total_flux)
            injections[str(total_flux)]={
                "target_epoch_aperture_njy":measured,
                "incremental_aperture_njy":shifts.tolist(),
                "target_epoch_recovery_fractions":(shifts/total_flux).tolist(),
                "weighted_aperture_increment_njy":float(weights@shifts),
                "weighted_aperture_recovery_fraction":effective,
                "weighted_centered_target_after_injection_njy":observed_centered+float(weights@shifts),
            }
        # 95 distinct, matched-sky source positions each receive the SAME
        # point-source injection in all four actual October science images.
        # This propagates native subpixel phase, local annulus, and sky
        # structure; spatial positions are correlated, NOT 95 independent noise draws.
        test_flux=5.5
        injected=[]
        for ra,dec,_,_ in coords:
            row=[measure_local(f,ra,dec,test_flux,model,params)["flux_njy"]
                 for f in frames]
            injected.append(row)
        injected=np.array(injected)
        response=injected-controls
        response_weighted=response@weights
        result_centered=(injected-median)@weights
        noiseless_expected=observed_centered+injections["5.5"]["weighted_aperture_increment_njy"]
        # Astrometry stress test: leave the extraction position fixed while
        # moving the injected source across 8 astronomical position angles.
        # PA measured eastward from north, NOT a measured motion vector.
        shift_trials=[]
        for radius in [.025,.05,.075,.10,.15,.20]:
            directional=[]
            for pa in range(0,360,45):
                a=math.radians(pa)
                offset=(radius*math.sin(a),radius*math.cos(a))
                shifted=[measure_local(f,RA,DEC,test_flux,model,params,
                                        injection_offset_en_arcsec=offset)["flux_njy"]
                         for f in frames]
                recovered=float(weights@(np.array(shifted)-baseline_array))/test_flux
                directional.append({"position_angle_deg_east_of_north":pa,
                                    "effective_recovery_fraction":recovered})
            shift_trials.append({
                "injected_offset_arcsec":radius,
                "min_recovery_fraction":min(x["effective_recovery_fraction"] for x in directional),
                "median_recovery_fraction":float(np.median([x["effective_recovery_fraction"] for x in directional])),
                "max_recovery_fraction":max(x["effective_recovery_fraction"] for x in directional),
                "all_eight_position_angles":directional,
            })
        # This counts observed source-like aperture deficits using controls
        # ONLY, not a calibrated frequentist p-value or proof of a transient.
        below=int(np.sum(result_centered<=observed_centered))
        profile={
            "psf_parameters":list(params),
            "official_anchor_ee":{"r0.051_arcsec":.5,"r0.134_arcsec":.8},
            "model_predicted_ee_at_aperture_radius":psf_ee(AP_ARCSEC,model,params),
            "target_injections_njy":injections,
            "astrometric_miscentering_stress_test_flux_5p5_njy":shift_trials,
            "matched_control_injections_5p5_njy":{
                "num_matched_positions":len(coords),
                "weighted_recovery_fraction_quantiles":percentile(response_weighted/test_flux),
                "injected_weighted_centered_njy_quantiles":percentile(result_centered),
                "injected_control_positions_at_or_below_uninjected_target":below,
                "injected_control_positions_exceed_uninjected_target":len(coords)-below,
                "warning":"Sky control positions are spatially correlated; fractions are descriptive, not calibrated p-values."
            },
            "nominal_flux_5p5_diagnostic":{
                "original_spectral_flux_inferred_f150w_total_njy":5.5,
                "spectral_model_bandpass_flux_1sigma_njy":.66,
                "observed_weighted_centered_aperture_njy":observed_centered,
                "counterfactual_persistent_source_weighted_aperture_njy":injections["5.5"]["weighted_aperture_increment_njy"],
                "counterfactual_source_plus_observed_target_njy":noiseless_expected,
                "model_difference_over_empirical_matched_sky_covariance_sd_NOT_formal_significance":
                    injections["5.5"]["weighted_aperture_increment_njy"]/cov_sigma,
                "model_difference_over_empirical_sd_and_spectral_model_uncertainty_NOT_formal_significance":
                    injections["5.5"]["weighted_aperture_increment_njy"]/
                    math.sqrt(cov_sigma**2+(.66*injections["5.5"]["weighted_aperture_recovery_fraction"])**2),
            }
        }
        if model == "two_gaussian":
            # Gaussian intrinsic extended source convolved with our two-
            # Gaussian PSF is analytic (variance adds). Unlike the point-
            # source test, this addresses the AGN-illuminated gas-cloud case.
            # No physical kpc size is inferred without an identified redshift.
            broadened=[]
            core,halo,weight=params
            for intrinsic_sigma in [0,.03,.06,.10,.15,.20]:
                bigger=(math.hypot(core,intrinsic_sigma),
                        math.hypot(halo,intrinsic_sigma),weight)
                summed=np.array([measure_local(f,RA,DEC,5.5,"two_gaussian",bigger)["flux_njy"]
                                 for f in frames])
                recovered=float(weights@(summed-baseline_array))
                broadened.append({
                    "intrinsic_angular_gaussian_sigma_arcsec":intrinsic_sigma,
                    "intrinsic_angular_gaussian_fwhm_arcsec":2*math.sqrt(2*math.log(2))*intrinsic_sigma,
                    "aperture_recovered_fraction":recovered/5.5,
                    "predicted_centered_aperture_njy_if_persistent":recovered,
                    "conditional_model_difference_over_matched_covariance_sd_not_formal_significance":recovered/cov_sigma,
                })
            profile["intrinsic_extended_emission_stress_test_5p5_njy"]=broadened
        profiles[model]=profile

    out={
        "schema":"janus-cosmos-actual-f150w-mast-pixels-analytic-psf-injection-v3",
        "status":"REAL_ARCHIVE_PIXEL_INJECTION_APPROXIMATE_PSF_NOT_CALIBRATED_PHYSICAL_IDENTIFICATION",
        "source_ra_dec_deg":[RA,DEC],
        "science_primary_source":"https://arxiv.org/html/2610.12553v1",
        "stscI_official_psf_source":"https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-point-spread-functions",
        "method":{
            "original_i2d_sci_pixel_data":True,
            "all_four_original_sci_sha256_verified":True,
            "observed_image_pixel_values_edited_on_disk":False,
            "formal_stpsf_or_extracted_field_psf":False,
            "analytic_psf_bracket":"Moffat and core+halo Gaussian mixture, each exactly calibrated at official i2d F150W EE50 and EE80 radii",
            "i2d_pixel_subsamples_per_axis":PSF_OVERSAMPLE,
            "aperture_arcsec":AP_ARCSEC,
            "background_annulus_arcsec":[SKY_IN_ARCSEC,SKY_OUT_ARCSEC],
            "target_subpixel_centers_x_y_in_sci_stripe": [[b["x"],b["y"]] for b in baseline],
            "sky_control_grid":"0.9 arcsec matched coordinates inside 5.4 arcsec, bright-sky positions jointly clipped, 95 positions",
            "source_light_profile_assumption":"unresolved point source, equal observed total flux at all four Oct2023 epochs",
            "injection_is_in_synthetic_local_pixel_stamps":"Equivalent to inserting 2D PSF pixel-integrals into original SCI images before existing aperture and local-background estimator; original images stay unchanged",
            "wcs_pixel_scales_arcsec": [b["target"]["pixel_scale_arcsec"] for b in frames],
            "unproven_assumptions":["Pointlike emission", "Absolute JWST NIRSpec line flux from paper is correct", "NIRCam synthetic photometry 5.5nJy from NIRSpec line alone is correct", "Official global F150W encircled-energy anchors describe these mosaics", "Sky-control sample approximates target noise and confusion despite spatial correlations"]
        },
        "input_epochs":[{"product":f["name"],"utc_date":f["date"],"original_sci_subset_sha256":f["hash"]} for f in frames],
        "matched_sky_controls":len(coords),
        "matched_control_aperture_central_njy":median.tolist(),
        "fixed_weights":weights.tolist(),
        "original_target_epoch_local_background_subtracted_njy":baseline_array.tolist(),
        "original_target_weighted_with_this_grid_centered_njy":observed_centered,
        "measured_matched_sky_covariance_sigma_njy":cov_sigma,
        "note_different_baseline":"0.371nJy target value here uses exactly the 95 clipped matched-grid control medians. Previously published 0.278nJy uses separate first-stage per-epoch control medians; do not compare/mix conventions.",
        "profiles":profiles,
        "scientific_limits":{
            "official_or_field_specific_instrument_psf_injection":False,
            "validated_confidence_or_false_negative_probability":False,
            "calibrated_spec2_nirspec_line_flux":False,
            "2023_october_source_intrinsic_variability_confirmed":False,
            "source_class":None,"redshift":None,"unique_host":None,"direction_of_physical_propagation":None,"directed_gas_trail":None,
            "interpretation":"If the reported NIRSpec line generated 5.5 nJy of persistent unresolved F150W emission in Oct2023, both instrument-EE-anchored approximate PSFs predict an aperture excess over the observed flux; this does not identify the line emitter or prove a transient."
        }
    }
    output=Path(args.output); output.parent.mkdir(exist_ok=True,parents=True)
    output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status":out["status"],"observed_weighted_njy":observed_centered,
                      "matched_sky_sigma_njy":cov_sigma,"models":{
                          name:{"target_recovery_frac":p["target_injections_njy"]["5.5"]["weighted_aperture_recovery_fraction"],
                                "conditional_diagnostic":p["nominal_flux_5p5_diagnostic"]["model_difference_over_empirical_sd_and_spectral_model_uncertainty_NOT_formal_significance"],
                                "below_target":p["matched_control_injections_5p5_njy"]["injected_control_positions_at_or_below_uninjected_target"]}
                          for name,p in profiles.items()},"output":str(output)},indent=2))

if __name__ == "__main__":
    main()
