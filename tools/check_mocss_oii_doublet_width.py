#!/usr/bin/env python3
"""Sensitivity calculation for an [O II] interpretation of MOCSS 16004 Å.

This is NOT a spectral extraction or an [O II] rejection test. A two-Gaussian
mixture's second central moment is converted to a Gaussian-equivalent FWHM;
the published single-Gaussian fitted FWHM need not equal that moment. In
particular, the effective slitless MOCSS line-spread function is not measured
by the nominal grating R. The result is a preregistered LSF calibration target.
"""

import argparse
import json
import math
from pathlib import Path


LAMBDA_OBS_A = 16004.0
FWHM_OBS_A = 13.2
FWHM_ERR_A = 2.0
EW_OBS_LOWER_A = 7620.0
F_OII_ASSUMED_ERG_S_CM2 = 2.5e-18
LAMBDA_OII1_VAC_A = 3727.092
LAMBDA_OII2_VAC_A = 3729.875
RATIO_GRID = (0.35, 1.0, 1.5)
RESOLUTION_GRID = (1000, 1200, 1500, 1800, 2200, 3000)
FWHM_TO_SIGMA = 2.0 * math.sqrt(2.0 * math.log(2.0))


def compute():
    tests = []
    for r in RATIO_GRID:
        p = r / (1.0 + r)
        rest_centroid = (LAMBDA_OII1_VAC_A + r * LAMBDA_OII2_VAC_A) / (1.0 + r)
        z = LAMBDA_OBS_A / rest_centroid - 1.0
        sep = (LAMBDA_OII2_VAC_A - LAMBDA_OII1_VAC_A) * (1.0 + z)
        blend_fwhm_equivalent = FWHM_TO_SIGMA * sep * math.sqrt(p * (1.0 - p))
        predictions = []
        for resolving_power in RESOLUTION_GRID:
            lsf_fwhm = LAMBDA_OBS_A / resolving_power
            width = math.hypot(lsf_fwhm, blend_fwhm_equivalent)
            predictions.append({
                "assumed_effective_resolving_power": resolving_power,
                "assumed_single_component_lsf_fwhm_A": round(lsf_fwhm, 4),
                "predicted_second_moment_equivalent_fwhm_A": round(width, 4),
                "heuristic_difference_over_published_fwhm_error": round(
                    (width - FWHM_OBS_A) / FWHM_ERR_A, 3
                ),
            })

        min_r_by_upper_bound = {}
        for n_sigma in (0, 1, 2):
            bound = FWHM_OBS_A + n_sigma * FWHM_ERR_A
            if blend_fwhm_equivalent >= bound:
                min_r_by_upper_bound[f"published_plus_{n_sigma}sigma"] = None
            else:
                max_lsf_fwhm = math.sqrt(bound**2 - blend_fwhm_equivalent**2)
                min_r_by_upper_bound[f"published_plus_{n_sigma}sigma"] = round(
                    LAMBDA_OBS_A / max_lsf_fwhm, 1
                )

        tests.append({
            "assumed_oii_3729_to_3726_intensity_ratio": r,
            "vacuum_rest_flux_weighted_centroid_A": round(rest_centroid, 5),
            "redshift_if_oii": round(z, 7),
            "expected_observed_doublet_separation_A": round(sep, 5),
            "doublet_only_second_moment_equivalent_fwhm_A": round(
                blend_fwhm_equivalent, 4
            ),
            "rest_ew_lower_bound_if_persistent_A": round(
                EW_OBS_LOWER_A / (1.0 + z), 1
            ),
            "minimum_effective_resolving_power_required": min_r_by_upper_bound,
            "conditional_gaussian_lsf_predictions": predictions,
        })

    return {
        "schema": "janus-mocss-oii-profile-second-moment-conditional-v1",
        "data_origin": "Published Huberty et al. 2026 single-line fit, not a remeasured profile",
        "published_line": {
            "observed_wavelength_A": LAMBDA_OBS_A,
            "single_gaussian_fwhm_A": FWHM_OBS_A,
            "fwhm_1sigma_A": FWHM_ERR_A,
            "reported_flux_erg_s_cm2": F_OII_ASSUMED_ERG_S_CM2,
            "reported_observed_ew_bound_if_persistent_A": EW_OBS_LOWER_A,
        },
        "oii_rest_vacuum_A": [LAMBDA_OII1_VAC_A, LAMBDA_OII2_VAC_A],
        "density_limited_oii_3729_to_3726_ratio_range": [0.35, 1.5],
        "formula": "FWHM_equiv^2=(lambda_obs/R_eff)^2+[2sqrt(2ln2)*delta_lambda_obs*sqrt(p*(1-p))]^2; p=r/(1+r)",
        "assumptions": [
            "Two intrinsically infinitesimal [O II] components and Gaussian instrumental LSF",
            "Exact photon-flux-weighted centroid, same-width component LSF",
            "Published Gaussian-fit FWHM approximated by Gaussian-equivalent second moment",
            "The true local MOCSS line-spread function and slitless source morphology are UNMEASURED",
            "Ratio endpoints are idealized density limits, not an astrophysical likelihood prior",
            "EW bound combines non-contemporaneous continuum imaging and line spectroscopy: not a transient-EW bound",
            "No noise covariance or profile-fitting likelihood enters heuristic differences",
        ],
        "ratio_cases": tests,
        "science_gate": {
            "oii_ruled_out": False,
            "oii_identified": False,
            "redshift_or_distance_measured": False,
            "spatial_trail_measured": False,
            "reported_ew_bound_applies_to_transient": False,
            "decisive_next_test": "Measure effective MOCSS LSF on a compact known narrow line at the same detector wavelength and fit convolved [O II] doublet directly to all original A-D spectra with injection-calibrated uncertainties",
        },
        "primary_references": [
            "https://arxiv.org/abs/2610.12553",
            "https://arxiv.org/abs/astro-ph/0510099",
            "https://physics.nist.gov/PhysRefData/ASD/lines_form.html",
            "https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-instrumentation/nirspec-dispersers-and-filters",
        ],
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    document = json.dumps(compute(), indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(document)
        print("MOCSS_OII_CONDITIONAL_LSF_SENSITIVITY_RECORDED")
    else:
        print(document, end="")


if __name__ == "__main__":
    main()
