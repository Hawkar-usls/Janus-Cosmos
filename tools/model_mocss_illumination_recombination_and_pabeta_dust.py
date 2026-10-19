#!/usr/bin/env python3
"""Conditional physics constraints on the 16004-A MOCSS orphan.

This is a *hypothesis stress test*, not a measurement of source identity:
 - source-to-cloud sky vectors and minimum projected one-way light path;
 - Pa beta / Pa gamma extinction needed to evade the published 3-sigma range;
 - hydrogen case-B recombination density needed for a fast October dip.

The October 2023 F150W non-detection is aperture/covariance diagnostic, so
an actual dip is NOT an input fact.  All geometry assumes a shared z=0.24856.
"""

import argparse
from datetime import date
import json
import math
from pathlib import Path


PARSEC_CM = 3.0856775814913673e18
KPC_LY = 3261.5637771674333
ALPHA_B_CM3_S_AT_10000K = 2.6e-13
SECONDS_PER_DAY = 86400
RV = 3.1


def ccm89_ir_a_over_av(wavelength_um, rv=RV):
    """Cardelli, Clayton, Mathis 1989; x=1/lambda, 0.3<=x<=1.1."""
    x = 1.0 / wavelength_um
    assert 0.3 <= x <= 1.1
    return (0.574 - 0.527 / rv) * x ** 1.61


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--spatial", type=Path,
                        default=Path("domains/space/JWST_MOCSS_HUDF_ORIGIN_SPATIAL_TEST_2026-10-16.json"))
    parser.add_argument("--photometry", type=Path,
                        default=Path("domains/space/NIRCAM_F150W_MATCHED_SKY_COVARIANCE_OBSERVED_2026-10-18.json"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    geo = json.loads(args.spatial.read_text())
    phot = json.loads(args.photometry.read_text())
    assert geo["sky_direction"]["right_ascension_deg"] == 53.13127
    assert phot["source_ra_dec_deg"] == [53.13127, -27.78908]
    pa = next(v for v in geo["conditional_distance_scenarios"]["cases"]
              if v["candidate"] == "Pa_beta / Type_IIn_possible")
    assert abs(pa["z"] - (16004 / 12818 - 1)) < 0.0001
    z = pa["z"]
    output = {
      "schema": "janus-mocss-conditional-source-cloud-path-dust-time-v1",
      "observed_line_wavelength_A": 16004,
      "source_ra_dec_deg": [53.13127, -27.78908],
      "scientific_status": "CONDITIONAL_PHYSICS_BOUNDS_NOT_SOURCE_OR_TRAIL_IDENTIFICATION",
      "hypothesis": "Pa beta at z about 0.249, possibly dense-gas Type IIn or off-nuclear AGN illuminated cloud",
      "line_model": {
        "Pa_beta_vacuum_rest_A_approx_paper": 12818,
        "z": z,
        "luminosity_distance_gpc_conditional": pa["luminosity_distance_gpc"],
        "proper_kpc_per_arcsec_conditional": pa["proper_kpc_per_arcsec"],
        "observed_line_flux_erg_s_cm2_from_paper": 2.5e-18,
        "observed_line_luminosity_erg_s_conditional":
          4 * math.pi * (pa["luminosity_distance_gpc"] * 1e9 * PARSEC_CM)**2 * 2.5e-18,
      },
      "geometry": {},
      "dust_model": {},
      "time_model": {},
      "separation_vs_time_delay_warning": "Projected physical separation constrains straight source-to-cloud path below, NOT directly the difference in observed light arrival times: 3D echo-delay geometry is unknown.",
      "evidence_inputs": {
        "source_paper": "https://arxiv.org/html/2610.12553v1",
        "dust_extinction_law": "Cardelli, Clayton, Mathis (1989) ApJ 345,245; https://adsabs.harvard.edu/pdf/1989ApJ...345..245C",
        "case_B_recombination_coefficient": "Hydrogen T=10^4 K alpha_B=2.6e-13 cm^3/s",
        "nircam_contemporaneous_diagnostic": "domains/space/NIRCAM_F150W_MATCHED_SKY_COVARIANCE_OBSERVED_2026-10-18.json"
      },
      "not_demonstrated": [
        "Association between MOCSS emitter and any neighboring galaxy",
        "Measured emission beaming axis, shock direction, or gas filament",
        "Uniform source and gas redshift, physical distance, or 3D orientation",
        "Valid science-grade MOCSS line flux or Pa-gamma line limit re-extracted independently",
        "True 2023 October turn-off and on, pending NIRCam PSF injection and covariance calibration",
        "Case-B recombination ratio or CCM89 dust screen necessarily applies to shocked, dense supernova gas"
      ]
    }
    target_ra, target_dec = output["source_ra_dec_deg"]
    for neighbor in geo["neighbor_positions_from_paper_table_3"]:
        dra = (target_ra - neighbor["ra_deg"]) * 3600 * math.cos(math.radians(target_dec))
        ddec = (target_dec - neighbor["dec_deg"]) * 3600
        separation = math.hypot(dra, ddec)
        bearing_from_galaxy_to_line = math.degrees(math.atan2(dra, ddec)) % 360
        kpc = separation * pa["proper_kpc_per_arcsec"]
        output["geometry"][neighbor["label"]] = {
          "separation_arcsec": separation,
          "galaxy_to_line_east_and_north_arcsec": [dra, ddec],
          "candidate_illumination_bearing_deg_east_of_north": bearing_from_galaxy_to_line,
          "projected_distance_kpc_if_shared_z": kpc,
          "minimum_direct_light_travel_years_if_shared_z": kpc*KPC_LY,
          "galaxy_photometric_redshift": neighbor["photometric_redshift"],
          "path_is_measured": False,
          "association_is_measured": False,
        }

    pa_beta_um, pa_gamma_um = 1.2818, 1.0938
    a_beta = ccm89_ir_a_over_av(pa_beta_um)
    a_gamma = ccm89_ir_a_over_av(pa_gamma_um)
    differential = a_gamma - a_beta
    intrinsic_ratio_beta_over_gamma = 1.6
    predicted_gamma_no_dust = 2.5e-18 / intrinsic_ratio_beta_over_gamma
    output["dust_model"] = {
      "Pa_beta_rest_um_approx": pa_beta_um,
      "Pa_gamma_rest_um_approx": pa_gamma_um,
      "assumed_caseB_intrinsic_Pa_beta_over_Pa_gamma": intrinsic_ratio_beta_over_gamma,
      "expected_Pa_gamma_erg_s_cm2_if_no_dust": predicted_gamma_no_dust,
      "published_global_three_sigma_line_limit_range_erg_s_cm2_NOT_VERIFIED_AT_Pagamma": [8.2e-19, 10.6e-19],
      "law": "CCM89 near-IR extinction at Rv3.1: A(lambda)/A(V)=(0.574-0.527/Rv)*(1/lambda_um)^1.61",
      "A_Pabeta_divided_by_Av": a_beta,
      "A_Pagamma_divided_by_Av": a_gamma,
      "delta_A_Pagamma_minus_Pabeta_divided_by_Av": differential,
      "limits": []
    }
    for limit in (10.6e-19, 8.2e-19):
        differential_mag = 2.5 * math.log10(predicted_gamma_no_dust / limit)
        av = differential_mag / differential
        a_beta_mag = av * a_beta
        output["dust_model"]["limits"].append({
          "assumed_Pagamma_3sigma_limit_erg_s_cm2": limit,
          "required_minimal_differential_attenuation_mag": differential_mag,
          "required_Av_mag_IF_CaseB_CCM_and_limit_hold": av,
          "intrinsic_Pabeta_luminosity_lower_erg_s_IF_dust_screen":
             output["line_model"]["observed_line_luminosity_erg_s_conditional"]*10**(0.4*a_beta_mag),
        })

    august = date(2023,8,23)
    november = date(2023,11,19)
    earliest_oct = date(2023,10,1)
    latest_oct = date(2023,10,9)
    gap_aug_oct = (earliest_oct-august).days
    gap_oct_nov = (november-latest_oct).days
    rest_days = min(gap_aug_oct,gap_oct_nov)/(1+z)
    density_bound = 1/(ALPHA_B_CM3_S_AT_10000K*rest_days*SECONDS_PER_DAY)
    output["time_model"] = {
      "spectral_epoch_dates_from_paper": [str(august),str(november)],
      "intermediate_imaging_dates_from_archived_F150W": phot["epochs_utc"],
      "non_detection_F150W_measured_matched_aperture_njy": phot["source_aperture_flux_njy"],
      "non_detection_F150W_covariance_sigma_njy_diagnostic": phot["matched_sky_covariance_sigma_njy"],
      "line_only_F150W_expected_total_njy_if_persistent_reported_paper": 5.5,
      "spectroscopy_Aug_to_first_Oct_imaging_days": gap_aug_oct,
      "last_Oct_imaging_to_spectroscopy_Nov_days": gap_oct_nov,
      "shortest_rest_frame_fade_or_rebrighten_timescale_days_if_real": rest_days,
      "assumed_H_caseB_recombination_alpha_cm3_s_T10000K": ALPHA_B_CM3_S_AT_10000K,
      "electron_density_at_least_cm_minus3_for_recombination_response_on_this_rest_timescale": density_bound,
      "comparison_recombination_times_years": {
        str(n): 1/(ALPHA_B_CM3_S_AT_10000K*n*SECONDS_PER_DAY*365.25)
        for n in (100,10000,1000000)
      },
      "critical_warning": "This density follows ONLY IF the October aperture non-detection is a genuine large, rapid fade at constant sky coordinates, and IF Case-B recombination dominates the response. It is not a measured density nor a proven variability event."
    }
    assert output["geometry"]["NORTHWEST"]["minimum_direct_light_travel_years_if_shared_z"] > 8200
    assert output["dust_model"]["limits"][0]["required_Av_mag_IF_CaseB_CCM_and_limit_hold"] > 5
    assert output["time_model"]["electron_density_at_least_cm_minus3_for_recombination_response_on_this_rest_timescale"] > 1e6
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    print("IF_PABETA Dust Av:",
          [round(x["required_Av_mag_IF_CaseB_CCM_and_limit_hold"],2)
           for x in output["dust_model"]["limits"]])
    print("IF_REAL_OCT_DIP Electron density >=", round(density_bound), "cm^-3")
    for name, data in output["geometry"].items():
        print(name, "bearing", round(data["candidate_illumination_bearing_deg_east_of_north"],1),
              "light-time", round(data["minimum_direct_light_travel_years_if_shared_z"]), "yr")
    print("NO IDENTIFICATION OR BEAM MEASUREMENT MADE")


if __name__ == "__main__":
    main()
