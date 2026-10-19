#!/usr/bin/env python3
"""Search the public, *uncalibrated* four MOCSS rate-stack images for directions.

Inputs are immutable original-RATE WCS-aligned cutouts frozen by the COSMOS
workflow, and the published-coordinate local JWST WCS Jacobians.  We sample
the same observed wavelength at three nearby galaxies and matched-radius sky
directions, rather than interpreting detector dispersion as a celestial trail.

This does not reproduce the special MOCSS row-illumination background model,
propagate interpolation covariance, deblend the central optical PSF, or assign
science-grade significance / calibrated flux.  Circular controls overlap;
their ranks are descriptive, not probabilities.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


SOURCE_KEY = "primary"
GROUPS = "ABCD"
RING_ANGLES_DEG = tuple(range(0, 360, 5))
DIRECTED_RADII_ARCSEC = (0.35, 0.6, 0.85, 1.1)
HALF_X, HALF_Y = 70, 15
X, Y = np.meshgrid(np.arange(-HALF_X, HALF_X + 1),
                   np.arange(-HALF_Y, HALF_Y + 1))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def matched_pixel_response(frame, dx, dy):
    """Aperture matched to a diagnostic Gaussian, not to NIRSpec's real PSF."""
    weights = np.exp(-0.5 * ((X - dx) / 1.2) ** 2
                     -0.5 * ((Y - dy) / 1.5) ** 2)
    support = weights > 0.0001
    usable = np.isfinite(frame) & support
    if np.count_nonzero(usable) < 15:
        return None
    return float(np.sum(frame[usable] * weights[usable]) /
                 np.sum(weights[usable] ** 2))


def map_sky(jacobian, east_arcsec, north_arcsec):
    return (jacobian[0][0] * east_arcsec + jacobian[0][1] * north_arcsec,
            jacobian[1][0] * east_arcsec + jacobian[1][1] * north_arcsec)


def azimuth_control(frame, jacobian, separation_arcsec):
    values = []
    for degrees in RING_ANGLES_DEG:
        az = math.radians(degrees)
        dx, dy = map_sky(jacobian, separation_arcsec * math.sin(az),
                         separation_arcsec * math.cos(az))
        value = matched_pixel_response(frame, dx, dy)
        if value is not None:
            values.append(value)
    return values


def summarize_control(candidate_value, controls):
    median = float(np.median(controls))
    mad = float(1.4826 * np.median(np.abs(np.asarray(controls) - median)))
    return {"relative_value_dn_per_s_not_calibrated": candidate_value,
            "ring_control_n_correlated": len(controls),
            "ring_median_dn_per_s": median,
            "ring_robust_mad_dn_per_s": mad,
            "offset_from_ring_median_in_mad_not_snr":
                ((candidate_value - median) / mad if mad > 0 else None),
            "empirical_fraction_ring_directions_at_least_as_bright_not_pvalue":
                float(np.mean(np.asarray(controls) >= candidate_value))}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--arrays", type=Path, required=True)
    p.add_argument("--geometry", type=Path,
                   default=Path("domains/space/MOCSS_LOCAL_SKY_WCS_PROJECTION_2026-10-16.json"))
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    geo = json.loads(args.geometry.read_text())
    positions = geo["neighbors"]
    assert set(geo["groups"]) == set(GROUPS)
    out = {
      "schema": "janus-mocss-four-rate-stack-same-wavelength-neighbor-azimuth-v1",
      "case_id": "JWST-MOCSS-HUDF-16004A",
      "inputs": {
        "original_Rate_stack_npz_sha256": sha256(args.arrays),
        "sky_geometry_json_sha256": sha256(args.geometry),
        "expected_stack_npz_sha256": "0b8974e9f761e111abce46d16e0db322f368fd67306b7e8920cfc92fefbe8190",
        "16_original_RATE_and_2_MSA_verified_upstream": True,
        "upstream_run": 39230522837,
        "upstream_artifact": 12176970014,
        "archive_validation": "No new raw FITS re-reduction; checks upstream artifact hash"
      },
      "line_lambda_A": 16004,
      "positions": positions,
      "group_results": {},
      "coordinate_meaning": "Offsets east/north are tangent-plane angular separations from line-emitter; 2x2 WCS Jacobians map same observed lambda to a nearby detector location.",
      "statistics_meaning": "Annular 5-degree angles are correlated and central source PSF may contaminate small radii; ring fractions are NOT independent probabilities or calibrated SNR.",
      "a_priori_directional_track_radii_arcsec": DIRECTED_RADII_ARCSEC,
      "radius_coverage_warning": "Tests extend only to 1.1 arcsec along three preselected bearings; cutout edge, approximate WCS, incomplete optics and raw-sky defects prohibit any physical length upper limit.",
      "science_limits": {
        "original_row_shutter_illumination_sky_model_reproduced": False,
        "independent_line_flux_detection_or_upper_limit": False,
        "NIRSpec_dispersed_source_PSF_and_neighbor_deblending": False,
        "trials_corrected_spatial_detection": False,
        "line_spatial_extent_established": False,
        "direction_or_physical_trail_established": False,
        "emitting_source_identity_and_distance": None
      }
    }
    with np.load(args.arrays) as archive:
        if sha256(args.arrays) != out["inputs"]["expected_stack_npz_sha256"]:
            raise RuntimeError("Frozen original-stack NPZ SHA mismatch")
        assert len(archive.files) == 12, archive.files
        for group in GROUPS:
            im = archive[group + "_" + SOURCE_KEY]
            if im.shape != (31, 141):
                raise ValueError((group, im.shape))
            jac = geo["groups"][group]["jacobian_detector_pixels_per_arcsec"]
            result = {
                "image_finite_fraction": float(np.isfinite(im).mean()),
                "position_tests": {},
                "directed_radii_tests": {},
            }
            for neighbor, location in positions.items():
                e, n = location["east_arcsec"], location["north_arcsec"]
                r = math.hypot(e, n)
                dx, dy = map_sky(jac, e, n)
                v = matched_pixel_response(im, dx, dy)
                assert v is not None
                controls = azimuth_control(im, jac, r)
                result["position_tests"][neighbor] = {
                   "radius_arcsec": r,
                   "position_angle_degrees_east_of_north": math.degrees(math.atan2(e,n)) % 360,
                   "predicted_detector_offset_x_y_px_at_16004A": [dx, dy],
                   **summarize_control(v, controls)
                }
                result["directed_radii_tests"][neighbor] = []
                for radius in DIRECTED_RADII_ARCSEC:
                    radial_e, radial_n = radius * e / r, radius * n / r
                    r_dx, r_dy = map_sky(jac, radial_e, radial_n)
                    radial_value = matched_pixel_response(im, r_dx, r_dy)
                    assert radial_value is not None
                    radial_controls = azimuth_control(im, jac, radius)
                    result["directed_radii_tests"][neighbor].append({
                       "radius_arcsec": radius,
                       "detector_offset_x_y_px": [r_dx, r_dy],
                       **summarize_control(radial_value, radial_controls)
                    })
            out["group_results"][group] = result
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    for g, info in out["group_results"].items():
        print(g, {k: round(v["offset_from_ring_median_in_mad_not_snr"], 2)
                  for k, v in info["position_tests"].items()})
    print("WARNING: RAW_RATE_SKY_PROXY_ONLY; NO PHYSICAL_TRAIL_OR_CALIBRATED_DETECTION")


if __name__ == "__main__":
    main()
