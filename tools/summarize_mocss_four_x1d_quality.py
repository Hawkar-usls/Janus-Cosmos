#!/usr/bin/env python3
"""Freeze quality-gated diagnostics from four SHA-validated MOCSS x1d FITS.

Consumes the independent existing real-FITS stress receipt, which SHA-checks
all four original diagnostic X1D files and measures off-line controls. This
summary deliberately never turns formal ratios into detection significances.
"""

import argparse
import hashlib
import json
from pathlib import Path


ORIGINAL_RUN = 39172038085
ARCHIVE_ZIP_DIGESTS = {
    "A_C": {
        "artifact_id": 12153067969,
        "sha256": "b9dcc599a13094729e61e3dfdf5a8db3413c3051a1441f0e93a7a2d68269f941",
    },
    "B_D": {
        "artifact_id": 12153078118,
        "sha256": "37907dde426f898452de43a3e996d43db71d5ac189474feb90843ab04e071fe5",
    },
}


def freeze(receipt):
    assert receipt["upstream_action_run_id"] == ORIGINAL_RUN
    assert receipt["test_source_info_invented"] is True
    assert receipt["valid_mocss_row_background_subtraction"] is False
    assert receipt["line_detection_from_this_run"] is None
    assert set(receipt["results"]) == set("ABCD")
    labels = {"16004A", "OIII_4960", "Pa_gamma"}
    groups = {}
    rejected = 0
    misleading = []
    nulls_total = 0
    nulls_formal_ge3 = 0
    off_line_max = 0.0
    for group in "ABCD":
        fields = receipt["results"][group]["line_windows"]
        assert set(fields) == labels
        windows = {}
        for label in sorted(labels):
            src = fields[label]
            chi = float(src["reduced_chi2"])
            formal = float(src["formal_fit_ratio_NOT_DETECTION_SIGNIFICANCE"])
            off = src["off_line"]
            assert off["status"] == "EMPIRICAL_DIAGNOSTIC_CONTROLS_NOT_SCIENCE_NULL"
            n = int(off["n"])
            count = round(n * float(off["fraction_control_abs_formal_ratio_ge_3"]))
            assert n >= 12 and 0 <= count <= n
            bad = chi > 3
            assert bad is src["chi2_model_rejected_for_unit_variance_if_gt_3"]
            rejected += int(bad)
            nulls_total += n
            nulls_formal_ge3 += count
            off_line_max = max(off_line_max, float(off["maximum_absolute_formal_template_ratio"]))
            if abs(formal) >= 5 and bad:
                misleading.append({"group": group, "wavelength_window": label,
                                   "formal_amplitude_over_error_NOT_SIGNIFICANCE": round(formal, 3),
                                   "reduced_chi2": round(chi, 3)})
            windows[label] = {
                "reduced_chi2": round(chi, 3),
                "formal_amplitude_over_error_NOT_SIGNIFICANCE": round(formal, 3),
                "quality_model_rejected_at_reduced_chi2_gt_3": bad,
                "off_line_control_window_count": n,
                "off_line_controls_formal_ratio_ge_3_count": count,
                "off_line_largest_absolute_formal_ratio_NOT_SIGNIFICANCE": round(
                    float(off["maximum_absolute_formal_template_ratio"]), 3
                ),
            }
        groups[group] = {"sha256_of_actual_real_rate_derived_diagnostic_x1d":
                         receipt["results"][group]["x1d_sha256"], "line_windows": windows}
    return {
        "schema": "janus-mocss-four-x1d-calibration-failure-and-local-control-audit-v1",
        "date_utc": "2026-10-19",
        "upstream_action_run_id": ORIGINAL_RUN,
        "upstream_action_url": f"https://github.com/Hawkar-usls/Janus-Cosmos/actions/runs/{ORIGINAL_RUN}",
        "upstream_full_stress_reproducer": "tools/stress_mocss_four_x1d_offline_and_injection.py",
        "replay_environment": "Python 3.12, numpy 2.5.4, astropy 8.0.2",
        "source_archive_zip_artifact_sha256": ARCHIVE_ZIP_DIGESTS,
        "source": "Four actual MAST RATE-derived NIRSpec x1d FITS from synthetic-source three-shutter Stage2 compatibility experiment; NOT original MOCSS pipeline output",
        "source_is_original_MPT_target": False,
        "pre_fix_three_shutter_axis": "CROSS_DISPERSION_COLUMNS__WRONG_FOR_MOCSS_DISPERSION_ROW_BAND",
        "groups": groups,
        "aggregate": {
            "groups": 4,
            "fixed_window_fits": 12,
            "fixed_window_reduced_chi2_above_3": rejected,
            "fixed_windows_with_formal_absolute_ratio_ge_5_AND_bad_fit": misleading,
            "local_control_window_measurements_not_independent": nulls_total,
            "local_control_formal_absolute_ratio_ge_3_not_independent": nulls_formal_ge3,
            "largest_offline_control_absolute_formal_ratio_NOT_SIGNIFICANCE": round(off_line_max, 3),
        },
        "interpretation": "Control-window formal ratios and severe fit residuals demonstrate that nominal per-bin x1d errors do not calibrate detection significance in this synthetic-MSA extraction. Off-line regions can contain real lines and share wavelength windows, so they are NOT an independent blank-sky false-positive distribution. Any nominal x1d line flux, significance or limit is invalid.",
        "science_gate": {
            "independent_valid_mocss_line_detection": False,
            "oiii_companion_at_15850A_verified": False,
            "pagamma_at_13656A_verified": False,
            "source_class_or_distance_verified": False,
            "directed_spatial_trail_verified": False,
            "next_test": "Reconstruct original open-MSA dispersion-row geometry, four-dither background and PSF/bar throughput before calibrated profile/line-ratio inference",
        },
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    raw = args.input.read_bytes()
    result = freeze(json.loads(raw))
    result["full_stress_receipt_sha256"] = hashlib.sha256(raw).hexdigest()
    data = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(data)
    else:
        print(data, end="")
    print("FOUR_X1D_CONTROL_QUALITY_GATES_APPLIED_NO_SCIENCE_LINE_DETECTION")


if __name__ == "__main__":
    main()
