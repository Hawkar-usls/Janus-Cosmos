#!/usr/bin/env python3
"""Integrity and epistemic checks for the MOCSS 2026-10-19 direction study."""

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def data(path):
    return json.loads((ROOT / path).read_text())


def main():
    scan = data("domains/space/MOCSS_FOUR_STACK_NEIGHBOR_AZIMUTH_NULLS_2026-10-19.json")
    physical = data("domains/space/MOCSS_NEIGHBOR_ILLUMINATION_DUST_RECOMBINATION_2026-10-19.json")
    assert scan["inputs"]["expected_stack_npz_sha256"] == scan["inputs"]["original_Rate_stack_npz_sha256"]
    assert len(scan["group_results"]) == 4
    assert len(scan["positions"]) == 3
    assert scan["science_limits"]["direction_or_physical_trail_established"] is False
    assert scan["science_limits"]["independent_line_flux_detection_or_upper_limit"] is False
    assert physical["scientific_status"].startswith("CONDITIONAL_")
    assert len(physical["not_demonstrated"]) >= 5
    assert physical["line_model"]["z"] < .25
    expected_names = set(scan["positions"])
    for group in "ABCD":
        tests = scan["group_results"][group]["position_tests"]
        profiles = scan["group_results"][group]["directed_radii_tests"]
        assert set(tests) == set(profiles) == expected_names
        for name in expected_names:
            assert tests[name]["ring_control_n_correlated"] == 72
            assert len(profiles[name]) == 4
            assert all(math.isfinite(t["relative_value_dn_per_s_not_calibrated"])
                       for t in profiles[name])
    # No single pre-selected galaxy direction rises consistently above its
    # own matched-radius ring median in all four orientation/dither stacks.
    assert all(not all(scan["group_results"][g]["position_tests"][name][
                    "offset_from_ring_median_in_mad_not_snr"] > 0
                    for g in "ABCD") for name in expected_names)
    north_west = physical["geometry"]["NORTHWEST"]
    south = physical["geometry"]["SOUTH"]
    assert 145 < north_west["candidate_illumination_bearing_deg_east_of_north"] < 150
    assert 352 < south["candidate_illumination_bearing_deg_east_of_north"] < 357
    assert 8000 < north_west["minimum_direct_light_travel_years_if_shared_z"] < 9000
    assert 13000 < south["minimum_direct_light_travel_years_if_shared_z"] < 14000
    assert 5 < physical["dust_model"]["limits"][0]["required_Av_mag_IF_CaseB_CCM_and_limit_hold"] < 6
    assert 8 < physical["dust_model"]["limits"][1]["required_Av_mag_IF_CaseB_CCM_and_limit_hold"] < 10
    assert 1e6 < physical["time_model"]["electron_density_at_least_cm_minus3_for_recombination_response_on_this_rest_timescale"] < 2e6
    print("PASS: 12 pre-selected galaxy/stack positions and 48 directional radii; no all-stack positive neighbor; physical assumptions and science gates valid")


if __name__ == "__main__":
    main()
