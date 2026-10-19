#!/usr/bin/env python3
"""Check the frozen [O II] width sensitivity calculation and its claim ceiling."""

import json
import math
from pathlib import Path

from check_mocss_oii_doublet_width import compute


RECEIPT = Path(__file__).resolve().parents[1] / "domains/space/MOCSS_OII_PROFILE_LSF_SENSITIVITY_2026-10-19.json"


def main():
    frozen = json.loads(RECEIPT.read_text())
    assert frozen == compute(), "Frozen receipt differs from reproducible calculation"
    cases = frozen["ratio_cases"]
    assert len(cases) == 3
    for case in cases:
        assert 3.29 < case["redshift_if_oii"] < 3.30
        assert math.isclose(case["expected_observed_doublet_separation_A"], 11.945, abs_tol=0.004)
        assert case["rest_equivalent_width_lower_bound_A"] > 1700
        previous = float("inf")
        for item in case["conditional_gaussian_lsf_predictions"]:
            width = item["predicted_second_moment_equivalent_fwhm_A"]
            assert 0 < width < previous
            previous = width
    assert frozen["science_gate"]["oii_ruled_out"] is False
    assert frozen["science_gate"]["redshift_or_distance_measured"] is False
    print("MOCSS_OII_LSF_MOMENT_SENSITIVITY_REPRODUCED_NOT_AN_IDENTIFICATION")


if __name__ == "__main__":
    main()
