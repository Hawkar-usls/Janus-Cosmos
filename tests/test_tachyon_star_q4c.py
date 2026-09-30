import unittest

import numpy as np

from experiments.luci.run_tachyon_star_q4c import (
    EXPECTED_CONTROLS,
    GRID_X,
    RANK_MAX,
    RESIDUAL_Z_MIN,
    TARGET_X,
    add_fixed_gaussian,
    empirical_rank_fraction,
    robust_center_scale,
    unique_gate,
)


class T(unittest.TestCase):
    def test_dense_grid_is_frozen_and_nonoverlapping(self):
        self.assertEqual(len(GRID_X), EXPECTED_CONTROLS)
        self.assertEqual(EXPECTED_CONTROLS, 55)
        self.assertTrue(all(abs(x - TARGET_X) >= 128.0 for x in GRID_X))
        self.assertTrue(all((b - a) >= 32.0 for a, b in zip(GRID_X, GRID_X[1:])))

    def test_rank_resolution(self):
        controls = [float(i) for i in range(EXPECTED_CONTROLS)]
        q = empirical_rank_fraction(100.0, controls)
        self.assertAlmostEqual(q["minimum_possible_rank_fraction"], 1.0 / 56.0)
        self.assertEqual(q["rank_descending_1_is_largest"], 1)
        self.assertAlmostEqual(q["empirical_one_sided_rank_fraction"], 1.0 / 56.0)

    def test_unique_gate_requires_both_rank_and_residual(self):
        controls = [float(i) for i in range(EXPECTED_CONTROLS)]
        center, scale, _ = robust_center_scale(controls)
        strong = unique_gate(100.0, center, scale, controls)
        self.assertTrue(strong["target_residual_z"] >= RESIDUAL_Z_MIN)
        self.assertTrue(strong["rank"]["empirical_one_sided_rank_fraction"] <= RANK_MAX)
        self.assertTrue(strong["passed"])
        weak = unique_gate(center, center, scale, controls)
        self.assertFalse(weak["passed"])

    def test_fixed_gaussian_injection_is_positive_at_target(self):
        img = np.zeros((2048, 2048), dtype=float)
        out = add_fixed_gaussian(img, TARGET_X, 2038.7227543099102, 6.0313418764, 10.0)
        self.assertGreater(float(np.max(out)), 9.0)
        self.assertEqual(out.shape, img.shape)


if __name__ == "__main__":
    unittest.main()
