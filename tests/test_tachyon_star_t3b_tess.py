import unittest

import numpy as np

from experiments.luci.run_tachyon_star_t3b_tess import (
    aperture_series,
    candidate_metrics,
    classify_trial,
    isolated_contrast,
    robust_center_sigma,
)


class T(unittest.TestCase):
    def test_aperture_background_subtraction(self):
        cube = np.full((3, 11, 11), 100.0, dtype=float)
        bkg = np.full_like(cube, 10.0)
        cube[:, 4:7, 4:7] += 5.0
        s, meta = aperture_series(cube, bkg)
        self.assertTrue(np.allclose(s, 45.0))
        self.assertEqual(meta["target_pixels"], 9)
        self.assertGreaterEqual(meta["background_pixels"], 20)

    def test_isolated_contrast(self):
        s = np.array([1.0, 2.0, 7.0, 2.0, 1.0])
        self.assertEqual(isolated_contrast(s, 2), 5.0)

    def test_candidate_metrics(self):
        s = np.zeros(9, dtype=float)
        s[4] = 10.0
        m = candidate_metrics(s, 4, 0.0, 1.0)
        self.assertTrue(m["passed"])
        self.assertGreaterEqual(m["contrast_z"], 8.0)
        self.assertGreaterEqual(m["b_minus_a_sigma"], 4.0)
        self.assertGreaterEqual(m["b_minus_c_sigma"], 4.0)

    def test_robust_sigma(self):
        x = [(-1.0 if i % 2 else 1.0) * (1.0 + (i % 5) * 0.1) for i in range(40)]
        center, sigma, mode = robust_center_sigma(x)
        self.assertTrue(np.isfinite(center))
        self.assertGreater(sigma, 0)
        self.assertIn(mode, ("MAD", "STD_FALLBACK"))

    def test_null_requires_injection_recovery(self):
        rng = np.random.default_rng(12345)
        s = rng.normal(0.0, 1.0, 401)
        q = np.zeros(401, dtype=np.int64)
        b = 200
        # Keep the frozen raw B unremarkable.
        s[b-1:b+2] = 0.0
        r = classify_trial(s, q, b)
        self.assertFalse(r["raw"]["passed"])
        recovered = {float(x["injected_snr"]): bool(x["candidate_recovered"]) for x in r["injections"]}
        self.assertTrue(recovered[10.0])
        self.assertTrue(recovered[12.0])
        self.assertEqual(r["classification"], "QUALIFIED_NO_ISOLATED_B_EVENT")


if __name__ == "__main__":
    unittest.main()
