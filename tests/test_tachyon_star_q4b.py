import unittest

import numpy as np

from experiments.luci.run_tachyon_star_q4b import (
    B_Z_MIN,
    COMPAT_Z_MAX,
    classify,
    empirical_rank_p,
    forced_gaussian_plane,
    temporal_contrast,
)


class T(unittest.TestCase):
    def test_forced_fit_recovers_fixed_gaussian_amplitude(self):
        h = w = 41
        y0 = x0 = 20.35
        yy, xx = np.indices((h, w), dtype=float)
        fwhm = 6.0
        sig = fwhm / 2.354820045
        g = np.exp(-0.5 * ((xx - x0) ** 2 + (yy - y0) ** 2) / (sig * sig))
        image = 100.0 + 0.2 * (xx - x0) - 0.1 * (yy - y0) + 50.0 * g
        rng = np.random.default_rng(1234)
        image += rng.normal(0.0, 1.0, image.shape)
        q = forced_gaussian_plane(image, x0, y0, fwhm)
        self.assertAlmostEqual(q["amplitude"], 50.0, delta=1.5)
        self.assertGreater(q["forced_z"], 10.0)

    def test_temporal_contrast(self):
        a = {"amplitude": 1.0}
        b = {"amplitude": 7.0}
        c = {"amplitude": 3.0}
        self.assertEqual(temporal_contrast(a, b, c), 5.0)

    def test_rank_p_has_frozen_resolution(self):
        controls = [float(i) for i in range(13)]
        q = empirical_rank_p(20.0, controls)
        self.assertEqual(q["one_sided_rank_p"], 1.0 / 14.0)
        self.assertEqual(q["rank_descending_1_is_largest"], 1)

    def test_classify_sky_persistence(self):
        b = {"amplitude": 10.0, "amplitude_sigma": 1.0, "forced_z": 10.0}
        sky = {
            "A_BEFORE": {"amplitude": 9.5, "amplitude_sigma": 1.0},
            "B_CANDIDATE": b,
            "C_AFTER": {"amplitude": 10.5, "amplitude_sigma": 1.0},
        }
        det = {
            "A_BEFORE": {"amplitude": 0.0, "amplitude_sigma": 1.0},
            "B_CANDIDATE": b,
            "C_AFTER": {"amplitude": 0.0, "amplitude_sigma": 1.0},
        }
        q = classify(sky, det)
        self.assertEqual(q["classification"], "EVIDENCE_FAVORS_SKY_FIXED_PERSISTENCE")
        self.assertTrue(q["sky_persistence_compatible"])
        self.assertFalse(q["detector_persistence_compatible"])

    def test_classify_one_frame_indeterminate(self):
        b = {"amplitude": 10.0, "amplitude_sigma": 1.0, "forced_z": 10.0}
        blank = {"amplitude": 0.0, "amplitude_sigma": 1.0}
        sky = {"A_BEFORE": blank, "B_CANDIDATE": b, "C_AFTER": blank}
        det = {"A_BEFORE": blank, "B_CANDIDATE": b, "C_AFTER": blank}
        q = classify(sky, det)
        self.assertEqual(q["classification"], "ONE_FRAME_EXCESS_ORIGIN_INDETERMINATE")

    def test_b_not_confirmed(self):
        b = {"amplitude": 1.0, "amplitude_sigma": 1.0, "forced_z": B_Z_MIN - 0.1}
        x = {"amplitude": 0.0, "amplitude_sigma": 1.0}
        q = classify({"A_BEFORE": x, "B_CANDIDATE": b, "C_AFTER": x}, {"A_BEFORE": x, "B_CANDIDATE": b, "C_AFTER": x})
        self.assertEqual(q["classification"], "B_FIXED_EXCESS_NOT_CONFIRMED")


if __name__ == "__main__":
    unittest.main()
