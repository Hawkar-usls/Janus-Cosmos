from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import unittest

import numpy as np


MODULE_PATH = Path(__file__).resolve().parents[1] / "experiments" / "luci" / "run_palomar_star_morphology_transfer.py"
spec = importlib.util.spec_from_file_location("luci_palomar_transfer", MODULE_PATH)
mod = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(mod)


class LuciPalomarTransferTests(unittest.TestCase):
    def test_angular_separation_identity_and_degree(self):
        self.assertAlmostEqual(mod.angular_sep_deg(10.0, 20.0, 10.0, 20.0), 0.0, places=10)
        self.assertAlmostEqual(mod.angular_sep_deg(0.0, 0.0, 1.0, 0.0), 1.0, places=8)

    def test_detect_and_measure_synthetic_star_field(self):
        rng = np.random.default_rng(20260815)
        h = w = 384
        image = rng.normal(1000.0, 2.0, size=(h, w)).astype(float)
        yy, xx = np.indices((h, w), dtype=float)
        stars = []
        for iy in range(4):
            for ix in range(5):
                y = 45 + iy * 80
                x = 42 + ix * 70
                amp = 90.0 + 4.0 * (iy + ix)
                sigma = 1.7 + 0.04 * ((iy + ix) % 3)
                image += amp * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2.0 * sigma * sigma))
                stars.append((y, x))

        measured, qa = mod.detect_and_measure(image)
        self.assertGreaterEqual(len(measured), 12)
        self.assertGreaterEqual(qa["finite_fwhm_elongation_fraction"], 0.95)
        self.assertGreaterEqual(qa["local_reference_ge3_fraction"], 0.50)
        self.assertTrue(qa["frame_transfer_gate_pass"])
        med_elong = float(np.median([m["elongation"] for m in measured]))
        self.assertTrue(math.isfinite(med_elong))
        self.assertLess(med_elong, 1.5)

    def test_inverted_tail_polarity_is_recorded(self):
        rng = np.random.default_rng(7)
        h = w = 256
        image = rng.normal(500.0, 1.0, size=(h, w)).astype(float)
        yy, xx = np.indices((h, w), dtype=float)
        for iy in range(4):
            for ix in range(4):
                y = 35 + iy * 60
                x = 35 + ix * 60
                image -= 70.0 * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2.0 * 1.6**2))
        measured, qa = mod.detect_and_measure(image)
        self.assertEqual(qa["polarity"], "inverted_by_tail_test")
        self.assertGreaterEqual(len(measured), 12)


if __name__ == "__main__":
    unittest.main()
