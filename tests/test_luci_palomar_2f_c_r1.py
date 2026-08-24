from __future__ import annotations

import unittest

import numpy as np

from janus_cosmos.luci_psf import _inject_gaussian, detect_psf_sources
from janus_cosmos.luci_psf_r1 import (
    estimate_native_psf_fwhm,
    measure_psf_at,
    psf_relative_injection_recovery_gate,
)


class LuciPalomar2FCR1Tests(unittest.TestCase):
    def test_local_measurement_is_not_censored_by_global_top_256(self):
        rng = np.random.default_rng(101)
        a = rng.normal(1000.0, 1.0, (512, 512))
        # 289 resolved sources: a global catalogue capped at 256 must truncate them.
        for y in range(32, 481, 28):
            for x in range(32, 481, 28):
                _inject_gaussian(a, float(y), float(x), 3.0, 18.0)
        global_cat = detect_psf_sources(a, max_sources=256)
        self.assertEqual(len(global_cat), 256)

        # Add another lower-SNR resolved source on a clean local patch. The R1
        # recovery measurement must classify it at the known coordinate even
        # though a global top-N catalogue is already saturated.
        y, x = 18.0, 250.0
        _inject_gaussian(a, y, x, 3.0, 8.0)
        local = measure_psf_at(a, y, x)
        self.assertIsNotNone(local)
        self.assertGreaterEqual(local.area_px, 5)

    def test_single_pixel_impulse_is_rejected_by_local_classifier(self):
        rng = np.random.default_rng(102)
        a = rng.normal(0.0, 1.0, (160, 160))
        a[80, 80] += 30.0
        self.assertIsNone(measure_psf_at(a, 80.0, 80.0))

    def test_psf_relative_gate_passes_resolved_synthetic_background(self):
        rng = np.random.default_rng(103)
        a = rng.normal(1000.0, 1.5, (320, 320))
        for y, x, peak in [
            (55, 60, 30), (65, 245, 28), (130, 100, 35), (155, 250, 26),
            (235, 70, 31), (250, 210, 34), (105, 190, 29), (210, 155, 32),
        ]:
            _inject_gaussian(a, float(y), float(x), 3.2, float(peak))
        fwhm, n = estimate_native_psf_fwhm(a)
        self.assertIsNotNone(fwhm)
        self.assertGreaterEqual(n, 2)
        gate = psf_relative_injection_recovery_gate(a, seed=104)
        self.assertTrue(gate["passed"], gate)
        self.assertGreaterEqual(gate["star_recovery_fraction_all"], 0.80)
        self.assertGreaterEqual(gate["star_recovery_fraction_snr_ge_8"], 0.90)
        self.assertLessEqual(gate["hot_pixel_acceptance_fraction"], 0.05)


if __name__ == "__main__":
    unittest.main()
