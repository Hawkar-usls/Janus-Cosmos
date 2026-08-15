from __future__ import annotations

import unittest

import numpy as np

from janus_cosmos.luci_psf import detect_psf_sources, injection_recovery_gate


def add_gaussian(a: np.ndarray, y: float, x: float, fwhm: float, peak: float) -> None:
    sig = fwhm / 2.354820045
    yy, xx = np.indices(a.shape, dtype=float)
    a += peak * np.exp(-0.5 * (((xx-x)/sig)**2 + ((yy-y)/sig)**2))


class LuciPalomar2FCTests(unittest.TestCase):
    def test_psf_detector_accepts_multipixel_stars_and_rejects_hot_pixels(self):
        rng = np.random.default_rng(7)
        a = rng.normal(0.0, 1.0, (160, 160))
        add_gaussian(a, 50, 50, 3.0, 20.0)
        add_gaussian(a, 105, 90, 2.5, 16.0)
        a[40, 120] += 35.0
        a[125, 35] += 30.0
        src = detect_psf_sources(a)
        self.assertTrue(any((q.x-50)**2 + (q.y-50)**2 < 4 for q in src))
        self.assertTrue(any((q.x-90)**2 + (q.y-105)**2 < 4 for q in src))
        self.assertFalse(any((q.x-120)**2 + (q.y-40)**2 < 4 for q in src))
        self.assertFalse(any((q.x-35)**2 + (q.y-125)**2 < 4 for q in src))
        self.assertTrue(all(q.area_px >= 5 and q.fwhm_minor_px >= 0.8 for q in src))

    def test_injection_recovery_gate_on_resolved_psf_grid(self):
        rng = np.random.default_rng(11)
        a = rng.normal(1000.0, 2.0, (256, 256))
        for y, x, p in [(60,60,35),(70,180,28),(150,90,31),(190,190,40),(130,210,25)]:
            add_gaussian(a, y, x, 3.0, p)
        g = injection_recovery_gate(
            a, seed=123, fwhm_grid_px=(2.0, 2.5, 3.0), snr_grid=(8.0, 12.0), replicates=3,
            min_all_recovery=0.80, min_high_snr_recovery=0.80, max_hot_pixel_acceptance=0.05,
        )
        self.assertTrue(g["passed"], g)
        self.assertLessEqual(g["hot_pixel_acceptance_fraction"], 0.05)

    def test_hot_pixel_only_field_does_not_create_psf_catalogue(self):
        rng = np.random.default_rng(17)
        a = rng.normal(0.0, 1.0, (128, 128))
        for y, x in [(30,30),(45,90),(96,55),(80,80)]:
            a[y, x] += 50.0
        self.assertEqual(detect_psf_sources(a), [])


if __name__ == "__main__":
    unittest.main()
