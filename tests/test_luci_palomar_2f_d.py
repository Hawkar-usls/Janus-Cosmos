from __future__ import annotations

import math
import unittest

import numpy as np

from experiments.luci.run_palomar_2f_d import (
    archive_first_spatial_crossmatch,
    counterpart_with_matched_controls,
    metadata_wcs_contains,
)
from janus_cosmos.luci_psf import _inject_gaussian


class LuciPalomar2FDTests(unittest.TestCase):
    def test_archive_first_crossmatch_handles_ra_wrap(self):
        sources = [
            {"src_id": "A", "ra_deg": 359.99, "dec_deg": 0.0},
            {"src_id": "B", "ra_deg": 20.0, "dec_deg": 20.0},
        ]
        inventory = [{
            "file_name": "f.fits",
            "file_url": "https://example.invalid/f.fits",
            "instrument": "LUCI2",
            "target": "x",
            "filters": "J",
            "date_obs": "2020-01-01",
            "crval1": 0.0,
            "crval2": 0.0,
            "naxis1": 2048.0,
            "naxis2": 2048.0,
            "pixscale": 0.25,
            "cd1_1": float("nan"),
            "cd1_2": float("nan"),
            "cd2_1": float("nan"),
            "cd2_2": float("nan"),
        }]
        pairs, meta = archive_first_spatial_crossmatch(sources, inventory)
        self.assertEqual([p["src_id"] for p in pairs], ["A"])
        self.assertEqual(meta["coarse_pair_count"], 1)

    def test_metadata_wcs_contains_tan_field(self):
        frame = {
            "crval1": 10.0, "crval2": 20.0,
            "crpix1": 50.5, "crpix2": 50.5,
            "cd1_1": -1.0/3600.0, "cd1_2": 0.0,
            "cd2_1": 0.0, "cd2_2": 1.0/3600.0,
            "ctype1": "RA---TAN", "ctype2": "DEC--TAN",
            "naxis1": 100.0, "naxis2": 100.0,
        }
        inside, wm = metadata_wcs_contains(frame, 10.0, 20.0)
        self.assertTrue(inside)
        self.assertTrue(math.isfinite(wm["x"]))
        outside, _ = metadata_wcs_contains(frame, 11.0, 20.0)
        self.assertFalse(outside)

    def test_local_r1_counterpart_rejects_hot_pixel(self):
        rng = np.random.default_rng(7)
        image = rng.normal(0.0, 1.0, (128, 128))
        image[64, 64] += 20.0
        result = counterpart_with_matched_controls(image, 64.0, 64.0)
        self.assertFalse(result["counterpart_present"])

    def test_local_r1_counterpart_accepts_resolved_psf(self):
        rng = np.random.default_rng(11)
        image = rng.normal(0.0, 1.0, (256, 256))
        _inject_gaussian(image, 128.0, 128.0, 4.0, 20.0)
        for y in range(30, 231, 35):
            for x in range(30, 231, 35):
                if abs(x-128) < 15 and abs(y-128) < 15:
                    continue
                _inject_gaussian(image, float(y), float(x), 4.1, 15.0)
        result = counterpart_with_matched_controls(image, 128.0, 128.0)
        self.assertTrue(result["counterpart_present"])
        self.assertIn(result["morphology_status"], {
            "MATCHED_LOCAL_CONTROL_COMPARISON_AVAILABLE",
            "INSUFFICIENT_MATCHED_LOCAL_CONTROLS",
        })


if __name__ == "__main__":
    unittest.main()
