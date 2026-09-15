from __future__ import annotations

import unittest

import numpy as np

from experiments.luci.run_palomar_2f_e import (
    deterministic_header_shards,
    exposure_admissible,
    local_coordinate_injection,
    standard_broadband,
    temporal_representatives,
)
from janus_cosmos.luci_psf import _inject_gaussian


class LuciPalomar2FETests(unittest.TestCase):
    def test_exposure_admissibility_is_literal_and_conservative(self):
        self.assertFalse(exposure_admissible({"filters": "blind blind"}))
        self.assertFalse(exposure_admissible({"filters": "PV lens J"}))
        self.assertFalse(exposure_admissible({"filters": "clear pv LENS H"}))
        self.assertTrue(exposure_admissible({"filters": "clear Ks"}))
        self.assertTrue(exposure_admissible({"filters": "Br_gam clear"}))

    def test_standard_broadband_contract(self):
        for s in ("clear z", "clear J", "clear H", "clear K", "clear Ks"):
            self.assertTrue(standard_broadband(s))
        self.assertFalse(standard_broadband("Br_gam Ks"))
        self.assertFalse(standard_broadband("J_high J"))

    def test_header_shards_preserve_all_files_and_old_cap(self):
        rows = [
            {"file_name": f"f{i:03d}.fits.gz"}
            for i in range(403)
        ]
        shards = deterministic_header_shards(rows)
        self.assertEqual([len(x) for x in shards], [225, 178])
        self.assertLessEqual(max(map(len, shards)), 250)
        flat = [x for shard in shards for x in shard]
        self.assertEqual(flat, sorted({r["file_name"] for r in rows}))

    def test_temporal_selection_prefers_broadband_then_earliest_median_latest(self):
        rows = []
        for i in range(7):
            rows.append({
                "src_id": "S1", "file_name": f"b{i}.fits.gz", "date_obs": f"2020-01-{i+1:02d}",
                "filters": "clear J", "ra_deg": 1.0, "dec_deg": 2.0,
            })
        rows.append({
            "src_id": "S1", "file_name": "narrow.fits.gz", "date_obs": "2010-01-01",
            "filters": "Br_gam clear", "ra_deg": 1.0, "dec_deg": 2.0,
        })
        got = temporal_representatives(rows)
        self.assertEqual({r["file_name"] for r in got}, {"b0.fits.gz", "b3.fits.gz", "b6.fits.gz"})
        self.assertTrue(all(r["selection_pool"] == "STANDARD_BROADBAND" for r in got))

    def test_local_coordinate_injection_recovers_resolved_psf(self):
        rng = np.random.default_rng(12)
        image = rng.normal(0.0, 1.0, (160, 160))
        # Add reference stars elsewhere so the local test resembles a valid frame.
        for y, x in ((30,30),(30,120),(120,30),(120,120)):
            _inject_gaussian(image, float(y), float(x), 4.0, 15.0)
        got = local_coordinate_injection(image, 80.0, 80.0, 4.0)
        self.assertTrue(got["passed"])
        self.assertEqual([t["snr"] for t in got["trials"]], [8.0, 12.0])

    def test_local_coordinate_injection_blocks_edge(self):
        image = np.zeros((100, 100), dtype=float)
        image += np.random.default_rng(1).normal(0.0, 1.0, image.shape)
        got = local_coordinate_injection(image, 3.0, 3.0, 4.0)
        self.assertFalse(got["passed"])
        self.assertEqual(got["reason"], "COORDINATE_TOO_CLOSE_TO_EDGE")


if __name__ == "__main__":
    unittest.main()
