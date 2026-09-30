import inspect
import unittest

import numpy as np

from experiments.luci.run_tachyon_star_t3a_tess import (
    FORBIDDEN_VALUE_COLUMNS,
    choose_prepointed_triple,
    eligible_triples,
    read_time_quality_only,
    tesscut_url,
)


class T(unittest.TestCase):
    def test_eligible_200s(self):
        step = 200.0 / 86400.0
        t = np.array([100.0 + i * step for i in range(7)], dtype=float)
        q = np.zeros(7, dtype=np.int64)
        x = eligible_triples(t, q, 200.0)
        self.assertEqual(len(x), 5)
        self.assertAlmostEqual(x[0]["delta_pre_s"], 200.0, places=5)
        self.assertAlmostEqual(x[0]["delta_post_s"], 200.0, places=5)

    def test_quality_break_rejects_triples_touching_bad_row(self):
        step = 600.0 / 86400.0
        t = np.array([10.0 + i * step for i in range(7)], dtype=float)
        q = np.zeros(7, dtype=np.int64)
        q[3] = 1
        x = eligible_triples(t, q, 600.0)
        self.assertTrue(all(3 not in (r["a_row"], r["b_row"], r["c_row"]) for r in x))

    def test_midpoint_rule(self):
        step = 200.0 / 86400.0
        t = np.array([50.0 + i * step for i in range(9)], dtype=float)
        q = np.zeros(9, dtype=np.int64)
        c = np.arange(1000, 1009, dtype=np.int64)
        x = choose_prepointed_triple(t, q, c, 200.0)
        self.assertEqual(x["b_row"], 4)
        self.assertEqual(x["a_row"], 3)
        self.assertEqual(x["c_row"], 5)
        self.assertEqual(x["b_cadenceno"], 1004)

    def test_tesscut_url_freezes_sector_and_size(self):
        u = tesscut_url({"ra_deg": "12.3", "dec_deg": "45.6", "sector": "71"})
        self.assertIn("sector=71", u)
        self.assertIn("x=11", u)
        self.assertIn("y=11", u)

    def test_time_reader_source_has_no_flux_value_dereference(self):
        src = inspect.getsource(read_time_quality_only)
        for name in FORBIDDEN_VALUE_COLUMNS:
            self.assertNotIn(f'upper_map["{name}"]', src)
            self.assertNotIn(f"upper_map['{name}']", src)


if __name__ == "__main__":
    unittest.main()
