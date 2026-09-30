import inspect
import unittest

from experiments.luci.run_tachyon_star_q4a import (
    CANDIDATE_X,
    CANDIDATE_Y,
    CONTROL_X,
    PATCH_HALF,
    inside_patch,
    header_wcs,
)


class T(unittest.TestCase):
    def test_control_grid_is_fixed_and_avoids_candidate_neighborhood(self):
        self.assertEqual(len(CONTROL_X), 13)
        self.assertTrue(all(abs(x - CANDIDATE_X) >= 128.0 for x in CONTROL_X))

    def test_candidate_patch_fits_2048_geometry(self):
        self.assertTrue(inside_patch((2048, 2048), CANDIDATE_X, CANDIDATE_Y, PATCH_HALF))

    def test_all_controls_fit_2048_geometry(self):
        self.assertTrue(all(inside_patch((2048, 2048), x, CANDIDATE_Y, PATCH_HALF) for x in CONTROL_X))

    def test_header_reader_does_not_dereference_hdu_data(self):
        src = inspect.getsource(header_wcs)
        executable = "\n".join(
            line for line in src.splitlines()
            if not line.lstrip().startswith("#")
        )
        self.assertNotIn(".data", executable)
        self.assertIn("fits.getheader", executable)

    def test_no_target_recentering_constant_is_mutable(self):
        self.assertAlmostEqual(CANDIDATE_X, 1437.216755721343, places=12)
        self.assertAlmostEqual(CANDIDATE_Y, 2038.7227543099102, places=12)


if __name__ == "__main__":
    unittest.main()
