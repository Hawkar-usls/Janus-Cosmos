import unittest

from experiments.luci.run_tachyon_star_t2c_ztf import (
    RETAINED_FRACTION_MIN,
    UNRESOLVED_IDS,
    _product_urls,
    classify_trial,
    gaussian_retained_fraction,
)


class T(unittest.TestCase):
    def test_frozen_unresolved_set(self):
        self.assertEqual(len(UNRESOLVED_IDS), 8)
        self.assertEqual(len(set(UNRESOLVED_IDS)), 8)
        self.assertEqual(UNRESOLVED_IDS[0], "TS-T2B-ZTF-03")
        self.assertEqual(UNRESOLVED_IDS[-1], "TS-T2B-ZTF-40")

    def test_gaussian_retained_center(self):
        self.assertGreater(gaussian_retained_fraction((100, 100), 50.0, 50.0, 3.0), 0.999999)

    def test_gaussian_retained_edge_blocks(self):
        self.assertLess(gaussian_retained_fraction((100, 100), 0.0, 50.0, 3.0), RETAINED_FRACTION_MIN)

    def test_same_product_fallback_never_changes_epoch(self):
        row = {
            "a_filefracday": "20240101512345",
            "a_imgtypecode": "o",
            "field": "123",
            "ccdid": "4",
            "qid": "2",
            "filtercode": "zr",
            "ra_deg": "10.25",
            "dec_deg": "20.5",
        }
        urls = _product_urls(row, "a", "sciimg.fits")
        self.assertEqual([q[0] for q in urls], ["CUTOUT_360_ARCSEC", "SAME_PRODUCT_FULL_FRAME"])
        cut = urls[0][1].split("?", 1)[0]
        full = urls[1][1]
        self.assertEqual(cut, full)
        self.assertIn("ztf_20240101512345_000123_zr_c04_o_q2_sciimg.fits", full)

    def test_trial_classification(self):
        q = {"status": "QUALIFIED_ABSENCE"}
        p = {"status": "SOURCE_PRESENT"}
        x = {"status": "BLOCKED_NATIVE_PSF_REFERENCE"}
        self.assertEqual(classify_trial(q, q, q), "NO_ISOLATED_B_EVENT")
        self.assertEqual(classify_trial(q, p, q), "ISOLATED_B_L0")
        self.assertEqual(classify_trial(p, q, q), "NON_ISOLATED_SOURCE_PATTERN")
        self.assertEqual(classify_trial(q, x, q), "UNRESOLVED_TRIAL")


if __name__ == "__main__":
    unittest.main()
