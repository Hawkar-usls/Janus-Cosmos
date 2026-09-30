import unittest

from experiments.luci.run_tachyon_star_t3c_tess import predicted_candidate, required_injected_snr


class T(unittest.TestCase):
    def test_trial_15_exact_floor(self):
        raw = {
            "contrast_z": -3.0210452396422274,
            "b_minus_a_sigma": -2.548671064219578,
            "b_minus_c_sigma": -3.540625574654188,
        }
        self.assertAlmostEqual(required_injected_snr(raw), 11.021045239642227, places=12)
        self.assertFalse(predicted_candidate(raw, 10.0))
        self.assertTrue(predicted_candidate(raw, 12.0))

    def test_trial_39_exact_floor(self):
        raw = {
            "contrast_z": -2.09740353666039,
            "b_minus_a_sigma": -1.481312020972129,
            "b_minus_c_sigma": -2.650206332926144,
        }
        self.assertAlmostEqual(required_injected_snr(raw), 10.09740353666039, places=12)
        self.assertFalse(predicted_candidate(raw, 10.0))
        self.assertTrue(predicted_candidate(raw, 12.0))

    def test_already_at_candidate_gate_requires_zero_extra(self):
        raw = {
            "contrast_z": 8.0,
            "b_minus_a_sigma": 4.0,
            "b_minus_c_sigma": 4.0,
        }
        self.assertEqual(required_injected_snr(raw), 0.0)
        self.assertTrue(predicted_candidate(raw, 0.0))

    def test_neighbor_requirement_can_dominate(self):
        raw = {
            "contrast_z": 7.5,
            "b_minus_a_sigma": -10.0,
            "b_minus_c_sigma": 3.0,
        }
        self.assertEqual(required_injected_snr(raw), 14.0)
        self.assertFalse(predicted_candidate(raw, 13.99))
        self.assertTrue(predicted_candidate(raw, 14.0))


if __name__ == "__main__":
    unittest.main()
