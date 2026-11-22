import unittest
import numpy as np
from retina_5d_mosaic import analyze_spectrogram, inspect_retina_video_tiles


class TestRetina5DMosaic(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(20261123)
        self.noise = self.rng.normal(10, 1, (32, 64))

    def line(self, data, center=30, amplitude=6):
        x = data.copy()
        x[:, center] += amplitude
        return x

    def test_native_signal_ranking_and_nonmutating(self):
        on = self.line(self.noise)
        before = on.copy()
        a = analyze_spectrogram(on, source_id="on", baseline_alert=True)
        n = analyze_spectrogram(self.noise)
        self.assertGreater(a["evidence"]["native_research_score"],
                           n["evidence"]["native_research_score"])
        self.assertTrue(np.array_equal(on, before))
        self.assertTrue(a["baseline_alert_after"])
        self.assertTrue(a["evidence"]["mirror_parity"])

    def test_shared_rfi_never_hard_vetoes(self):
        on = self.line(self.noise, amplitude=8)
        off = self.line(self.noise * 1.01, amplitude=8)
        a = analyze_spectrogram(on, off_power=off, source_id="src", reference_id="ref",
                                native_grid_verified=True, baseline_alert=True)
        self.assertTrue(a["evidence"]["shared_pattern_advisory_only"])
        self.assertTrue(a["baseline_alert_after"])
        self.assertFalse(a["production_adoption"])

    def test_missing_provenance_fails_open(self):
        a = analyze_spectrogram(self.line(self.noise),
                                off_power=self.line(self.noise), baseline_alert=True)
        self.assertFalse(a["evidence"]["caller_attested_distinct_native_scan_pair"])
        self.assertTrue(a["baseline_alert_after"])

    def test_invalid_grid_and_missing_data_rejected(self):
        with self.assertRaises(ValueError):
            analyze_spectrogram(self.noise, off_power=np.zeros((16, 64)))
        broken = self.noise.copy()
        broken[0, 0] = np.nan
        with self.assertRaises(ValueError):
            analyze_spectrogram(broken)

    def test_video_separate_from_radio(self):
        frames = np.zeros((3, 16, 16))
        frames[:, :4, :4] = 255
        out = inspect_retina_video_tiles(frames)
        self.assertEqual(out["class"], "VISUAL_PATTERN_QC_ONLY_NOT_SIGNAL_ORIGIN")
        self.assertTrue(out["mirror_inverse_map_parity"])
        self.assertFalse(out["mosaic_votes_independent"])

    def test_determinism_and_native_source_hash(self):
        one = analyze_spectrogram(self.noise)
        two = analyze_spectrogram(self.noise)
        self.assertEqual(one["origin_prime_evidence_sha256"],
                         two["origin_prime_evidence_sha256"])
        self.assertEqual(analyze_spectrogram(self.noise, prior_evidence_sha256=one["origin_prime_evidence_sha256"])["spiral_return"], "NO_INFORMATION_GAIN__STOP")
        self.assertEqual(one["evidence"]["source_sha256"],
                         two["evidence"]["source_sha256"])


if __name__ == "__main__":
    unittest.main()
