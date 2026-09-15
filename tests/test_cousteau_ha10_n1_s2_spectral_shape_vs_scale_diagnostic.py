from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "workspace" / "cousteau_ha10_n1_s2_spectral_shape_vs_scale_diagnostic.py"
spec = importlib.util.spec_from_file_location("spectral_shape_gate", MODULE)
gate = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(gate)

SUBBANDS = ["low", "mid", "high"]


def pair(raw, corrected):
    return {
        "subband_ratios_db": {
            name: {
                "raw_s2_minus_n1_db": raw[index],
                "corrected_s2_minus_n1_db": corrected[index],
            }
            for index, name in enumerate(SUBBANDS)
        }
    }


class SpectralShapeTests(unittest.TestCase):
    def test_blocked_below_minimum_pairs(self):
        verdict, aggregate = gate.decide(
            [pair((10.0, 10.5, 11.0), (10.1, 10.6, 11.1))] * 14,
            subband_names=SUBBANDS,
            minimum_pairs=15,
            floor_db=6.0,
            max_spread_db=3.0,
        )
        self.assertEqual(verdict, "BLOCKED_SPECTRAL_SHAPE_DIAGNOSTIC_DATA_ACCESS")
        self.assertEqual(aggregate["complete_paired_windows"], 14)

    def test_broadband_scale_like(self):
        verdict, aggregate = gate.decide(
            [pair((12.0, 13.0, 14.0), (12.2, 13.2, 14.2))] * 15,
            subband_names=SUBBANDS,
            minimum_pairs=15,
            floor_db=6.0,
            max_spread_db=3.0,
        )
        self.assertEqual(verdict, "BROADBAND_SCALE_LIKE_H10S2_ASYMMETRY")
        self.assertLessEqual(aggregate["raw_median_ratio_spread_db"], 3.0)
        self.assertTrue(aggregate["all_raw_and_corrected_subband_medians_ge_floor"])

    def test_large_spread_is_nonuniform(self):
        verdict, aggregate = gate.decide(
            [pair((8.0, 12.0, 15.0), (8.2, 12.2, 15.2))] * 15,
            subband_names=SUBBANDS,
            minimum_pairs=15,
            floor_db=6.0,
            max_spread_db=3.0,
        )
        self.assertEqual(verdict, "FREQUENCY_SELECTIVE_OR_NONUNIFORM_H10S2_ASYMMETRY")
        self.assertGreater(aggregate["raw_median_ratio_spread_db"], 3.0)

    def test_subband_below_floor_is_nonuniform(self):
        verdict, aggregate = gate.decide(
            [pair((5.9, 12.0, 12.5), (6.2, 12.1, 12.6))] * 15,
            subband_names=SUBBANDS,
            minimum_pairs=15,
            floor_db=6.0,
            max_spread_db=7.0,
        )
        self.assertEqual(verdict, "FREQUENCY_SELECTIVE_OR_NONUNIFORM_H10S2_ASYMMETRY")
        self.assertFalse(aggregate["all_raw_and_corrected_subband_medians_ge_floor"])

    def test_runner_refuses_canonical_output(self):
        with self.assertRaisesRegex(RuntimeError, "CANONICAL_DATA_WRITE_FORBIDDEN"):
            gate.ensure_noncanonical_output(gate.DATA / "forbidden.json")

    def test_runner_allows_ephemeral_output(self):
        with tempfile.TemporaryDirectory() as td:
            gate.ensure_noncanonical_output(Path(td) / "run.json")

    def test_frozen_contracts(self):
        protocol, replication, frozen = gate.verify_frozen_contracts()
        self.assertEqual(protocol["gate_id"], "COUSTEAU_HA10_N1_S2_SPECTRAL_SHAPE_VS_SCALE_DIAGNOSTIC_V1")
        self.assertEqual(replication["result"]["verdict"], gate.EXPECTED_REPLICATION_VERDICT)
        self.assertEqual(frozen["summary"]["verdict"], gate.EXPECTED_FROZEN_119_VERDICT)
        self.assertEqual(protocol["epistemic_position"]["authority_delta_for_119hz"], 0)


if __name__ == "__main__":
    unittest.main()
