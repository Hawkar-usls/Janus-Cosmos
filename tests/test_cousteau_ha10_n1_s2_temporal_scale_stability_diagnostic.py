from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "workspace" / "cousteau_ha10_n1_s2_temporal_scale_stability_diagnostic.py"
spec = importlib.util.spec_from_file_location("temporal_scale_gate", MODULE)
gate = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(gate)


def row(raw: float, corrected: float):
    return {
        "raw_s2_minus_n1_db": raw,
        "corrected_s2_minus_n1_db": corrected,
    }


class TemporalScaleStabilityTests(unittest.TestCase):
    def test_blocked_below_minimum_pairs(self):
        verdict, stats = gate.decide(
            [row(12.0, 12.1)] * 14,
            minimum_pairs=15,
            floor_db=6.0,
            max_iqr_db=3.0,
            minimum_fraction=0.8,
        )
        self.assertEqual(verdict, "BLOCKED_TEMPORAL_SCALE_STABILITY_DATA_ACCESS")
        self.assertEqual(stats["complete_paired_windows"], 14)

    def test_stable_large_offset(self):
        rows = [row(12.0 + (i % 3) * 0.4, 12.1 + (i % 3) * 0.4) for i in range(15)]
        verdict, stats = gate.decide(
            rows,
            minimum_pairs=15,
            floor_db=6.0,
            max_iqr_db=3.0,
            minimum_fraction=0.8,
        )
        self.assertEqual(verdict, "TEMPORALLY_STABLE_BROADBAND_H10S2_SCALE_OFFSET")
        self.assertTrue(stats["medians_ge_floor"])
        self.assertTrue(stats["iqrs_le_maximum"])
        self.assertTrue(stats["fractions_ge_minimum"])

    def test_large_iqr_is_variable(self):
        rows = [row(v, v + 0.1) for v in [7, 7, 7, 7, 7, 12, 12, 12, 15, 15, 15, 18, 18, 18, 18]]
        verdict, stats = gate.decide(
            rows,
            minimum_pairs=15,
            floor_db=6.0,
            max_iqr_db=3.0,
            minimum_fraction=0.8,
        )
        self.assertEqual(verdict, "TEMPORALLY_VARIABLE_H10S2_SCALE_OFFSET")
        self.assertFalse(stats["iqrs_le_maximum"])

    def test_too_few_windows_above_floor_is_variable(self):
        rows = [row(12.0, 12.1)] * 11 + [row(2.0, 2.1)] * 4
        verdict, stats = gate.decide(
            rows,
            minimum_pairs=15,
            floor_db=6.0,
            max_iqr_db=20.0,
            minimum_fraction=0.8,
        )
        self.assertEqual(verdict, "TEMPORALLY_VARIABLE_H10S2_SCALE_OFFSET")
        self.assertFalse(stats["fractions_ge_minimum"])

    def test_runner_refuses_canonical_output(self):
        with self.assertRaisesRegex(RuntimeError, "CANONICAL_DATA_WRITE_FORBIDDEN"):
            gate.ensure_noncanonical_output(gate.DATA / "forbidden.json")

    def test_runner_allows_ephemeral_output(self):
        with tempfile.TemporaryDirectory() as td:
            gate.ensure_noncanonical_output(Path(td) / "run.json")

    def test_frozen_contracts(self):
        protocol, shape, frozen = gate.verify_frozen_contracts()
        self.assertEqual(protocol["gate_id"], "COUSTEAU_HA10_N1_S2_TEMPORAL_SCALE_STABILITY_DIAGNOSTIC_V1")
        self.assertEqual(shape["result"]["verdict"], gate.EXPECTED_SHAPE_VERDICT)
        self.assertEqual(frozen["summary"]["verdict"], gate.EXPECTED_FROZEN_119_VERDICT)
        self.assertEqual(protocol["new_window_contract"]["window_start_time_utc_each_day"], "06:00:00Z")
        self.assertEqual(protocol["epistemic_position"]["authority_delta_for_119hz"], 0)


if __name__ == "__main__":
    unittest.main()
