from __future__ import annotations

import importlib.util
import math
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "workspace" / "cousteau_ha10_n1_s2_independent_noise_asymmetry_replication.py"
spec = importlib.util.spec_from_file_location(
    "cousteau_ha10_n1_s2_independent_noise_asymmetry_replication", MODULE
)
gate = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(gate)


def pair(raw_db: float, corrected_db: float):
    return {
        "s2_minus_n1_raw_power_db": raw_db,
        "s2_minus_n1_corrected_power_db": corrected_db,
    }


class CousteauHA10IndependentNoiseAsymmetryTests(unittest.TestCase):
    def test_ratio_db(self):
        self.assertAlmostEqual(gate.ratio_db(4.0, 1.0), 10.0 * math.log10(4.0))
        self.assertTrue(math.isnan(gate.ratio_db(0.0, 1.0)))

    def test_blocked_below_minimum_pairs(self):
        verdict, aggregate = gate.decide(
            [pair(12.0, 12.0)] * 14,
            minimum_pairs=15,
            threshold_db=6.0,
        )
        self.assertEqual(verdict, "BLOCKED_INDEPENDENT_NOISE_REPLICATION_DATA_ACCESS")
        self.assertEqual(aggregate["complete_paired_windows"], 14)

    def test_replication_requires_both_raw_and_corrected_effects(self):
        verdict, aggregate = gate.decide(
            [pair(7.0, 6.5)] * 15,
            minimum_pairs=15,
            threshold_db=6.0,
        )
        self.assertEqual(
            verdict, "REPLICATED_LARGE_H10S2_BASELINE_POWER_ASYMMETRY"
        )
        self.assertGreaterEqual(
            aggregate["median_per_window_raw_s2_minus_n1_db"], 6.0
        )
        self.assertGreaterEqual(
            aggregate["median_per_window_corrected_s2_minus_n1_db"], 6.0
        )

    def test_raw_only_effect_does_not_replicate(self):
        verdict, _ = gate.decide(
            [pair(9.0, 5.9)] * 15,
            minimum_pairs=15,
            threshold_db=6.0,
        )
        self.assertEqual(
            verdict, "NOT_REPLICATED_LARGE_H10S2_BASELINE_POWER_ASYMMETRY"
        )

    def test_runner_refuses_canonical_data_output(self):
        with self.assertRaisesRegex(
            RuntimeError,
            "CANONICAL_DATA_WRITE_FORBIDDEN_BY_NOISE_REPLICATION_RUNNER",
        ):
            gate.ensure_noncanonical_output(gate.DATA / "forbidden.json")

    def test_runner_allows_ephemeral_output(self):
        with tempfile.TemporaryDirectory() as td:
            gate.ensure_noncanonical_output(Path(td) / "run.json")

    def test_frozen_contracts_bind_failed_control_and_119_negative(self):
        protocol, summary, frozen = gate.verify_frozen_contracts()
        self.assertEqual(
            protocol["status"],
            "PREREGISTERED_BEFORE_INDEPENDENT_NOISE_WINDOWS_ARE_DOWNLOADED",
        )
        self.assertEqual(
            summary["control_summary"]["verdict"],
            gate.EXPECTED_POSITIVE_CONTROL_VERDICT,
        )
        self.assertEqual(
            frozen["summary"]["verdict"], gate.EXPECTED_FROZEN_119_VERDICT
        )
        self.assertEqual(
            protocol["epistemic_position"]["authority_delta_for_119hz"], 0
        )
        self.assertEqual(
            len(protocol["independent_window_contract"]["dates"]), 20
        )


if __name__ == "__main__":
    unittest.main()
