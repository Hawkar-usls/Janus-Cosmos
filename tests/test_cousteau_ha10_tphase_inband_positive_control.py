from __future__ import annotations

import importlib.util
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "workspace" / "cousteau_ha10_tphase_inband_positive_control.py"
spec = importlib.util.spec_from_file_location(
    "cousteau_ha10_tphase_inband_positive_control", MODULE
)
gate = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(gate)


class CousteauHA10InbandPositiveControlTests(unittest.TestCase):
    def test_git_blob_identity_is_content_bound(self):
        first = gate.git_blob_sha1_bytes(b"abc")
        second = gate.git_blob_sha1_bytes(b"abd")
        self.assertEqual(len(first), 40)
        self.assertNotEqual(first, second)

    def test_db_ratio(self):
        self.assertAlmostEqual(gate.db_ratio(2.0, 1.0), 10.0 * math.log10(2.0))
        self.assertTrue(math.isnan(gate.db_ratio(0.0, 1.0)))

    def test_integrated_band_power_uses_full_declared_band(self):
        frequencies = np.linspace(0.0, 100.0, 1001)
        psd = np.ones_like(frequencies)
        power, bins = gate.integrated_band_power(frequencies, psd, [10.0, 80.0])
        self.assertAlmostEqual(power, 70.0, places=6)
        self.assertGreater(bins, 2)

    def test_decision_blocks_when_too_few_complete_events(self):
        events = [
            {
                "stations": {
                    "a": {"data_status": "ANALYZED"},
                    "b": {"data_status": "BLOCKED_PAIR"},
                },
                "positive_control_replicated_both_stations": False,
            }
            for _ in range(8)
        ]
        verdict, complete, passed = gate.decide_control(events, 3)
        self.assertEqual(verdict, "BLOCKED_POSITIVE_CONTROL_DATA_ACCESS_OR_RESPONSE")
        self.assertEqual(complete, 0)
        self.assertEqual(passed, 0)

    def test_decision_fail_is_scientific_outcome_not_software_error(self):
        events = [
            {
                "stations": {
                    "a": {"data_status": "ANALYZED"},
                    "b": {"data_status": "ANALYZED"},
                },
                "positive_control_replicated_both_stations": index < 2,
            }
            for index in range(8)
        ]
        verdict, complete, passed = gate.decide_control(events, 3)
        self.assertEqual(verdict, "FAIL_HA10_INBAND_TPHASE_PIPELINE_CONTROL")
        self.assertEqual(complete, 8)
        self.assertEqual(passed, 2)

    def test_decision_pass_requires_three_cross_station_events(self):
        events = [
            {
                "stations": {
                    "a": {"data_status": "ANALYZED"},
                    "b": {"data_status": "ANALYZED"},
                },
                "positive_control_replicated_both_stations": index < 3,
            }
            for index in range(8)
        ]
        verdict, complete, passed = gate.decide_control(events, 3)
        self.assertEqual(verdict, "PASS_HA10_INBAND_TPHASE_PIPELINE_CONTROL")
        self.assertEqual(complete, 8)
        self.assertEqual(passed, 3)

    def test_runner_refuses_canonical_data_output(self):
        with self.assertRaisesRegex(
            RuntimeError, "CANONICAL_DATA_WRITE_FORBIDDEN_BY_POSITIVE_CONTROL_RUNNER"
        ):
            gate.ensure_noncanonical_output(gate.DATA / "forbidden.json")

    def test_runner_allows_ephemeral_artifact_output(self):
        with tempfile.TemporaryDirectory() as td:
            gate.ensure_noncanonical_output(Path(td) / "receipt.json")

    def test_frozen_contracts_bind_existing_negative_without_reinterpretation(self):
        protocol, windows, frozen = gate.verify_frozen_contracts()
        self.assertEqual(
            protocol["epistemic_position"]["authority_delta_for_119hz"], 0
        )
        self.assertFalse(
            protocol["control_frequency_contract"]["119hz_or_117_121hz_bins_may_be_used"]
        )
        self.assertEqual(
            windows["window_freeze_sha256"], gate.EXPECTED_WINDOW_FREEZE_SHA256
        )
        self.assertEqual(
            frozen["summary"]["verdict"], gate.EXPECTED_FROZEN_119_VERDICT
        )


if __name__ == "__main__":
    unittest.main()
