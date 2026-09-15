from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "workspace" / "cousteau_ha10_n1_s2_january_0600_same_scalar_holdout.py"
spec = importlib.util.spec_from_file_location(
    "cousteau_ha10_n1_s2_january_0600_same_scalar_holdout", MODULE
)
gate = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(gate)


class January0600SameScalarTests(unittest.TestCase):
    def test_verdict_mapping(self):
        self.assertEqual(
            gate.map_verdict("SCALAR_COLLAPSE_HOLDOUT_PASS"),
            "JANUARY_0600_SAME_SCALAR_PASS",
        )
        self.assertEqual(
            gate.map_verdict("SCALAR_COLLAPSE_HOLDOUT_FAIL"),
            "JANUARY_0600_SAME_SCALAR_FAIL",
        )
        self.assertEqual(
            gate.map_verdict("BLOCKED_SCALAR_COLLAPSE_HOLDOUT_DATA_ACCESS"),
            "BLOCKED_JANUARY_0600_SAME_SCALAR_DATA_ACCESS",
        )

    def test_unknown_verdict_fails_closed(self):
        with self.assertRaisesRegex(RuntimeError, "UNEXPECTED_PRIOR_DECISION"):
            gate.map_verdict("UNKNOWN")

    def test_frozen_contract_binds_prior_fail_and_exact_scalar(self):
        protocol, jan1800, decscalar, frozen = gate.verify_frozen_contracts()
        self.assertEqual(
            jan1800["result"]["verdict"], gate.EXPECTED_JAN1800_VERDICT
        )
        self.assertEqual(
            decscalar["result"]["verdict"], gate.EXPECTED_DEC_SCALAR_VERDICT
        )
        self.assertEqual(
            frozen["summary"]["verdict"], gate.EXPECTED_FROZEN_119_VERDICT
        )
        self.assertEqual(protocol["frozen_scalar"]["source_value_db"], gate.EXPECTED_SCALAR_DB)
        self.assertFalse(protocol["frozen_scalar"]["holdout_refit_allowed"])
        self.assertFalse(protocol["frozen_scalar"]["alternate_scalar_search_allowed"])
        self.assertFalse(protocol["epistemic_position"]["may_erase_prior_january_1800_fail"])

    def test_january_0600_block_is_exact(self):
        protocol, _, _, _ = gate.verify_frozen_contracts()
        wc = protocol["holdout_window_contract"]
        self.assertEqual(wc["window_start_time_utc_each_day"], "06:00:00Z")
        self.assertEqual(len(wc["dates"]), 20)
        self.assertEqual(wc["dates"][0], "2015-01-01")
        self.assertEqual(wc["dates"][-1], "2015-01-20")
        self.assertTrue(wc["no_replacement_for_missing_windows"])
        self.assertTrue(wc["no_posthoc_date_exclusion"])

    def test_thresholds_unchanged(self):
        protocol, _, _, _ = gate.verify_frozen_contracts()
        t = protocol["classification_thresholds"]
        self.assertEqual(t["maximum_absolute_broadband_residual_median_db"], 2.0)
        self.assertEqual(t["maximum_broadband_residual_iqr_db"], 3.0)
        self.assertEqual(t["window_absolute_residual_tolerance_db"], 3.0)
        self.assertEqual(t["minimum_fraction_windows_within_tolerance"], 0.8)
        self.assertEqual(t["maximum_absolute_fixed_subband_residual_median_db"], 2.5)

    def test_runner_refuses_canonical_output(self):
        with self.assertRaisesRegex(
            RuntimeError, "CANONICAL_DATA_WRITE_FORBIDDEN_BY_JANUARY_0600_SAME_SCALAR"
        ):
            gate.ensure_noncanonical_output(gate.DATA / "forbidden.json")

    def test_runner_allows_ephemeral_output(self):
        with tempfile.TemporaryDirectory() as td:
            gate.ensure_noncanonical_output(Path(td) / "run.json")


if __name__ == "__main__":
    unittest.main()
