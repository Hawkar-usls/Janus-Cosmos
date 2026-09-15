from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "workspace" / "cousteau_ha10_n1_s2_out_of_month_scalar_holdout.py"
spec = importlib.util.spec_from_file_location(
    "cousteau_ha10_n1_s2_out_of_month_scalar_holdout", MODULE
)
gate = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(gate)


class CousteauHA10OutOfMonthScalarHoldoutTests(unittest.TestCase):
    def test_verdict_mapping_preserves_semantics(self):
        self.assertEqual(
            gate.map_verdict("SCALAR_COLLAPSE_HOLDOUT_PASS"),
            "OUT_OF_MONTH_SCALAR_HOLDOUT_PASS",
        )
        self.assertEqual(
            gate.map_verdict("SCALAR_COLLAPSE_HOLDOUT_FAIL"),
            "OUT_OF_MONTH_SCALAR_HOLDOUT_FAIL",
        )
        self.assertEqual(
            gate.map_verdict("BLOCKED_SCALAR_COLLAPSE_HOLDOUT_DATA_ACCESS"),
            "BLOCKED_OUT_OF_MONTH_SCALAR_HOLDOUT_DATA_ACCESS",
        )

    def test_unknown_verdict_fails_closed(self):
        with self.assertRaisesRegex(RuntimeError, "UNEXPECTED_PRIOR_DECISION"):
            gate.map_verdict("MAYBE")

    def test_runner_refuses_canonical_output(self):
        with self.assertRaisesRegex(
            RuntimeError,
            "CANONICAL_DATA_WRITE_FORBIDDEN_BY_OUT_OF_MONTH_SCALAR_HOLDOUT",
        ):
            gate.ensure_noncanonical_output(gate.DATA / "forbidden.json")

    def test_runner_allows_ephemeral_output(self):
        with tempfile.TemporaryDirectory() as td:
            gate.ensure_noncanonical_output(Path(td) / "run.json")

    def test_frozen_contract_binds_exact_scalar_and_january_block(self):
        protocol, upstream, frozen = gate.verify_frozen_contracts()
        self.assertEqual(
            protocol["status"],
            "PREREGISTERED_BEFORE_JANUARY_HOLDOUT_WINDOWS_ARE_DOWNLOADED",
        )
        self.assertEqual(
            upstream["result"]["verdict"], gate.EXPECTED_PRIOR_VERDICT
        )
        self.assertEqual(
            frozen["summary"]["verdict"], gate.EXPECTED_FROZEN_119_VERDICT
        )
        self.assertEqual(
            protocol["frozen_scalar"]["source_value_db"], gate.EXPECTED_SCALAR_DB
        )
        self.assertFalse(protocol["frozen_scalar"]["holdout_refit_allowed"])
        self.assertFalse(protocol["frozen_scalar"]["alternate_scalar_search_allowed"])
        dates = protocol["holdout_window_contract"]["dates"]
        self.assertEqual(len(dates), 20)
        self.assertEqual(dates[0], "2015-01-01")
        self.assertEqual(dates[-1], "2015-01-20")
        self.assertTrue(all(date.startswith("2015-01-") for date in dates))

    def test_thresholds_are_identical_to_prior_holdout(self):
        protocol, _, _ = gate.verify_frozen_contracts()
        thresholds = protocol["classification_thresholds"]
        self.assertEqual(
            thresholds["maximum_absolute_broadband_residual_median_db"], 2.0
        )
        self.assertEqual(thresholds["maximum_broadband_residual_iqr_db"], 3.0)
        self.assertEqual(thresholds["window_absolute_residual_tolerance_db"], 3.0)
        self.assertEqual(thresholds["minimum_fraction_windows_within_tolerance"], 0.8)
        self.assertEqual(
            thresholds["maximum_absolute_fixed_subband_residual_median_db"], 2.5
        )

    def test_no_target_band_or_replacement_authority(self):
        protocol, _, _ = gate.verify_frozen_contracts()
        self.assertFalse(
            protocol["frequency_contract"]["119hz_or_117_121hz_bins_may_be_used"]
        )
        self.assertFalse(
            protocol["holdout_window_contract"]["no_replacement_for_missing_windows"]
            is False
        )
        self.assertTrue(protocol["holdout_window_contract"]["no_posthoc_date_exclusion"])
        self.assertEqual(protocol["epistemic_position"]["authority_delta_for_119hz"], 0)


if __name__ == "__main__":
    unittest.main()
