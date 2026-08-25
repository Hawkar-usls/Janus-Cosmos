from __future__ import annotations

import math
import tempfile
import unittest
from pathlib import Path

from workspace import cousteau_ha10_n1_s2_independent_scalar_collapse_holdout as holdout


class ScalarCollapseHoldoutTests(unittest.TestCase):
    def make_row(
        self,
        raw_residual: float,
        corrected_residual: float,
        *,
        subband_residual: float | None = None,
    ) -> dict:
        s = corrected_residual if subband_residual is None else subband_residual
        return {
            "broadband": {
                "raw_residual_db": raw_residual,
                "corrected_residual_db": corrected_residual,
            },
            "subbands": {
                "low": {
                    "band_hz": [10.0, 30.0],
                    "raw_residual_db": s,
                    "corrected_residual_db": s,
                },
                "mid": {
                    "band_hz": [30.0, 55.0],
                    "raw_residual_db": s,
                    "corrected_residual_db": s,
                },
                "high": {
                    "band_hz": [55.0, 80.0],
                    "raw_residual_db": s,
                    "corrected_residual_db": s,
                },
            },
        }

    def test_frozen_contracts_and_scalar_binding(self) -> None:
        protocol, temporal, frozen = holdout.verify_frozen_contracts()
        self.assertEqual(
            protocol["gate_id"],
            "COUSTEAU_HA10_N1_S2_INDEPENDENT_SCALAR_COLLAPSE_HOLDOUT_V1",
        )
        self.assertEqual(
            temporal["result"]["verdict"],
            "TEMPORALLY_STABLE_BROADBAND_H10S2_SCALE_OFFSET",
        )
        self.assertEqual(
            frozen["summary"]["verdict"],
            "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE",
        )
        self.assertEqual(
            protocol["frozen_scalar"]["source_value_db"],
            holdout.EXPECTED_SCALAR_DB,
        )
        expected_amp = 10.0 ** (-holdout.EXPECTED_SCALAR_DB / 20.0)
        expected_power = 10.0 ** (-holdout.EXPECTED_SCALAR_DB / 10.0)
        self.assertTrue(
            math.isclose(
                protocol["frozen_scalar"]["s2_waveform_amplitude_factor"],
                expected_amp,
                rel_tol=0.0,
                abs_tol=1e-15,
            )
        )
        self.assertTrue(
            math.isclose(
                protocol["frozen_scalar"]["s2_power_factor"],
                expected_power,
                rel_tol=0.0,
                abs_tol=1e-15,
            )
        )

    def test_pass_requires_broadband_time_and_subband_collapse(self) -> None:
        rows = [
            self.make_row(
                raw_residual=(-0.8 + 0.1 * (i % 10)),
                corrected_residual=(-0.7 + 0.1 * (i % 10)),
                subband_residual=0.4,
            )
            for i in range(18)
        ]
        verdict, broadband, subbands, gates = holdout.decide(
            rows,
            minimum_pairs=15,
            max_abs_median_db=2.0,
            max_iqr_db=3.0,
            tolerance_db=3.0,
            minimum_fraction=0.8,
            max_abs_subband_median_db=2.5,
        )
        self.assertEqual(verdict, "SCALAR_COLLAPSE_HOLDOUT_PASS")
        self.assertIsNotNone(broadband)
        self.assertIsNotNone(subbands)
        self.assertTrue(all(gates.values()))

    def test_subband_failure_prevents_pass(self) -> None:
        rows = [
            self.make_row(0.2, 0.3, subband_residual=3.1)
            for _ in range(18)
        ]
        verdict, _, _, gates = holdout.decide(
            rows,
            minimum_pairs=15,
            max_abs_median_db=2.0,
            max_iqr_db=3.0,
            tolerance_db=3.0,
            minimum_fraction=0.8,
            max_abs_subband_median_db=2.5,
        )
        self.assertEqual(verdict, "SCALAR_COLLAPSE_HOLDOUT_FAIL")
        self.assertFalse(gates["fixed_subbands"])

    def test_temporal_dispersion_failure_prevents_pass(self) -> None:
        residuals = [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0] * 2
        rows = [self.make_row(x, x, subband_residual=0.0) for x in residuals]
        verdict, _, _, gates = holdout.decide(
            rows,
            minimum_pairs=15,
            max_abs_median_db=2.0,
            max_iqr_db=3.0,
            tolerance_db=3.0,
            minimum_fraction=0.8,
            max_abs_subband_median_db=2.5,
        )
        self.assertEqual(verdict, "SCALAR_COLLAPSE_HOLDOUT_FAIL")
        self.assertFalse(gates["broadband_iqrs"])

    def test_blocked_below_minimum_pairs(self) -> None:
        rows = [self.make_row(0.0, 0.0) for _ in range(14)]
        verdict, broadband, subbands, gates = holdout.decide(
            rows,
            minimum_pairs=15,
            max_abs_median_db=2.0,
            max_iqr_db=3.0,
            tolerance_db=3.0,
            minimum_fraction=0.8,
            max_abs_subband_median_db=2.5,
        )
        self.assertEqual(verdict, "BLOCKED_SCALAR_COLLAPSE_HOLDOUT_DATA_ACCESS")
        self.assertIsNone(broadband)
        self.assertIsNone(subbands)
        self.assertFalse(gates["minimum_complete_pairs"])

    def test_runner_refuses_canonical_output(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "CANONICAL_DATA_WRITE_FORBIDDEN"):
            holdout.ensure_noncanonical_output(
                holdout.DATA / "forbidden-scalar-collapse-output.json"
            )

    def test_runner_allows_ephemeral_output(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            holdout.ensure_noncanonical_output(Path(td) / "receipt.json")


if __name__ == "__main__":
    unittest.main()
