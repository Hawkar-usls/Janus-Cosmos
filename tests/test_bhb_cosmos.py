"""Fail-closed frozen controls for JANUS Cosmos BHB R1."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from decimal import Decimal, localcontext
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workspace"))
from janus_bhb import BHBError, evaluate_batch, evaluate_case, replay_upstream  # noqa: E402


def case(r: str, second: str = "in", context: str = "collapsing_black_hole",
         origin: str = "infinity", e1: str = "1", e2: str = "1") -> dict:
    return {
        "r_over_rs": r,
        "physical_context": context,
        "particle_1": {"direction": "in", "origin": "infinity", "specific_energy": e1},
        "particle_2": {"direction": second, "origin": origin, "specific_energy": e2},
    }


class BHBCosmosTests(unittest.TestCase):
    def test_frozen_meta_registry_replay(self) -> None:
        receipt = replay_upstream()
        self.assertEqual(receipt["checked_scenarios"], 14)
        self.assertEqual(receipt["status"], "PASS_FROZEN_R0_INDEPENDENT_REPLAY")

    def test_tampered_upstream_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as path:
            mutated = Path(path) / "invalid.json"
            raw = (ROOT / "data/bhb/BHB-R0-UPSTREAM-2026-10-17.json").read_text()
            mutated.write_text(raw.replace('"r_over_rs": "10"', '"r_over_rs": "11"', 1))
            with self.assertRaisesRegex(BHBError, "SHA256_MISMATCH"):
                replay_upstream(mutated)

    def test_infall_same_e_is_finite_at_horizon_limit(self) -> None:
        for r in ("10", "1.1", "1.000001", "1.000000000000000000000001"):
            result = evaluate_case(case(r))
            self.assertEqual(Decimal(result["E_cm_GeV"]), Decimal("1.876544"))
            self.assertEqual(result["trajectory_source_gate"], "IDEALIZED_GEODESIC_KINEMATICS_ADMISSIBLE")
            self.assertGreater(Decimal(result["particle_1"]["one_minus_abs_beta"]), 0)
            self.assertEqual(result["kinematics_cross_check"], "PASS")

    def test_head_on_analytic_formula_and_white_hole_gate(self) -> None:
        result = evaluate_case(case("1.000001", "out", "eternal_schwarzschild", "white_hole"))
        with localcontext() as ctx:
            ctx.prec = 70
            x = Decimal("1.000001")
            expected = 2 * Decimal("0.938272") / ((x - 1) / x).sqrt()
            self.assertLess(abs(Decimal(result["E_cm_GeV"]) - expected), Decimal("1e-20"))
        self.assertEqual(result["trajectory_source_gate"], "CONDITIONAL_ETERNAL_WHITE_HOLE_ONLY")
        self.assertEqual(result["escaped_energy_gate"], "OPEN_NO_RAY_TRACING_OR_EMISSION_MODEL")

    def test_collapsing_black_hole_cannot_promote_white_hole(self) -> None:
        result = evaluate_case(case("1.000001", "out", "collapsing_black_hole", "white_hole"))
        self.assertEqual(result["trajectory_source_gate"], "NOT_ESTABLISHED_OR_INCOMPATIBLE")
        self.assertIn("WHITE_HOLE_WORLDLINE_NOT_IN_COLLAPSING_BH_SPACETIME", result["trajectory_source_reasons"])

    def test_outgoing_local_production_creates_energy_debt(self) -> None:
        result = evaluate_case(case("1.01", "out", "collapsing_black_hole", "local_production", e2="1.2"))
        self.assertEqual(result["trajectory_source_gate"], "NOT_ESTABLISHED_OR_INCOMPATIBLE")
        self.assertIn("LOCAL_PRODUCTION_REQUIRES_INDEPENDENT_ENERGY_BUDGET_AND_CAUSAL_SOURCE", result["trajectory_source_reasons"])

    def test_unequal_infall_remains_finite(self) -> None:
        r1 = evaluate_case(case("1.000001", e2="2"))
        r2 = evaluate_case(case("1.000000000001", e2="2"))
        self.assertTrue(Decimal("1.9") < Decimal(r1["E_cm_GeV"]) < Decimal("2.1"))
        self.assertLess(abs(Decimal(r1["E_cm_GeV"]) - Decimal(r2["E_cm_GeV"])), Decimal("1e-5"))

    def test_invalid_inputs_rejected(self) -> None:
        for bad in [case("1"), case("0.9"), case("NaN"), case("Infinity")]:
            with self.assertRaises(BHBError):
                evaluate_case(bad)
        bad = case("2")
        bad["particle_1"]["rest_energy_GeV"] = -1
        with self.assertRaises(BHBError):
            evaluate_case(bad)
        bad = case("2")
        bad["particle_1"]["specific_energy"] = "0.1"
        with self.assertRaises(BHBError):
            evaluate_case(bad)
        bad = case("2")
        bad["spacetime"] = "kerr"
        with self.assertRaises(BHBError):
            evaluate_case(bad)

    def test_deterministic_batch_and_no_scientific_promotion(self) -> None:
        document = json.loads((ROOT / "data/bhb/BHB-COSMOS-R1-CASES-2026-10-17.json").read_text())
        a = evaluate_batch(document)
        b = evaluate_batch(document)
        self.assertEqual(a["result_sha256"], b["result_sha256"])
        self.assertEqual(len(a["cases"]), 5)
        self.assertEqual(a["gates"]["ALL_TRAJECTORIES_ADMISSIBLE"], "NOT_ESTABLISHED")
        self.assertEqual(a["gates"]["REAL_BHB_DETECTION"], "NOT_ESTABLISHED")
        self.assertEqual(a["gates"]["KERR_ROTATION"], "NOT_MODELED")
        self.assertTrue(all(row["observation_gate"] == "NO_OBSERVATION_INPUT" for row in a["cases"]))


if __name__ == "__main__":
    unittest.main()
