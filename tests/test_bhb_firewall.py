"""Fail-closed AMPS quantum-firewall hypothesis registry, no detection claims."""

from __future__ import annotations

import copy
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workspace"))
from janus_bhb import BHBError  # noqa: E402
from janus_bhb_firewall import DEFAULT_HYPOTHESES, DEFAULT_R2, evaluate_overlay  # noqa: E402


class BHBFirewallTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.spec = json.loads(DEFAULT_HYPOTHESES.read_text())
        cls.r2_bytes = DEFAULT_R2.read_bytes()
        cls.r2 = json.loads(cls.r2_bytes)

    def test_frozen_overlay_and_seven_cases(self):
        result = evaluate_overlay(self.spec, self.r2, self.r2_bytes)
        frozen = json.loads((ROOT / "data/bhb/BHB-FIREWALL-FW1-RESULT-2026-10-18.json").read_text())
        self.assertEqual(result, frozen)
        self.assertEqual(len(result["cases"]), 7)
        self.assertEqual(result["gates"]["CLASSICAL_R2_RECEIPT_INTEGRITY"], "PASS")

    def test_firewall_not_inferred_from_high_local_energy(self):
        result = evaluate_overlay(self.spec, self.r2, self.r2_bytes)
        near = next(c for c in result["cases"] if c["scenario_id"] == "R2-K1-CRITICAL-NEAR-HORIZON-LIMIT")
        self.assertGreater(float(near["classical_R2_local_E_cm_GeV"]), 100)
        self.assertIsNone(near["quantum_firewall_collision_E_cm_GeV"])
        self.assertIsNone(near["quantum_firewall_photon_spectrum"])
        self.assertEqual(near["firewall_presence"], "UNDETERMINED_NOT_MEASURED")

    def test_refuse_hypothesis_promotion(self):
        for field, value in (("physical_firewall_presence", True), ("scientific_status", "CONFIRMED")):
            mutated = copy.deepcopy(self.spec)
            mutated[field] = value
            with self.subTest(field=field), self.assertRaises(BHBError):
                evaluate_overlay(mutated, self.r2, self.r2_bytes)

    def test_refuse_missing_quantum_premises(self):
        mutated = copy.deepcopy(self.spec)
        mutated["premises"].pop()
        with self.assertRaises(BHBError):
            evaluate_overlay(mutated, self.r2, self.r2_bytes)
        mutated = copy.deepcopy(self.spec)
        mutated["applicability_conditions"] = ["EVAPORATING_BLACK_HOLE"]
        with self.assertRaises(BHBError):
            evaluate_overlay(mutated, self.r2, self.r2_bytes)

    def test_tampered_r2_fails_before_adding_firewall(self):
        corrupted = self.r2_bytes.replace(b'"spin_a_over_M": "0.9"',b'"spin_a_over_M": "0.8"',1)
        self.assertNotEqual(corrupted, self.r2_bytes)
        with self.assertRaises(BHBError):
            evaluate_overlay(self.spec, json.loads(corrupted), corrupted)

    def test_r2_detection_gate_cannot_be_promoted(self):
        changed = copy.deepcopy(self.r2)
        changed["gates"]["REAL_BHB_DETECTION"] = "CONFIRMED"
        with self.assertRaises(BHBError):
            evaluate_overlay(self.spec, changed, self.r2_bytes)


if __name__ == "__main__":
    unittest.main()
