"""Video-inspired JANUS Cosmos tests: reproducible math, negative science gates."""

from __future__ import annotations

import copy
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workspace"))
from janus_cosmos_simulation_audit import (  # noqa: E402
    AuditError, C_SI, FROZEN_R3A_PATH, MODELS, evaluate_batch,
    evaluate_case, verify_frozen_r3a,
)

CASES = ROOT / "data/cosmos/COSMOS-ALI-SIMULATION-CLAIMS-2026-10-18.json"
RESULT = ROOT / "data/cosmos/COSMOS-ALI-SIMULATION-AUDIT-RESULT-2026-10-18.json"


def read_cases():
    return json.loads(CASES.read_text(encoding="utf-8"))


class ALIVideoSimulationAuditTests(unittest.TestCase):
    def test_six_models_and_18_claims_pin_to_receipt(self):
        result = evaluate_batch(read_cases())
        self.assertEqual(result, json.loads(RESULT.read_text(encoding="utf-8")))
        self.assertEqual(result["case_count"], 16)
        self.assertEqual(result["claims_count"], 18)
        self.assertEqual(set(result["models"]), MODELS)
        self.assertEqual(result["gates"]["SIMULATION_HYPOTHESIS_OBSERVED"], "NOT_ESTABLISHED")
        self.assertEqual(result["gates"]["AMPS_FIREWALL_PHYSICAL_EXISTENCE"], "OPEN_UNCONFIRMED")
        self.assertEqual(result["gates"]["BHB_EXPLOSION_OBSERVED"], "NOT_ESTABLISHED")

    def test_lorentz_gamma_and_length_at_half_c(self):
        row = evaluate_case({"case_id":"half", "model":"lorentz", "speed_over_c":"0.5"})["outputs"]
        self.assertAlmostEqual(float(row["gamma"]), 1.1547005383792515)
        self.assertAlmostEqual(float(row["contracted_length_over_rest_length"]), .8660254037844386)
        self.assertEqual(row["transformed_signal_speed_over_c"], "1")
        self.assertEqual(row["exact_si_light_speed_m_per_s"], C_SI)

    def test_signal_c_is_invariant_for_boosts_both_directions(self):
        for beta in ("0", "0.5", "-0.85", "0.999999999"):
            with self.subTest(beta=beta):
                row = evaluate_case({"case_id":"light", "model":"lorentz", "speed_over_c":beta})["outputs"]
                self.assertEqual(Decimal(row["transformed_signal_speed_over_c"]), Decimal(1))
                self.assertGreaterEqual(Decimal(row["gamma"]), Decimal(1))

    def test_superluminal_and_nan_and_boolean_rejected(self):
        for bad in ("1", "-1", "1.00000001", "NaN", "Infinity", True, 0.5):
            with self.subTest(bad=bad), self.assertRaises(AuditError):
                evaluate_case({"case_id":"bad", "model":"lorentz", "speed_over_c":bad})

    def test_horizon_massive_stays_subluminal(self):
        near = evaluate_case({"case_id":"n", "model":"schwarzschild", "r_over_rs":"1.000000000000000000000001"})["outputs"]
        self.assertLess(Decimal(near["infall_speed_relative_to_static_observer_over_c"]), 1)
        self.assertGreater(Decimal(near["one_minus_infall_speed_over_c"]), 0)
        self.assertEqual(near["outgoing_radial_photon_local_speed_over_c"], "1")
        self.assertEqual(near["firewall_existence"], "OPEN_UNCONFIRMED")

    def test_at_or_below_horizon_rejected(self):
        for x in ("1", "0.9999999", "0", "-1"):
            with self.subTest(x=x), self.assertRaises(AuditError):
                evaluate_case({"case_id":"bad", "model":"schwarzschild", "r_over_rs":x})

    def test_area_entropy_quadruples_when_radius_doubles(self):
        def entropy(radius):
            return Decimal(evaluate_case({"case_id":str(radius), "model":"bh_area_entropy",
                                          "schwarzschild_radius_over_planck_length":str(radius)})["outputs"]["bekenstein_hawking_entropy_over_kB"])
        self.assertEqual(entropy(20) / entropy(10), 4)

    def test_compact_spatially_flat_torus_is_counterexample(self):
        row = evaluate_case({"case_id":"top", "model":"flat_torus", "Lx":"2", "Ly":"3", "Lz":"5"})["outputs"]
        self.assertTrue(row["compact"])
        self.assertFalse(row["has_boundary"])
        self.assertEqual(row["volume"], "30")
        self.assertEqual(row["intrinsic_spatial_curvature"], "0")

    def test_lattice_dispersion_departs_from_continuum(self):
        def vg(theta):
            return Decimal(evaluate_case({"case_id":"k", "model":"lattice_wave", "wave_number_times_spacing":theta})["outputs"]["group_velocity_over_c"])
        low, mid = vg("0.001"), vg("0.5")
        self.assertLess(abs(low-1), Decimal("0.000001"))
        self.assertGreater(low, mid)
        self.assertLess(abs(vg("3.141592653589793238462643383279502884197169399375105820974944592307816406286")), Decimal("1e-65"))

    def test_lattice_outside_first_brillouin_half_zone_rejected(self):
        for x in ("0", "-0.1", "3.2"):
            with self.subTest(x=x), self.assertRaises(AuditError):
                evaluate_case({"case_id":"bad", "model":"lattice_wave", "wave_number_times_spacing":x})

    def test_bostrom_count_is_conditional_not_posterior(self):
        def count(sim):
            return evaluate_case({"case_id":"n", "model":"conditional_observer_count",
                                  "non_simulated_observers":"1", "simulated_observers":str(sim)})["outputs"]
        self.assertEqual(Decimal(count(0)["illustrative_simulation_fraction"]), Decimal(0))
        self.assertEqual(Decimal(count(1)["illustrative_simulation_fraction"]), Decimal("0.5"))
        self.assertEqual(Decimal(count(999)["illustrative_simulation_fraction"]), Decimal("0.999"))
        self.assertIsNone(count(999)["empirical_likelihood_or_prior"])

    def test_video_provenance_and_fabricated_confirmation_rejected(self):
        b = read_cases()
        b["video"]["transcript_status"] = "VERIFIED_FROM_ORIGINAL"
        with self.assertRaises(AuditError):
            evaluate_batch(b)
        b = read_cases()
        b["claims"][0]["assessment"] = "SIMULATION_CONFIRMED"
        with self.assertRaises(AuditError):
            evaluate_batch(b)
        b = read_cases()
        b["claims"][0]["evidence_status"] = "DIRECT_TRANSCRIPT"
        with self.assertRaises(AuditError):
            evaluate_batch(b)

    def test_missing_source_case_duplicate_and_model_rejected(self):
        for mutation in ("source", "case", "model", "reference"):
            b = read_cases()
            if mutation == "source":
                del b["sources"]["nist_si"]
            elif mutation == "case":
                b["cases"][1]["case_id"] = b["cases"][0]["case_id"]
            elif mutation == "model":
                b["cases"][0]["model"] = "proof_of_simulation"
            else:
                b["claims"][0]["reference_ids"] = ["fabricated_science"]
            with self.subTest(mutation=mutation), self.assertRaises(AuditError):
                evaluate_batch(b)

    def test_tampered_r3a_frozen_receipt_rejected(self):
        self.assertEqual(verify_frozen_r3a()["file_sha256"], hashlib.sha256(FROZEN_R3A_PATH.read_bytes()).hexdigest())
        with tempfile.TemporaryDirectory() as folder:
            altered = Path(folder) / "altered.json"
            altered.write_bytes(FROZEN_R3A_PATH.read_bytes() + b"\n")
            with self.assertRaisesRegex(AuditError, "SHA256_MISMATCH"):
                verify_frozen_r3a(altered)

    def test_cli_replays_frozen_json_byte_for_byte(self):
        out = subprocess.check_output([sys.executable, str(ROOT / "workspace/janus_cosmos_simulation_audit.py"),
                                       "--batch", str(CASES)], cwd=ROOT)
        self.assertEqual(out, RESULT.read_bytes())


if __name__ == "__main__":
    unittest.main()
