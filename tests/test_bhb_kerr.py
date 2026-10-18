"""BHB-R2 Kerr science gate: analytic potential minima and local invariants."""

from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workspace"))

from janus_bhb import BHBError  # noqa: E402
from janus_bhb_kerr import evaluate_case, evaluate_batch, minimum_massive_q, minimum_photon_n  # noqa: E402


def sample(r="5", a="0", l1="0", l2="0", nphi="0", direction2="in", origin2="infinity"):
    return {"spin_a_over_M": a, "r_over_M": r,
            "particle_1": {"direction": "in", "origin": "infinity", "specific_Lz_over_M": l1},
            "particle_2": {"direction": direction2, "origin": origin2, "specific_Lz_over_M": l2},
            "photon": {"source": "hypothetical_zamo_emission", "azimuth_sine": nphi}}


class BHBCosmosKerrTests(unittest.TestCase):
    def test_schwarzschild_consistency_and_redshift(self):
        for r in ("2.01", "3", "100"):
            result = evaluate_case(sample(r=r))
            self.assertEqual(Decimal(result["collision"]["E_cm_GeV"]), Decimal("1.876544"))
            self.assertEqual(result["collision"]["tetrad_cross_check"], "PASS")
            self.assertEqual(result["test_photon"]["escape_gate"], "EQUATORIAL_NULL_GEODESIC_REACHES_INFINITY")
            with localcontext() as ctx:
                ctx.prec = 50
                expected = (1 - Decimal(2) / Decimal(r)).sqrt()
                self.assertLess(abs(Decimal(result["test_photon"]["E_inf_over_E_zamo"]) - expected), Decimal("1e-23"))

    def test_near_extremal_critical_limit_grows_without_horizon_collision(self):
        low = evaluate_case(sample("1.01", "1", l2="2"))
        high = evaluate_case(sample("1.0001", "1", l2="2"))
        self.assertEqual(high["trajectory_pair_gate"], "PASS")
        self.assertGreater(Decimal(high["collision"]["E_cm_GeV"]), 5 * Decimal(low["collision"]["E_cm_GeV"]))
        self.assertEqual(high["escaped_collision_signal_gate"], "OPEN_NO_COLLISION_PRODUCT_FOUR_MOMENTUM")

    def test_nonextremal_critical_particle_not_admitted(self):
        row = evaluate_case(sample("1.5", "0.9", l2="2.7"))
        self.assertNotEqual(row["trajectory_pair_gate"], "PASS")
        self.assertEqual(row["particle_2"]["trajectory_gate"], "REJECTED_FORBIDDEN_RADIAL_BARRIER_FROM_INFINITY")
        self.assertEqual(row["collision"]["tetrad_cross_check"], "PASS")

    def test_escape_cone_barrier_and_outside_photon_sphere(self):
        blocked = evaluate_case(sample("2.5", nphi="0.99"))
        escaped = evaluate_case(sample("4", nphi="0.99"))
        self.assertEqual(blocked["test_photon"]["escape_gate"], "BLOCKED_OR_CRITICAL_RADIAL_POTENTIAL")
        self.assertEqual(escaped["test_photon"]["escape_gate"], "EQUATORIAL_NULL_GEODESIC_REACHES_INFINITY")

    def test_no_outgoing_infinity_provenance_and_no_collision_product_claim(self):
        row = evaluate_case(sample("1.6", "0.9", direction2="out", origin2="local_production"))
        self.assertEqual(row["trajectory_pair_gate"], "NOT_ESTABLISHED_OR_REJECTED")
        self.assertEqual(row["test_photon"]["collision_product_gate"], "NOT_MODELED_NO_EMISSION_SPECTRUM_OR_YIELD")

    def test_reflection_symmetry_and_mass_scaling(self):
        p = sample("2", "0.9", l1="0", l2="2", nphi="0.3")
        forward = evaluate_case(p)
        q = sample("2", "-0.9", l1="0", l2="-2", nphi="-0.3")
        reverse = evaluate_case(q)
        self.assertEqual(forward["collision"]["E_cm_GeV"], reverse["collision"]["E_cm_GeV"])
        self.assertEqual(forward["test_photon"]["E_inf_over_E_zamo"], reverse["test_photon"]["E_inf_over_E_zamo"])
        p["particle_1"]["rest_energy_GeV"] = "1.876544"
        p["particle_2"]["rest_energy_GeV"] = "1.876544"
        twice = evaluate_case(p)
        self.assertLess(abs(Decimal(twice["collision"]["E_cm_GeV"]) - 2 * Decimal(forward["collision"]["E_cm_GeV"])), Decimal("1e-20"))

    def test_exact_potential_minima_match_analytic_schwarzschild_limits(self):
        with localcontext() as ctx:
            ctx.prec = 80
            # Marginally bound E=1, L=4 has double root at r=4 in Schwarzschild.
            self.assertEqual(minimum_massive_q(Decimal("2"),Decimal("0"),Decimal("1"),Decimal("4")), Decimal(0))
            # Impact b=sqrt(27) gives double photon radial root at r=3.
            b = Decimal(27).sqrt()
            self.assertLess(abs(minimum_photon_n(Decimal("2.5"),Decimal(0),b)), Decimal("1e-70"))

    def test_double_root_separatrix_is_not_finite_arrival(self):
        row = evaluate_case(sample("2.5", l2="4"))
        self.assertEqual(row["particle_2"]["trajectory_gate"], "REJECTED_FORBIDDEN_RADIAL_BARRIER_FROM_INFINITY")
        self.assertEqual(row["trajectory_pair_gate"], "NOT_ESTABLISHED_OR_REJECTED")

    def test_invalid_geometry_or_photon_rejected(self):
        bad_cases = [sample("2", "0"), sample("1", "1"), sample("4", "1.01"), sample("5", "0", nphi="1"), sample("nan", "0")]
        for bad in bad_cases:
            with self.subTest(bad=bad), self.assertRaises(BHBError):
                evaluate_case(bad)

    def test_frozen_batch(self):
        inp = ROOT / "data/bhb/BHB-COSMOS-R2-CASES-2026-10-18.json"
        frozen = ROOT / "data/bhb/BHB-COSMOS-R2-RESULT-2026-10-18.json"
        result = evaluate_batch(json.loads(inp.read_text()))
        self.assertEqual(result, json.loads(frozen.read_text()))
        self.assertEqual(len(result["cases"]), 7)
        self.assertEqual(result["gates"]["REAL_BHB_DETECTION"], "NOT_ESTABLISHED")
        self.assertTrue(all(row["test_photon"]["null_frame_cross_check"] == "PASS" for row in result["cases"]))


if __name__ == "__main__":
    unittest.main()
