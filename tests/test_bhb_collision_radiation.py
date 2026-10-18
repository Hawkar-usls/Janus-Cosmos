"""BHB-R3A conservation-complete local annihilation and escape controls."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from decimal import Decimal, localcontext
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"workspace"))
from janus_bhb import BHBError  # noqa: E402
from janus_bhb_collision_radiation import evaluate_case,evaluate_batch,verify_fw1,FW1_PATH  # noqa: E402


def case(x="1.01",e="1"):
    return {"r_over_rs":x,"electron_specific_energy":e,
            "process":"electron_positron_to_two_radial_photons",
            "target":"externally_supported_static_positron"}


class BHBQEDRadiationTests(unittest.TestCase):
    def test_frozen_fw1_provenance(self):
        self.assertEqual(verify_fw1()["fw1_sha256"],"b16a3f086f63a2df06423c2f399dd847d40e4857b0c2bc04c5f7a6b0d9923384")
        with tempfile.TemporaryDirectory() as d:
            bad=Path(d)/"altered.json"
            bad.write_bytes(FW1_PATH.read_bytes()+b"\n")
            with self.assertRaisesRegex(BHBError,"SHA256_MISMATCH"):
                verify_fw1(bad)

    def test_frozen_seven_case_receipt(self):
        p=ROOT/"data/bhb/BHB-R3A-QED-COLLISION-CASES-2026-10-18.json"
        q=ROOT/"data/bhb/BHB-R3A-QED-COLLISION-RESULT-2026-10-18.json"
        out=evaluate_batch(json.loads(p.read_text()))
        self.assertEqual(out,json.loads(q.read_text()))
        self.assertEqual(len(out["cases"]),7)
        self.assertEqual(out["gates"]["OBSERVED_BHB_RADIATION"],"NOT_ESTABLISHED")

    def test_four_momentum_and_killing_energy_exact(self):
        with localcontext() as ctx:
            ctx.prec=80
            for x in ("100", "1.1", "1.000001", "1.000000000000000000000001"):
                row=evaluate_case(case(x))
                self.assertGreater(Decimal(row["r_over_rs"]),Decimal(1))
                E=Decimal(row["incoming_four_momentum_tetrad_GeV"]["energy"])
                P=Decimal(row["incoming_four_momentum_tetrad_GeV"]["radial_momentum"])
                out=Decimal(row["outgoing_photon"]["local_energy_GeV"])
                inward=Decimal(row["ingoing_photon"]["local_energy_GeV"])
                self.assertLess(abs(E-out-inward),Decimal("1e-19")*E)
                self.assertLess(abs(P-out+inward),Decimal("1e-19")*E)
                initial=Decimal(row["initial_Killing_energy_GeV"])
                outgoing=Decimal(row["outgoing_photon"]["E_infinity_GeV"])
                incoming=Decimal(row["ingoing_photon"]["conserved_Killing_energy_GeV"])
                self.assertLess(abs(initial-outgoing-incoming),Decimal("1e-19")*initial)
                self.assertLess(abs(Decimal(row["electron_local_speed_over_c_signed"])),Decimal(1))
                self.assertGreater(Decimal(row["one_minus_absolute_electron_speed_over_c"]),Decimal(0))

    def test_high_local_ecm_does_not_mean_high_escape_energy(self):
        r1=evaluate_case(case("1.01"))
        r2=evaluate_case(case("1.000000000000000000000001"))
        self.assertGreater(Decimal(r2["E_cm_local_GeV"]),100*Decimal(r1["E_cm_local_GeV"]))
        self.assertLess(Decimal(r2["outgoing_photon"]["E_infinity_GeV"]),Decimal("1e-9")*Decimal(r1["outgoing_photon"]["E_infinity_GeV"]))
        self.assertEqual(r2["null_transport_gate"],"PASS_R2_SCHWARZSCHILD_LIMIT")
        self.assertEqual(r2["target_source_gate"],"EXTERNAL_STATIC_SUPPORT_REQUIRED_NOT_ESTABLISHED")

    def test_holding_positron_more_costly_near_horizon(self):
        a=evaluate_case(case("1.01"))
        b=evaluate_case(case("1.000001"))
        self.assertGreater(Decimal(b["static_target_required_proper_acceleration_times_M"]),Decimal(a["static_target_required_proper_acceleration_times_M"]))
        self.assertIsNone(b["annihilation_event_rate"])
        self.assertIsNone(b["qed_differential_cross_section"])

    def test_invalid_or_false_firewall_source_rejected(self):
        for bad in (case("1"),case("0.99"),case("1.01","0.99"),case("nan")):
            with self.subTest(bad=bad),self.assertRaises(BHBError):
                evaluate_case(bad)
        bad=case()
        bad["target"]="confirmed_AMPS_firewall"
        with self.assertRaises(BHBError):
            evaluate_case(bad)
        bad=case()
        bad["spacetime"]="kerr"
        with self.assertRaises(BHBError):
            evaluate_case(bad)


if __name__=="__main__":
    unittest.main()
