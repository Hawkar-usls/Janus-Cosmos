"""Reproducible, provenance-preserving synthetic lattice photon TOF gates."""

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
from janus_cosmos_lattice_tof import (  # noqa: E402
    AuditError, V1_FILE, frozen_v1_gate, lattice_delay_seconds, prepare,
    reveal, score, validate_protocol,
)

BASE = ROOT / "data/cosmos"
PREFIX = "COSMOS-ALI-V2-LATTICE-TOF-"
SUFFIX = "-2026-10-18.json"
PROTOCOL_PATH = BASE / (PREFIX + "PROTOCOL" + SUFFIX)
TRUTH_PATH = BASE / (PREFIX + "SEALED-TRUTH" + SUFFIX)
OBS_PATH = BASE / (PREFIX + "OBSERVATIONS" + SUFFIX)
SCORE_PATH = BASE / (PREFIX + "BLIND-SCORE" + SUFFIX)
REVEAL_PATH = BASE / (PREFIX + "REVEAL" + SUFFIX)


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


class CosmosLatticeTofV2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.protocol = load(PROTOCOL_PATH)
        cls.truth = load(TRUTH_PATH)
        cls.observations = load(OBS_PATH)
        cls.scoring = load(SCORE_PATH)
        cls.revealed = load(REVEAL_PATH)
        cls.protocol_hash = hashlib.sha256(PROTOCOL_PATH.read_bytes()).hexdigest()
        cls.truth_hash = hashlib.sha256(TRUTH_PATH.read_bytes()).hexdigest()
        cls.by_id = {c["case_id"]:c for c in cls.protocol["cases"]}
        cls.scores_by_id = {c["case_id"]:c for c in cls.scoring["cases"]}

    def test_frozen_v1_boundaries_are_preserved(self):
        self.assertEqual(frozen_v1_gate()["simulation_gate"], "NOT_ESTABLISHED")
        with tempfile.TemporaryDirectory() as folder:
            altered=Path(folder)/"altered.json"
            altered.write_bytes(V1_FILE.read_bytes()+b"\n")
            with self.assertRaisesRegex(AuditError,"SHA256_MISMATCH"):
                frozen_v1_gate(altered)

    def test_protocol_rejects_missing_science_and_source_gates(self):
        for mutation in ("ceiling","origin","source","duplicate","negative_a","supercritical_q","negative_sigma","nonflat"):
            p=copy.deepcopy(self.protocol)
            if mutation=="ceiling": p["scientific_ceiling"]="LIV_OBSERVED"
            if mutation=="origin": p["all_observations"]="REAL"
            if mutation=="source": del p["sources"]["magic"]
            if mutation=="duplicate": p["cases"][1]["case_id"]=p["cases"][0]["case_id"]
            if mutation=="negative_a": p["cases"][0]["lattice_spacing_m"]="-1"
            if mutation=="supercritical_q": p["cases"][0]["lattice_spacing_m"]="1e-10"
            if mutation=="negative_sigma": p["cases"][0]["timing_sigma_s"]="0"
            if mutation=="nonflat": p["cosmology"]["Omega_lambda"]="0.9"
            with self.subTest(mutation=mutation), self.assertRaises(AuditError):
                validate_protocol(p)

    def test_prediction_monotonic_with_energy_spacing_and_redshift(self):
        baseline=self.by_id["V2-LATTICE-KNOWN"]
        value=lattice_delay_seconds(baseline,self.protocol["cosmology"])
        self.assertGreater(value,Decimal("15"))
        self.assertLess(value,Decimal("16"))
        for field,updated in (("lattice_spacing_m","2e-27"),
                              ("energy_high_eV","2e12"),("redshift","1")):
            changed=copy.deepcopy(baseline)
            changed[field]=updated
            new=lattice_delay_seconds(changed,self.protocol["cosmology"])
            with self.subTest(field=field):
                if field=="lattice_spacing_m": self.assertLess(new,value)
                else: self.assertGreater(new,value)

    def test_zero_distance_or_energy_separation_has_no_effect(self):
        for name in ("V2-EQUAL-ENERGY","V2-ZERO-BASELINE"):
            self.assertEqual(lattice_delay_seconds(self.by_id[name],self.protocol["cosmology"]),Decimal(0))

    def test_planck_scale_choice_unmeasurable_at_fixture_precision(self):
        item=self.scores_by_id["V2-PLANCK-NEGLIGIBLE"]
        self.assertLess(Decimal(item["predicted_lattice_minus_lorentz_lag_s"]),Decimal("1e-12"))
        self.assertEqual(item["model_preference"],"INCONCLUSIVE_SYNTHETIC")

    def test_fixed_known_zero_source_recovers_both_injected_hypotheses(self):
        for id_,expected in (("V2-LI-KNOWN","LORENTZ_PREFERRED_SYNTHETIC"),
                             ("V2-LATTICE-KNOWN","LATTICE_PREFERRED_SYNTHETIC")):
            self.assertEqual(self.scores_by_id[id_]["model_preference"],expected)
            self.assertGreater(abs(Decimal(self.scores_by_id[id_]["delta_chi2_lorentz_minus_lattice"])),Decimal(25))

    def test_free_intrinsic_lag_is_exactly_unidentifiable(self):
        for id_ in ("V2-LATTICE-FREE","V2-LI-MIMIC-FREE"):
            item=self.scores_by_id[id_]
            self.assertEqual(item["model_preference"],"UNIDENTIFIABLE_FREE_INTRINSIC_LAG")
            self.assertIsNone(item["delta_chi2_lorentz_minus_lattice"])
            self.assertIsNone(item["chi2_lorentz"])
            self.assertIsNone(item["chi2_lattice"])

    def test_source_mimic_generates_an_explicit_false_positive(self):
        item=self.scores_by_id["V2-LI-MIMIC-FALSEPOSITIVE"]
        self.assertEqual(item["model_preference"],"LATTICE_PREFERRED_SYNTHETIC")
        self.assertEqual(self.revealed["counters"]["known_source_lag_mismatch_false_positive_controls"],1)
        self.assertEqual(self.revealed["counters"]["incorrect_preferences"],1)

    def test_bounded_source_lag_and_uncertainty_behave_as_preregistered(self):
        self.assertEqual(self.scores_by_id["V2-LI-BOUNDED"]["model_preference"],"LORENTZ_PREFERRED_SYNTHETIC")
        self.assertEqual(self.scores_by_id["V2-LATTICE-BOUNDED"]["model_preference"],"LATTICE_PREFERRED_SYNTHETIC")
        for id_ in ("V2-LARGE-UNCERTAINTY","V2-SUBTHRESHOLD"):
            self.assertEqual(self.scores_by_id[id_]["model_preference"],"INCONCLUSIVE_SYNTHETIC")

    def test_prepare_does_not_leak_injection_model_or_true_source_lag(self):
        self.assertEqual(self.observations["case_count"],12)
        for obs in self.observations["observations"]:
            self.assertEqual(obs["origin"],"SYNTHETIC_ONLY")
            self.assertNotIn("injected_model",obs)
            self.assertNotIn("intrinsic_source_lag_s",obs)
        self.assertNotIn("injections",self.observations)

    def test_sha_tamper_controls_and_exact_reveal_reproduction(self):
        regenerated=prepare(self.protocol,self.truth,self.protocol_hash,self.truth_hash)
        self.assertEqual(regenerated,self.observations)
        self.assertEqual(score(self.protocol,self.observations,self.protocol_hash),self.scoring)
        self.assertEqual(reveal(self.protocol,self.observations,self.scoring,self.truth,self.protocol_hash,self.truth_hash),self.revealed)
        bad_obs=copy.deepcopy(self.observations)
        bad_obs["observations"][0]["observed_high_minus_low_arrival_s"]="42"
        with self.assertRaisesRegex(AuditError,"integrity"):
            score(self.protocol,bad_obs,self.protocol_hash)
        bad_scoring=copy.deepcopy(self.scoring)
        bad_scoring["cases"][0]["model_preference"]="PROOF_OF_SIMULATION"
        with self.assertRaisesRegex(AuditError,"tampered"):
            reveal(self.protocol,self.observations,bad_scoring,self.truth,self.protocol_hash,self.truth_hash)
        with self.assertRaisesRegex(AuditError,"TRUTH_SHA256"):
            reveal(self.protocol,self.observations,self.scoring,self.truth,self.protocol_hash,"badsha")

    def test_all_three_cli_phases_replay_fixture_bytes(self):
        script=ROOT/"workspace/janus_cosmos_lattice_tof.py"
        base=[sys.executable,str(script),"--protocol",str(PROTOCOL_PATH)]
        for phase,extra,expected in (
            ("--prepare",["--truth",str(TRUTH_PATH)],OBS_PATH),
            ("--score",["--observations",str(OBS_PATH)],SCORE_PATH),
            ("--reveal",["--observations",str(OBS_PATH),"--blind-score",str(SCORE_PATH),"--truth",str(TRUTH_PATH)],REVEAL_PATH),
        ):
            with self.subTest(phase=phase):
                result=subprocess.check_output(base+[phase]+extra,cwd=ROOT)
                self.assertEqual(result,expected.read_bytes())

    def test_science_claim_ceiling_still_explicitly_negative(self):
        self.assertEqual(self.scoring["gates"]["QUADRATIC_LIV_OBSERVED"],"NOT_ESTABLISHED")
        self.assertEqual(self.scoring["gates"]["BHB_QUANTUM_FIREWALL"],"OPEN_UNCONFIRMED")
        self.assertEqual(self.scoring["gates"]["SIMULATION_HYPOTHESIS_OBSERVED"],"NOT_ESTABLISHED")
        self.assertEqual(self.revealed["gates"]["ACTUAL_LIV_OBSERVED"],"NOT_ESTABLISHED")


if __name__=="__main__":
    unittest.main()
