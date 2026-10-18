"""JANUS Cosmos V3 guards against interpreting real Fermi GRB summaries as LIV."""

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

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"workspace"))
from janus_cosmos_fermi_observation_gate import (  # noqa: E402
    AuditError, CATALOG_PATH, CATALOG_SHA256, FROZEN_V2_PATH,
    evaluate_catalog, validate_catalog, verify_history,
)

RESULT=ROOT/"data/cosmos/COSMOS-ALI-V3-FERMI-PUBLISHED-TOF-METHOD-RESULT-2026-10-18.json"


class FermiPublishedObservationGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog=json.loads(CATALOG_PATH.read_text(encoding="utf-8"))
        cls.receipt=json.loads(RESULT.read_text(encoding="utf-8"))

    def test_source_transcription_exact_sha_and_v1_v2_continuity(self):
        self.assertEqual(hashlib.sha256(CATALOG_PATH.read_bytes()).hexdigest(),CATALOG_SHA256)
        frozen=verify_history()
        self.assertEqual(frozen["liv_status"],"NOT_ESTABLISHED")
        self.assertEqual(frozen["amps_status"],"OPEN_UNCONFIRMED")
        with tempfile.TemporaryDirectory() as tmp:
            bad=Path(tmp)/"v2.json"
            bad.write_bytes(FROZEN_V2_PATH.read_bytes()+b"\n")
            with self.assertRaisesRegex(AuditError,"V2_REVEAL_SHA256_MISMATCH"):
                verify_history(bad)

    def test_real_published_events_count_and_match_frozen_receipt(self):
        self.assertEqual(evaluate_catalog(self.catalog),self.receipt)
        self.assertEqual(self.receipt["counts"]["bursts"],11)
        self.assertEqual(self.receipt["counts"]["published_high_energy_photons"],38)
        self.assertEqual(self.receipt["counts"]["within_burst_energy_order_pairs"],162)
        self.assertEqual(self.receipt["counts"]["higher_energy_arrived_earlier"],61)
        self.assertEqual(self.receipt["counts"]["higher_energy_arrived_later"],101)

    def test_real_source_values_distinct_from_hypothetical_model_reference(self):
        self.assertEqual(self.receipt["observation_type"],"REAL_PUBLISHED_FERMI_LAT_2FLGC_ROUNDED_TABLE_7_VALUES_NOT_RAW_EVENT_FITS")
        self.assertEqual(self.receipt["reference_photon"],"HYPOTHETICAL_NOT_AN_OBSERVED_EVENT")
        self.assertEqual(self.receipt["counts"]["matched_low_energy_photons"],0)
        self.assertEqual(self.receipt["counts"]["events_with_verified_emission_time"],0)

    def test_photon_model_overlays_are_nonnegative_and_planck_negligible(self):
        for burst in self.receipt["grbs"]:
            for p in burst["photons"]:
                d=p["lattice_predicted_delay_vs_unobserved_1GeV_reference_s"]
                with self.subTest(grb=burst["grb"],energy=p["measured_energy_GeV"]):
                    self.assertGreater(Decimal(d["4e-27"]),Decimal(0))
                    self.assertLess(Decimal(d["1.616255e-35"]),Decimal("1e-12"))
                    self.assertIsNone(p["fitted_lattice_spacing_m"])
                    self.assertEqual(p["physical_speed_inference"],"BLOCKED_UNKNOWN_SOURCE_EMISSION_TIME")

    def test_known_photon_example_has_mixed_arrival_order(self):
        e=self.receipt["negative_control_example"]
        self.assertEqual(e["grb"],"090902B")
        self.assertEqual(Decimal(e["observed_pair_time_difference_high_minus_low_s"]),Decimal("-250.2"))
        self.assertEqual(self.receipt["gates"]["REAL_LIV_LIKELIHOOD_FIT"],"BLOCKED_NOT_IDENTIFIABLE")

    def test_source_forgery_or_missing_event_rejected(self):
        for mutation in ("citation","raw_data","precision","time_zero","selection","truncated","duplicate","type"):
            c=copy.deepcopy(self.catalog)
            if mutation=="citation": c["source"]["versioned_html_table_url"]="https://example.org/some-paper"
            if mutation=="raw_data": c["source"]["event_catalog_fits"]="VERIFIED"
            if mutation=="precision": c["source"]["precision"]="NANOSECOND_EXACT"
            if mutation=="time_zero": c["source"]["time_zero"]="COMMON_EMISSION_TIME"
            if mutation=="selection": c["source"]["association_selection"]="P=1_EACH"
            if mutation=="truncated": c["bursts"][5]["events"].pop()
            if mutation=="duplicate": c["bursts"][1]["grb"]=c["bursts"][0]["grb"]
            if mutation=="type": c["bursts"][1]["class"]="LONG"
            with self.subTest(mutation=mutation),self.assertRaises(AuditError):
                validate_catalog(c)

    def test_no_simultaneous_emission_or_liv_promotion(self):
        for mutation in ("assumed_pair","calibrated_source","instrument","raw","calibrated_trials","analysis_ceiling","fake_ref"):
            c=copy.deepcopy(self.catalog)
            a=c["analysis_protocol"]
            if mutation=="assumed_pair": a["true_event_pair_lag_available"]=True
            if mutation=="calibrated_source": a["source_intrinsic_lag_calibrated"]=True
            if mutation=="instrument": a["instrument_response_reconstructed"]=True
            if mutation=="raw": a["raw_event_files_ingested"]=True
            if mutation=="calibrated_trials": a["statistical_trials_calibrated"]=True
            if mutation=="analysis_ceiling": a["inference_ceiling"]="LIV_OBSERVED"
            if mutation=="fake_ref": a["reference_energy_status"]="OBSERVED_MATCHED_PHOTON"
            with self.subTest(mutation=mutation),self.assertRaises(AuditError):
                validate_catalog(c)

    def test_invalid_numeric_energy_time_and_redshift_fail_closed(self):
        for mutation in ("negative_time","low_energy","not_sorted","bad_z","NaN_energy","probability"):
            c=copy.deepcopy(self.catalog)
            b=c["bursts"][2]
            if mutation=="negative_time": b["events"][0]["time_s"]="-1"
            if mutation=="low_energy": b["events"][-1]["energy_GeV"]="9.9"
            if mutation=="not_sorted": b["events"][0]["energy_GeV"]="11"
            if mutation=="bad_z": b["redshift"]="NaN"
            if mutation=="NaN_energy": b["events"][0]["energy_GeV"]="NaN"
            if mutation=="probability": b["events"][0]["probability"]="1.0"
            with self.subTest(mutation=mutation),self.assertRaises(AuditError):
                validate_catalog(c)

    def test_no_observational_claim_or_bhb_promotion(self):
        g=self.receipt["gates"]
        self.assertEqual(g["LATTICE_SPACING_CONSTRAINT"],"NOT_DERIVED")
        self.assertEqual(g["LIV_OBSERVED"],"NOT_ESTABLISHED")
        self.assertEqual(g["SIMULATION_HYPOTHESIS"],"NOT_ESTABLISHED")
        self.assertEqual(g["AMPS_FIREWALL"],"OPEN_UNCONFIRMED")
        self.assertEqual(g["BHB_BOOM"],"NOT_ESTABLISHED")
        self.assertEqual(g["OSIRIS_PROOF_AUTHORITY"],"NO_CHANGE")

    def test_cli_replays_exact_published_observation_gate_receipt(self):
        script=ROOT/"workspace/janus_cosmos_fermi_observation_gate.py"
        result=subprocess.check_output([sys.executable,str(script),"--batch",str(CATALOG_PATH)],cwd=ROOT)
        self.assertEqual(result,RESULT.read_bytes())

    def test_cli_rejects_even_small_transcription_sha_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad=Path(tmp)/"bad.json"
            bad.write_bytes(CATALOG_PATH.read_bytes()+b"\n")
            r=subprocess.run([sys.executable,str(ROOT/"workspace/janus_cosmos_fermi_observation_gate.py"),"--batch",str(bad)],cwd=ROOT,capture_output=True,text=True)
            self.assertEqual(r.returncode,2)
            self.assertIn("TRANSCRIPTION_SHA256_MISMATCH",r.stderr)


if __name__=="__main__":
    unittest.main()
