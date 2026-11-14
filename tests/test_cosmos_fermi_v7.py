#!/usr/bin/env python3
"""Independent control tests on archived V4 JSON, not independent astrophysical replication."""
import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"tools"))
from build_cosmos_fermi_v7 import (
    V4_MANIFEST,V4_EVENTS,V6_GATES,EvidenceError,
    audited_snapshot,canonical_event_ids,load,main
)

class FermiV7AuditTests(unittest.TestCase):
    def setUp(self):
        self.manifest=load(V4_MANIFEST)
        self.events=load(V4_EVENTS)
        self.v6=load(V6_GATES)

    def build(self,manifest=None,events=None,v6=None,raw_dir=None):
        return audited_snapshot(
            self.manifest if manifest is None else manifest,
            self.events if events is None else events,
            self.v6 if v6 is None else v6,
            raw_dir
        )

    def test_real_archived_events_and_claim_ceiling(self):
        d=self.build()
        self.assertEqual(d["counts"]["archive_query_roi_events"],6520)
        self.assertEqual(d["counts"]["diagnostic_selected"],175)
        self.assertEqual(d["counts"]["grbs"],3)
        self.assertEqual([x["grb"] for x in d["sources"]],["090510","090902B","130427A"])
        self.assertEqual(len(d["events"]),175)
        self.assertEqual(len({(e["grb"],e["run_id"],e["event_id"]) for e in d["events"]}),175)
        self.assertEqual(d["gates"]["source_byte_recheck_this_build"],"RAW_FITS_NOT_RECHECKED_THIS_BUILD")
        self.assertEqual(d["gates"]["multi_burst_LIV_inference"],"BLOCKED")
        self.assertIsNone(d["gates"]["observational_LIV_p_value"])
        self.assertIsNone(d["gates"]["lattice_spacing_fit"])

    def test_provenance_sha_change_is_rejected(self):
        m=copy.deepcopy(self.manifest)
        m["bursts"][0]["event"]["sha256"]="0"*64
        with self.assertRaises(EvidenceError):self.build(manifest=m)

    def test_missing_event_and_inflated_count_are_rejected(self):
        r=copy.deepcopy(self.events)
        r["sources"][0]["events"].pop()
        with self.assertRaises(EvidenceError):self.build(events=r)

    def test_duplicate_run_event_identity_is_rejected(self):
        r=copy.deepcopy(self.events)
        r["sources"][1]["events"][1]["run_id"]=r["sources"][1]["events"][0]["run_id"]
        r["sources"][1]["events"][1]["event_id"]=r["sources"][1]["events"][0]["event_id"]
        with self.assertRaises(EvidenceError):self.build(events=r)

    def test_wrong_trigger_relative_time_is_rejected(self):
        r=copy.deepcopy(self.events)
        r["sources"][2]["events"][0]["after_trigger_s"]+=1
        with self.assertRaises(EvidenceError):self.build(events=r)

    def test_nonfinite_energy_is_rejected(self):
        r=copy.deepcopy(self.events)
        r["sources"][0]["events"][0]["energy_MeV"]=float("nan")
        with self.assertRaises(EvidenceError):self.build(events=r)

    def test_removed_source_class_bit_is_rejected(self):
        r=copy.deepcopy(self.events)
        r["sources"][0]["events"][0]["event_class_mask"]=0
        with self.assertRaises(EvidenceError):self.build(events=r)

    def test_order_invariance_of_identity_listing_but_not_time_order_of_archive(self):
        src=self.events["sources"][1]
        a=canonical_event_ids(src)
        b=copy.deepcopy(src)
        b["events"]=list(reversed(b["events"]))
        self.assertEqual(a,canonical_event_ids(b))
        mutated=copy.deepcopy(self.events);mutated["sources"][1]=b
        with self.assertRaises(EvidenceError):self.build(events=mutated)

    def test_cannot_upgrade_historical_v6_physics_gate(self):
        d=copy.deepcopy(self.v6)
        d["scientific_gates"]["OBSERVED_MULTIBURST_LIV_LRT"]="PASS"
        with self.assertRaises(EvidenceError):self.build(v6=d)

    def test_raw_mode_fails_closed_without_exact_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(EvidenceError):self.build(raw_dir=Path(directory))

if __name__=="__main__":
    unittest.main(verbosity=2)
