#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import unittest

import cousteau_synesthetic_research_core as research


class CousteauSynestheticResearchCoreTests(unittest.TestCase):
    def setUp(self):
        self.base = {
            "em122_depth": 3421.4,
            "ea600_depth": 3419.1,
            "em122_minus_ea600": 2.3,
            "latitude": -7.845673,
            "longitude": -14.48023,
            "heading": 359.0,
            "speed": 4.2,
            "cadence_jitter": 0.04,
            "missing_fraction": 0.02,
            "rolling_depth_mad": 2.1,
            "depth_local_range": 15.0,
            "depth_local_slope": 0.4,
            "outlier_score": 0.2,
            "verdict": "H1_REAL_MORPHOLOGY",
            "hypothesis": "PREFERRED_STORY_MUST_NOT_LEAK",
            "target_label": "0012",
        }
        self.raw_binding = {
            "source_raw_sha256": hashlib.sha256(b"hannah-fixture").hexdigest(),
            "parser_sha256": hashlib.sha256(b"parser-v1").hexdigest(),
            "window_id": "JR15001:0012:60s:fixture",
        }

    def build_synth(self, **kwargs):
        options = {
            "payload": self.base,
            "event_id": "SYNTH-0012",
            "direction": "HEAD_FORWARD",
            "scale": "60s",
            "profile": "SYNTHETIC_TEST",
        }
        options.update(kwargs)
        return research.build_research_passport(**options)

    def build_hannah(self, **kwargs):
        options = {
            "payload": self.base,
            "event_id": "JR15001:0012:60s",
            "direction": "HEAD_FORWARD",
            "scale": "60s",
            "profile": "HANNAH_BODC",
            "raw_binding": self.raw_binding,
        }
        options.update(kwargs)
        return research.build_research_passport(**options)

    def test_self_test(self):
        self.assertEqual(research.self_test()["status"], "PASS")

    def test_contract_hash_frozen(self):
        contract = research.load_and_verify_contract()
        self.assertEqual(research.digest(contract), research.PROTOCOL_CONTRACT_SHA256)

    def test_hannah_without_raw_provenance_is_fail_closed(self):
        p = research.build_research_passport(
            self.base,
            event_id="JR15001:0012:60s",
            direction="HEAD_FORWARD",
            scale="60s",
            profile="HANNAH_BODC",
        )
        self.assertEqual(p["status"], "BLOCKED")
        self.assertFalse(p["scientific_measurement_use_allowed"])
        packet = research.export_handshake_packet(p)
        self.assertIsNone(packet["measurement_fingerprint"])
        self.assertEqual(packet["epistemic_state"]["overall_state"], "BLOCKED")

    def test_hannah_raw_bound_window_is_ready(self):
        p = self.build_hannah()
        self.assertEqual(p["status"], "READY")
        self.assertTrue(p["scientific_measurement_use_allowed"])
        self.assertEqual(p["source_binding"]["status"], "RAW_SHA_BOUND_DERIVED_WINDOW")
        self.assertTrue(research.verify_handshake_packet(research.export_handshake_packet(p)))

    def test_raw_bytes_hash_bound_this_call(self):
        p = research.build_research_passport(
            self.base,
            event_id="JR15001:0012:raw",
            direction="HEAD_FORWARD",
            scale="60s",
            profile="HANNAH_BODC",
            raw_bytes=b"real bytes fixture",
        )
        self.assertEqual(p["source_binding"]["status"], "RAW_BYTES_VERIFIED_THIS_CALL")
        self.assertEqual(
            p["source_binding"]["source_raw_sha256"],
            hashlib.sha256(b"real bytes fixture").hexdigest(),
        )

    def test_story_labels_cannot_change_exported_measurement_fingerprint(self):
        a = self.build_hannah()
        mutated = dict(self.base)
        mutated.update(
            verdict="H0_INSTRUMENT_ONLY",
            hypothesis="TOTALLY_DIFFERENT_STORY",
            target_label="0037",
            pyramid="YES_PLEASE_BUT_MUST_NOT_LEAK",
        )
        b = self.build_hannah(payload=mutated, event_id="JR15001:0037:60s", direction="TAIL_REVERSE")
        af = research.export_handshake_packet(a)["measurement_fingerprint"]
        bf = research.export_handshake_packet(b)["measurement_fingerprint"]
        self.assertEqual(af["sha256"], bf["sha256"])
        self.assertEqual(af["blake2b_256"], bf["blake2b_256"])
        self.assertEqual(af["embedding"], bf["embedding"])

    def test_epistemic_overlay_does_not_change_measurement_fingerprint(self):
        a = self.build_hannah()
        b = self.build_hannah(field_states={"em122_depth": "STALE", "ea600_depth": "UNKNOWN"})
        af = research.export_handshake_packet(a)["measurement_fingerprint"]
        bf = research.export_handshake_packet(b)["measurement_fingerprint"]
        self.assertEqual(af, bf)
        self.assertNotEqual(a["epistemic"]["overall_state"], b["epistemic"]["overall_state"])

    def test_forbidden_epistemic_key_rejected(self):
        with self.assertRaises(ValueError):
            self.build_hannah(field_states={"target_label": "STALE"})

    def test_handshake_tamper_rejected(self):
        packet = research.export_handshake_packet(self.build_hannah())
        tampered = json.loads(json.dumps(packet))
        tampered["measurement_fingerprint"]["embedding"][0] += 0.01
        self.assertFalse(research.verify_handshake_packet(tampered))

    def test_quality_adjustment_never_becomes_scientific_claim(self):
        a = self.build_hannah()
        b = self.build_hannah(event_id="JR15001:0012:TAIL", direction="TAIL_REVERSE")
        cmp = research.compare_research_passports(a, b)
        self.assertGreater(cmp["quality_adjusted_review_score"], 0.0)
        self.assertFalse(cmp["scientific_convergence_claim"])
        self.assertEqual(cmp["authority"], "REVIEW_PRIORITY_ONLY")

    def test_multiscale_bundle_is_chain_bound(self):
        p60 = self.build_hannah(event_id="JR15001:0012:60s", scale="60s")
        p300 = self.build_hannah(event_id="JR15001:0012:300s", scale="300s")
        bundle = research.build_multiscale_bundle([p60, p300], bundle_id="JR15001:0012:HEAD")
        self.assertEqual(len(bundle["entries"]), 2)
        self.assertTrue(bundle["all_scientific_convergence_claims_false"])
        self.assertEqual(len(bundle["packet_chain_sha256"]), 64)
        self.assertEqual(len(bundle["bundle_sha256"]), 64)

    def test_cross_front_ranking_is_review_only(self):
        head = [self.build_hannah(event_id="HEAD-60", direction="HEAD_FORWARD")]
        tail = [self.build_hannah(event_id="TAIL-60", direction="TAIL_REVERSE")]
        rows = research.rank_cross_front_research_pairs(head, tail)
        self.assertEqual(len(rows), 1)
        self.assertFalse(rows[0]["scientific_convergence_claim"])


if __name__ == "__main__":
    unittest.main()
