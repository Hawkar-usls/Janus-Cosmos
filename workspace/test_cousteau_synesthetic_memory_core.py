#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "cousteau_synesthetic_memory_core",
    HERE / "cousteau_synesthetic_memory_core.py",
)
mod = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = mod
SPEC.loader.exec_module(mod)


class SynestheticMemoryCoreTests(unittest.TestCase):
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
            "verdict": "H1_REAL_MORPHOLOGY",
            "hypothesis": "PREFERRED_STORY_MUST_NOT_LEAK",
            "target_label": "0012",
        }

    def test_self_test(self):
        self.assertEqual(mod.self_test()["status"], "PASS")

    def test_forbidden_labels_do_not_influence_measurement_fingerprint(self):
        a = mod.build_passport(self.base, direction="HEAD_FORWARD", scale="60s")
        mutated = dict(self.base)
        mutated["verdict"] = "H0_INSTRUMENT_ONLY"
        mutated["hypothesis"] = "DIFFERENT_STORY"
        mutated["target_label"] = "0037"
        b = mod.build_passport(mutated, direction="TAIL_REVERSE", scale="60s")
        self.assertEqual(
            a["measurement_fingerprint"]["sha256"],
            b["measurement_fingerprint"]["sha256"],
        )
        self.assertEqual(
            a["measurement_fingerprint"]["embedding"],
            b["measurement_fingerprint"]["embedding"],
        )

    def test_measurement_change_does_influence_fingerprint(self):
        a = mod.build_passport(self.base)
        changed = dict(self.base)
        changed["em122_depth"] += 900
        b = mod.build_passport(changed)
        self.assertNotEqual(
            a["measurement_fingerprint"]["sha256"],
            b["measurement_fingerprint"]["sha256"],
        )

    def test_direction_is_overlay_not_measurement(self):
        a = mod.build_passport(self.base, direction="HEAD_FORWARD", scale="60s")
        b = mod.build_passport(self.base, direction="TAIL_REVERSE", scale="60s")
        self.assertEqual(
            a["measurement_fingerprint"]["sha256"],
            b["measurement_fingerprint"]["sha256"],
        )
        self.assertNotEqual(a["passport_id"], b["passport_id"])

    def test_blocked_never_synthesizes_measurement(self):
        p = mod.build_blocked_passport(
            source_receipt={"status": "BLOCKED_RAW_BYTES_NOT_MOUNTED"}
        )
        self.assertEqual(p["status"], "BLOCKED_NULL")
        self.assertIsNone(p["measurement_fingerprint"])
        self.assertFalse(p["measurement_claims_allowed"])
        self.assertEqual(p["sensory_channels"]["audio"]["mode"], "SILENCE")

    def test_similarity_identical(self):
        a = mod.build_passport(self.base)
        cmp = mod.compare_passports(a, a)
        self.assertEqual(cmp["common_measurement_similarity"], 1.0)
        self.assertFalse(cmp["scientific_convergence_claim"])

    def test_heading_wrap_is_close(self):
        a = mod.build_passport({**self.base, "heading": 359.0})
        b = mod.build_passport({**self.base, "heading": 1.0})
        cmp = mod.compare_passports(a, b)
        self.assertGreater(cmp["common_measurement_similarity"], 0.99)

    def test_raw_byte_hash_is_named_as_raw(self):
        p = mod.build_passport(self.base, raw_bytes=b"actual bytes")
        self.assertEqual(p["source_identity"]["source_hash_kind"], "RAW_BYTES_SHA256")
        self.assertEqual(
            p["source_identity"]["source_sha256"],
            __import__("hashlib").sha256(b"actual bytes").hexdigest(),
        )

    def test_cross_front_rank_is_only_review_priority(self):
        h = [mod.build_passport(self.base, direction="HEAD_FORWARD", scale="60s")]
        t = [mod.build_passport(self.base, direction="TAIL_REVERSE", scale="60s")]
        ranked = mod.rank_cross_front_pairs(h, t)
        self.assertEqual(len(ranked), 1)
        self.assertFalse(ranked[0]["scientific_convergence_claim"])


if __name__ == "__main__":
    unittest.main()
