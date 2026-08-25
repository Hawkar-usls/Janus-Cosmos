#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


core = load("cousteau_synesthetic_memory_core", "cousteau_synesthetic_memory_core.py")
overlay = load("cousteau_synesthetic_semantic_overlay", "cousteau_synesthetic_semantic_overlay.py")


class CousteauSemanticOverlayTests(unittest.TestCase):
    def passport(self, **changes):
        payload = {
            "em122_depth": 3000.0,
            "ea600_depth": 2998.0,
            "em122_minus_ea600": 2.0,
            "rolling_depth_mad": 4.0,
            "depth_local_range": 12.0,
            "depth_local_slope": 0.4,
            "outlier_score": 0.2,
            "cadence_jitter": 0.03,
            "missing_fraction": 0.01,
        }
        payload.update(changes)
        return core.build_passport(payload, direction="HEAD_FORWARD", scale="60s")

    def test_deeper_depth_maps_to_lower_register(self):
        shallow = overlay.enrich_passport(self.passport(em122_depth=100.0))
        deep = overlay.enrich_passport(self.passport(em122_depth=8000.0))
        self.assertGreater(
            shallow["cousteau_semantic_overlay"]["depth_register"]["frequency_hz"],
            deep["cousteau_semantic_overlay"]["depth_register"]["frequency_hz"],
        )

    def test_larger_cross_sensor_disagreement_maps_to_stronger_beating(self):
        small = overlay.enrich_passport(self.passport(em122_minus_ea600=0.1))
        large = overlay.enrich_passport(self.passport(em122_minus_ea600=500.0))
        self.assertLess(
            small["cousteau_semantic_overlay"]["cross_sensor_beating"]["beat_hz"],
            large["cousteau_semantic_overlay"]["cross_sensor_beating"]["beat_hz"],
        )

    def test_direction_changes_pan_not_measurement(self):
        payload = {"em122_depth": 3000.0, "missing_fraction": 0.0}
        a = core.build_passport(payload, direction="HEAD_FORWARD", scale="60s")
        b = core.build_passport(payload, direction="TAIL_REVERSE", scale="60s")
        ea = overlay.enrich_passport(a)
        eb = overlay.enrich_passport(b)
        self.assertEqual(
            ea["measurement_fingerprint"]["sha256"],
            eb["measurement_fingerprint"]["sha256"],
        )
        self.assertLess(ea["cousteau_semantic_overlay"]["track_pan"], 0)
        self.assertGreater(eb["cousteau_semantic_overlay"]["track_pan"], 0)

    def test_missingness_becomes_fog(self):
        clear = overlay.enrich_passport(self.passport(missing_fraction=0.0))
        foggy = overlay.enrich_passport(self.passport(missing_fraction=0.8))
        self.assertLess(
            clear["cousteau_semantic_overlay"]["missingness"]["fog"],
            foggy["cousteau_semantic_overlay"]["missingness"]["fog"],
        )

    def test_blocked_is_silence_and_fog(self):
        blocked = core.build_blocked_passport(source_receipt={"status": "BLOCKED"})
        enriched = overlay.enrich_passport(blocked)
        sem = enriched["cousteau_semantic_overlay"]
        self.assertEqual(sem["status"], "BLOCKED_OR_NO_MEASUREMENT")
        self.assertEqual(sem["mnemonic_sentence"], "SILENCE | FOG | NO MEASUREMENT")
        self.assertFalse(blocked["measurement_claims_allowed"])

    def test_overlay_never_mutates_measurement_fingerprint(self):
        p = self.passport()
        original = p["measurement_fingerprint"]["sha256"]
        enriched = overlay.enrich_passport(p)
        self.assertEqual(original, enriched["measurement_fingerprint"]["sha256"])
        self.assertFalse(enriched["cousteau_semantic_overlay"]["scientific_claim"])


if __name__ == "__main__":
    unittest.main()
