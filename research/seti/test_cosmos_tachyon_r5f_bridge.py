"""COSMOS × TACHYON R5F bridge: deterministic metadata-only negative controls.

The real end-to-end 1,024-spectrum native/frozen scorer is executed in CI by
cosmos_tachyon_r5f_bridge.py; these unit tests do not call a remote model.
"""
import copy
import unittest
import numpy as np

from cosmos_tachyon_r5f_bridge import (
    ABSTAIN, ARMS, BUDGET, CATEGORIES, original_manifest,
    tamper_metadata, overlay_assignments, mirror_case, mirrored_exact,
    mirrored_replay, spectrum_sha, path_sha, score_sha, _metadata_check,
    cosmos_only, mirror_graph_only, mirror_replay, strong_direct,
)


class CosmosTachyonBridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw=np.zeros((16,128),dtype=np.float32)
        cls.path=np.full(16,48,dtype=np.int32)
        cls.score=0.315
        cls.lock="ab"*32
        cls.valid=original_manifest(0,cls.raw,cls.path,cls.score,cls.lock)

    def call(self,method,record):
        return method(record,self.raw,self.path,self.score,self.lock)

    def test_five_arms(self):
        self.assertEqual(len(ARMS),5)

    def test_stock_operation_budget(self):
        self.assertEqual(BUDGET,500)

    def test_native_evidence_hashes(self):
        self.assertEqual(self.valid["raw_sha256"],spectrum_sha(self.raw))
        self.assertEqual(self.valid["candidate_path_sha256"],path_sha(self.path))
        self.assertEqual(self.valid["frozen_score_sha256"],score_sha(self.score))

    def test_case_no_generator_label_leak(self):
        rec=self.valid
        for name in ("label","truth","attack_category","actual_signal","source_type"):
            self.assertNotIn(name,rec)
            self.assertNotIn(name,rec["graph"])

    def test_stages_are_computational_not_physical_clocks(self):
        self.assertEqual(self.valid["graph"]["witness_reports"],[])
        self.assertNotIn("capture_utc",self.valid)

    def test_full_reverse_graph_mirror(self):
        m=mirror_case(self.valid["graph"])
        self.assertEqual(m["graph_direction"],"C_TO_B_TO_A")
        self.assertEqual(m["anchor_source"],"score:0000")
        self.assertEqual(mirrored_exact(m)["prediction"],"raw:0000|path:0000|score:0000")

    def test_intact_preserved_by_all_nonnull_arms(self):
        for name,fn in ARMS.items():
            ok,ops=self.call(fn,self.valid)
            self.assertEqual(ok,name!="E_ALWAYS_ABSTAIN",name)
            self.assertGreaterEqual(ops,0)

    def test_seven_corruption_categories_detected_by_both_strong_arms(self):
        for category in CATEGORIES:
            with self.subTest(category=category):
                changed=tamper_metadata(self.valid,category)
                expected=category=="INTACT"
                self.assertEqual(self.call(mirror_replay,changed)[0],expected)
                self.assertEqual(self.call(strong_direct,changed)[0],expected)

    def test_weak_cosmos_baseline_does_not_check_metadata(self):
        for category in CATEGORIES:
            rec=tamper_metadata(self.valid,category)
            self.assertTrue(self.call(cosmos_only,rec)[0])

    def test_graph_only_mirror_fails_to_detect_hash_change(self):
        rec=tamper_metadata(self.valid,"BAD_RAW_DIGEST")
        self.assertTrue(self.call(mirror_graph_only,rec)[0])
        self.assertFalse(self.call(mirror_replay,rec)[0])

    def test_reject_ghost_reverse_link(self):
        rec=tamper_metadata(self.valid,"GHOST_REVERSE_PREDECESSOR")
        self.assertFalse(_metadata_check(rec,self.raw,self.path,self.score,self.lock))
        self.assertFalse(self.call(mirror_replay,rec)[0])

    def test_reject_untrusted_missing_provenance(self):
        rec=tamper_metadata(self.valid,"MISSING_RAW_DIGEST")
        self.assertFalse(self.call(strong_direct,rec)[0])

    def test_incomplete_reverse_index_abstains(self):
        rec=tamper_metadata(self.valid,"INCOMPLETE_PREDECESSOR_INDEX")
        self.assertEqual(mirrored_replay(mirror_case(rec["graph"]))["prediction"],ABSTAIN)
        self.assertFalse(self.call(mirror_replay,rec)[0])

    def test_path_change_detected_no_gold_needed(self):
        rec=copy.deepcopy(self.valid)
        other=self.path.copy()
        other[2]+=1
        self.assertFalse(_metadata_check(rec,self.raw,other,self.score,self.lock))

    def test_score_change_detected_no_gold_needed(self):
        self.assertFalse(_metadata_check(self.valid,self.raw,self.path,self.score+0.001,self.lock))

    def test_raw_array_change_detected_no_gold_needed(self):
        other=self.raw.copy()
        other[4,40]=0.125
        self.assertFalse(_metadata_check(self.valid,other,self.path,self.score,self.lock))

    def test_tamper_has_no_mutation_of_original(self):
        old=copy.deepcopy(self.valid)
        for category in CATEGORIES:
            tamper_metadata(self.valid,category)
        self.assertEqual(old,self.valid)

    def test_overlay_balanced(self):
        cats=overlay_assignments()
        self.assertEqual(len(cats),1024)
        self.assertEqual({x:cats.count(x) for x in CATEGORIES},dict.fromkeys(CATEGORIES,128))

    def test_overlay_reproducible(self):
        self.assertEqual(overlay_assignments(),overlay_assignments())

    def test_overlay_seed_no_test_truth_dependency(self):
        # No label or truth input is part of overlay_assignments signature.
        self.assertEqual(overlay_assignments(16).count("INTACT"),2)

    def test_no_model_updates(self):
        from cosmos_tachyon_r5f_bridge import PROTOCOL
        self.assertFalse(PROTOCOL["safety"]["model_training"])
        self.assertFalse(PROTOCOL["safety"]["automatic_RFI_veto"])
        self.assertFalse(PROTOCOL["safety"]["production_detector_modified"])

    def test_no_real_causality_claim(self):
        from cosmos_tachyon_r5f_bridge import PROTOCOL
        self.assertFalse(PROTOCOL["safety"]["real_telescope_ONOFF_used"])
        self.assertFalse(PROTOCOL["safety"]["external_physical_causality_or_retrocausality_proved"])


if __name__=="__main__":
    unittest.main()
