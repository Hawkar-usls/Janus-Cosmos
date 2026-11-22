"""JANUS SIT2B unit checks; all local arrays are synthetic fixtures.

The real data experiment is required to re-download all six checksum-locked
GBT scans. Synthetic fixtures do not qualify as observational measurements.
"""
import copy
import unittest
import numpy as np
import janus_sit2b_temporal_recovery as s

class SIT2BControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.panel=np.random.default_rng(17).normal(3,0.5,(16,128)).astype(np.float32)
        cls.source={"scan":"0061","role":"ON","panel_idx":0,
                    "panel":cls.panel,"source_file_sha256":"a"*64,
                    "measured_panel_sha256":s.SIT1.hash_array(cls.panel,dtype="<f4")}

    def test_frozen_hypotheses_672(self):
        self.assertEqual(s.PATHS.shape,(672,16))
        self.assertEqual(len(s.SLOPES),7)

    def test_same_native_path_grid(self):
        npath,nval=s.SIT1.cosmos.blind_line_candidate(self.panel[None,:,:])
        dpath,dval,_=s.direct_matched(self.panel,False)
        self.assertTrue(np.array_equal(npath[0],dpath))
        self.assertAlmostEqual(float(nval[0]),dval,places=9)

    def test_reverse_and_forward_same_sideband_path(self):
        dpath,dscore,_=s.direct_matched(self.panel,True)
        rev=s.scores_for_paths(self.panel,True,True)
        rpath,rscore,_=s.best_path(rev)
        self.assertTrue(np.array_equal(dpath,rpath))
        self.assertAlmostEqual(dscore,rscore,places=9)

    def test_reverse_replay_actual_pinned_graph(self):
        row=copy.deepcopy(self.source)
        a,b,info=s.reverse_sideband_verified(row,"a"*64)
        self.assertEqual(len(a),16)
        self.assertTrue(np.isfinite(b))
        self.assertGreaterEqual(info["graph_operations"],2)

    def test_source_not_mutated_by_reverse_search(self):
        row=copy.deepcopy(self.source)
        s.reverse_sideband_verified(row,"a"*64)
        self.assertTrue(np.array_equal(row["panel"],self.panel))

    def test_strong_control_amplitude_exceeds_weak(self):
        p=s.Q["preregistered_dataset"]
        self.assertGreater(p["strong_amplitude_range_robust_row_mad"][0],
                           p["weak_amplitude_range_robust_row_mad"][1])

    def test_strong_positive_clear_test_signal(self):
        st,truth=s.generate_overlay(self.panel,s.CATEGORIES[3],np.random.default_rng(42))
        self.assertIsNotNone(truth)
        rpath,score,_=s.direct_matched(st,True)
        self.assertTrue(s.localizes(rpath,truth))

    def test_weak_injection_known_truth_only_in_evaluator(self):
        p,truth=s.generate_overlay(self.panel,s.CATEGORIES[2],np.random.default_rng(12))
        self.assertEqual(len(truth),16)
        self.assertFalse(np.array_equal(p,self.panel))
        self.assertTrue(np.array_equal(self.panel,self.source["panel"]))

    def test_no_injection_preserves_original_bytes(self):
        p,truth=s.generate_overlay(self.panel,s.CATEGORIES[0],np.random.default_rng(12))
        self.assertIsNone(truth)
        self.assertTrue(np.array_equal(p,self.panel))

    def test_engineered_RFI_is_not_astronomical_truth(self):
        p,truth=s.generate_overlay(self.panel,s.CATEGORIES[1],np.random.default_rng(12))
        self.assertIsNone(truth)
        self.assertFalse(np.array_equal(p,self.panel))

    def test_seed_reproducibility(self):
        p1,t1=s.generate_overlay(self.panel,s.CATEGORIES[2],np.random.default_rng(2026))
        p2,t2=s.generate_overlay(self.panel,s.CATEGORIES[2],np.random.default_rng(2026))
        self.assertTrue(np.array_equal(p1,p2))
        self.assertTrue(np.array_equal(t1,t2))

    def test_all_cases_fixed_384(self):
        rows=[]
        for scan in ("0061","0062","0063","0064","0065","0066"):
            for j in range(4):
                row=copy.copy(self.source)
                row["scan"]=scan;row["panel_idx"]=j
                rows.append(row)
        x=s.all_cases(rows)
        self.assertEqual(len(x),384)
        self.assertEqual(sum(z["kind"]==s.CATEGORIES[2] for z in x),96)
        self.assertEqual(sum(z["kind"]==s.CATEGORIES[3] for z in x),96)

    def test_threshold_only_negative_64_dev(self):
        self.assertEqual(s.decide_threshold(np.arange(64)),60.0)
        with self.assertRaises(ValueError):s.decide_threshold(np.arange(63))

    def test_scientific_gate_requires_strong_positive(self):
        # Mock completely frozen case-score output for five arms.
        arms=("A_NATIVE_COSMOS_RIDGE","B_FROZEN_COSMOS_SHADOW",
              "C_STRONG_DIRECT_MATCHED_FILTER","D_STRONG_DIRECT_SIDEBAND",
              "E_TACHYON_REVERSE_SIDEBAND")
        x=[]
        for i in range(384):
            kind=s.CATEGORIES[i%4]
            scan=("0061","0062","0063","0064","0065","0066")[i//64]
            x.append({"case_id":i,"scan":scan,"split":"development" if scan in s.DEVELOPMENT else "evaluation",
                      "kind":kind,"scores":{k:0. for k in arms},
                      "located":{k:False for k in arms}})
        d,g=s.evaluate(x)
        self.assertFalse(g["all_arms_positive_guards_pass"])
        self.assertFalse(g["tachyon_specific_benchmark_gain_gate"])

    def test_eval_label_not_in_native_row(self):
        row=s.prev.matched_row(self.panel,self.source)
        self.assertNotIn("ground_truth_for_evaluation_only",row)
        self.assertNotIn("kind",row)

    def test_no_physical_reverse_time_inference(self):
        self.assertFalse(s.Q["physical_science_gate"]["physical_retrocausality_claim"])
        self.assertFalse(s.Q["physical_science_gate"]["original_real_emission_ground_truth"])

    def test_no_production_model_mutation(self):
        self.assertFalse(s.Q["physical_science_gate"]["production_updates"])

    def test_one_observational_cadence_dependent(self):
        self.assertEqual(s.Q["immutable_inputs"]["backgrounds"],24)
        self.assertEqual(s.Q["immutable_inputs"]["real_frequency_resolution_hz"],2861.02294921875)

if __name__=="__main__":
    unittest.main()
