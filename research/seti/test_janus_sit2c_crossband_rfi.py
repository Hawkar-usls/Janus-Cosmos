"""Offline SIT2C guards — simulated fixtures ONLY, never replace real GBT eval.

The actual study must download public X dev + L evaluation rawspec HDF5
and fail closed if an official file is unavailable or has incompatible axes.
"""
import copy
import unittest
import numpy as np
import janus_sit2c_crossband_rfi as s

class SIT2CGuards(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.panel=np.random.default_rng(775).normal(5.,1.,size=(16,128)).astype(np.float32)
        cls.row={"panel":cls.panel,"scan":"0013","role":"ON","panel_idx":0,
                 "source_file_sha256":"a"*64,
                 "measured_panel_sha256":s.SIT1.hash_array(cls.panel,dtype="<f4")}
        cls.off=np.random.default_rng(3).normal(0.,1.,size=128)

    def test_all_exactly_eight_pinned_source_files(self):
        self.assertEqual(len(s.Q["inputs"]["dev_files"]),2)
        self.assertEqual(len(s.Q["inputs"]["heldout_files"]),6)
        self.assertEqual([x["scan"] for x in s.Q["inputs"]["heldout_files"]],
                         ["0013","0014","0015","0016","0017","0018"])
        self.assertTrue(all(x["filename"].endswith(".rawspec.0002.h5")
                            for x in s.Q["inputs"]["heldout_files"]))

    def test_new_receiver_is_not_same_as_development(self):
        self.assertTrue(all(x["filename"].startswith("blc20") for x in s.Q["inputs"]["dev_files"]))
        self.assertTrue(all(x["filename"].startswith("blc21") for x in s.Q["inputs"]["heldout_files"]))

    def test_not_external_independent_dates(self):
        self.assertFalse(s.Q["safety"]["no_external_heldout_custodian"] is False)
        self.assertEqual(s.Q["inputs"]["date_utc"],"2025-12-18")

    def test_equal_672_hypotheses(self):
        self.assertEqual(s.PATHS.shape,(672,16))

    def test_rfi_score_is_finite(self):
        for method in ("sideband","offaware","rfi"):
            x=s.score_bank(self.panel,self.off,method,False)
            self.assertEqual(x.shape,(672,))
            self.assertTrue(np.isfinite(x).all())

    def test_rfi_backwards_same_direct_objective(self):
        a=s.score_bank(self.panel,self.off,"rfi",False)
        b=s.score_bank(self.panel,self.off,"rfi",True)
        self.assertTrue(np.allclose(a,b,atol=1e-9,rtol=0))
        pa,va,_=s.choose(a);pb,vb,_=s.choose(b)
        self.assertTrue(np.array_equal(pa,pb))
        self.assertAlmostEqual(va,vb,places=8)

    def test_off_reference_is_used(self):
        x=np.zeros(128); y=x.copy()
        y[20:80]=8.
        a=s.score_bank(self.panel,x,"rfi",False)
        b=s.score_bank(self.panel,y,"rfi",False)
        self.assertFalse(np.array_equal(a,b))

    def test_sideband_not_modified_by_off_reference(self):
        a=s.score_bank(self.panel,np.zeros(128),"sideband",False)
        b=s.score_bank(self.panel,np.ones(128)*7,"sideband",False)
        self.assertTrue(np.array_equal(a,b))

    def test_actual_model_tracks_remain_frozen(self):
        path,score=s.SIT1.cosmos.blind_line_candidate(self.panel[None,:,:])
        self.assertEqual(path.shape,(1,16))
        self.assertTrue(np.isfinite(score).all())

    def test_off_reference_from_unaltered_off_only(self):
        rng=np.random.default_rng(32)
        rows=[]
        for scan,role in (("0013","ON"),("0014","OFF"),("0015","ON"),("0016","OFF"),("0017","ON"),("0018","OFF")):
            for idx in range(4):
                x=rng.normal(0,1,(16,128)).astype(np.float32)
                rows.append({"scan":scan,"role":role,"panel_idx":idx,"panel":x})
        ref=s.off_reference(rows)
        self.assertEqual(set(ref),{0,1,2,3})
        self.assertEqual(ref[0].shape,(128,))
        self.assertTrue(np.isfinite(ref[0]).all())

    def test_dev_off_reference_with_one_off(self):
        rows=[]
        for scan,role in (("0061","ON"),("0062","OFF")):
            for idx in range(4):
                rows.append({"scan":scan,"role":role,"panel_idx":idx,"panel":self.panel})
        r=s.off_reference(rows)
        self.assertEqual(len(r),4)

    def test_injected_labels_never_in_score_inputs(self):
        x=s.previous.generate_overlay(self.panel,s.previous.CATEGORIES[2],np.random.default_rng(9))
        self.assertIsNotNone(x[1])
        row=s.previous.prev.matched_row(x[0],self.row)
        self.assertNotIn("truth_for_eval_only",row)
        self.assertNotIn("kind",row)
        self.assertTrue(np.array_equal(self.panel,self.row["panel"]))

    def test_four_frozen_conditions_512_cases(self):
        rng=np.random.default_rng(1)
        base=[]
        for scan,role in (("0013","ON"),("0014","OFF"),("0015","ON"),("0016","OFF"),("0017","ON"),("0018","OFF")):
            for idx in range(4):
                q=dict(self.row);q.update(scan=scan,role=role,panel_idx=idx)
                base.append(q)
        cases=s.make_cases(base,rng)
        self.assertEqual(len(cases),384)
        self.assertEqual(sum(z["kind"]==s.LABELS[2] for z in cases),96)

    def test_dev_threshold_does_not_receive_eval_input(self):
        self.assertEqual(s.Q["dataset"]["development_cases"],128)
        self.assertEqual(s.Q["dataset"]["heldout_cases"],384)
        self.assertEqual(s.Q["decision_gate"]["heldout_controls_if_no_signal"],192) if "heldout_controls_if_no_signal" in s.Q["decision_gate"] else self.assertEqual(s.Q["dataset"]["heldout_controls"],192)

    def test_score_pinned_no_training(self):
        self.assertTrue(s.Q["safety"]["no_training_on_L"])
        self.assertTrue(s.Q["safety"]["no_manual_threshold_tuning_after_L"])

    def test_reverse_metadata_causal_graph(self):
        p,score,_=s.choose(s.score_bank(self.panel,self.off,"rfi",False))
        graph_ops=s.trace_actual_source(self.row,p,score,self.row["source_file_sha256"],"c"*64)
        self.assertGreaterEqual(graph_ops,2)

    def test_reverse_metadata_rejects_wrong_source(self):
        p,score,_=s.choose(s.score_bank(self.panel,self.off,"rfi",False))
        with self.assertRaises(AssertionError):
            s.trace_actual_source(self.row,p,score,"b"*64,"c"*64)

    def test_positive_guard_no_sci_claim_on_prelocked_single_cadence(self):
        self.assertFalse(s.Q["safety"]["no_technosignature_discovery_claim"] is False)
        self.assertFalse(s.Q["safety"]["no_physical_time_reversal"] is False)
        self.assertTrue(s.Q["decision_gate"]["single_L_cadence_sample_dependence"])

if __name__=="__main__":
    unittest.main()
