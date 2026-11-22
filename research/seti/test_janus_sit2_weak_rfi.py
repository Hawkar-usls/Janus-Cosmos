"""SIT2 offline contracts; real evidence is fetched ONLY in the workflow audit."""
import copy
import unittest
import numpy as np

import janus_sit2_weak_rfi as s

class SIT2WeakRfiContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base=np.random.default_rng(223).normal(2,0.3,size=(16,128)).astype(np.float32)
        cls.row={"scan":"0061","role":"ON","panel_idx":0,
                 "panel":cls.base,"source_file_sha256":"ab"*32,
                 "measured_panel_sha256":s.sit1.hash_array(cls.base,dtype="<f4")}

    def test_frozen_protocol_uses_authentic_sources(self):
        self.assertEqual(s.P["benchmark"]["base_panels"],24)
        self.assertEqual(s.P["benchmark"]["total_controlled_cases"],288)

    def test_exact_injection_conditions(self):
        self.assertEqual(len(s.TYPE),3)
        self.assertEqual(s.TYPE[1],"WEAK_ARTIFICIAL_DRIFT_INJECTION")

    def test_injected_weak_line_not_native_astronomy(self):
        inp=self.base.copy()
        changed,truth=s.inject_from_real(inp,s.TYPE[1],np.random.default_rng(12))
        self.assertIsNotNone(truth)
        self.assertEqual(len(truth),16)
        self.assertTrue(np.array_equal(inp,self.base))
        self.assertFalse(np.array_equal(changed,self.base))

    def test_control_has_no_known_truth(self):
        changed,truth=s.inject_from_real(self.base,s.TYPE[0],np.random.default_rng(12))
        self.assertIsNone(truth)
        self.assertTrue(np.array_equal(changed,self.base))

    def test_rfi_is_derived_synthetic_control(self):
        changed,truth=s.inject_from_real(self.base,s.TYPE[2],np.random.default_rng(12))
        self.assertIsNone(truth)
        self.assertFalse(np.array_equal(changed,self.base))

    def test_fixed_random_seed_reproducible(self):
        a,ta=s.inject_from_real(self.base,s.TYPE[1],np.random.default_rng(7))
        b,tb=s.inject_from_real(self.base,s.TYPE[1],np.random.default_rng(7))
        self.assertTrue(np.array_equal(a,b))
        self.assertTrue(np.array_equal(ta,tb))

    def test_no_data_aliasing_in_derived_row(self):
        m=s.matched_row(self.base+0.3,self.row)
        self.assertNotEqual(m["measured_panel_sha256"],self.row["measured_panel_sha256"])
        self.assertEqual(m["source_file_sha256"],self.row["source_file_sha256"])
        self.assertTrue(np.array_equal(self.row["panel"],self.base))

    def test_frozen_threshold_is_negative_only(self):
        vals=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
        self.assertEqual(s.choose_threshold(vals),20)

    def test_scoring_truth_is_evaluation_only(self):
        dummy=[]
        for i in range(100):
            dummy.append({"kind":s.TYPE[1] if i%3==0 else s.TYPE[0],
                          "split":"development" if i<40 else "evaluation",
                          "localized":True,
                          "blind_ridge_statistic":float(i),
                          "frozen_shadow_score":float(i)})
        res=s.assess(dummy)
        for metric in ("blind_ridge_statistic","frozen_shadow_score"):
            self.assertIn("threshold_fixed_on_development",res[metric])
            self.assertEqual(res[metric]["evaluation_cases"],60)

    def test_four_replicates_and_six_scans(self):
        d=[]
        for scan in ("0061","0062","0063","0064","0065","0066"):
            for k in range(4):
                row=copy.copy(self.row)
                row["scan"]=scan;row["panel_idx"]=k
                d.append(row)
        x=s.build_test_cases(d)
        self.assertEqual(len(x),288)
        self.assertEqual(sum(z["kind"]==s.TYPE[1] for z in x),96)

    def test_injection_truth_not_in_radio_row(self):
        d=[]
        for scan in ("0061","0062","0063","0064","0065","0066"):
            for k in range(4):
                row=copy.copy(self.row)
                row["scan"]=scan;row["panel_idx"]=k
                d.append(row)
        x=s.build_test_cases(d)
        for z in x:
            self.assertNotIn("truth_for_evaluation_only",z["row"])
            self.assertNotIn("kind",z["row"])

    def test_no_real_rfi_or_sky_truth_claim(self):
        self.assertFalse(s.P["gates"]["retrocausality_claim"])
        self.assertFalse(s.P["gates"]["external_custodian"])
        self.assertFalse(s.P["gates"]["promotion_allowed"])
        self.assertEqual(s.P["data"]["independent_observational_groups"],1)

    def test_production_unchanged(self):
        self.assertFalse(s.P["gates"]["promotion_allowed"])
        self.assertFalse(s.P["gates"]["retrocausality_claim"])

if __name__=="__main__":
    unittest.main()
