"""SIT2D scientific guards, OFFLINE SYNTHETIC FIXTURES ONLY.

These cannot replace the source-verified actual GBT HDF5 examination.
Exact historical R41S1H1 donor source is imported unmodified.
"""
import unittest
import numpy as np
import janus_sit2d_archive_donors as s

class SIT2DContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng=np.random.default_rng(211)
        cls.mats=[rng.normal(10.,1.,(16,128)).astype(np.float32) for _ in range(4)]
        cls.group={
            "scans":["0061","0062","0063","0064"],
            "window":0,
            "four_original_sha256":["0"*64,"1"*64,"2"*64,"3"*64],
            "matrix":cls.mats
        }

    def test_donor_is_real_archive_module(self):
        self.assertTrue(hasattr(s.donor,"templates"))
        self.assertTrue(hasattr(s.donor,"score"))
        self.assertEqual(s.Q["archive_audit"]["exact_donor_blob_sha1"],
                         "43683dee60433ea2612bb1a433d8e7afc9704494")

    def test_offline_template_bank_original_physics(self):
        b,p=s.donor.templates(2861.02294921875,18.725)
        self.assertEqual(b.shape[1:],(16,49))
        self.assertEqual(b.shape[0],270)
        self.assertEqual(p.shape,(270,4))

    def test_all_valid_paths_have_49_support(self):
        self.assertTrue(np.all(s.BANK>=24))
        self.assertTrue(np.all(s.BANK<=103))

    def test_blind_selector_invariant_to_truth_label(self):
        p,v=s.select_blind_path(self.mats)
        self.assertEqual(p.shape,(16,))
        self.assertTrue(np.isfinite(v))

    def test_extracted_four_roles_have_exact_native_49_channels(self):
        p,_=s.select_blind_path(self.mats)
        cube=s.donor_cutout(self.mats,p)
        self.assertEqual(cube.shape,(4,16,49))

    def test_h1_original_physical_and_fixed_outputs_finite(self):
        p,_=s.select_blind_path(self.mats)
        template=s.donor.templates(2861.02294921875,18.725)
        output=s.donor_evidence(self.mats,p,template,2861.02294921875)
        self.assertTrue(np.isfinite(output["physical"]))
        self.assertTrue(np.isfinite(output["fixed"]))

    def test_direct_four_role_same_observations(self):
        p,_=s.select_blind_path(self.mats)
        score=s.direct_four_role(self.mats,p)
        self.assertTrue(np.isfinite(score))

    def test_real_background_unchanged_on_no_injection(self):
        rng=np.random.default_rng(20)
        x,t=s.overlay(self.mats,s.KIND[0],rng)
        self.assertIsNone(t)
        self.assertTrue(all(np.array_equal(a,b) for a,b in zip(x,self.mats)))

    def test_all_modified_overlay_cases_are_explicitly_synthetic(self):
        for label in s.KIND[1:]:
            x,truth=s.overlay(self.mats,label,np.random.default_rng(31))
            self.assertFalse(all(np.array_equal(a,b) for a,b in zip(x,self.mats)))
            if label in s.KIND[2:]:
                self.assertEqual(len(truth),16)

    def test_injection_only_changes_both_ON_roles(self):
        x,t=s.overlay(self.mats,s.KIND[2],np.random.default_rng(34))
        self.assertIsNotNone(t)
        self.assertTrue(np.array_equal(x[1],self.mats[1]))
        self.assertTrue(np.array_equal(x[3],self.mats[3]))
        self.assertFalse(np.array_equal(x[0],self.mats[0]))
        self.assertFalse(np.array_equal(x[2],self.mats[2]))

    def test_new_receiver_files_not_previously_examined(self):
        q=s.Q["real_radio"]
        self.assertEqual([x["scan"] for x in q["evaluation_files"]],
                         ["0013","0014","0015","0016","0017","0018"])
        self.assertTrue(all(x["filename"].startswith("blc22") for x in q["evaluation_files"]))

    def test_dev_source_requires_original_sha256(self):
        self.assertEqual(len(s.Q["real_radio"]["dev_files"]),6)
        self.assertTrue(all(x["filename"].startswith("blc20") for x in s.Q["real_radio"]["dev_files"]))

    def test_role_groups_from_real_like_scan_labels(self):
        rows=[]
        for scan,role in (("0013","ON"),("0014","OFF"),("0015","ON"),
                          ("0016","OFF"),("0017","ON"),("0018","OFF")):
            for win in range(4):
                rows.append({"scan":scan,"role":role,"panel_idx":win,
                             "panel":self.mats[0],
                             "source_file_sha256":scan*9})
        g=s.role_groups(rows)
        self.assertEqual(len(g),12)
        self.assertEqual(g[0]["scans"],["0013","0014","0015","0016"])

    def test_generation_144_with_72_controls(self):
        x=s.create_cases([self.group]*12,np.random.default_rng(35))
        self.assertEqual(len(x),144)
        self.assertEqual(sum(z["kind_for_evaluation_only"]==s.KIND[2] for z in x),36)
        self.assertEqual(sum(z["kind_for_evaluation_only"] in s.KIND[:2] for z in x),72)

    def test_tachyon_replay_no_extra_information(self):
        p,_=s.select_blind_path(self.mats)
        op=s.tachyon_evidence("",p,0.123,self.group["four_original_sha256"])
        self.assertGreaterEqual(op,2)

    def test_no_production_update_and_no_physical_retrocausal_claim(self):
        self.assertTrue(s.Q["gates"]["no_prod_changes"])
        self.assertTrue(s.Q["gates"]["no_train_on_eval"])
        self.assertFalse(s.Q["gates"]["external_replication"])

if __name__=="__main__":
    unittest.main()
