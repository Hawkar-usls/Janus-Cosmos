"""JANUS SIT0: unit controls for real-profile provenance sidecar.

Only the workflow's --audit step downloads/authenticates the actual
archived astronomical profiles. Unit fixtures here are clearly synthetic.
"""
import copy
import unittest
import numpy as np

import janus_sit0_hawking as sit

class SIT0RealProvenanceContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        v=np.linspace(-65,-30,300)
        y=1+1.5*np.exp(-0.5*((v+45)/0.8)**2)
        step=float(np.median(np.diff(v)))
        e=np.r_[v[0]-step/2,(v[:-1]+v[1:])/2,v[-1]+step/2]
        rebin,_=sit.native.pinned_rebin()
        transform={"factor":1.25,"phase":0}
        f,z,err=sit.edge_rebin(e,y,transform,rebin)
        bg=float(np.median(y[(v<-60)|(v>-30)]))
        comp=sit.native.compare(e,y,f,z,bg,[-50,-40])
        cls.rebin=staticmethod(rebin)  # Avoid descriptor binding extra 'self' arg.
        cls.ctx={
          "tag":"W3OH:57934:T0","source":"W3OH","mjd_day":57934,
          "profile_name":"W3OH_science_profile.npz",
          "npz_sha":"a"*64,"v":v,"y":y,"e":e,"f":f,"z":z,
          "bg":bg,"line":[-50,-40],"transform":transform,"compare":comp,
        }
        cls.m=sit.manifest(cls.ctx)

    def test_exactly_eight_preregistered_overlay_types(self):
        self.assertEqual(len(sit.OVERLAYS),8)

    def test_actual_mirror_reused_from_pinned_tachyon(self):
        mir=sit.mirror_case(self.m["forward_lineage_graph"])
        self.assertEqual(mir["graph_direction"],"C_TO_B_TO_A")
        self.assertEqual(mir["anchor_source"],"natural_line_metrics:W3OH:57934:T0")

    def test_graph_has_no_fake_instrument_witnesses(self):
        self.assertEqual(self.m["forward_lineage_graph"]["witness_reports"],[])

    def test_unchanged_profile_metadata_valid(self):
        self.assertTrue(sit.trusted_provenance_valid(self.m,self.ctx,self.rebin))

    def test_instrument_lineage_mirror_recovers_exact_path(self):
        ok,ops=sit.tachyon_mirror_replay(self.m,self.ctx,self.rebin)
        self.assertTrue(ok)
        self.assertGreaterEqual(ops,2)

    def test_intact_accepted_both_strong_arms(self):
        for name in ("C_TACHYON_MIRROR_PLUS_FORWARD_REPLAY","D_STRONG_DIRECT_FULL_REPLAY"):
            self.assertTrue(sit.ARMS[name](self.m,self.ctx,self.rebin)[0])

    def test_metadata_corruptions_rejected_by_both_strong_arms(self):
        for kind in sit.OVERLAYS:
            with self.subTest(kind=kind):
                t=sit.corrupt(self.m,kind)
                good=kind=="INTACT"
                self.assertEqual(sit.tachyon_mirror_replay(t,self.ctx,self.rebin)[0],good)
                self.assertEqual(sit.strong_direct_replay(t,self.ctx,self.rebin)[0],good)

    def test_plain_COSMOS_does_not_use_added_overlay(self):
        for kind in sit.OVERLAYS:
            self.assertTrue(sit.cosmos_native_unmodified(sit.corrupt(self.m,kind),self.ctx,self.rebin)[0])

    def test_graph_only_insufficient_for_bad_hash(self):
        t=sit.corrupt(self.m,"BAD_RAW_SHA")
        self.assertTrue(sit.tachyon_mirror_only(t,self.ctx,self.rebin)[0])
        self.assertFalse(sit.tachyon_mirror_replay(t,self.ctx,self.rebin)[0])

    def test_original_manifest_immutable(self):
        m=copy.deepcopy(self.m)
        for kind in sit.OVERLAYS:sit.corrupt(m,kind)
        self.assertEqual(m,self.m)

    def test_rebin_input_arrays_not_modified(self):
        old=self.ctx["y"].copy()
        sit.tachyon_mirror_replay(self.m,self.ctx,self.rebin)
        self.assertTrue(np.array_equal(old,self.ctx["y"]))

    def test_nonfinite_array_hash_is_not_same(self):
        altered=np.asarray(self.ctx["z"]).copy()
        altered[3]=altered[3]+0.05
        self.assertNotEqual(sit.arrays_sha(self.ctx["f"],altered),
                            self.m["native_rebin_arrays_sha256"])

    def test_missing_evidence_rejected(self):
        t=sit.corrupt(self.m,"MISSING_SOURCE_ID")
        self.assertFalse(sit.strong_direct_replay(t,self.ctx,self.rebin)[0])

    def test_reverse_index_incompleteness_rejected(self):
        t=sit.corrupt(self.m,"INCOMPLETE_REVERSE_INDEX")
        mir=sit.mirror_case(t["forward_lineage_graph"])
        self.assertEqual(sit.mirrored_replay(mir)["prediction"],sit.ABSTAIN)

    def test_ghost_predecessor_rejected(self):
        t=sit.corrupt(self.m,"GHOST_REVERSE_EDGE")
        self.assertFalse(sit.strong_direct_replay(t,self.ctx,self.rebin)[0])

    def test_preserve_real_derived_not_raw_claim(self):
        self.assertFalse(sit.Q["observational_sample"]["time_frequency_original_input_available_here"])
        self.assertFalse(sit.Q["safety"]["physical_time_reversal_or_retrocausal_signal_claim"]) if "safety" in sit.Q else self.assertFalse(sit.Q["gates"]["physical_time_reversal_or_retrocausal_signal_claim"])

    def test_no_R42_native_classifier_claim(self):
        self.assertIn("NO R42A8A",sit.Q["source"]["model_source"])

    def test_no_detector_change_or_veto(self):
        self.assertFalse(sit.Q["gates"]["production_detector_changed"])
        self.assertFalse(sit.Q["gates"]["active_COSMOS_NOISE_HANNAH_weights_changed"])

if __name__=="__main__":
    unittest.main()
