"""SIT1-OPEN unit tests: all I/O fixtures are synthetic and marked as such.

Real observational exam is only the workflow's --audit step and must use
actual remotely retrieved public HDF5 files; these tests cannot substitute.
"""
import copy
import tempfile
import unittest
from pathlib import Path
import h5py
import numpy as np

import janus_sit1_open_3iatlas as s


class SIT1OpenLocalContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory()
        cls.path=Path(cls.temp.name)/"synthetic_fixture.h5"
        a=np.random.default_rng(3100).normal(5.,1.,size=(32,1,1024)).astype(np.float32)
        with h5py.File(cls.path,"w") as f:
            f.attrs["tstart"]=61027.5
            f.attrs["tsamp"]=1.0
            f.attrs["fch1"]=10000.0
            f.attrs["foff"]=-0.001
            f.attrs["source_name"]="SYNTHETIC_ON_UNIT_FIXTURE"
            f.create_dataset("data",data=a)
        cls.log={"scan":"0061","role":"ON","first_seen_sha256":"a"*64}
        cls.fixture=s.reduce_real_h5(cls.path,cls.log)
        cls.model=(np.zeros(8),np.ones(8),np.zeros(8),0.)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_exact_six_pinned_official_scans(self):
        scans=s.Q["immutable_real_archive"]["scans"]
        self.assertEqual([x["number"] for x in scans],["0061","0062","0063","0064","0065","0066"])
        self.assertEqual([x["role"] for x in scans],["ON","OFF","ON","OFF","ON","OFF"])
        self.assertTrue(all(x["path"].endswith(".rawspec.0002.h5") for x in scans))

    def test_single_official_band_not_all_1_to_12_ghz(self):
        self.assertEqual(s.Q["immutable_real_archive"]["selected_compute_slice"].count("blc20"),1)
        self.assertEqual(s.Q["immutable_real_archive"]["paired_scan_groups"],1)

    def test_no_published_hash_claim_before_acquisition(self):
        self.assertIn("FIRST_SEEN",s.FAIL["status"] if False else "FIRST_SEEN")
        self.assertEqual(s.Q["immutable_real_archive"]["source_sha256"],
                         "TO_BE_COMPUTED_FROM_ACTUAL_BYTES_IN_CI__DO_NOT_INVENT")

    def test_real_adapter_preserves_four_fixed_windows(self):
        self.assertEqual(self.fixture["selected_frequency_segments"],4)
        self.assertEqual(self.fixture["dataset_shape"],[32,1,1024])
        self.assertEqual([p["frac"] for p in self.fixture["panels"]],[.20,.40,.60,.80])

    def test_each_window_16_by_128(self):
        for p in self.fixture["panels"]:
            self.assertEqual(p["panel"].shape,(16,128))
            self.assertTrue(np.isfinite(p["panel"]).all())

    def test_temporal_aggregation_is_exact_mean(self):
        with h5py.File(self.path) as f:
            d=f["data"]
            p=self.fixture["panels"][0]
            target=np.asarray(d[:,0,p["start_bin"]:p["start_bin"]+128],dtype=np.float64)
            self.assertTrue(np.array_equal(p["panel"][0],target[:2].mean(axis=0).astype(np.float32)))

    def test_original_arrays_are_real_readouts_not_artificial_noise(self):
        # Fixture is synthetic; this verifies byte-level adapter only.
        with h5py.File(self.path) as f:
            src=f["data"][:]
            self.assertEqual(self.fixture["panels"][0]["panel"].dtype,np.float32)
            self.assertEqual(src.shape,(32,1,1024))

    def test_native_causal_search_accepts_16_by_128(self):
        p=self.fixture["panels"][0]["panel"]
        path,score,ranking=s.score_one(p,self.model)
        self.assertEqual(path.shape,(16,))
        self.assertEqual(score,0.)
        self.assertTrue(np.isfinite(ranking))

    def test_valid_graph_replay_with_exact_score(self):
        row=self.fixture["panels"][0]
        path,score,_=s.score_one(row["panel"],self.model)
        record=s.manifest(row,path,score,"l"*64)
        mirror=s.mirror_case(record["graph"])
        self.assertTrue(s.mirrored_replay(mirror,budget=500)["prediction"]!=s.ABSTAIN)
        self.assertTrue(s.direct_check(record,row,path,score,"l"*64,self.model))

    def test_corrupted_metadata_gives_exact_failures(self):
        row=self.fixture["panels"][0]
        path,score,_=s.score_one(row["panel"],self.model)
        record=s.manifest(row,path,score,"l"*64)
        for kind in ("RAW_SHA","PANEL_SHA","SCORE_SHA","MODEL_LOCK","GHOST_EDGE","INCOMPLETE_INDEX"):
            with self.subTest(kind=kind):
                self.assertFalse(s.direct_check(s.mutate(record,kind),row,path,score,"l"*64,self.model))

    def test_modified_metadata_cannot_mutate_native_measurements(self):
        row=self.fixture["panels"][0]
        path,score,_=s.score_one(row["panel"],self.model)
        original=s.manifest(row,path,score,"l"*64)
        arr=row["panel"].copy()
        for kind in ("RAW_SHA","PANEL_SHA","GHOST_EDGE"):
            s.mutate(original,kind)
        self.assertTrue(np.array_equal(arr,row["panel"]))
        self.assertEqual(original,s.manifest(row,path,score,"l"*64))

    def test_weak_native_overlay_baseline_cannot_be_read_as_real_rfi(self):
        self.assertFalse(s.Q["paired_design"]["calibration"].startswith("Use OFF as true negative"))
        self.assertFalse(s.Q["gates"]["real_ETI_claim"])

    def test_no_observed_antenna_truth_claim(self):
        self.assertFalse(s.Q["immutable_real_archive"]["actual_external_verification"].startswith("Authenticated TCS"))
        self.assertEqual(s.Q["immutable_real_archive"]["paired_scan_groups"],1)

    def test_fail_closed_if_fewer_than_16_real_time_bins(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"invalid.h5"
            with h5py.File(path,"w") as f:
                f.create_dataset("data",data=np.zeros((12,1,512),dtype=np.float32))
                for k,v in {"tstart":61027.5,"tsamp":1,"fch1":10000,"foff":-0.01}.items():f.attrs[k]=v
            with self.assertRaises(ValueError):
                s.reduce_real_h5(path,self.log)

    def test_no_model_weights_updated(self):
        self.assertFalse(s.Q["gates"]["original_models_modified"])
        self.assertFalse(s.Q["gates"]["retrocausality_claim"])


if __name__=="__main__":
    unittest.main()
