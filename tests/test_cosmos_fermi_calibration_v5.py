"""Offline V5 integrity, calibration boundaries and source-lag null gates."""

import copy
import json
import math
import struct
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workspace"))
sys.path.insert(0, str(ROOT / "scripts"))

from janus_cosmos_fermi_calibration_v5 import (
    MANIFEST_PATH, V4_MANIFEST_PATH, V4_PATH, SCHEMA, AuditError,
    _tte_events, area_lookup, audit_caldb, checked_file, fits_large_tte,
    run, sha256,
)
from janus_cosmos_fermi_joint_v5 import run as joint_v5
from janus_cosmos_fermi_raw_v4 import HDU, BLOCK

OBS = ROOT / "data/cosmos/COSMOS-ALI-V5-CALDB-GBM-LAT-INSTRUMENT-PULSE-2026-10-19.json"
JOINT = ROOT / "data/cosmos/COSMOS-ALI-V5-SOURCE-LAG-TOY-SENSITIVITY-2026-10-19.json"


def fits_card(key, value):
    raw = "'" + value + "'" if isinstance(value, str) else str(value).rjust(20)
    return f"{key:<8}= {raw:<20}".ljust(80).encode("ascii")


def fits_hdu(cards, data=b""):
    header = b"".join(fits_card(k, v) for k, v in cards) + b"END".ljust(80, b" ")
    return header.ljust(math.ceil(len(header) / BLOCK) * BLOCK, b" ") + data.ljust(math.ceil(len(data) / BLOCK) * BLOCK, b"\0")


class TestV5Calibration(unittest.TestCase):
    def test_manifest_complete_and_source_correct(self):
        m = json.loads(MANIFEST_PATH.read_text())
        self.assertEqual(m["lat_caldb"]["irf_name"], "P8R3_SOURCE_V3")
        self.assertEqual(len(m["lat_caldb"]["caldb_files"]), 9)
        self.assertEqual([b["grb"] for b in m["gbm_bursts"]], ["090510", "090902B", "130427A"])
        self.assertEqual(sum(len(b["detectors"]) for b in m["gbm_bursts"]), 9)
        self.assertEqual(m["v4_frozen_observed_json_sha256"], sha256(V4_PATH))
        self.assertEqual(m["v4_frozen_manifest_sha256"], sha256(V4_MANIFEST_PATH))
        for burst in m["gbm_bursts"]:
            for detector in burst["detectors"]:
                self.assertEqual(set(detector["event_and_response_files"]), {"tte", "cspec", "rsp2"})
                for kind, f in detector["event_and_response_files"].items():
                    self.assertIn("/fermi/data/gbm/triggers/", f["origin_url"])
                    self.assertEqual(len(f["sha256"]), 64)
                    self.assertGreater(f["bytes"], 10000)

    def test_pinned_binary_rejects_path_escape_and_tampering(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            source = base / "fixture.fit"
            source.write_bytes(b"GBM")
            from hashlib import sha256 as digest
            item = {"name": source.name, "bytes": 3, "sha256": digest(b"GBM").hexdigest()}
            self.assertEqual(checked_file(base, item), source)
            source.write_bytes(b"GBN")
            with self.assertRaises(AuditError):
                checked_file(base, item)
            with self.assertRaises(AuditError):
                checked_file(base, {**item, "name": "../outside.fit"})

    def test_v2_caldb_masquerade_fails_closed(self):
        manifest = json.loads(MANIFEST_PATH.read_text())
        manifest["lat_caldb"]["irf_name"] = "P8R3_SOURCE_V2"
        with self.assertRaisesRegex(AuditError, "V2"):
            audit_caldb(manifest, Path("/tmp/nonexistent-caldb"))

    def test_big_tte_reader_scalar_channel_and_monotonicity(self):
        primary = fits_hdu([("SIMPLE", True), ("BITPIX", 8), ("NAXIS", 0), ("EXTEND", True)])
        bounds = fits_hdu([("XTENSION", "BINTABLE"), ("BITPIX", 8), ("NAXIS", 2),
            ("NAXIS1", 10), ("NAXIS2", 1), ("PCOUNT", 0), ("TFIELDS", 1),
            ("TTYPE1", "DUMMY"), ("TFORM1", "10B"), ("EXTNAME", "EBOUNDS")], b"\0" * 10)
        def make_events(vals):
            data = b"".join(struct.pack(">dh", t, p) for t, p in vals)
            return fits_hdu([("XTENSION", "BINTABLE"), ("BITPIX", 8), ("NAXIS", 2),
                ("NAXIS1", 10), ("NAXIS2", len(vals)), ("PCOUNT", 0),
                ("TFIELDS", 2), ("TTYPE1", "TIME"), ("TFORM1", "1D"),
                ("TTYPE2", "PHA"), ("TFORM2", "1I"), ("EXTNAME", "EVENTS")], data)
        gti = fits_hdu([("XTENSION", "BINTABLE"), ("BITPIX", 8), ("NAXIS", 2),
            ("NAXIS1", 16), ("NAXIS2", 1), ("PCOUNT", 0), ("TFIELDS", 2),
            ("TTYPE1", "START"), ("TFORM1", "1D"), ("TTYPE2", "STOP"),
            ("TFORM2", "1D"), ("EXTNAME", "GTI")], struct.pack(">dd", -20, 10))
        with tempfile.TemporaryDirectory() as temp:
            file = Path(temp) / "tte.fit"
            file.write_bytes(primary + bounds + make_events([(2., 7), (3., 9)]) + gti)
            hdus = fits_large_tte(file)
            t, p = _tte_events(hdus[2])
            self.assertEqual(t, [2., 3.])
            self.assertEqual(p, [7, 9])
            file.write_bytes(primary + bounds + make_events([(3., 7), (2., 9)]) + gti)
            with self.assertRaisesRegex(AuditError, "unsorted"):
                _tte_events(fits_large_tte(file)[2])

    def test_observed_v5_frozen_result_and_nonclaims(self):
        obs = json.loads(OBS.read_text())
        self.assertEqual(obs["schema"], SCHEMA)
        self.assertEqual(obs["counts"]["gbm_tte_events_all_selected_detectors"], 6962621)
        self.assertEqual(obs["counts"]["lat_v4_selected_roi_candidates"], 175)
        self.assertEqual(obs["counts"]["lat_irf_files"], 9)
        self.assertEqual([b["exact_gbm_trigger_met_s"] for b in obs["gbm_bursts"]],
                         [263607781.97109, 273582310.312714, 388741629.420158])
        self.assertGreater(obs["gbm_bursts"][2]["detectors"][0]["cspec_vs_tte_ebounds_max_relative_discrepancy"], .02)
        self.assertEqual(obs["gbm_bursts"][1]["detectors"][0]["rsp2_zero_duration_terminal_matrices"], 1)
        self.assertEqual(obs["scientific_gates"]["joint_observational_liv_hypothesis_test"], "BLOCKED_NO_VALIDATED_LIKELIHOOD")
        self.assertIsNone(obs["observational_joint_liv_p_value"])
        self.assertFalse(obs["black_hole_firewall_detection"])
        for b in obs["gbm_bursts"]:
            for detector in b["detectors"]:
                self.assertTrue(detector["rsp2_covers_prompt"])
                for band in detector["pulse_bands"].values():
                    self.assertEqual(len(band["time_bin_counts"]), round((2 if b["grb"] == "090510" else 40 if b["grb"] == "090902B" else 140) / detector["time_bin_width_s"]))
                    self.assertGreater(band["prompt_counts"], 0)

    def test_source_lag_never_becomes_observed_claim(self):
        obs = json.loads(OBS.read_text())
        result = joint_v5(obs, seed=20261019)
        self.assertEqual(result, json.loads(JOINT.read_text()))
        self.assertEqual(result["unconstrained_source_lag_true_identifiability"]["twice_delta_log_likelihood"], 0)
        self.assertIsNone(result["unconstrained_source_lag_true_identifiability"]["eta_mle"])
        self.assertIsNone(result["observational_liv_p_value"])
        self.assertEqual(result["actual_lat_event_count_used_in_likelihood"], 0)
        self.assertEqual(len(result["source_lag_prior_sensitivity_synthetic_only"]), 3)
        broken = copy.deepcopy(obs)
        broken["scientific_gates"]["joint_observational_liv_hypothesis_test"] = "DISCOVERED"
        with self.assertRaises(AuditError):
            joint_v5(broken)


if __name__ == "__main__":
    unittest.main()
