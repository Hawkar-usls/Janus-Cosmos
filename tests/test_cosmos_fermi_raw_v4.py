"""Offline V4 source integrity, FITS parsing and non-detection controls."""

import hashlib
import json
import math
import struct
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"workspace"))
sys.path.insert(0,str(ROOT/"scripts"))
from janus_cosmos_fermi_raw_v4 import (AuditError, BLOCK, MANIFEST, MANIFEST_SCHEMA,
    _column_info, angular_sep_deg, read_fits, rows, run, sha256)
from janus_cosmos_fermi_joint_v4 import fit_grid, run as joint_run, toy_bursts
from janus_cosmos_simulation_audit import digest
from fetch_fermi_lat_v4 import trusted_url, verified

RAW_RESULT=ROOT/"data/cosmos/COSMOS-ALI-V4-FERMI-RAW-EVENT-DIAGNOSTICS-2026-10-18.json"
JOINT_RESULT=ROOT/"data/cosmos/COSMOS-ALI-V4-JOINT-TOY-LIKELIHOOD-2026-10-18.json"


def _card(key,value):
    if isinstance(value,str):v="'"+value+"'"
    else:v=str(value).rjust(20)
    return f"{key:<8}= {v:<20}".ljust(80).encode("ascii")


def _hdu(cards,data=b""):
    head=b"".join(_card(k,v) for k,v in cards)+b"END".ljust(80,b" ")
    return head.ljust(math.ceil(len(head)/BLOCK)*BLOCK,b" ")+data.ljust(math.ceil(len(data)/BLOCK)*BLOCK,b"\0")


class TestV4Raw(unittest.TestCase):
    def test_manifest_has_only_retrievable_shapinned_extended(self):
        m=json.loads(MANIFEST.read_text())
        self.assertEqual(m["schema"],MANIFEST_SCHEMA)
        self.assertEqual(len(m["bursts"]),3)
        self.assertEqual([x["grb"] for x in m["bursts"]],["090510","090902B","130427A"])
        for burst in m["bursts"]:
            for role,ext in (("event","EV00"),("spacecraft","SC00")):
                f=burst[role]
                self.assertTrue(f["url"].startswith("https://fermi.gsfc.nasa.gov/FTP/fermi/data/lat/queries/"))
                self.assertTrue(f["url"].endswith("_"+ext+".fits"))
                self.assertEqual(len(f["sha256"]),64)
                self.assertEqual(f["bytes"]%2880,0)
            self.assertEqual(burst["selection"],m["common_selection"])
        self.assertEqual(m["source_lag_constraint"],"NONE")

    def test_frozen_actual_raw_counts_and_claim_ceiling(self):
        r=json.loads(RAW_RESULT.read_text())
        m=json.loads(MANIFEST.read_text())
        self.assertEqual(r["counts"],{"bursts":3,"diagnostic_source_region_candidates":175,"raw_events_in_query_rois":6520})
        self.assertEqual([x["counts"]["after_diagnostic_selection"] for x in r["sources"]],[32,67,76])
        self.assertEqual(r["gates"]["RAW_EXTENDED_FITS"],"PASS_SHA256_VERIFIED")
        self.assertEqual(r["gates"]["MULTIBURST_PHYSICAL_LIKELIHOOD"],"BLOCKED_SOURCE_AND_INSTRUMENT_DEGENERACY")
        self.assertEqual(r["gates"]["LIV_OBSERVED"],"NOT_ESTABLISHED")
        self.assertEqual(r["gates"]["AMPS_FIREWALL"],"OPEN_UNCONFIRMED")
        self.assertEqual(r["manifest_sha256"],hashlib.sha256(json.dumps(m,sort_keys=True,separators=(',',':')).encode()).hexdigest())
        copied=dict(r);result_sha=copied.pop("result_sha256")
        self.assertEqual(result_sha,digest(copied))
        for burst in r["sources"]:
            self.assertEqual(burst["counts"]["after_diagnostic_selection"],len(burst["events"]))
            t=[x["met_s"] for x in burst["events"]]
            self.assertEqual(t,sorted(t))
            self.assertTrue(all(x["event_class_mask"]&128 for x in burst["events"]))
            self.assertTrue(all(1000<=x["energy_MeV"]<=300000 for x in burst["events"]))
            self.assertTrue(all(x["angular_sep_deg"]<=5 for x in burst["events"]))
            self.assertTrue(all(x["zenith_deg"]<=90 for x in burst["events"]))

    def test_read_small_binary_fits_and_bit_mask(self):
        primary=_hdu([("SIMPLE",True),("BITPIX",8),("NAXIS",0),("EXTEND",True)])
        table=_hdu([("XTENSION","BINTABLE"),("BITPIX",8),("NAXIS",2),
            ("NAXIS1",8),("NAXIS2",1),("PCOUNT",0),("GCOUNT",1),
            ("TFIELDS",2),("TTYPE1","ENERGY"),("TFORM1","E"),
            ("TTYPE2","EVENT_CLASS"),("TFORM2","32X"),("EXTNAME","EVENTS")],
            struct.pack(">fI",1234.5,128))
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"tiny.fits";p.write_bytes(primary+table)
            h=read_fits(p)
            self.assertEqual(len(h),2)
            rec=rows(h[1],{"ENERGY","EVENT_CLASS"})[0]
            self.assertAlmostEqual(rec["ENERGY"],1234.5)
            self.assertEqual(rec["EVENT_CLASS"],128)
            p.write_bytes((primary+table)[:-5])
            with self.assertRaises(AuditError):read_fits(p)

    def test_reject_non_fits_and_malformed_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"bogus.fits";p.write_bytes(b"garbage"*300)
            with self.assertRaises(AuditError):read_fits(p)
        with self.assertRaises(AuditError):_column_info({"TFIELDS":1,"TTYPE1":"BAD","TFORM1":"1P","NAXIS1":8})
        with self.assertRaises(AuditError):_column_info({"TFIELDS":1,"TTYPE1":"BAD","TFORM1":"1E","NAXIS1":6})

    def test_spherical_angle(self):
        self.assertAlmostEqual(angular_sep_deg(5,20,5,20),0,places=6)
        self.assertAlmostEqual(angular_sep_deg(0,0,180,0),180,places=6)

    def test_fail_closed_missing_binary(self):
        m=json.loads(MANIFEST.read_text())
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(AuditError):run(m,Path(tmp))

    def test_download_allowlist_and_bad_checksum(self):
        self.assertTrue(trusted_url("https://fermi.gsfc.nasa.gov/FTP/fermi/data/lat/queries/L261018113415F357373F80_EV00.fits"))
        self.assertFalse(trusted_url("https://fermi.gsfc.nasa.gov.evil.org/FTP/fermi/data/lat/queries/L1_EV00.fits"))
        self.assertFalse(trusted_url("http://fermi.gsfc.nasa.gov/FTP/fermi/data/lat/queries/L1_EV00.fits"))
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"bad.fits";p.write_bytes(b"bad")
            self.assertFalse(verified(p,{"bytes":3,"sha256":"0"*64}))


class TestV4Joint(unittest.TestCase):
    def test_recovery_only_given_intrinsic_lag_known_zero(self):
        positive=toy_bursts(seed=20261018,eta_true=0.8)
        fixed=fit_grid(positive,"known_zero")
        self.assertAlmostEqual(fixed["eta_mle"],0.8,delta=.02)
        free=fit_grid(positive,"unbounded")
        self.assertIsNone(free["eta_mle"])
        self.assertEqual(free["twice_delta_log_likelihood"],0.)

    def test_intrinsic_delay_can_mimic_identical_signal(self):
        intrinsic=toy_bursts(seed=20261018,eta_true=0.,source_beta=0.8)
        naive=fit_grid(intrinsic,"known_zero")
        self.assertAlmostEqual(naive["eta_mle"],0.8,delta=.02)
        self.assertIsNone(fit_grid(intrinsic,"unbounded")["eta_mle"])

    def test_per_burst_prior_not_empirically_established(self):
        synthetic=toy_bursts()
        with self.assertRaises(AuditError):fit_grid(synthetic,"gaussian_prior",beta_prior_sigma=None)
        self.assertGreater(fit_grid(synthetic,"gaussian_prior",beta_prior_sigma=.5)["twice_delta_log_likelihood"],0)

    def test_frozen_joint_null_only_toy_and_no_real_p(self):
        frozen=json.loads(JOINT_RESULT.read_text())
        self.assertEqual(frozen,joint_run())
        self.assertTrue(frozen["synthetic_null_calibration"]["not_an_observed_p_value"])
        self.assertEqual(frozen["gates"]["ASTROPHYSICAL_MULTIBURST_LRT"],"BLOCKED_NO_FORWARD_FOLDED_INSTRUMENT_AND_SOURCE_MODEL")
        self.assertIsNone(frozen["real_lat_likelihood_result"])
        self.assertIsNone(frozen["unbounded_intrinsic_lag"]["fit"]["eta_mle"])


if __name__=="__main__":unittest.main()
