import unittest
from workspace.cousteau_polar_observability import profile_gate, lag_diagnostic, polar_audit, latitude_band, longitude_sector

def make(qc="1", **opts):
    p={"wmo":"1900001","position_qc":qc,"data_mode":"D",
       "latitude_deg":-67.2,"longitude_deg":40.2,
       "source_url":"https://data-argo.ifremer.fr/dac/abc/1900001/profiles/D1900001_010.nc",
       "source_sha256":"a"*64,"index_date_utc":"20230620000000",
       "sample_levels":[{"pressure_dbar":100,"temp_c":-1.2,"psal_psu":34.3}]}
    p.update(opts)
    return p

class PolarObservabilityTest(unittest.TestCase):
    def test_qc1_admitted(self):
        self.assertEqual(profile_gate(make())["tier"],"ADMIT_GEOGRAPHIC_QC1")
    def test_qc8_physical_not_exact_map(self):
        x=profile_gate(make("8"))
        self.assertTrue(x["physical_samples_present"])
        self.assertFalse(x["geography_usable"])
        self.assertEqual(x["ice_presence"],"UNKNOWN_WITHOUT_INDEPENDENT_SEA_ICE_PRODUCT")
    def test_qc9_no_map(self):
        self.assertFalse(profile_gate(make("9"))["geography_usable"])
    def test_qc8_not_assumed_ice(self):
        self.assertFalse(profile_gate(make("8"))["geography_usable"])
    def test_no_provenance_no_science(self):
        self.assertEqual(profile_gate(make(source_sha256="fake"))["tier"],"REJECT_PHYSICAL_QC_OR_PROVENANCE")
    def test_bad_physics(self):
        self.assertFalse(profile_gate(make(sample_levels=[{"temp_c":None}]))["physical_samples_present"])
    def test_missing_telemetry_delay_not_zero(self):
        x=lag_diagnostic(make())
        self.assertIsNone(x["lag_days"])
        self.assertEqual(x["state"],"TELEMETRY_DELAY_UNKNOWN__NOT_ZERO")
    def test_negative_telemetry_lag(self):
        x=lag_diagnostic(make(observed_at_utc="2025-01-02T00:00:00Z",received_at_utc="2025-01-01T00:00:00Z"))
        self.assertEqual(x["state"],"CLOCK_OR_PROVENANCE_CONFLICT")
    def test_lag_correct(self):
        x=lag_diagnostic(make(observed_at_utc="2025-01-01T00:00:00Z",received_at_utc="2025-01-11T12:00:00Z"))
        self.assertAlmostEqual(x["lag_days"],10.5)
    def test_no_claim_from_empty(self):
        x=polar_audit({"profiles":[]})
        self.assertEqual(x["profile_count"],0)
        self.assertIn("NOT_PROOF",x["absence_interpretation"])
    def test_qc1_coverage_only(self):
        x=polar_audit({"profiles":[make("1"),make("8"),make("9")]})
        self.assertEqual(x["observed_geographic_cells"]["[-70,-60)_[0,60)"],1)
        self.assertEqual(x["admission_tiers"]["ADMIT_GEOGRAPHIC_QC1"],1)
        self.assertEqual(x["geographic_cells_with_no_qc1_sample_in_subset"],23)
    def test_band_boundary(self):
        self.assertEqual(latitude_band(-50),"[-60,-50]")
        self.assertEqual(longitude_sector(180),"[-180,-120)")
    def test_missing_wmo_does_not_generate_zeros(self):
        x=polar_audit({"profiles":[make(wmo="")]})
        self.assertEqual(x["profile_count"],1)
        self.assertEqual(x["distinct_wmo"],0)

if __name__=="__main__":
    unittest.main(verbosity=2)
