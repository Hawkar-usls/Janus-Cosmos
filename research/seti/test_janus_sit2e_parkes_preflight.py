"""Offline SIT2E preflight contracts. No synthetic files can qualify archive."""
import unittest
import janus_sit2e_parkes_preflight as p

class SIT2EPreflightContracts(unittest.TestCase):
    def test_source_dates_are_real_distinct(self):
        self.assertEqual([x["date"] for x in p.P["archive_sessions"]],
                         ["2025-07-31","2025-09-12","2025-10-05"])
    def test_october_is_future_holdout(self):
        self.assertEqual(p.P["archive_sessions"][-1]["role"],"future_blinded_evaluation")
    def test_only_official_host_qualifies(self):
        self.assertTrue(p.safe_url("https://bldata.berkeley.edu/ATLAS/PKS_Oct_05_2025/test.h5",p.BASE))
        self.assertFalse(p.safe_url("https://badsite.test/ATLAS/PKS_Oct_05_2025/test.h5",p.BASE))
    def test_escaped_path_cannot_leak(self):
        self.assertFalse(p.safe_url("https://bldata.berkeley.edu/other/file.h5",
                                    "https://bldata.berkeley.edu/ATLAS/"))
    def test_HDF5_file_classes(self):
        self.assertEqual(p.product_hint("Parkes_60900_12345_3IATLAS_S_mid.h5"),"mid")
        self.assertEqual(p.product_hint("Parkes_60900_12345_3IATLAS_R_fine.h5"),"fine")
        self.assertEqual(p.product_hint("Parkes_60900_12345_3IATLAS_S_time.h5"),"time")
    def test_role_hint_requires_explicit_label(self):
        self.assertEqual(p.role_hint("Parkes_60900_12345_3IATLAS_S_mid.h5"),"SOURCE_NAME_ONLY")
        self.assertEqual(p.role_hint("Parkes_60900_12345_3IATLAS_R_mid.h5"),"REFERENCE_NAME_ONLY")
        self.assertEqual(p.role_hint("Parkes_60900_12345_3IATLAS_mid.h5"),"UNKNOWN")
    def test_html_parser(self):
        x=p.LinkParser();x.feed('<a href="file_S_mid.h5">a</a><a href="file_R_mid.h5">b</a>')
        self.assertEqual(x.hrefs,["file_S_mid.h5","file_R_mid.h5"])
    def test_only_metadata_and_sha(self):
        r=p.P["scientific_gate"]
        self.assertFalse(r["raw_power_used_in_preflight"])
        self.assertFalse(r["October_scores_read_before_threshold_lock"])
        self.assertFalse(r["thresholds_fitted"])
        self.assertFalse(r["science_qualified"])
    def test_download_bounded(self):
        self.assertLessEqual(p.MAXFILE,90_000_000)
        self.assertEqual(p.P["resources"]["selected_files_per_date_max"],2)
    def test_original_source_and_reference_not_pointing_certificate(self):
        self.assertFalse(p.P["scientific_gate"]["external_telescope_pointing_authenticated"])

if __name__=="__main__":unittest.main()
