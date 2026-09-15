import unittest

from experiments.luci.run_tachyon_star_t2d_ztf import (
    PROVENANCE_TARGETS,
    SENSITIVITY_TARGETS,
    exact_metadata_match,
    sensitivity_floor,
)


class T(unittest.TestCase):
    def test_target_cardinalities(self):
        self.assertEqual(sum(len(v) for v in SENSITIVITY_TARGETS.values()), 6)
        self.assertEqual(sum(len(v) for v in PROVENANCE_TARGETS.values()), 3)
        self.assertEqual(set(SENSITIVITY_TARGETS), {"TS-T2B-ZTF-03", "TS-T2B-ZTF-09", "TS-T2B-ZTF-10"})
        self.assertEqual(set(PROVENANCE_TARGETS), {"TS-T2B-ZTF-16", "TS-T2B-ZTF-39"})

    def test_sensitivity_floor_monotonic_tail(self):
        self.assertEqual(sensitivity_floor({8.0: False, 10.0: True, 12.0: True, 15.0: True, 20.0: True}), 10.0)
        self.assertEqual(sensitivity_floor({8.0: True, 10.0: False, 12.0: True, 15.0: True, 20.0: True}), 12.0)
        self.assertIsNone(sensitivity_floor({8.0: False, 10.0: False, 12.0: False, 15.0: False, 20.0: False}))

    def test_exact_metadata_match(self):
        target = {"filefracday": "20180413268183", "pid": "467268180115", "field": "755", "ccdid": "1", "qid": "2", "filtercode": "zr"}
        row = {"filefracday": "20180413268183", "pid": "467268180115", "field": "000755", "ccdid": "01", "qid": "2", "filtercode": "zr"}
        self.assertTrue(exact_metadata_match(row, target))
        self.assertFalse(exact_metadata_match({**row, "qid": "3"}, target))
        self.assertFalse(exact_metadata_match({**row, "filefracday": "20180413268184"}, target))

    def test_pid_optional_only_if_service_omits_it(self):
        target = {"filefracday": "20191026441655", "pid": "1028441703315", "field": "605", "ccdid": "9", "qid": "2", "filtercode": "zr"}
        row = {"filefracday": "20191026441655", "field": "605", "ccdid": "9", "qid": "2", "filtercode": "zr"}
        self.assertTrue(exact_metadata_match(row, target))


if __name__ == "__main__":
    unittest.main()
