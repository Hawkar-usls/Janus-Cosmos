"""No-network tamper tests for historical Palomar archive reuse."""
import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import verify_palomar_archived_ledger as m


class VerifyPalomarArchive(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.valid = json.loads((ROOT / "domains/space/PALOMAR_OBSERVATION_DAY_LEDGER_v0.1.json").read_text())

    def test_existing_archive_is_not_represented_as_fresh_remote_materialization(self):
        d = m.verify(self.valid)
        self.assertEqual(d["status"], "ARCHIVE_VERIFIED_OFFLINE_NOT_REMOTE_REFRESH")
        self.assertFalse(d["upstream_current_state_verified"])
        self.assertEqual(d["authority_delta"], 0)

    def test_tampered_archive_never_becomes_a_daily_success(self):
        variants = []
        v = copy.deepcopy(self.valid); v["source_freeze"]["upstream_commit"] = "f"*40; variants.append(v)
        v = copy.deepcopy(self.valid); v["source_freeze"]["stage_sha256_uncompressed"] = "0"*64; variants.append(v)
        v = copy.deepcopy(self.valid); v["rows"][0]["candidate_count"] += 1; variants.append(v)
        v = copy.deepcopy(self.valid); v["rows"][1]["plate_ids"] = list(v["rows"][0]["plate_ids"]); variants.append(v)
        v = copy.deepcopy(self.valid); v["rows"][0]["observed_flag"] = False; variants.append(v)
        v = copy.deepcopy(self.valid); v["rows"][0]["sky_position"]["plate_centers_ra_dec"] = []; variants.append(v)
        v = copy.deepcopy(self.valid); v["rows"] = v["rows"][:-1]; variants.append(v)
        v = copy.deepcopy(self.valid); v["rows"][0]["date"] = "not-a-date"; variants.append(v)
        for row in variants:
            with self.subTest(index=variants.index(row)), self.assertRaises(m.ArchiveHold):
                m.verify(row)

    def test_daily_schedule_is_offline_and_remote_refresh_explicit(self):
        flow = (ROOT / ".github/workflows/palomar-day-materialization.yml").read_text()
        self.assertIn("schedule:", flow)
        self.assertIn("if: github.event_name == 'workflow_dispatch'", flow)
        self.assertIn("tools/verify_palomar_archived_ledger.py", flow)
        self.assertIn("ARCHIVE_VERIFIED_OFFLINE_NOT_REMOTE_REFRESH", (ROOT / "tools/verify_palomar_archived_ledger.py").read_text())


if __name__ == "__main__":
    unittest.main()
