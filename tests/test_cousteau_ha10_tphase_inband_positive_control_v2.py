from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "workspace" / "cousteau_ha10_tphase_inband_positive_control_v2.py"
spec = importlib.util.spec_from_file_location(
    "cousteau_ha10_tphase_inband_positive_control_v2", MODULE
)
gate = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(gate)


def event(code: str, start_s: int, *, passed: bool = True, analyzed: bool = True):
    def iso(second: int) -> str:
        return f"2015-01-01T00:{second // 60:02d}:{second % 60:02d}Z"

    station = {
        "event_window": {
            "start_utc": iso(start_s),
            "end_utc": iso(start_s + 60),
        },
        "data_status": "ANALYZED" if analyzed else "BLOCKED_PAIR",
        "per_station_event_pass": passed,
    }
    return {
        "source_time_code": code,
        "stations": {
            "IM.H10N1..EDH": dict(station),
            "IM.H10S2..EDH": dict(station),
        },
        "positive_control_replicated_both_stations": passed,
    }


class CousteauHA10PositiveControlV2Tests(unittest.TestCase):
    def test_overlap_is_not_independent(self):
        first = event("A", 0)
        second = event("B", 30)
        self.assertTrue(gate.events_conflict(first, second))
        selected = gate.maximum_nonoverlapping_passing_subset([first, second])
        self.assertEqual(len(selected), 1)

    def test_touching_endpoints_are_independent(self):
        first = event("A", 0)
        second = event("B", 60)
        self.assertFalse(gate.events_conflict(first, second))
        selected = gate.maximum_nonoverlapping_passing_subset([first, second])
        self.assertEqual(selected, ["A", "B"])

    def test_three_raw_passes_with_overlap_cannot_promote(self):
        events = [event("A", 0), event("B", 30), event("C", 45)]
        verdict, diagnostic = gate.decide_control(events, 3)
        self.assertEqual(
            diagnostic["raw_candidate_events_passing_both_stations"], 3
        )
        self.assertEqual(
            diagnostic["independent_events_passing_both_stations"], 1
        )
        self.assertEqual(verdict, "FAIL_HA10_INBAND_TPHASE_PIPELINE_CONTROL")

    def test_three_nonoverlapping_passes_can_pass(self):
        events = [event("A", 0), event("B", 60), event("C", 120)]
        verdict, diagnostic = gate.decide_control(events, 3)
        self.assertEqual(
            diagnostic["independent_events_passing_both_stations"], 3
        )
        self.assertEqual(verdict, "PASS_HA10_INBAND_TPHASE_PIPELINE_CONTROL")

    def test_data_access_block_is_preserved(self):
        events = [
            event("A", 0, analyzed=True),
            event("B", 60, analyzed=True),
            event("C", 120, analyzed=False),
        ]
        verdict, diagnostic = gate.decide_control(events, 3)
        self.assertEqual(
            diagnostic["complete_events_analyzed_on_both_stations"], 2
        )
        self.assertEqual(
            verdict, "BLOCKED_POSITIVE_CONTROL_DATA_ACCESS_OR_RESPONSE"
        )

    def test_network_retry_budget_is_bounded(self):
        self.assertEqual(gate.HTTP_ATTEMPTS, 2)
        self.assertEqual(gate.HTTP_TIMEOUT_S, 15)
        self.assertLessEqual(
            gate.HTTP_ATTEMPTS * gate.HTTP_TIMEOUT_S + gate.HTTP_BACKOFF_S,
            31,
        )

    def test_runner_refuses_canonical_output(self):
        with self.assertRaisesRegex(
            RuntimeError, "CANONICAL_DATA_WRITE_FORBIDDEN_BY_POSITIVE_CONTROL_V2"
        ):
            gate.ensure_noncanonical_output(gate.DATA / "forbidden.json")

    def test_runner_allows_ephemeral_output(self):
        with tempfile.TemporaryDirectory() as td:
            gate.ensure_noncanonical_output(Path(td) / "run.json")

    def test_frozen_contracts_bind_original_protocol_addendum_and_negative(self):
        protocol, windows, frozen = gate.verify_frozen_contracts()
        self.assertEqual(
            protocol["gate_id"],
            "COUSTEAU_HA10_TPHASE_INBAND_ARRIVAL_POSITIVE_CONTROL_V1",
        )
        self.assertEqual(windows["status"], "WINDOWS_FROZEN_READY_FOR_FFT")
        self.assertEqual(
            frozen["summary"]["verdict"], gate.EXPECTED_FROZEN_119_VERDICT
        )


if __name__ == "__main__":
    unittest.main()
