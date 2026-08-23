#!/usr/bin/env python3
"""JANUS COSMOS — TESLA SWEEP metadata triage.

This first-stage tool does not claim source identification from frequency alone.
It annotates human-active bands, searches for same-target/same-drift pairs,
and tests near-integer MHz spacing as an RFI fingerprint candidate.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

HUMAN_ACTIVE_BANDS = [
    (1164.0, 1215.0, "RNSS/GNSS / aeronautical radionavigation"),
    (1300.0, 1350.0, "aeronautical radionavigation / radiolocation"),
    (1350.0, 1390.0, "fixed/mobile/radiolocation; RFI-heavy region"),
    (1435.0, 1525.0, "aeronautical telemetry/telecommand"),
    (1710.0, 1755.0, "AWS/mobile + federal fixed/mobile use"),
]


def load_rows(path: Path):
    rows = []
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            row["freq_mhz"] = float(row["freq_mhz"])
            row["mjd"] = float(row["mjd"])
            row["drift_hz_s"] = float(row["drift_hz_s"])
            row["snr"] = float(row["snr"])
            rows.append(row)
    return rows


def tags_for_frequency(freq_mhz: float):
    return [label for lo, hi, label in HUMAN_ACTIVE_BANDS if lo <= freq_mhz <= hi]


def analyze(rows, drift_tolerance_hz_s=0.02, integer_spacing_tolerance_hz=2000.0):
    annotated = []
    for row in rows:
        item = dict(row)
        item["human_active_band_tags"] = tags_for_frequency(row["freq_mhz"])
        annotated.append(item)

    pairs = []
    for i, a in enumerate(rows):
        for b in rows[i + 1 :]:
            if a["target"] != b["target"]:
                continue
            drift_delta = abs(a["drift_hz_s"] - b["drift_hz_s"])
            spacing_mhz = abs(a["freq_mhz"] - b["freq_mhz"])
            nearest_integer_mhz = round(spacing_mhz)
            residual_hz = abs(spacing_mhz - nearest_integer_mhz) * 1e6
            dt_s = abs(a["mjd"] - b["mjd"]) * 86400.0
            pair = {
                "a": a["id"],
                "b": b["id"],
                "target": a["target"],
                "drift_delta_hz_s": drift_delta,
                "frequency_spacing_mhz": spacing_mhz,
                "nearest_integer_spacing_mhz": nearest_integer_mhz,
                "integer_spacing_residual_hz": residual_hz,
                "time_separation_s": dt_s,
                "same_drift_gate": drift_delta <= drift_tolerance_hz_s,
                "integer_spacing_gate": nearest_integer_mhz > 0 and residual_hz <= integer_spacing_tolerance_hz,
            }
            pair["fingerprint_gate"] = pair["same_drift_gate"] and pair["integer_spacing_gate"]
            pairs.append(pair)

    return {
        "schema": "JANUS_TESLA_SWEEP_METADATA_TRIAGE",
        "principle": "BAND_MEMBERSHIP_IS_PRIOR_NOT_PROOF",
        "rows": annotated,
        "summary": {
            "row_count": len(rows),
            "rows_in_known_human_active_bands": sum(bool(x["human_active_band_tags"]) for x in annotated),
            "same_target_pairs": len(pairs),
            "fingerprint_pairs": [p for p in pairs if p["fingerprint_gate"]],
        },
        "same_target_pairs": pairs,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", type=Path)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    report = analyze(load_rows(args.csv))
    text = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
