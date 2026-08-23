#!/usr/bin/env python3
"""TOPA bathymetry candidate engine v1.1.

Repair frozen after TOPA-TRAINING-GROUND-003A synthetic adversarial QA and before
numeric FKt230303 field-grid scoring. It reuses v1.0 feature extraction but adds a
lineation-confounder penalty. This is still a morphology candidate generator only.
"""
from __future__ import annotations

import numpy as np
import topa_bathymetry_benchmark as base

LINEATION_PENALTY = 2.0
LINEATION_FLAG_THRESHOLD = 0.55


def score_rows_v1_1(rows):
    keys = [
        "local_relief_p95_p05",
        "slope_p95",
        "abs_curvature_p95",
        "rugosity_residual_std",
        "pinnacle_density",
        "multiscale_positive_persistence",
    ]
    weights = np.array([0.15, 0.10, 0.15, 0.15, 0.25, 0.20], dtype=float)
    mat = np.array([[float(r["metrics"][k]) for k in keys] for r in rows], dtype=float)
    z = np.column_stack([base.robust_z(mat[:, i]) for i in range(mat.shape[1])])
    score_v1_0 = z @ weights

    anis = np.array([float(r["metrics"]["orientation_anisotropy"]) for r in rows], dtype=float)
    anis_z = base.robust_z(anis)
    score = score_v1_0 - LINEATION_PENALTY * np.clip(anis_z, 0.0, None)

    for i, r in enumerate(rows):
        r["morphology_score_v1_0"] = float(score_v1_0[i])
        r["orientation_anisotropy_robust_z"] = float(anis_z[i])
        r["lineation_penalty_applied"] = float(LINEATION_PENALTY * max(float(anis_z[i]), 0.0))
        r["morphology_score"] = float(score[i])
        r["lineation_confounder_flag"] = bool(anis[i] >= LINEATION_FLAG_THRESHOLD)
        r["classification_ceiling"] = "MORPHOLOGICAL_CANDIDATE_ONLY"

    rows.sort(key=lambda r: float(r["morphology_score"]), reverse=True)
    for rank, r in enumerate(rows, 1):
        r["rank"] = rank


# Monkey-patch only the scoring stage. Feature extraction and I/O remain v1.0.
base.score_rows = score_rows_v1_1

# Ensure receipts can identify the frozen scoring revision.
_original_run_benchmark = base.run_benchmark

def run_benchmark_v1_1(a, tile_size, stride, top_k):
    out = _original_run_benchmark(a, tile_size, stride, top_k)
    out["metric_rule"] = "FROZEN_GENERIC_MORPHOLOGY_V1_1_AFTER_SYNTHETIC_QA"
    out["lineation_penalty"] = LINEATION_PENALTY
    out["repair_parent"] = "TOPA-TRAINING-GROUND-003A"
    return out

base.run_benchmark = run_benchmark_v1_1

if __name__ == "__main__":
    base.main()
