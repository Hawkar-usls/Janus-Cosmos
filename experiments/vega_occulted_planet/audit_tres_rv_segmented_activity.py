#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
TSV = OUT / "data" / "vega_tres_rv.tsv"
REPORT = OUT / "vega_tres_rv_segmented_activity_audit.json"

ROTATION_D = 0.676
CANDIDATE_D = 2.43
LONG_CONTROL_D = 196.4
CHUNK_DAYS = [90.0, 180.0, 360.0]
MIN_POINTS_PER_ACTIVITY_SEGMENT = 8


def parse_rows(path: Path) -> list[tuple[float, float, float]]:
    rows: list[tuple[float, float, float]] = []
    header_seen = False
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip("\ufeff\r\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if not header_seen:
            if [p.strip() for p in parts[:4]] == ["Seq", "BJD", "RVel", "e_RVel"]:
                header_seen = True
            continue
        if len(parts) < 4:
            continue
        try:
            bjd = float(parts[1].strip())
            rv = float(parts[2].strip())
            erv = float(parts[3].strip())
        except ValueError:
            continue
        if erv > 0 and all(math.isfinite(x) for x in (bjd, rv, erv)):
            rows.append((bjd, rv, erv))
    rows.sort(key=lambda r: r[0])
    return rows


def solve_linear(a: list[list[float]], b: list[float]) -> list[float]:
    n = len(b)
    m = [row[:] + [rhs] for row, rhs in zip(a, b)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[pivot][col]) < 1e-12:
            raise RuntimeError(f"singular normal matrix at column {col}")
        if pivot != col:
            m[col], m[pivot] = m[pivot], m[col]
        div = m[col][col]
        for j in range(col, n + 1):
            m[col][j] /= div
        for r in range(n):
            if r == col:
                continue
            f = m[r][col]
            if abs(f) < 1e-30:
                continue
            for j in range(col, n + 1):
                m[r][j] -= f * m[col][j]
    return [m[i][n] for i in range(n)]


def base_bin(t: float, t0: float, chunk_days: float) -> int:
    return int(math.floor((t - t0) / chunk_days))


def build_segments(
    rows: list[tuple[float, float, float]],
    chunk_days: float,
    min_points: int = MIN_POINTS_PER_ACTIVITY_SEGMENT,
) -> tuple[dict[int, int], list[dict]]:
    """Merge sparse adjacent nominal time bins before fitting local activity.

    A free offset + sin/cos activity model needs enough observations in each
    segment. The raw TRES cadence contains a few 90/180-day bins with only
    1-5 points, so treating every occupied bin as three free parameters makes
    the normal matrix singular. This deterministic pre-registered merge keeps
    the nominal activity timescale while forbidding underdetermined segments.
    """
    t0 = rows[0][0]
    counts: dict[int, int] = {}
    for t, _, _ in rows:
        b = base_bin(t, t0, chunk_days)
        counts[b] = counts.get(b, 0) + 1

    groups: list[list[int]] = []
    current: list[int] = []
    current_n = 0
    for b in sorted(counts):
        current.append(b)
        current_n += counts[b]
        if current_n >= min_points:
            groups.append(current)
            current = []
            current_n = 0
    if current:
        if groups:
            groups[-1].extend(current)
        else:
            groups.append(current)

    mapping: dict[int, int] = {}
    diagnostics: list[dict] = []
    for slot, bins in enumerate(groups):
        n = sum(counts[b] for b in bins)
        if n < min_points:
            raise RuntimeError(f"activity segment {slot} has only {n} points")
        for b in bins:
            mapping[b] = slot
        diagnostics.append(
            {
                "slot": slot,
                "base_bins": bins,
                "points": n,
                "approx_start_days_from_first_epoch": bins[0] * chunk_days,
                "approx_end_days_from_first_epoch": (bins[-1] + 1) * chunk_days,
            }
        )
    return mapping, diagnostics


def design_row(
    t: float,
    t0: float,
    chunk_days: float,
    bin_to_segment: dict[int, int],
    n_segments: int,
    global_periods: list[float],
) -> list[float]:
    n_local = 3 * n_segments
    x = [0.0] * (n_local + 2 * len(global_periods))
    b = base_bin(t, t0, chunk_days)
    slot = bin_to_segment[b]
    local = 3 * slot
    dt = t - t0
    x[local] = 1.0
    rot_phase = 2.0 * math.pi * dt / ROTATION_D
    x[local + 1] = math.sin(rot_phase)
    x[local + 2] = math.cos(rot_phase)
    pos = n_local
    for p in global_periods:
        phase = 2.0 * math.pi * dt / p
        x[pos] = math.sin(phase)
        x[pos + 1] = math.cos(phase)
        pos += 2
    return x


def fit(
    rows: list[tuple[float, float, float]],
    chunk_days: float,
    global_periods: list[float],
) -> dict:
    t0 = rows[0][0]
    bin_to_segment, segment_diag = build_segments(rows, chunk_days)
    n_segments = len(segment_diag)
    k = 3 * n_segments + 2 * len(global_periods)
    ata = [[0.0 for _ in range(k)] for _ in range(k)]
    aty = [0.0 for _ in range(k)]
    design: list[tuple[list[float], float, float]] = []

    for t, y, err in rows:
        x = design_row(t, t0, chunk_days, bin_to_segment, n_segments, global_periods)
        design.append((x, y, err))
        w = 1.0 / (err * err)
        for i in range(k):
            xi = x[i]
            if xi == 0.0:
                continue
            aty[i] += w * xi * y
            for j in range(i, k):
                xj = x[j]
                if xj == 0.0:
                    continue
                ata[i][j] += w * xi * xj

    for i in range(k):
        for j in range(i):
            ata[i][j] = ata[j][i]

    beta = solve_linear(ata, aty)
    chi2 = 0.0
    rss = 0.0
    for x, y, err in design:
        yhat = sum(c * q for c, q in zip(beta, x))
        resid = y - yhat
        chi2 += (resid / err) ** 2
        rss += resid * resid

    n = len(rows)
    bic = chi2 + k * math.log(n)
    n_local = 3 * n_segments
    amps: dict[str, float] = {}
    pos = n_local
    for p in global_periods:
        amps[str(p)] = math.hypot(beta[pos], beta[pos + 1])
        pos += 2

    return {
        "chunk_days": chunk_days,
        "activity_segments": n_segments,
        "minimum_points_per_segment": MIN_POINTS_PER_ACTIVITY_SEGMENT,
        "segment_diagnostics": segment_diag,
        "n": n,
        "k": k,
        "chi2": chi2,
        "bic": bic,
        "rms_m_s": math.sqrt(rss / n),
        "global_period_amplitudes_m_s": amps,
    }


def classify(results: list[dict]) -> dict:
    cand_deltas = [r["delta_bic_candidate"] for r in results]
    ctrl_deltas = [r["delta_bic_long_control"] for r in results]
    cand_strong = sum(v >= 10.0 for v in cand_deltas)
    ctrl_strong = sum(v >= 10.0 for v in ctrl_deltas)
    candidate_beats_control = sum((c - l) >= 10.0 for c, l in zip(cand_deltas, ctrl_deltas))
    candidate_higher_all = all(c > l for c, l in zip(cand_deltas, ctrl_deltas))

    if cand_strong == 0:
        status = "ACTIVITY_FLEXIBILITY_ABSORBS_2P43_CANDIDATE"
    elif cand_strong == len(results) and ctrl_strong == 0:
        status = "2P43_SURVIVES_WHILE_196P4_CONTROL_COLLAPSES"
    elif cand_strong >= 2 and candidate_beats_control >= 2:
        status = "2P43_PREFERRED_OVER_196P4_CONTROL_BUT_NOT_CONFIRMED"
    elif cand_strong >= 2 and ctrl_strong >= 2:
        status = "NON_SPECIFIC_BOTH_SIGNALS_SURVIVE_ACTIVITY_MODEL"
    else:
        status = "MIXED_ACTIVITY_SPECIFICITY_RESULT"

    return {
        "status": status,
        "candidate_strong_chunk_scales": cand_strong,
        "long_control_strong_chunk_scales": ctrl_strong,
        "candidate_beats_control_by_delta_bic_10_chunk_scales": candidate_beats_control,
        "candidate_delta_bic_higher_than_control_at_all_chunk_scales": candidate_higher_all,
        "rule": "This is a specificity diagnostic, not a planet confirmation. The 196.4-day signal is used as a preregistered negative-control-like comparator because Hurt et al. 2021 explicitly judged it not good evidence for a planet after more complete analysis."
    }


def main() -> int:
    if not TSV.exists():
        raise SystemExit("missing output/data/vega_tres_rv.tsv; run fetch_tres_rv.py first")
    rows = parse_rows(TSV)
    if len(rows) != 1524:
        raise SystemExit(f"expected 1524 RV rows, got {len(rows)}")

    scales: list[dict] = []
    for chunk_days in CHUNK_DAYS:
        base = fit(rows, chunk_days, [])
        candidate = fit(rows, chunk_days, [CANDIDATE_D])
        control = fit(rows, chunk_days, [LONG_CONTROL_D])
        both = fit(rows, chunk_days, [CANDIDATE_D, LONG_CONTROL_D])
        scales.append(
            {
                "chunk_days": chunk_days,
                "activity_only": base,
                "activity_plus_2p43": candidate,
                "activity_plus_196p4": control,
                "activity_plus_both": both,
                "delta_bic_candidate": base["bic"] - candidate["bic"],
                "delta_bic_long_control": base["bic"] - control["bic"],
                "candidate_amplitude_m_s": candidate["global_period_amplitudes_m_s"][str(CANDIDATE_D)],
                "long_control_amplitude_m_s": control["global_period_amplitudes_m_s"][str(LONG_CONTROL_D)],
            }
        )

    specificity = classify(scales)
    report = {
        "schema": "janus.cosmos.vega.tres_rv_segmented_activity_audit.v1.2.1",
        "source": "VizieR J/AJ/161/157/table2 / Hurt et al. 2021",
        "row_count": len(rows),
        "stellar_rotation_period_days": ROTATION_D,
        "candidate_period_days": CANDIDATE_D,
        "negative_control_period_days": LONG_CONTROL_D,
        "activity_model": {
            "type": "piecewise_local_offset_plus_rotation_sinusoid_with_sparse_bin_merge",
            "nominal_chunk_days": CHUNK_DAYS,
            "minimum_points_per_activity_segment": MIN_POINTS_PER_ACTIVITY_SEGMENT,
            "description": "Nominal fixed time bins are deterministically merged forward when sparse; each resulting segment has its own offset and independent sin/cos coefficients at Vega's 0.676-day stellar rotation period. Candidate/control periods remain global across all segments.",
            "not_equivalent_to": "The quasi-periodic Gaussian-process activity model used in Hurt et al. 2021."
        },
        "scales": scales,
        "specificity": specificity,
        "claim_ceiling": "ACTIVITY_SPECIFICITY_DIAGNOSTIC_ONLY_NO_PLANET_CONFIRMATION"
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("VEGA SEGMENTED ACTIVITY AUDIT PASS")
    print("specificity_status =", specificity["status"])
    for r in scales:
        print(
            "chunk_days =", r["chunk_days"],
            "segments =", r["activity_only"]["activity_segments"],
            "delta_BIC_2p43 =", r["delta_bic_candidate"],
            "delta_BIC_196p4 =", r["delta_bic_long_control"],
            "amp_2p43 =", r["candidate_amplitude_m_s"],
            "amp_196p4 =", r["long_control_amplitude_m_s"],
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
