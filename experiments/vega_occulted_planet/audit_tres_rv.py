#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
DATA = OUT / "data"
TSV = DATA / "vega_tres_rv.tsv"
REPORT = OUT / "vega_tres_rv_audit.json"

PERIOD_ROT_D = 0.676
PERIOD_CANDIDATE_D = 2.43
PERIOD_LONG_D = 196.4


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
        if erv > 0 and math.isfinite(bjd) and math.isfinite(rv) and math.isfinite(erv):
            rows.append((bjd, rv, erv))
    return rows


def solve_linear(a: list[list[float]], b: list[float]) -> list[float]:
    n = len(b)
    m = [row[:] + [rhs] for row, rhs in zip(a, b)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[pivot][col]) < 1e-14:
            raise RuntimeError("singular normal matrix")
        if pivot != col:
            m[col], m[pivot] = m[pivot], m[col]
        div = m[col][col]
        for j in range(col, n + 1):
            m[col][j] /= div
        for r in range(n):
            if r == col:
                continue
            f = m[r][col]
            if f == 0:
                continue
            for j in range(col, n + 1):
                m[r][j] -= f * m[col][j]
    return [m[i][n] for i in range(n)]


def design_row(t: float, t0: float, periods: list[float]) -> list[float]:
    x = [1.0]
    dt = t - t0
    for p in periods:
        phase = 2.0 * math.pi * dt / p
        x.extend([math.sin(phase), math.cos(phase)])
    return x


def fit(rows: list[tuple[float, float, float]], periods: list[float]) -> dict:
    t0 = min(r[0] for r in rows)
    k = 1 + 2 * len(periods)
    ata = [[0.0 for _ in range(k)] for _ in range(k)]
    aty = [0.0 for _ in range(k)]
    for t, y, err in rows:
        x = design_row(t, t0, periods)
        w = 1.0 / (err * err)
        for i in range(k):
            aty[i] += w * x[i] * y
            for j in range(k):
                ata[i][j] += w * x[i] * x[j]
    beta = solve_linear(ata, aty)
    chi2 = 0.0
    wrss = 0.0
    for t, y, err in rows:
        x = design_row(t, t0, periods)
        yhat = sum(c * q for c, q in zip(beta, x))
        resid = y - yhat
        chi2 += (resid / err) ** 2
        wrss += resid * resid
    n = len(rows)
    bic = chi2 + k * math.log(n)
    amplitudes = {}
    idx = 1
    for p in periods:
        s = beta[idx]
        c = beta[idx + 1]
        amplitudes[str(p)] = math.hypot(s, c)
        idx += 2
    return {
        "periods_days": periods,
        "n": n,
        "k": k,
        "chi2": chi2,
        "bic": bic,
        "rms_m_s": math.sqrt(wrss / n),
        "amplitudes_m_s": amplitudes,
        "coefficients": beta,
    }


def subset(rows: list[tuple[float, float, float]], lo: float, hi: float) -> list[tuple[float, float, float]]:
    return [r for r in rows if lo <= r[0] <= hi]


def main() -> int:
    if not TSV.exists():
        raise SystemExit("missing output/data/vega_tres_rv.tsv; run fetch_tres_rv.py first")
    rows = parse_rows(TSV)
    if len(rows) != 1524:
        raise SystemExit(f"expected 1524 parsed RV rows, got {len(rows)}")
    rows.sort(key=lambda r: r[0])
    tmin = rows[0][0]
    tmax = rows[-1][0]
    tmid = 0.5 * (tmin + tmax)

    models = {
        "constant": fit(rows, []),
        "rotation_only": fit(rows, [PERIOD_ROT_D]),
        "rotation_plus_2p43": fit(rows, [PERIOD_ROT_D, PERIOD_CANDIDATE_D]),
        "rotation_plus_196p4": fit(rows, [PERIOD_ROT_D, PERIOD_LONG_D]),
        "rotation_plus_both": fit(rows, [PERIOD_ROT_D, PERIOD_CANDIDATE_D, PERIOD_LONG_D]),
    }
    early = subset(rows, tmin, tmid)
    late = subset(rows, tmid, tmax)
    early_fit = fit(early, [PERIOD_ROT_D, PERIOD_CANDIDATE_D])
    late_fit = fit(late, [PERIOD_ROT_D, PERIOD_CANDIDATE_D])

    delta_bic_candidate = models["rotation_only"]["bic"] - models["rotation_plus_2p43"]["bic"]
    delta_bic_long = models["rotation_only"]["bic"] - models["rotation_plus_196p4"]["bic"]
    report = {
        "schema": "janus.cosmos.vega.tres_rv_audit.v1.1",
        "source": "VizieR J/AJ/161/157/table2 / Hurt et al. 2021",
        "row_count": len(rows),
        "bjd_min": tmin,
        "bjd_max": tmax,
        "timespan_days": tmax - tmin,
        "fixed_periods_days": {
            "stellar_rotation": PERIOD_ROT_D,
            "published_candidate": PERIOD_CANDIDATE_D,
            "published_long_signal_not_good_evidence": PERIOD_LONG_D,
        },
        "models": models,
        "diagnostics": {
            "delta_bic_rotation_to_rotation_plus_2p43": delta_bic_candidate,
            "delta_bic_rotation_to_rotation_plus_196p4": delta_bic_long,
            "candidate_amplitude_full_m_s": models["rotation_plus_2p43"]["amplitudes_m_s"][str(PERIOD_CANDIDATE_D)],
            "candidate_amplitude_early_m_s": early_fit["amplitudes_m_s"][str(PERIOD_CANDIDATE_D)],
            "candidate_amplitude_late_m_s": late_fit["amplitudes_m_s"][str(PERIOD_CANDIDATE_D)],
            "early_rows": len(early),
            "late_rows": len(late),
        },
        "interpretation": {
            "status": "DIAGNOSTIC_ONLY_NO_CONFIRMATION",
            "warning": "This fixed-period weighted-sinusoid audit is intentionally simpler than the activity-aware analysis in Hurt et al. 2021. It cannot confirm a planet and must not be used as an independent detection claim.",
            "purpose": "Verify live ingestion, quantify how much a fixed 2.43-day component changes a simple rotation model, and expose time-split stability for TOPA prioritization."
        }
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("VEGA TRES RV AUDIT PASS")
    print("rows =", len(rows))
    print("delta_bic_2p43 =", delta_bic_candidate)
    print("candidate_amp_full_m_s =", report["diagnostics"]["candidate_amplitude_full_m_s"])
    print("candidate_amp_early_m_s =", report["diagnostics"]["candidate_amplitude_early_m_s"])
    print("candidate_amp_late_m_s =", report["diagnostics"]["candidate_amplitude_late_m_s"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
