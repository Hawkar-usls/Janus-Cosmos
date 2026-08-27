#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import minimize
from celerite2 import GaussianProcess, terms

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
TSV = OUT / "data" / "vega_tres_rv.tsv"
REPORT = OUT / "vega_tres_rv_quasiperiodic_activity_audit.json"

ROTATION_D = 0.676
CANDIDATE_D = 2.43
LONG_CONTROL_D = 196.4
EXPECTED_ROWS = 1524

# Scalable Hurt-like approximation:
# Hurt et al. 2021 used a quasi-periodic Gaussian process and inferred a
# characteristic activity evolution timescale near 180 d. Here we keep the
# published rotation period fixed and use celerite2 RotationTerm, a two-SHO
# quasi-periodic kernel. This is intentionally not represented as an exact
# reimplementation of the paper's kernel/hyperparameterization.


def parse_rows(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    if len(rows) != EXPECTED_ROWS:
        raise RuntimeError(f"expected {EXPECTED_ROWS} RV rows, got {len(rows)}")
    a = np.asarray(rows, dtype=float)
    t = a[:, 0] - a[0, 0]
    return t, a[:, 1], a[:, 2]


def design_matrix(t: np.ndarray, periods: list[float]) -> np.ndarray:
    cols = [np.ones_like(t)]
    for p in periods:
        phase = 2.0 * np.pi * t / p
        cols.extend([np.sin(phase), np.cos(phase)])
    return np.column_stack(cols)


def weighted_linear_start(
    t: np.ndarray,
    y: np.ndarray,
    yerr: np.ndarray,
    periods: list[float],
) -> np.ndarray:
    x = design_matrix(t, periods)
    w = 1.0 / np.maximum(yerr, 1e-9)
    beta, *_ = np.linalg.lstsq(x * w[:, None], y * w, rcond=None)
    return beta


def unpack(
    theta: np.ndarray,
    n_beta: int,
) -> tuple[float, float, float, float, float, np.ndarray]:
    log_sigma, log_q0, log_dq, logit_f, log_jitter = theta[:5]
    sigma = float(np.exp(log_sigma))
    q0 = float(np.exp(log_q0))
    dq = float(np.exp(log_dq))
    f = float(1.0 / (1.0 + np.exp(-logit_f)))
    jitter = float(np.exp(log_jitter))
    beta = np.asarray(theta[5 : 5 + n_beta], dtype=float)
    return sigma, q0, dq, f, jitter, beta


def fit_model(
    t: np.ndarray,
    y: np.ndarray,
    yerr: np.ndarray,
    periods: list[float],
) -> dict:
    x = design_matrix(t, periods)
    beta0 = weighted_linear_start(t, y, yerr, periods)
    n_beta = x.shape[1]
    ystd = max(float(np.std(y)), 1.0)
    emed = max(float(np.median(yerr)), 0.1)

    beta_bounds = [(-1000.0, 1000.0)] + [(-100.0, 100.0)] * (n_beta - 1)
    bounds = [
        (math.log(0.1), math.log(200.0)),
        (math.log(0.05), math.log(2.0e4)),
        (math.log(0.05), math.log(2.0e4)),
        (-6.0, 6.0),
        (math.log(0.01), math.log(200.0)),
        *beta_bounds,
    ]

    def objective(theta: np.ndarray) -> float:
        sigma, q0, dq, f, jitter, beta = unpack(theta, n_beta)
        try:
            kernel = terms.RotationTerm(
                sigma=sigma,
                period=ROTATION_D,
                Q0=q0,
                dQ=dq,
                f=f,
            )
            gp = GaussianProcess(kernel, mean=0.0)
            diag = np.sqrt(yerr * yerr + jitter * jitter)
            gp.compute(t, yerr=diag)
            resid = y - x @ beta
            ll = float(gp.log_likelihood(resid))
            if not math.isfinite(ll):
                return 1e100
            return -ll
        except Exception:
            return 1e100

    starts: list[np.ndarray] = []
    for tau_days in (30.0, 180.0, 1000.0):
        q_guess = max(math.pi * tau_days / ROTATION_D, 0.1)
        theta0 = np.concatenate(
            [
                np.array(
                    [
                        math.log(min(max(ystd * 0.7, 1.0), 100.0)),
                        math.log(min(q_guess, 1.0e4)),
                        math.log(max(q_guess * 0.3, 0.1)),
                        0.0,
                        math.log(min(max(emed, 0.5), 100.0)),
                    ],
                    dtype=float,
                ),
                beta0.astype(float),
            ]
        )
        starts.append(theta0)

    best = None
    for theta0 in starts:
        result = minimize(
            objective,
            theta0,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 500, "ftol": 1e-10, "gtol": 1e-7, "maxls": 40},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result

    assert best is not None
    sigma, q0, dq, f, jitter, beta = unpack(np.asarray(best.x, dtype=float), n_beta)
    loglike = -float(best.fun)
    k = 5 + n_beta
    n = len(y)
    bic = k * math.log(n) - 2.0 * loglike
    q_eff = q0 + dq + 0.5
    coherence_days_proxy = q_eff * ROTATION_D / math.pi

    amps: dict[str, float] = {}
    pos = 1
    for p in periods:
        amps[str(p)] = float(math.hypot(beta[pos], beta[pos + 1]))
        pos += 2

    return {
        "periods_days": periods,
        "n": n,
        "k": k,
        "log_likelihood": loglike,
        "bic": bic,
        "optimizer_success": bool(best.success),
        "optimizer_status": int(best.status),
        "optimizer_message": str(best.message),
        "optimizer_nit": int(getattr(best, "nit", -1)),
        "activity_hyperparameters": {
            "rotation_period_days_fixed": ROTATION_D,
            "sigma_m_s": sigma,
            "Q0": q0,
            "dQ": dq,
            "mix_f": f,
            "extra_jitter_m_s": jitter,
            "coherence_timescale_proxy_days": coherence_days_proxy,
        },
        "mean_coefficients": beta.tolist(),
        "global_period_amplitudes_m_s": amps,
    }


def classify(activity: dict, candidate: dict, control: dict) -> dict:
    dbic_cand = activity["bic"] - candidate["bic"]
    dbic_ctrl = activity["bic"] - control["bic"]
    cand_over_ctrl = control["bic"] - candidate["bic"]

    if dbic_cand >= 10.0 and dbic_ctrl < 10.0 and cand_over_ctrl >= 10.0:
        status = "2P43_SURVIVES_HURT_LIKE_QP_ACTIVITY_WITH_SPECIFICITY"
    elif dbic_cand >= 10.0 and dbic_ctrl >= 10.0:
        status = "NON_SPECIFIC_BOTH_PERIODS_SURVIVE_HURT_LIKE_QP_ACTIVITY"
    elif dbic_cand > 0.0 and dbic_ctrl <= 0.0:
        status = "WEAK_2P43_PREFERENCE_AFTER_HURT_LIKE_QP_ACTIVITY"
    elif dbic_cand <= 0.0:
        status = "HURT_LIKE_QP_ACTIVITY_ABSORBS_2P43"
    else:
        status = "MIXED_HURT_LIKE_QP_ACTIVITY_SPECIFICITY_RESULT"

    return {
        "status": status,
        "delta_bic_candidate_vs_activity": dbic_cand,
        "delta_bic_control_vs_activity": dbic_ctrl,
        "delta_bic_candidate_over_control": cand_over_ctrl,
        "threshold_convention": "Delta BIC >= 10 is treated as strong model-selection preference for this diagnostic only.",
    }


def main() -> int:
    if not TSV.exists():
        raise SystemExit("missing output/data/vega_tres_rv.tsv; run fetch_tres_rv.py first")

    t, y, yerr = parse_rows(TSV)
    models = {
        "quasiperiodic_activity_only": fit_model(t, y, yerr, []),
        "quasiperiodic_activity_plus_2p43": fit_model(t, y, yerr, [CANDIDATE_D]),
        "quasiperiodic_activity_plus_196p4": fit_model(t, y, yerr, [LONG_CONTROL_D]),
        "quasiperiodic_activity_plus_both": fit_model(t, y, yerr, [CANDIDATE_D, LONG_CONTROL_D]),
    }
    specificity = classify(
        models["quasiperiodic_activity_only"],
        models["quasiperiodic_activity_plus_2p43"],
        models["quasiperiodic_activity_plus_196p4"],
    )

    import celerite2

    report = {
        "schema": "janus.cosmos.vega.tres_rv_hurt_like_quasiperiodic_activity_audit.v1.3",
        "source": "VizieR J/AJ/161/157/table2 / Hurt et al. 2021",
        "row_count": int(len(y)),
        "fixed_periods_days": {
            "stellar_rotation": ROTATION_D,
            "published_candidate": CANDIDATE_D,
            "published_long_signal_negative_control_like": LONG_CONTROL_D,
        },
        "model": {
            "class": "celerite2_RotationTerm_two_SHO_quasiperiodic_GP",
            "relationship_to_hurt_2021": "SCALABLE_HURT_LIKE_APPROXIMATION_NOT_EXACT_REIMPLEMENTATION",
            "reason": "Hurt et al. used a quasi-periodic Gaussian process and inferred a characteristic activity evolution timescale near 180 days. This audit fixes the published rotation period and uses a scalable two-SHO quasi-periodic kernel so all 1524 RV points can be fit directly.",
            "important_difference": "The kernel family and hyperparameterization are not identical to Hurt et al. 2021, so agreement or disagreement is a diagnostic reproduction, not an independent confirmation.",
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "celerite2": celerite2.__version__,
        },
        "models": models,
        "specificity": specificity,
        "claim_ceiling": "ACTIVITY_AWARE_DIAGNOSTIC_ONLY_NO_PLANET_CONFIRMATION",
        "claim_firewall": [
            "A positive 2.43-day model preference does not establish planetary origin.",
            "The 196.4-day signal is retained as a negative-control-like comparator because Hurt et al. judged it not good evidence for a planet.",
            "Independent/new RV epochs or an independent instrument are still required for confirmation.",
        ],
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print("VEGA HURT-LIKE QUASIPERIODIC ACTIVITY AUDIT PASS")
    print("specificity_status =", specificity["status"])
    print("delta_BIC_2p43 =", specificity["delta_bic_candidate_vs_activity"])
    print("delta_BIC_196p4 =", specificity["delta_bic_control_vs_activity"])
    print("delta_BIC_2p43_over_196p4 =", specificity["delta_bic_candidate_over_control"])
    print(
        "amp_2p43_m_s =",
        models["quasiperiodic_activity_plus_2p43"]["global_period_amplitudes_m_s"][str(CANDIDATE_D)],
    )
    print(
        "coherence_proxy_days =",
        models["quasiperiodic_activity_only"]["activity_hyperparameters"]["coherence_timescale_proxy_days"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
