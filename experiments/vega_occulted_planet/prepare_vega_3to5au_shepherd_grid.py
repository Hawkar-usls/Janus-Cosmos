#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
REPORT = OUT / "vega_3to5au_shepherd_analytic_grid.json"
MANIFEST = OUT / "vega_3to5au_injection_recovery_manifest.json"

MSTAR_SOLAR = 2.15
M_EARTH_PER_JUPITER = 317.83
M_SUN_PER_EARTH = 332946.0
TARGET_INNER_EDGE_AU = (3.0, 5.0)

PLANET_A_AU = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
PLANET_M_EARTH = [10.0, 17.0, 30.0, 50.0, 95.0]
INCLINATION_DEG = [3.0, 6.5, 10.0, 30.0, 60.0, 90.0]
ECCENTRICITY = [0.0, 0.05, 0.10, 0.20]


def canonical_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def orbital_period_years(a_au: float) -> float:
    return math.sqrt(a_au**3 / MSTAR_SOLAR)


def rv_semiamplitude_m_s(m_earth: float, a_au: float, inc_deg: float, ecc: float) -> float:
    p_years = orbital_period_years(a_au)
    m_jup = m_earth / M_EARTH_PER_JUPITER
    return (
        28.4329
        * m_jup
        * math.sin(math.radians(inc_deg))
        * (p_years ** (-1.0 / 3.0))
        * (MSTAR_SOLAR ** (-2.0 / 3.0))
        / math.sqrt(max(1.0 - ecc * ecc, 1e-12))
    )


def chaotic_zone_half_width_au(m_earth: float, a_au: float) -> float:
    # Low-eccentricity Wisdom-style screening estimate:
    # Delta a ~= 1.5 * mu^(2/7) * a_p.
    # This is a pre-injection proxy, not an N-body result.
    mu = (m_earth / M_SUN_PER_EARTH) / MSTAR_SOLAR
    return 1.5 * (mu ** (2.0 / 7.0)) * a_au


def build() -> tuple[dict, dict]:
    rows: list[dict] = []
    for a in PLANET_A_AU:
        period = orbital_period_years(a)
        for m in PLANET_M_EARTH:
            dz = chaotic_zone_half_width_au(m, a)
            outer_edge = a + dz
            edge_match = TARGET_INNER_EDGE_AU[0] <= outer_edge <= TARGET_INNER_EDGE_AU[1]
            for inc in INCLINATION_DEG:
                for ecc in ECCENTRICITY:
                    k = rv_semiamplitude_m_s(m, a, inc, ecc)
                    rows.append(
                        {
                            "a_au": a,
                            "mass_earth": m,
                            "inclination_deg": inc,
                            "eccentricity": ecc,
                            "period_years": period,
                            "chaotic_zone_half_width_au": dz,
                            "analytic_outer_edge_au": outer_edge,
                            "warm_disk_inner_edge_match": edge_match,
                            "rv_semiamplitude_m_s": k,
                            "near_pole_on_geometry": inc <= 10.0,
                            "neptune_class_mass": 10.0 <= m <= 30.0,
                        }
                    )

    feasible = [
        r
        for r in rows
        if r["warm_disk_inner_edge_match"]
        and r["near_pole_on_geometry"]
        and r["neptune_class_mass"]
    ]
    report = {
        "schema": "janus.cosmos.vega.shepherd_analytic_screen.v1.3",
        "target": "Vega warm-disk inner edge",
        "source_constraints": {
            "warm_disk_inner_edge_au": list(TARGET_INNER_EDGE_AU),
            "miri_interpretation": "Su et al. 2024 explicitly allows a possible modest/Neptune-size shepherd near the 3-5 au warm-disk inner edge.",
            "stellar_mass_solar": MSTAR_SOLAR,
            "stellar_mass_source": "Monnier et al. 2012 rotating-star estimate, approximately 2.15 Msun.",
        },
        "grid": {
            "semi_major_axis_au": PLANET_A_AU,
            "mass_earth": PLANET_M_EARTH,
            "inclination_deg": INCLINATION_DEG,
            "eccentricity": ECCENTRICITY,
            "row_count": len(rows),
        },
        "analytic_model": {
            "disk_proxy": "Wisdom-style low-eccentricity chaotic-zone outer boundary a_p + 1.5 mu^(2/7) a_p",
            "rv_proxy": "Keplerian semi-amplitude from mass, semimajor axis, eccentricity, and inclination",
            "not_included": [
                "No N-body dust integration.",
                "No MIRI image-domain forward model.",
                "No PSF/convolution/noise injection.",
                "No recovered-likelihood comparison against the real JWST image.",
            ],
        },
        "summary": {
            "analytic_feasible_neptune_like_pole_on_rows": len(feasible),
            "minimum_rv_k_m_s_among_feasible": min(
                (r["rv_semiamplitude_m_s"] for r in feasible), default=None
            ),
            "maximum_rv_k_m_s_among_feasible": max(
                (r["rv_semiamplitude_m_s"] for r in feasible), default=None
            ),
            "status": (
                "ANALYTIC_FEASIBILITY_EXISTS_REAL_IMAGE_RECOVERY_PENDING"
                if feasible
                else "NO_ANALYTIC_NEPTUNE_LIKE_EDGE_MATCH_IN_FROZEN_GRID"
            ),
        },
        "rows": rows,
        "claim_ceiling": "ANALYTIC_SCREEN_ONLY_NOT_A_DISK_PLANET_DETECTION",
    }
    report["freeze_sha256"] = canonical_hash(report)

    manifest = {
        "schema": "janus.cosmos.vega.disk_injection_recovery_manifest.v1.3",
        "parent_analytic_hash": report["freeze_sha256"],
        "target": "H1B_3_TO_5_AU_NEPTUNE_SHEPHERD",
        "status": "READY_FOR_REAL_MIRI_PRODUCT",
        "required_input": [
            "Calibrated or science-ready MIRI image/product resolving the warm-disk inner edge.",
            "PSF or empirically justified convolution kernel.",
            "Pixel scale/WCS and uncertainty/noise model.",
            "Mask/support map for coronagraphic or invalid regions if applicable.",
        ],
        "frozen_test": {
            "planet_grid": {
                "semi_major_axis_au": PLANET_A_AU,
                "mass_earth": PLANET_M_EARTH,
                "inclination_deg": INCLINATION_DEG,
                "eccentricity": ECCENTRICITY,
            },
            "controls": [
                "NO_PLANET_PR_DRAG_ONLY",
                "MASS_SHUFFLED_PLANET_GRID",
                "AZIMUTH_PHASE_PLACEBO",
            ],
            "primary_recovery_metric": "Improvement in held-out image-domain residual likelihood after identical PSF/noise processing.",
            "success_rule": "A predeclared low-mass family must improve held-out recovery over no-planet and placebo controls without violating independent RV/direct-imaging constraints.",
            "failure_rule": "No-planet or placebo controls match as well or better, or the required mass/orbit conflicts with independent constraints.",
        },
        "claim_firewall": "This manifest pre-registers a future real-data injection/recovery. The analytic grid above is not itself a recovered planet.",
    }
    manifest["freeze_sha256"] = canonical_hash(manifest)
    return report, manifest


def main() -> int:
    OUT.mkdir(exist_ok=True)
    report, manifest = build()
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print("VEGA 3-5 AU SHEPHERD ANALYTIC GRID PASS")
    print("grid_rows =", report["grid"]["row_count"])
    print(
        "feasible_neptune_like_pole_on_rows =",
        report["summary"]["analytic_feasible_neptune_like_pole_on_rows"],
    )
    print("status =", report["summary"]["status"])
    print("analytic_freeze =", report["freeze_sha256"])
    print("injection_manifest_freeze =", manifest["freeze_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
