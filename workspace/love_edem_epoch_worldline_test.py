#!/usr/bin/env python3
"""
LOVE/EDEM epoch-aware astrometric worldline gate.

Uses the frozen source-classification receipt to identify exact Gaia stellar groups,
then optionally refreshes only astrometric error/correlation fields for those source
IDs from Gaia DR3 through VizieR. It never moves the frozen target centers and never
substitutes catalog release years for observation epochs.
"""
from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict

from janus_astrometric_worldlines import (
    AstrometricState,
    GAIA_DR3_REF_EPOCH_JYEAR,
    evaluate_worldline,
)

SCHEMA = "janus.cosmos.love-edem.epoch-worldlines.v1"
DEFAULT_SAMPLES = 8192
DEFAULT_SEED = 26016
EPOCH_WINDOW = (1900.0, 2100.0)
TARGET_CENTERS = {
    "LOVE": {"ra_deg": 204.30267916666668, "dec_deg": -36.78240527777778},
    "EDEM_SEARCH_CENTER_ZP": {"ra_deg": 139.22409686590188, "dec_deg": 30.26038779947318},
}
GROUPS = {
    "LOVE_GAIA_WISE_2MASS": {"target": "LOVE", "gaia_source_id": "6163586620213012352"},
    "EDEM_GAIA_WISE_SDSS_PS1": {"target": "EDEM_SEARCH_CENTER_ZP", "gaia_source_id": "699051350998534656"},
}
VIZIER_CORR_MAP = {
    "RADEcor": "ra_dec",
    "RAPlxcor": "ra_parallax",
    "RApmRAcor": "ra_pmra",
    "RApmDEcor": "ra_pmdec",
    "DEPlxcor": "dec_parallax",
    "DEpmRAcor": "dec_pmra",
    "DEpmDEcor": "dec_pmdec",
    "PlxpmRAcor": "parallax_pmra",
    "PlxpmDEcor": "parallax_pmdec",
    "pmRApmDEcor": "pmra_pmdec",
}


def _jsonable(v: Any) -> Any:
    if v is None:
        return None
    try:
        if getattr(v, "mask", False):
            return None
    except Exception:
        pass
    if isinstance(v, (str, int, float, bool)):
        return v
    try:
        return v.item()
    except Exception:
        return str(v)


def _case_get(d: Dict[str, Any], name: str) -> Any:
    low = {k.lower(): k for k in d}
    k = low.get(name.lower())
    return d.get(k) if k is not None else None


def _finite(d: Dict[str, Any], *names: str) -> float | None:
    for n in names:
        v = _case_get(d, n)
        if v is None:
            continue
        try:
            x = float(v)
            if math.isfinite(x):
                return x
        except Exception:
            pass
    return None


def _extract_frozen_gaia(group: Dict[str, Any]) -> Dict[str, Any] | None:
    selected = group.get("catalog_queries", {}).get("GAIA_DR3", {}).get("selected")
    if not selected:
        return None
    features = dict(selected.get("features", {}))
    features.setdefault("RA_ICRS", selected.get("ra_deg"))
    features.setdefault("DE_ICRS", selected.get("dec_deg"))
    features.setdefault("Source", selected.get("source_id"))
    return features


def query_gaia_live(source_id: str) -> Dict[str, Any]:
    from astroquery.vizier import Vizier
    viz = Vizier(columns=["**"], row_limit=20)
    tables = viz.query_constraints(catalog="I/355/gaiadr3", Source=str(source_id))
    rows = []
    for table in tables:
        for row in table:
            snap = {name: _jsonable(row[name]) for name in table.colnames}
            if str(_case_get(snap, "Source")).strip() == str(source_id):
                rows.append(snap)
    if not rows:
        raise RuntimeError(f"GAIA_SOURCE_NOT_RETURNED:{source_id}")
    return rows[0]


def state_from_features(source_id: str, features: Dict[str, Any]) -> AstrometricState:
    required = {
        "ra_deg": _finite(features, "RA_ICRS"),
        "dec_deg": _finite(features, "DE_ICRS"),
        "parallax_mas": _finite(features, "Plx"),
        "pmra_masyr": _finite(features, "pmRA"),
        "pmdec_masyr": _finite(features, "pmDE"),
        "ra_error_mas": _finite(features, "e_RA_ICRS"),
        "dec_error_mas": _finite(features, "e_DE_ICRS"),
        "parallax_error_mas": _finite(features, "e_Plx"),
        "pmra_error_masyr": _finite(features, "e_pmRA"),
        "pmdec_error_masyr": _finite(features, "e_pmDE"),
    }
    missing = [k for k, v in required.items() if v is None]
    if missing:
        raise ValueError("ASTROMETRY_INCOMPLETE:" + ",".join(missing))
    corrs: Dict[str, float] = {}
    for viz_name, canonical in VIZIER_CORR_MAP.items():
        v = _finite(features, viz_name)
        if v is not None:
            corrs[canonical] = v
    return AstrometricState(
        source_id=str(source_id),
        catalog="GAIA_DR3",
        ref_epoch_jyear=GAIA_DR3_REF_EPOCH_JYEAR,
        correlations=corrs,
        **{k: float(v) for k, v in required.items() if v is not None},
    )


def _cross_catalog_epoch_audit(group: Dict[str, Any]) -> Dict[str, Any]:
    entries = {}
    for catalog, query in group.get("catalog_queries", {}).items():
        selected = query.get("selected")
        if not selected:
            continue
        features = selected.get("features", {})
        epoch_fields = {
            k: v for k, v in features.items()
            if any(token in k.lower() for token in ("epoch", "mjd", "jd")) and v is not None
        }
        entries[catalog] = {
            "source_id": selected.get("source_id"),
            "epoch_fields_in_frozen_receipt": epoch_fields,
            "has_exact_observation_epoch": bool(epoch_fields),
        }
    complete = bool(entries) and all(x["has_exact_observation_epoch"] for x in entries.values())
    return {
        "catalogs": entries,
        "status": "READY" if complete else "BLOCKED_MISSING_PER_DETECTION_EPOCH",
        "rule": "CATALOG_RELEASE_YEAR_MUST_NOT_BE_SUBSTITUTED_FOR_OBSERVATION_EPOCH",
    }


def run(classification: Dict[str, Any], *, samples: int, seed: int, live_vizier: bool) -> Dict[str, Any]:
    results: Dict[str, Any] = {}
    for i, (group_name, cfg) in enumerate(GROUPS.items()):
        group = classification.get("groups", {}).get(group_name)
        if group is None:
            results[group_name] = {"status": "I_DO_NOT_KNOW", "reason": "GROUP_MISSING"}
            continue
        frozen = _extract_frozen_gaia(group)
        if frozen is None:
            results[group_name] = {"status": "I_DO_NOT_KNOW", "reason": "GAIA_STATE_MISSING"}
            continue

        source_id = cfg["gaia_source_id"]
        astrometry_source = "FROZEN_CLASSIFICATION_RECEIPT"
        features = frozen
        live_error = None
        if live_vizier:
            try:
                features = query_gaia_live(source_id)
                astrometry_source = "LIVE_VIZIER_GAIA_DR3_EXACT_SOURCE_ID"
            except Exception as exc:
                live_error = f"{type(exc).__name__}: {exc}"

        try:
            state = state_from_features(source_id, features)
        except Exception as exc:
            results[group_name] = {
                "status": "ASTROMETRY_INCOMPLETE",
                "reason": f"{type(exc).__name__}: {exc}",
                "astrometry_source": astrometry_source,
                "live_query_error": live_error,
                "frozen_available_fields": sorted(frozen.keys()),
                "cross_catalog_epoch_audit": _cross_catalog_epoch_audit(group),
            }
            continue

        target_name = cfg["target"]
        target = TARGET_CENTERS[target_name]
        worldline = evaluate_worldline(
            state,
            target["ra_deg"], target["dec_deg"],
            samples=samples,
            seed=seed + i,
            epoch_start=EPOCH_WINDOW[0],
            epoch_end=EPOCH_WINDOW[1],
        )
        results[group_name] = {
            "status": "WORLDLINE_EVALUATED",
            "target": target_name,
            "astrometry_source": astrometry_source,
            "live_query_error": live_error,
            "worldline": worldline,
            "cross_catalog_epoch_audit": _cross_catalog_epoch_audit(group),
        }

    nearest_non_gaia = {}
    for group_name in ("LOVE_NEAREST_WISE_ONLY", "EDEM_NEAREST_WISE_ONLY"):
        g = classification.get("groups", {}).get(group_name, {})
        nearest_non_gaia[group_name] = {
            "target": g.get("target"),
            "offset_from_frozen_target_center_arcsec": g.get("offset_from_frozen_target_center_arcsec"),
            "status": "NO_GAIA_5P_WORLDLINE_AVAILABLE",
            "reason": "Nearest source group has no Gaia counterpart within the frozen classification radius; do not invent proper motion/parallax.",
        }

    evaluated = [x for x in results.values() if x.get("status") == "WORLDLINE_EVALUATED"]
    all_zero_3 = bool(evaluated) and all(
        x["worldline"]["gate_counts_after_conservative_parallax_envelope"].get("3.0") == 0
        for x in evaluated
    )
    full_cov = bool(evaluated) and all(
        x["worldline"]["covariance_status"] == "FULL_CATALOG_COVARIANCE"
        for x in evaluated
    )
    return {
        "schema": SCHEMA,
        "experiment": "LOVE_EDEM_EPOCH_AWARE_WORLDLINES",
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "formula": "RESPICIENS_ET_PROSPICIENS",
        "input_experiment_id": classification.get("experiment_id"),
        "live_vizier_requested": live_vizier,
        "samples_per_gaia_group": samples,
        "seed_base": seed,
        "epoch_window_jyear": list(EPOCH_WINDOW),
        "gaia_reference_epoch_jyear": GAIA_DR3_REF_EPOCH_JYEAR,
        "gaia_groups": results,
        "nearest_non_gaia_groups": nearest_non_gaia,
        "summary": {
            "evaluated_gaia_worldline_count": len(evaluated),
            "all_evaluated_gaia_worldlines_zero_hits_at_3arcsec": all_zero_3,
            "all_evaluated_gaia_worldlines_have_full_catalog_covariance": full_cov,
            "independent_catalog_epoch_crossmatch": "BLOCKED_UNTIL_PER_DETECTION_EPOCHS_ARE_FROZEN",
        },
        "orbital_bridge": {
            "stellar_catalog_sources": "USE_ASTROMETRIC_WORLDLINE_LAYER",
            "solar_system_or_cartesian_moving_object_state": "USE_JANUS_CELESTIAL_DYNAMICS",
            "do_not_mix": True,
        },
        "epistemic_firewall": {
            "simulation_is_evidence": False,
            "propagated_worldline_is_observation": False,
            "release_year_is_observation_epoch": False,
            "missing_covariance_may_be_replaced_with_zero": False,
            "nearest_wise_only_source_may_be_assigned_gaia_motion": False,
            "negative_result_is_valid": True,
            "i_do_not_know_is_valid": True,
        },
        "claim_ceiling": "GAIA_STELLAR_WORLDLINE_ROBUSTNESS_IN_1900_2100_ONLY__NOT_LOVE_EDEM_IDENTITY",
    }


def self_test() -> None:
    fake = {
        "experiment_id": "fake",
        "groups": {
            "LOVE_GAIA_WISE_2MASS": {
                "target": "LOVE",
                "catalog_queries": {"GAIA_DR3": {"selected": {
                    "source_id": "6163586620213012352",
                    "ra_deg": 204.29573215827, "dec_deg": -36.77911588007,
                    "features": {"Plx": .7101, "e_Plx": .1093, "pmRA": -14.675,
                                 "e_pmRA": .123, "pmDE": -4.695, "e_pmDE": .078},
                }}},
            },
            "EDEM_GAIA_WISE_SDSS_PS1": {
                "target": "EDEM_SEARCH_CENTER_ZP",
                "catalog_queries": {"GAIA_DR3": {"selected": {
                    "source_id": "699051350998534656",
                    "ra_deg": 139.21802643661, "dec_deg": 30.26276387854,
                    "features": {"Plx": .7325, "e_Plx": .3615, "pmRA": 2.124,
                                 "e_pmRA": .337, "pmDE": -8.028, "e_pmDE": .239},
                }}},
            },
        },
    }
    out = run(fake, samples=64, seed=10, live_vizier=False)
    assert out["gaia_groups"]["LOVE_GAIA_WISE_2MASS"]["status"] == "ASTROMETRY_INCOMPLETE"
    assert out["gaia_groups"]["EDEM_GAIA_WISE_SDSS_PS1"]["status"] == "ASTROMETRY_INCOMPLETE"
    print("LOVE_EDEM_EPOCH_WORLDLINES_SELF_TEST=PASS")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input", default="data/love/LOVE-EDEM-SOURCE-CLASSIFICATION-v1-LATEST-RECEIPT.json")
    p.add_argument("--output", default=None)
    p.add_argument("--samples", type=int, default=DEFAULT_SAMPLES)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p.add_argument("--live-vizier", action="store_true")
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        self_test()
        return 0
    classification = json.loads(Path(args.input).read_text(encoding="utf-8"))
    out = run(classification, samples=args.samples, seed=args.seed, live_vizier=args.live_vizier)
    text = json.dumps(out, ensure_ascii=False, indent=2, sort_keys=True)
    if args.output:
        Path(args.output).write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    if args.live_vizier:
        statuses = [x.get("status") for x in out["gaia_groups"].values()]
        if statuses != ["WORLDLINE_EVALUATED", "WORLDLINE_EVALUATED"]:
            return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
