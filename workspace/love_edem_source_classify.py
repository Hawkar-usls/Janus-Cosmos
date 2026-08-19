from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

CATALOGS = {
    "GAIA_DR3": {
        "vizier": "I/355/gaiadr3",
        "id_candidates": ["Source"],
    },
    "ALLWISE": {
        "vizier": "II/328/allwise",
        "id_candidates": ["AllWISE", "WISE"],
    },
    "2MASS_PSC": {
        "vizier": "II/246/out",
        "id_candidates": ["_2MASS", "2MASS"],
    },
    "PANSTARRS_DR1": {
        "vizier": "II/349/ps1",
        "id_candidates": ["objID", "objName"],
    },
    "SDSS_DR16": {
        "vizier": "V/154/sdss16",
        "id_candidates": ["SDSS16", "objID"],
    },
}

GROUPS = {
    "LOVE_NEAREST_WISE_ONLY": {
        "target": "LOVE",
        "center": [204.2960668, -36.7804774],
        "expected": {"ALLWISE": "J133711.05-364649.7"},
        "catalogs": ["ALLWISE", "GAIA_DR3", "2MASS_PSC", "PANSTARRS_DR1", "SDSS_DR16"],
        "frozen_center_offset_arcsec": 20.289587649207107,
    },
    "LOVE_GAIA_WISE_2MASS": {
        "target": "LOVE",
        "center": [204.29573215827, -36.77911588007],
        "expected": {
            "GAIA_DR3": "6163586620213012352",
            "ALLWISE": "J133710.98-364644.9",
            "2MASS_PSC": "13371100-3646446",
        },
        "catalogs": ["GAIA_DR3", "ALLWISE", "2MASS_PSC", "PANSTARRS_DR1", "SDSS_DR16"],
        "frozen_center_offset_arcsec": 23.269252366704684,
    },
    "EDEM_NEAREST_WISE_ONLY": {
        "target": "EDEM_SEARCH_CENTER_ZP",
        "center": [139.2249571, 30.262781],
        "expected": {"ALLWISE": "J091653.98+301546.0"},
        "catalogs": ["ALLWISE", "GAIA_DR3", "2MASS_PSC", "PANSTARRS_DR1", "SDSS_DR16"],
        "frozen_center_offset_arcsec": 9.021198661597598,
    },
    "EDEM_OPTICAL_SDSS_PS1": {
        "target": "EDEM_SEARCH_CENTER_ZP",
        "center": [139.222681, 30.256918],
        "expected": {
            "PANSTARRS_DR1": "144301392226548847",
            "SDSS_DR16": "SDSS J091653.44+301524.9",
        },
        "catalogs": ["SDSS_DR16", "PANSTARRS_DR1", "GAIA_DR3", "ALLWISE", "2MASS_PSC"],
        "frozen_center_offset_arcsec": 13.244457192554428,
    },
    "EDEM_GAIA_WISE_SDSS_PS1": {
        "target": "EDEM_SEARCH_CENTER_ZP",
        "center": [139.21802643661, 30.26276387854],
        "expected": {
            "GAIA_DR3": "699051350998534656",
            "ALLWISE": "J091652.31+301545.8",
            "PANSTARRS_DR1": "144311392180445928",
            "SDSS_DR16": "SDSS J091652.33+301546.1",
        },
        "catalogs": ["GAIA_DR3", "ALLWISE", "PANSTARRS_DR1", "SDSS_DR16", "2MASS_PSC"],
        "frozen_center_offset_arcsec": 20.72338948797579,
    },
}

RA_CANDIDATES = ["RA_ICRS", "RAJ2000", "RAdeg", "RA"]
DEC_CANDIDATES = ["DE_ICRS", "DEJ2000", "DEdeg", "DEC", "DE"]

FEATURE_CANDIDATES = {
    "GAIA_DR3": [
        "Source", "RA_ICRS", "DE_ICRS", "Plx", "e_Plx", "pmRA", "e_pmRA", "pmDE", "e_pmDE",
        "Gmag", "BPmag", "RPmag", "BP-RP", "RUWE", "Teff", "logg", "[Fe/H]",
    ],
    "ALLWISE": [
        "AllWISE", "RAJ2000", "DEJ2000", "W1mag", "e_W1mag", "W2mag", "e_W2mag",
        "W3mag", "e_W3mag", "W4mag", "e_W4mag", "Jmag", "Hmag", "Kmag", "qph", "ccf", "ex",
    ],
    "2MASS_PSC": [
        "_2MASS", "RAJ2000", "DEJ2000", "Jmag", "e_Jmag", "Hmag", "e_Hmag", "Kmag", "e_Kmag", "Qflg", "Rflg", "Bflg", "Cflg",
    ],
    "PANSTARRS_DR1": [
        "objID", "RAJ2000", "DEJ2000", "gmag", "e_gmag", "rmag", "e_rmag", "imag", "e_imag", "zmag", "e_zmag", "ymag", "e_ymag", "Qual", "Nd",
    ],
    "SDSS_DR16": [
        "SDSS16", "objID", "RA_ICRS", "DE_ICRS", "umag", "e_umag", "gmag", "e_gmag", "rmag", "e_rmag", "imag", "e_imag", "zmag", "e_zmag", "class", "q_mode",
    ],
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
        x = v.item()
        if isinstance(x, (str, int, float, bool)):
            return x
        return str(x)
    except Exception:
        return str(v)


def _find_col(names: list[str], candidates: list[str]) -> str | None:
    low = {x.lower(): x for x in names}
    for c in candidates:
        if c.lower() in low:
            return low[c.lower()]
    return None


def _id_from_row(row, names: list[str], candidates: list[str]) -> str | None:
    for c in candidates:
        if c in names:
            v = _jsonable(row[c])
            if v not in (None, ""):
                return str(v)
    return None


def _float(d: dict[str, Any], *names: str) -> float | None:
    for name in names:
        if name in d and d[name] is not None:
            try:
                x = float(d[name])
                if math.isfinite(x):
                    return x
            except Exception:
                pass
    return None


def _query_catalog(group_name: str, group: dict[str, Any], catalog_key: str, radius: u.Quantity = 2.0 * u.arcsec) -> dict[str, Any]:
    cfg = CATALOGS[catalog_key]
    center = SkyCoord(group["center"][0] * u.deg, group["center"][1] * u.deg, frame="icrs")
    result: dict[str, Any] = {
        "catalog": catalog_key,
        "catalog_id": cfg["vizier"],
        "radius_arcsec": float(radius.to_value(u.arcsec)),
        "status": "UNKNOWN",
        "count": 0,
        "selected": None,
    }
    try:
        viz = Vizier(columns=["*"], row_limit=100)
        tables = viz.query_region(center, radius=radius, catalog=cfg["vizier"])
        rows = []
        table_for_row = []
        for table in tables:
            for row in table:
                rows.append(row)
                table_for_row.append(table)
        result["count"] = len(rows)
        result["status"] = "OK"
        if not rows:
            return result

        expected = group.get("expected", {}).get(catalog_key)
        candidates = []
        for row, table in zip(rows, table_for_row):
            names = list(table.colnames)
            ra_col = _find_col(names, RA_CANDIDATES)
            dec_col = _find_col(names, DEC_CANDIDATES)
            sep = None
            ra = dec = None
            if ra_col and dec_col:
                try:
                    ra = float(row[ra_col])
                    dec = float(row[dec_col])
                    sep = float(center.separation(SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")).to_value(u.arcsec))
                except Exception:
                    pass
            sid = _id_from_row(row, names, cfg["id_candidates"])
            snapshot: dict[str, Any] = {}
            for wanted in FEATURE_CANDIDATES[catalog_key]:
                if wanted in names:
                    snapshot[wanted] = _jsonable(row[wanted])
            # Preserve additional columns that are explicitly morphology/quality/astrometry relevant.
            for name in names:
                lname = name.lower()
                if any(token in lname for token in ["psf", "model", "type", "class", "prob", "ruwe", "plx", "pmra", "pmde", "snr", "qual"]):
                    if name not in snapshot:
                        snapshot[name] = _jsonable(row[name])
            candidates.append({
                "source_id": sid,
                "ra_deg": ra,
                "dec_deg": dec,
                "separation_from_group_arcsec": sep,
                "features": snapshot,
                "expected_id_match": bool(expected is not None and sid is not None and str(sid).strip() == str(expected).strip()),
            })

        chosen = None
        matches = [x for x in candidates if x["expected_id_match"]]
        if matches:
            chosen = min(matches, key=lambda x: x["separation_from_group_arcsec"] if x["separation_from_group_arcsec"] is not None else 1e99)
        else:
            chosen = min(candidates, key=lambda x: x["separation_from_group_arcsec"] if x["separation_from_group_arcsec"] is not None else 1e99)
        result["selected"] = chosen
        result["candidate_count"] = len(candidates)
        return result
    except Exception as exc:
        result["status"] = "QUERY_ERROR"
        result["error"] = f"{type(exc).__name__}: {exc}"
        return result


def _gaia_derived(selected: dict[str, Any] | None) -> dict[str, Any] | None:
    if not selected:
        return None
    f = selected.get("features", {})
    plx = _float(f, "Plx", "parallax")
    eplx = _float(f, "e_Plx", "parallax_error")
    pmra = _float(f, "pmRA", "pmra")
    epmra = _float(f, "e_pmRA", "pmra_error")
    pmde = _float(f, "pmDE", "pmdec")
    epmde = _float(f, "e_pmDE", "pmdec_error")
    gmag = _float(f, "Gmag", "phot_g_mean_mag")
    bpmag = _float(f, "BPmag", "phot_bp_mean_mag")
    rpmag = _float(f, "RPmag", "phot_rp_mean_mag")
    bp_rp = _float(f, "BP-RP", "bp_rp")
    ruwe = _float(f, "RUWE", "ruwe")

    plx_snr = plx / eplx if plx is not None and eplx not in (None, 0) else None
    pmra_snr = abs(pmra / epmra) if pmra is not None and epmra not in (None, 0) else None
    pmde_snr = abs(pmde / epmde) if pmde is not None and epmde not in (None, 0) else None
    pm_vector_snr = None
    if pmra_snr is not None and pmde_snr is not None:
        pm_vector_snr = math.sqrt(pmra_snr ** 2 + pmde_snr ** 2)

    distance_pc = None
    abs_g_no_extinction = None
    if plx is not None and plx > 0 and plx_snr is not None and plx_snr >= 5:
        distance_pc = 1000.0 / plx
        if gmag is not None:
            abs_g_no_extinction = gmag - 5.0 * math.log10(distance_pc / 10.0)

    return {
        "parallax_mas": plx,
        "parallax_error_mas": eplx,
        "parallax_snr": plx_snr,
        "pmra_masyr": pmra,
        "pmra_error_masyr": epmra,
        "pmdec_masyr": pmde,
        "pmdec_error_masyr": epmde,
        "proper_motion_vector_snr": pm_vector_snr,
        "g_mag": gmag,
        "bp_mag": bpmag,
        "rp_mag": rpmag,
        "bp_rp": bp_rp,
        "ruwe": ruwe,
        "naive_inverse_parallax_distance_pc_if_snr_ge_5": distance_pc,
        "absolute_g_no_extinction_if_distance_available": abs_g_no_extinction,
    }


def _wise_colors(selected: dict[str, Any] | None) -> dict[str, Any] | None:
    if not selected:
        return None
    f = selected.get("features", {})
    w1 = _float(f, "W1mag")
    w2 = _float(f, "W2mag")
    w3 = _float(f, "W3mag")
    w4 = _float(f, "W4mag")
    return {
        "W1": w1,
        "W2": w2,
        "W3": w3,
        "W4": w4,
        "W1_minus_W2": (w1 - w2) if w1 is not None and w2 is not None else None,
        "W2_minus_W3": (w2 - w3) if w2 is not None and w3 is not None else None,
        "quality": f.get("qph"),
        "contamination_flags": f.get("ccf"),
        "extended_flag": f.get("ex"),
    }


def _classify(group_name: str, query: dict[str, Any]) -> dict[str, Any]:
    gaia = query.get("GAIA_DR3", {}).get("selected")
    gd = _gaia_derived(gaia)
    wise = query.get("ALLWISE", {}).get("selected")
    twomass = query.get("2MASS_PSC", {}).get("selected")
    ps1 = query.get("PANSTARRS_DR1", {}).get("selected")
    sdss = query.get("SDSS_DR16", {}).get("selected")

    roots = [k for k, v in query.items() if v.get("status") == "OK" and v.get("selected") is not None]
    reasons: list[str] = []

    if gd:
        plx_snr = gd.get("parallax_snr")
        pm_snr = gd.get("proper_motion_vector_snr")
        if (plx_snr is not None and plx_snr >= 5) or (pm_snr is not None and pm_snr >= 5):
            reasons.append("Gaia astrometry shows statistically significant parallax and/or proper motion, which is strong evidence for a Galactic star.")
            if gd.get("ruwe") is not None and gd["ruwe"] < 1.4:
                reasons.append("Gaia RUWE is below 1.4, consistent with a well-behaved single-source astrometric solution.")
            return {
                "class": "GALACTIC_STAR",
                "confidence": "HIGH",
                "reasons": reasons,
                "independent_catalog_roots": roots,
                "gaia_derived": gd,
                "wise_derived": _wise_colors(wise),
            }
        reasons.append("A Gaia counterpart exists, but current parallax/proper-motion significance is not strong enough by the frozen >=5 sigma rule.")
        if len(roots) >= 3:
            reasons.append("The same compact position is detected by multiple independent surveys.")
        return {
            "class": "LIKELY_STAR_OR_COMPACT_SOURCE",
            "confidence": "MEDIUM",
            "reasons": reasons,
            "independent_catalog_roots": roots,
            "gaia_derived": gd,
            "wise_derived": _wise_colors(wise),
        }

    if sdss is not None and ps1 is not None:
        reasons.append("The source is independently detected in SDSS and Pan-STARRS at the same optical position.")
        reasons.append("No Gaia counterpart was recovered within the frozen 2 arcsec classification radius, so stellar versus compact-galaxy identity is not established by astrometry.")
        return {
            "class": "OPTICAL_SOURCE_STAR_OR_GALAXY_UNRESOLVED",
            "confidence": "MEDIUM",
            "reasons": reasons,
            "independent_catalog_roots": roots,
            "gaia_derived": None,
            "wise_derived": _wise_colors(wise),
        }

    if wise is not None and len(roots) == 1:
        reasons.append("Only an AllWISE counterpart is recovered within the frozen 2 arcsec classification radius.")
        reasons.append("A single-band or weak mid-infrared detection without optical/Gaia/2MASS support is insufficient to distinguish a faint star, galaxy, or low-S/N infrared source.")
        return {
            "class": "FAINT_INFRARED_SOURCE_UNCLASSIFIED",
            "confidence": "LOW_TO_MEDIUM",
            "reasons": reasons,
            "independent_catalog_roots": roots,
            "gaia_derived": None,
            "wise_derived": _wise_colors(wise),
        }

    reasons.append("The available cross-survey evidence does not satisfy any stronger frozen classification rule.")
    return {
        "class": "UNCLASSIFIED_CATALOG_SOURCE",
        "confidence": "LOW",
        "reasons": reasons,
        "independent_catalog_roots": roots,
        "gaia_derived": gd,
        "wise_derived": _wise_colors(wise),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    groups_out: dict[str, Any] = {}
    failures: list[dict[str, str]] = []
    for name, group in GROUPS.items():
        query: dict[str, Any] = {}
        for catalog in group["catalogs"]:
            r = _query_catalog(name, group, catalog)
            query[catalog] = r
            if r.get("status") == "QUERY_ERROR":
                failures.append({"group": name, "catalog": catalog, "error": r.get("error", "unknown")})
        groups_out[name] = {
            "target": group["target"],
            "representative_icrs": {"ra_deg": group["center"][0], "dec_deg": group["center"][1]},
            "offset_from_frozen_target_center_arcsec": group["frozen_center_offset_arcsec"],
            "catalog_queries": query,
            "classification": _classify(name, query),
        }

    payload = {
        "schema": "janus.cosmos.love_edem.source_classification.receipt.v1",
        "experiment_id": "LOVE-EDEM-SOURCE-CLASSIFICATION-v1",
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETE" if not failures else "COMPLETE_WITH_SERVICE_FAILURES",
        "mode": "ASTROPHYSICAL_CLASSIFICATION_OF_FROZEN_SOURCE_GROUPS",
        "anomaly_scoring_used": False,
        "classification_radius_arcsec": 2.0,
        "frozen_rules": {
            "galactic_star_high_confidence": "Gaia parallax S/N >= 5 OR combined proper-motion significance >= 5 sigma",
            "gaia_compact_medium": "Gaia counterpart exists but does not pass high-confidence astrometric star rule",
            "optical_unresolved_medium": "SDSS + Pan-STARRS counterpart, no Gaia within 2 arcsec",
            "wise_only": "AllWISE only within 2 arcsec => faint infrared source, class unresolved",
            "no_post_hoc_threshold_changes": True,
        },
        "groups": groups_out,
        "service_failures": failures,
        "interpretation_contract": {
            "classification_is_of_catalog_source_group_not_of_semantic_name": True,
            "group_near_LOVE_is_not_the_identity_LOVE": True,
            "group_near_EDEM_is_not_the_identity_EDEM": True,
            "catalog_classification_is_not_anomaly_detection": True,
            "physical_LOVE_EDEM_Orion_association_not_established": True,
            "claim_ceiling": "ASTROPHYSICAL_CLASSIFICATION_OF_NEARBY_CATALOG_SOURCE_GROUPS_ONLY",
        },
        "source_lineage": [
            "data/love/LOVE-EDEM-CATALOG-CROSSCHECK-v1-LATEST-RECEIPT.json",
            "research/love-g1-007-multieye-v2:data/love/LOVE-EDEM-GAIA-MIRROR-CLUSTER-v1-RUN-001-RECEIPT.json",
        ],
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "output": str(args.output),
        "classifications": {k: v["classification"] for k, v in groups_out.items()},
        "service_failures": failures,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
