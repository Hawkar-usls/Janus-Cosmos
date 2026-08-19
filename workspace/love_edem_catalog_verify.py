from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.mast import Observations
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier

TARGETS = {
    "LOVE": SkyCoord(ra=204.30267916666668 * u.deg, dec=-36.78240527777778 * u.deg, frame="icrs"),
    "EDEM_SEARCH_CENTER_ZP": SkyCoord(ra=139.22409686590188 * u.deg, dec=30.26038779947318 * u.deg, frame="icrs"),
}

CATALOGS = [
    {
        "key": "GAIA_DR3",
        "vizier": "I/355/gaiadr3",
        "independent_root": "ESA_GAIA_DR3",
        "id_columns": ["Source"],
    },
    {
        "key": "ALLWISE",
        "vizier": "II/328/allwise",
        "independent_root": "NASA_IPAC_ALLWISE",
        "id_columns": ["AllWISE", "WISE"],
    },
    {
        "key": "PANSTARRS_DR1",
        "vizier": "II/349/ps1",
        "independent_root": "PANSTARRS_DR1",
        "id_columns": ["objID", "objName"],
    },
    {
        "key": "SDSS_DR16",
        "vizier": "V/154/sdss16",
        "independent_root": "SDSS_DR16",
        "id_columns": ["SDSS16", "objID"],
    },
    {
        "key": "2MASS_PSC",
        "vizier": "II/246/out",
        "independent_root": "2MASS_PSC",
        "id_columns": ["_2MASS", "2MASS"],
    },
]

RA_CANDIDATES = ["RA_ICRS", "RAJ2000", "RAdeg", "RA"]
DEC_CANDIDATES = ["DE_ICRS", "DEJ2000", "DEdeg", "DEC", "DE"]


def _jsonable(value: Any) -> Any:
    if value is None:
        return None
    try:
        if getattr(value, "mask", False):
            return None
    except Exception:
        pass
    if isinstance(value, (str, int, float, bool)):
        return value
    try:
        return value.item()
    except Exception:
        return str(value)


def _find_col(names: list[str], candidates: list[str]) -> str | None:
    lookup = {name.lower(): name for name in names}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def _row_id(row, names: list[str], candidates: list[str]) -> str | None:
    for name in candidates:
        if name in names:
            value = _jsonable(row[name])
            if value not in (None, ""):
                return str(value)
    return None


def _query_vizier(target_name: str, center: SkyCoord, radius: u.Quantity = 30 * u.arcsec) -> dict[str, Any]:
    out: dict[str, Any] = {"target": target_name, "radius_arcsec": float(radius.to_value(u.arcsec)), "catalogs": []}
    detections: list[dict[str, Any]] = []

    for cfg in CATALOGS:
        item: dict[str, Any] = {
            "catalog": cfg["key"],
            "catalog_id": cfg["vizier"],
            "independent_root": cfg["independent_root"],
            "status": "UNKNOWN",
            "count": 0,
            "nearest": None,
        }
        try:
            viz = Vizier(columns=["*"], row_limit=500)
            tables = viz.query_region(center, radius=radius, catalog=cfg["vizier"])
            rows = []
            for table in tables:
                rows.extend(list(table))
            item["count"] = len(rows)
            item["status"] = "OK"
            if not rows:
                out["catalogs"].append(item)
                continue

            table = tables[0]
            names = list(table.colnames)
            ra_col = _find_col(names, RA_CANDIDATES)
            dec_col = _find_col(names, DEC_CANDIDATES)
            if not ra_col or not dec_col:
                item["status"] = "POSITION_COLUMNS_NOT_FOUND"
                item["columns"] = names
                out["catalogs"].append(item)
                continue

            cat_dets = []
            for row in rows:
                try:
                    ra = float(row[ra_col])
                    dec = float(row[dec_col])
                except Exception:
                    continue
                coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
                sep = float(center.separation(coord).to_value(u.arcsec))
                det = {
                    "catalog": cfg["key"],
                    "independent_root": cfg["independent_root"],
                    "source_id": _row_id(row, names, cfg["id_columns"]),
                    "ra_deg": ra,
                    "dec_deg": dec,
                    "separation_arcsec": sep,
                }
                detections.append(det)
                cat_dets.append(det)
            if cat_dets:
                item["nearest"] = min(cat_dets, key=lambda d: d["separation_arcsec"])
            out["catalogs"].append(item)
        except Exception as exc:
            item["status"] = "QUERY_ERROR"
            item["error"] = f"{type(exc).__name__}: {exc}"
            out["catalogs"].append(item)

    out["detections"] = sorted(detections, key=lambda d: d["separation_arcsec"])
    return out


def _cluster(detections: list[dict[str, Any]], threshold_arcsec: float = 1.5) -> list[dict[str, Any]]:
    n = len(detections)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    coords = [SkyCoord(d["ra_deg"] * u.deg, d["dec_deg"] * u.deg, frame="icrs") for d in detections]
    for i in range(n):
        for j in range(i + 1, n):
            if coords[i].separation(coords[j]).to_value(u.arcsec) <= threshold_arcsec:
                union(i, j)

    groups: dict[int, list[dict[str, Any]]] = {}
    for i, det in enumerate(detections):
        groups.setdefault(find(i), []).append(det)

    result = []
    for members in groups.values():
        roots = sorted({m["independent_root"] for m in members})
        catalogs = sorted({m["catalog"] for m in members})
        nearest = min(members, key=lambda d: d["separation_arcsec"])
        result.append({
            "nearest_to_search_center_arcsec": nearest["separation_arcsec"],
            "representative_ra_deg": nearest["ra_deg"],
            "representative_dec_deg": nearest["dec_deg"],
            "catalogs": catalogs,
            "independent_roots": roots,
            "independent_root_count": len(roots),
            "detection_count": len(members),
            "members": sorted(members, key=lambda d: d["separation_arcsec"]),
        })
    return sorted(result, key=lambda g: g["nearest_to_search_center_arcsec"])


def _query_simbad(target_name: str, center: SkyCoord, radius: u.Quantity = 30 * u.arcsec) -> dict[str, Any]:
    result: dict[str, Any] = {"target": target_name, "radius_arcsec": float(radius.to_value(u.arcsec)), "status": "UNKNOWN", "count": 0, "objects": []}
    try:
        sim = Simbad()
        try:
            sim.add_votable_fields("otype")
        except Exception:
            pass
        table = sim.query_region(center, radius=radius)
        if table is None:
            result["status"] = "OK"
            return result
        result["status"] = "OK"
        result["count"] = len(table)
        names = list(table.colnames)
        ra_col = _find_col(names, ["ra", "RA"])
        dec_col = _find_col(names, ["dec", "DEC"])
        id_col = _find_col(names, ["main_id", "MAIN_ID"])
        otype_col = _find_col(names, ["otype", "OTYPE"])
        for row in table:
            obj: dict[str, Any] = {
                "main_id": _jsonable(row[id_col]) if id_col else None,
                "otype": _jsonable(row[otype_col]) if otype_col else None,
            }
            if ra_col and dec_col:
                try:
                    ra = float(row[ra_col])
                    dec = float(row[dec_col])
                    c = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
                    obj.update({"ra_deg": ra, "dec_deg": dec, "separation_arcsec": float(center.separation(c).to_value(u.arcsec))})
                except Exception:
                    pass
            result["objects"].append(obj)
        result["objects"].sort(key=lambda x: x.get("separation_arcsec", 1e99))
        return result
    except Exception as exc:
        result["status"] = "QUERY_ERROR"
        result["error"] = f"{type(exc).__name__}: {exc}"
        return result


def _query_mast_hst(target_name: str, center: SkyCoord, radius: u.Quantity = 6 * u.arcmin) -> dict[str, Any]:
    result: dict[str, Any] = {"target": target_name, "radius_arcmin": float(radius.to_value(u.arcmin)), "status": "UNKNOWN", "hst_observation_count": None, "observations": []}
    try:
        obs = Observations.query_region(center, radius=radius)
        if obs is None or len(obs) == 0:
            result["status"] = "OK"
            result["hst_observation_count"] = 0
            return result
        hst = obs[obs["obs_collection"] == "HST"]
        result["status"] = "OK"
        result["hst_observation_count"] = int(len(hst))
        for row in hst[:50]:
            result["observations"].append({
                "obs_id": _jsonable(row["obs_id"]) if "obs_id" in hst.colnames else None,
                "instrument_name": _jsonable(row["instrument_name"]) if "instrument_name" in hst.colnames else None,
                "filters": _jsonable(row["filters"]) if "filters" in hst.colnames else None,
                "distance_arcmin": float(row["distance"]) if "distance" in hst.colnames and row["distance"] is not None else None,
            })
        return result
    except Exception as exc:
        result["status"] = "QUERY_ERROR"
        result["error"] = f"{type(exc).__name__}: {exc}"
        return result


def verify_target(name: str, center: SkyCoord) -> dict[str, Any]:
    viz = _query_vizier(name, center)
    groups = _cluster(viz["detections"], threshold_arcsec=1.5)
    multi = [g for g in groups if g["independent_root_count"] >= 2]
    return {
        "center": {"ra_deg": float(center.ra.deg), "dec_deg": float(center.dec.deg), "frame": "ICRS"},
        "vizier_inventory": viz,
        "source_groups_1p5arcsec": {
            "group_count": len(groups),
            "multi_independent_root_group_count": len(multi),
            "groups": groups,
        },
        "simbad": _query_simbad(name, center),
        "mast_hst_6arcmin": _query_mast_hst(name, center),
    }


def _service_failures(targets: dict[str, Any]) -> list[dict[str, Any]]:
    failures: list[dict[str, Any]] = []
    for name, item in targets.items():
        for catalog in item["vizier_inventory"]["catalogs"]:
            if catalog["status"] != "OK":
                failures.append({"target": name, "service": catalog["catalog"], "status": catalog["status"], "error": catalog.get("error")})
        if item["simbad"]["status"] != "OK":
            failures.append({"target": name, "service": "SIMBAD", "status": item["simbad"]["status"], "error": item["simbad"].get("error")})
        if item["mast_hst_6arcmin"]["status"] != "OK":
            failures.append({"target": name, "service": "MAST_HST", "status": item["mast_hst_6arcmin"]["status"], "error": item["mast_hst_6arcmin"].get("error")})
    return failures


def _compact_summary(payload: dict[str, Any]) -> dict[str, Any]:
    targets: dict[str, Any] = {}
    for name, item in payload["targets"].items():
        detections = item["vizier_inventory"]["detections"]
        groups = item["source_groups_1p5arcsec"]["groups"]
        targets[name] = {
            "center": item["center"],
            "catalog_counts": {
                c["catalog"]: {"status": c["status"], "count": c["count"], "nearest": c.get("nearest")}
                for c in item["vizier_inventory"]["catalogs"]
            },
            "catalog_detection_count": len(detections),
            "nearest_catalog_detection": detections[0] if detections else None,
            "source_group_count": len(groups),
            "multi_independent_root_group_count": item["source_groups_1p5arcsec"]["multi_independent_root_group_count"],
            "nearest_five_source_groups": groups[:5],
            "strongest_multi_root_groups": sorted(
                [g for g in groups if g["independent_root_count"] >= 2],
                key=lambda g: (-g["independent_root_count"], g["nearest_to_search_center_arcsec"]),
            )[:5],
            "simbad": {
                "status": item["simbad"]["status"],
                "count": item["simbad"]["count"],
                "nearest_objects": item["simbad"]["objects"][:10],
            },
            "mast_hst_6arcmin": item["mast_hst_6arcmin"],
            "plain_status": "CATALOG_SOURCES_PRESENT" if detections else "NO_CATALOG_SOURCES_WITHIN_FROZEN_RADIUS",
        }
    return {
        "schema": "janus.cosmos.love_edem.catalog_crosscheck.summary.v1",
        "experiment_id": payload["experiment_id"],
        "run_time_utc": payload["run_time_utc"],
        "status": payload["status"],
        "mode": payload["mode"],
        "anomaly_scoring_used": payload["anomaly_scoring_used"],
        "search_radius_arcsec": payload["search_radius_arcsec"],
        "cluster_threshold_arcsec": payload["cluster_threshold_arcsec"],
        "service_failures": payload["service_failures"],
        "targets": targets,
        "claim_ceiling": payload["interpretation_contract"]["claim_ceiling"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="JANUS-Cosmos catalog-only WHAT_IS_AT_THESE_COORDINATES verifier for LOVE/EDEM.")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    targets = {name: verify_target(name, center) for name, center in TARGETS.items()}
    failures = _service_failures(targets)
    payload = {
        "schema": "janus.cosmos.love_edem.catalog_crosscheck.receipt.v1",
        "experiment_id": "LOVE-EDEM-CATALOG-CROSSCHECK-v1",
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETE" if not failures else "PARTIAL_SERVICE_FAILURE",
        "mode": "WHAT_IS_AT_THESE_COORDINATES",
        "anomaly_scoring_used": False,
        "anomaly_gate_required": False,
        "search_radius_arcsec": 30.0,
        "cluster_threshold_arcsec": 1.5,
        "targets": targets,
        "service_failures": failures,
        "interpretation_contract": {
            "catalog_detection_means_cataloged_source_near_direction": True,
            "catalog_detection_is_not_anomaly": True,
            "catalog_detection_is_not_planet_claim": True,
            "catalog_detection_is_not_physical_love_edem_association": True,
            "simbad_nonmatch_does_not_mean_empty_sky": True,
            "hst_no_coverage_does_not_mean_empty_sky": True,
            "claim_ceiling": "CATALOG_SOURCE_INVENTORY_AROUND_FROZEN_DIRECTIONS_ONLY",
        },
        "provenance": {
            "catalog_frontend": "CDS VizieR via astroquery.vizier",
            "simbad": "CDS SIMBAD via astroquery.simbad",
            "hst": "MAST/STScI via astroquery.mast",
            "software_lineage": "Janus-Cosmos legacy astronomy stack preserved on legacy-cosmos-pre-osiris-2026-08-19",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    summary = _compact_summary(payload)
    if args.output.name.endswith("-RECEIPT.json"):
        summary_name = args.output.name[:-len("-RECEIPT.json")] + "-SUMMARY.json"
    else:
        summary_name = args.output.stem + "-SUMMARY.json"
    summary_path = args.output.with_name(summary_name)
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": payload["status"],
        "output": str(args.output),
        "summary_output": str(summary_path),
        "summary": {
            name: {
                "catalog_detections": len(item["vizier_inventory"]["detections"]),
                "source_groups": item["source_groups_1p5arcsec"]["group_count"],
                "multi_root_groups": item["source_groups_1p5arcsec"]["multi_independent_root_group_count"],
                "simbad_count": item["simbad"]["count"],
                "hst_observations": item["mast_hst_6arcmin"]["hst_observation_count"],
            }
            for name, item in payload["targets"].items()
        },
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
