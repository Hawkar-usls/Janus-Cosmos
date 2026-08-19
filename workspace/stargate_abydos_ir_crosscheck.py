from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

TARGET = SkyCoord(223.415064157 * u.deg, 33.979315670 * u.deg, frame="icrs")
FIELD_RADIUS_ARCSEC = 30.0
CENTER_RADIUS_ARCSEC = 2.0
GROUP_RADIUS_ARCSEC = 2.0

CATALOGS = [
    {"key": "ALLWISE", "vizier": "II/328/allwise", "root": "NASA_IPAC_ALLWISE", "ir": True, "ids": ["AllWISE", "WISE"]},
    {"key": "CATWISE2020", "vizier": "II/365/catwise", "root": "NASA_IPAC_CATWISE2020", "ir": True, "ids": ["Name"]},
    {"key": "UNWISE", "vizier": "II/363/unwise", "root": "UNWISE_CATALOG", "ir": True, "ids": ["objid", "unwise_objid", "ID"]},
    {"key": "2MASS_PSC", "vizier": "II/246/out", "root": "NASA_IPAC_2MASS_PSC", "ir": True, "ids": ["_2MASS", "2MASS"]},
    {"key": "GAIA_DR3", "vizier": "I/355/gaiadr3", "root": "ESA_GAIA_DR3", "ir": False, "ids": ["Source"]},
    {"key": "PANSTARRS_DR1", "vizier": "II/349/ps1", "root": "PANSTARRS_DR1", "ir": False, "ids": ["objID", "objName"]},
    {"key": "SDSS_DR16", "vizier": "V/154/sdss16", "root": "SDSS_DR16", "ir": False, "ids": ["SDSS16", "objID"]},
]

RA_CANDIDATES = ["RA_ICRS", "RAJ2000", "RAICRS", "RAdeg", "RA"]
DEC_CANDIDATES = ["DE_ICRS", "DEJ2000", "DEICRS", "DEdeg", "DEC", "DE"]
FEATURE_NAMES = [
    "AllWISE", "Name", "Source", "_2MASS", "2MASS", "objID", "SDSS16",
    "W1mag", "e_W1mag", "W2mag", "e_W2mag", "W3mag", "e_W3mag", "W4mag", "e_W4mag", "qph", "ccf", "ex",
    "W1mproPM", "e_W1mproPM", "W2mproPM", "e_W2mproPM", "snrW1pm", "snrW2pm", "pmRA", "pmDE", "e_pmRA", "e_pmDE", "pmQual", "abf",
    "FW1", "e_FW1", "FW2", "e_FW2", "q_W1", "q_W2", "fFW1", "fFW2",
    "Jmag", "e_Jmag", "Hmag", "e_Hmag", "Kmag", "e_Kmag", "Qflg", "Cflg", "Xflg",
    "Gmag", "BPmag", "RPmag", "BP-RP", "Plx", "e_Plx", "RUWE",
    "gmag", "e_gmag", "rmag", "e_rmag", "imag", "e_imag", "zmag", "e_zmag", "ymag", "e_ymag", "umag", "e_umag", "class",
]


def jsonable(value: Any) -> Any:
    if value is None:
        return None
    try:
        if getattr(value, "mask", False):
            return None
    except Exception:
        pass
    try:
        value = value.item()
    except Exception:
        pass
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def find_col(names: list[str], candidates: list[str]) -> str | None:
    lookup = {n.lower(): n for n in names}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def source_id(row, names: list[str], candidates: list[str]) -> str | None:
    for candidate in candidates:
        if candidate in names:
            value = jsonable(row[candidate])
            if value not in (None, ""):
                return str(value)
    return None


def query_catalog(cfg: dict[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {
        "catalog": cfg["key"],
        "catalog_id": cfg["vizier"],
        "independent_root": cfg["root"],
        "infrared": cfg["ir"],
        "radius_arcsec": FIELD_RADIUS_ARCSEC,
        "status": "UNKNOWN",
        "count": 0,
        "nearest": None,
        "detections": [],
    }
    try:
        viz = Vizier(columns=["*"], row_limit=500)
        tables = viz.query_region(TARGET, radius=FIELD_RADIUS_ARCSEC * u.arcsec, catalog=cfg["vizier"])
        detections: list[dict[str, Any]] = []
        position_table_seen = False
        for table in tables:
            names = list(table.colnames)
            ra_col = find_col(names, RA_CANDIDATES)
            dec_col = find_col(names, DEC_CANDIDATES)
            if not ra_col or not dec_col:
                continue
            position_table_seen = True
            for row in table:
                try:
                    ra = float(row[ra_col])
                    dec = float(row[dec_col])
                except Exception:
                    continue
                coord = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
                sep = float(TARGET.separation(coord).to_value(u.arcsec))
                features = {name: jsonable(row[name]) for name in FEATURE_NAMES if name in names}
                detections.append({
                    "catalog": cfg["key"],
                    "independent_root": cfg["root"],
                    "infrared": cfg["ir"],
                    "source_id": source_id(row, names, cfg["ids"]),
                    "ra_deg": ra,
                    "dec_deg": dec,
                    "separation_from_frozen_center_arcsec": sep,
                    "features": features,
                })
        detections.sort(key=lambda d: d["separation_from_frozen_center_arcsec"])
        result["count"] = len(detections)
        result["detections"] = detections[:50]
        result["nearest"] = detections[0] if detections else None
        if len(tables) == 0 or position_table_seen:
            result["status"] = "OK"
        else:
            result["status"] = "POSITION_COLUMNS_NOT_FOUND"
        return result
    except Exception as exc:
        result["status"] = "QUERY_ERROR"
        result["error"] = f"{type(exc).__name__}: {exc}"
        return result


def cluster(detections: list[dict[str, Any]]) -> list[dict[str, Any]]:
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
            if coords[i].separation(coords[j]).to_value(u.arcsec) <= GROUP_RADIUS_ARCSEC:
                union(i, j)

    grouped: dict[int, list[dict[str, Any]]] = {}
    for i, det in enumerate(detections):
        grouped.setdefault(find(i), []).append(det)

    out = []
    for members in grouped.values():
        members = sorted(members, key=lambda d: d["separation_from_frozen_center_arcsec"])
        ir_members = [m for m in members if m["infrared"]]
        roots = sorted({m["independent_root"] for m in members})
        ir_roots = sorted({m["independent_root"] for m in ir_members})
        out.append({
            "nearest_to_frozen_center_arcsec": members[0]["separation_from_frozen_center_arcsec"],
            "representative_ra_deg": members[0]["ra_deg"],
            "representative_dec_deg": members[0]["dec_deg"],
            "catalogs": sorted({m["catalog"] for m in members}),
            "independent_roots": roots,
            "independent_root_count": len(roots),
            "infrared_catalogs": sorted({m["catalog"] for m in ir_members}),
            "infrared_independent_roots": ir_roots,
            "infrared_independent_root_count": len(ir_roots),
            "members": members,
        })
    return sorted(out, key=lambda g: g["nearest_to_frozen_center_arcsec"])


def main() -> None:
    catalog_results = [query_catalog(cfg) for cfg in CATALOGS]
    detections = []
    for item in catalog_results:
        detections.extend(item["detections"])
    detections.sort(key=lambda d: d["separation_from_frozen_center_arcsec"])
    groups = cluster(detections)

    ir_catalog_keys = {cfg["key"] for cfg in CATALOGS if cfg["ir"]}
    ir_results = [item for item in catalog_results if item["catalog"] in ir_catalog_keys]
    ir_failures = [item for item in ir_results if item["status"] != "OK"]
    ir_detections = [d for d in detections if d["infrared"]]
    center_ir_groups = [
        g for g in groups
        if g["nearest_to_frozen_center_arcsec"] <= CENTER_RADIUS_ARCSEC
        and g["infrared_independent_root_count"] >= 1
    ]
    multi_center_ir = [g for g in center_ir_groups if g["infrared_independent_root_count"] >= 2]

    if multi_center_ir:
        verdict = "MULTI_CATALOG_IR_SOURCE_AT_FROZEN_CENTER"
        plain_ru = "В пределах замороженного радиуса 2 arcsec есть ИК-группа, подтверждённая минимум двумя независимыми ИК-каталогами. Это подтверждает каталогируемый ИК-источник у фиксированного центра, но не его природу."
    elif center_ir_groups:
        verdict = "SINGLE_CATALOG_IR_CANDIDATE_AT_FROZEN_CENTER"
        plain_ru = "В пределах 2 arcsec от замороженного центра есть ИК-каталожная детекция, но независимое многокаталожное подтверждение в этом проходе не получено."
    elif ir_failures:
        verdict = "INCONCLUSIVE_DUE_TO_IR_SERVICE_FAILURE"
        plain_ru = "Надёжно исключить каталогируемый ИК-источник у центра нельзя, потому что хотя бы один из замороженных ИК-каталогов не отработал корректно."
    elif ir_detections:
        verdict = "NO_IR_CATALOG_SOURCE_AT_FROZEN_CENTER__NEARBY_FIELD_SOURCES_EXIST"
        plain_ru = "В пределах 2 arcsec от замороженного центра ИК-каталожного источника нет, но в поле 30 arcsec присутствуют соседние ИК-источники."
    else:
        verdict = "NO_IR_CATALOG_SOURCES_WITHIN_30ARCSEC"
        plain_ru = "В замороженном поле радиусом 30 arcsec ИК-каталожных источников в запрошенных каталогах не найдено. Это не доказывает отсутствие слабого пиксельного ИК-потока ниже порога каталогизации."

    payload = {
        "schema": "janus.cosmos.stargate_abydos.what_is_at_coordinates_ir.receipt.v1",
        "experiment_id": "STARGATE-ABYDOS-WHAT-IS-AT-THESE-COORDINATES-IR-v1",
        "status": "COMPLETE",
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "WHAT_IS_AT_THESE_COORDINATES_WITH_INFRARED_PRIORITY",
        "anomaly_scoring_used": False,
        "upstream_geometry_receipt": "data/stargate/STARGATE-ABYDOS-GEOMETRY-v1-FROZEN-RECEIPT.json",
        "gate_receipt": "data/stargate/STARGATE-ABYDOS-WHAT-IS-AT-THESE-COORDINATES-IR-GATE-v1.json",
        "frozen_target": {"frame": "ICRS", "ra_deg": float(TARGET.ra.deg), "dec_deg": float(TARGET.dec.deg)},
        "frozen_rules": {
            "field_search_radius_arcsec": FIELD_RADIUS_ARCSEC,
            "center_counterpart_radius_arcsec": CENTER_RADIUS_ARCSEC,
            "source_group_link_radius_arcsec": GROUP_RADIUS_ARCSEC,
            "no_post_hoc_radius_changes": True,
            "raw_image_signal_is_distinct_from_catalog_detection": True,
        },
        "catalog_results": catalog_results,
        "source_groups": groups,
        "summary": {
            "catalog_detection_count": len(detections),
            "source_group_count": len(groups),
            "infrared_detection_count": len(ir_detections),
            "nearest_catalog_detection": detections[0] if detections else None,
            "nearest_infrared_detection": ir_detections[0] if ir_detections else None,
            "center_ir_groups": center_ir_groups,
            "multi_catalog_center_ir_groups": multi_center_ir,
            "infrared_service_failures": [
                {"catalog": x["catalog"], "status": x["status"], "error": x.get("error")}
                for x in ir_failures
            ],
            "verdict": verdict,
            "plain_ru": plain_ru,
        },
        "firewall": {
            "catalogued_ir_source_is_real_Abydos": False,
            "catalogued_ir_source_is_planet": False,
            "catalogued_ir_source_is_anomaly": False,
            "absence_from_catalog_proves_no_raw_ir_flux": False,
            "post_hoc_geometry_promoted_to_discovery": False,
            "claim_ceiling": "CATALOG_INVENTORY_AND_COUNTERPART_CLASSIFICATION_AROUND_FROZEN_COORDINATE_ONLY",
        },
        "next_gate": "RAW_WISE_UNWISE_IMAGE_CUTOUT_OR_FORCED_PHOTOMETRY_IF_SCIENTIFICALLY_WARRANTED",
    }

    out = Path("data/stargate/STARGATE-ABYDOS-WHAT-IS-AT-THESE-COORDINATES-IR-v1-LATEST-RECEIPT.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload["summary"], indent=2, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
