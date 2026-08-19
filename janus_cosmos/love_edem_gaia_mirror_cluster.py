from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

from janus_cosmos import love_edem_center_object_probe as primary_probe


TARGETS = primary_probe.TARGETS


@dataclass
class Detection:
    catalog: str
    source_id: str
    ra: float
    dec: float
    payload: dict


def _clean(v):
    try:
        if getattr(v, "mask", False):
            return None
    except Exception:
        pass
    try:
        v = v.item()
    except Exception:
        pass
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    try:
        if math.isnan(float(v)):
            return None
    except Exception:
        pass
    if str(v) in {"--", "nan", "masked"}:
        return None
    return v


def _sep_arcsec(a: Detection, b: Detection) -> float:
    c1 = SkyCoord(a.ra * u.deg, a.dec * u.deg, frame="icrs")
    c2 = SkyCoord(b.ra * u.deg, b.dec * u.deg, frame="icrs")
    return float(c1.separation(c2).arcsec)


def _sep_center_arcsec(ra0: float, dec0: float, ra: float, dec: float) -> float:
    c1 = SkyCoord(ra0 * u.deg, dec0 * u.deg, frame="icrs")
    c2 = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    return float(c1.separation(c2).arcsec)


def gaia_vizier(ra: float, dec: float, radius_arcsec: float) -> dict:
    coord = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    viz = Vizier(columns=["*", "+_r"], row_limit=200)
    tables = viz.query_region(coord, radius=radius_arcsec * u.arcsec, catalog="I/355/gaiadr3")
    if not tables:
        return {"status": "OK", "count": 0, "rows": []}
    t = tables[0]
    wanted = [
        "Source", "RA_ICRS", "DE_ICRS", "Plx", "e_Plx", "pmRA", "pmDE",
        "Gmag", "BPmag", "RPmag", "BP-RP", "RUWE", "_r",
    ]
    rows = []
    for row in t[:200]:
        d = {c: _clean(row[c]) for c in wanted if c in t.colnames}
        r = d.get("RA_ICRS")
        de = d.get("DE_ICRS")
        if r is None or de is None:
            continue
        d["separation_from_center_arcsec"] = _sep_center_arcsec(ra, dec, float(r), float(de))
        plx = d.get("Plx")
        eplx = d.get("e_Plx")
        try:
            if plx is not None and eplx not in (None, 0) and float(plx) > 0:
                snr = float(plx) / float(eplx)
                d["parallax_snr"] = snr
                if snr >= 5:
                    d["naive_inverse_parallax_distance_pc"] = 1000.0 / float(plx)
        except Exception:
            pass
        rows.append(d)
    rows.sort(key=lambda x: x.get("separation_from_center_arcsec", 1e99))
    return {"status": "OK", "count": len(rows), "rows": rows}


def _position(row: dict, options: list[tuple[str, str]]) -> tuple[float, float] | None:
    for rk, dk in options:
        r = row.get(rk)
        d = row.get(dk)
        try:
            if r is not None and d is not None and math.isfinite(float(r)) and math.isfinite(float(d)):
                return float(r), float(d)
        except Exception:
            pass
    return None


def _collect(primary: dict, gaia: dict) -> list[Detection]:
    out: list[Detection] = []

    for row in gaia.get("rows", []):
        pos = _position(row, [("RA_ICRS", "DE_ICRS")])
        if pos:
            out.append(Detection("gaia_dr3_vizier", str(row.get("Source", "")), pos[0], pos[1], row))

    specs = [
        ("allwise", "allwise", [("ra", "dec")], "designation"),
        ("2mass_psc", "2mass_psc", [("ra", "dec")], "designation"),
        ("sdss_dr17", "sdss_dr17", [("ra", "dec")], "objid"),
        ("panstarrs_vizier", "panstarrs_dr2", [("RAJ2000", "DEJ2000"), ("RAdeg", "DEdeg")], "objID"),
        ("simbad", "simbad", [("ra", "dec")], "main_id"),
    ]
    for key, catalog, poskeys, idkey in specs:
        block = primary.get(key, {})
        for row in block.get("rows", []) or []:
            pos = _position(row, poskeys)
            if not pos:
                continue
            sid = row.get(idkey)
            if sid is None:
                sid = f"{catalog}:{pos[0]:.8f},{pos[1]:.8f}"
            out.append(Detection(catalog, str(sid), pos[0], pos[1], row))
    return out


def _pair_tolerance(a: Detection, b: Detection, cfg: dict) -> float:
    cats = {a.catalog, b.catalog}
    if "allwise" in cats:
        return float(cfg["with_allwise"])
    if "2mass_psc" in cats:
        return float(cfg["with_2mass"])
    return float(cfg["default"])


def _spherical_mean(members: list[Detection]) -> tuple[float, float]:
    c = SkyCoord([m.ra for m in members] * u.deg, [m.dec for m in members] * u.deg, frame="icrs")
    cart = c.cartesian
    x = float(cart.x.value.mean())
    y = float(cart.y.value.mean())
    z = float(cart.z.value.mean())
    r = math.sqrt(x*x + y*y + z*z)
    x, y, z = x/r, y/r, z/r
    dec = math.degrees(math.asin(z))
    ra = math.degrees(math.atan2(y, x)) % 360.0
    return ra, dec


def _cluster(detections: list[Detection], cfg: dict, center_ra: float, center_dec: float) -> list[dict]:
    n = len(detections)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            if detections[i].catalog == detections[j].catalog:
                continue
            sep = _sep_arcsec(detections[i], detections[j])
            tol = _pair_tolerance(detections[i], detections[j], cfg)
            if sep <= tol:
                union(i, j)
                edges.append({
                    "a": f"{detections[i].catalog}:{detections[i].source_id}",
                    "b": f"{detections[j].catalog}:{detections[j].source_id}",
                    "separation_arcsec": sep,
                    "tolerance_arcsec": tol,
                })

    groups: dict[int, list[int]] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)

    clusters = []
    for cluster_index, idxs in enumerate(groups.values(), 1):
        members = [detections[i] for i in idxs]
        cra, cdec = _spherical_mean(members)
        catalogs = sorted(set(m.catalog for m in members))
        if "gaia_dr3_vizier" in catalogs:
            label = "GAIA_LINKED_SOURCE_GROUP"
        elif "sdss_dr17" in catalogs and "panstarrs_dr2" in catalogs:
            label = "OPTICAL_CROSS_SURVEY_SOURCE_GROUP"
        elif "allwise" in catalogs and len(catalogs) > 1:
            label = "IR_LINKED_SOURCE_GROUP"
        elif catalogs == ["allwise"]:
            label = "WISE_ONLY_CURRENT_MATCH_RADIUS"
        else:
            label = "SINGLE_CATALOG_SOURCE_GROUP" if len(catalogs) == 1 else "CATALOG_SOURCE_GROUP"

        cluster_edges = [
            e for e in edges
            if any(e["a"] == f"{m.catalog}:{m.source_id}" or e["b"] == f"{m.catalog}:{m.source_id}" for m in members)
        ]
        clusters.append({
            "cluster_id": f"C{cluster_index:03d}",
            "label": label,
            "catalogs": catalogs,
            "catalog_count": len(catalogs),
            "member_count": len(members),
            "center_ra_deg": cra,
            "center_dec_deg": cdec,
            "separation_from_requested_center_arcsec": _sep_center_arcsec(center_ra, center_dec, cra, cdec),
            "members": [
                {
                    "catalog": m.catalog,
                    "source_id": m.source_id,
                    "ra_deg": m.ra,
                    "dec_deg": m.dec,
                    "separation_from_requested_center_arcsec": _sep_center_arcsec(center_ra, center_dec, m.ra, m.dec),
                    "payload": m.payload,
                }
                for m in members
            ],
            "match_edges": cluster_edges,
        })

    clusters.sort(key=lambda c: c["separation_from_requested_center_arcsec"])
    return clusters


def run(prereg_path: Path, output_dir: Path) -> dict:
    prereg = json.loads(prereg_path.read_text(encoding="utf-8"))
    radius = float(prereg["search_radius_arcsec"])
    output_dir.mkdir(parents=True, exist_ok=True)

    primary_path = output_dir / "primary-catalog-result.json"
    primary = primary_probe.run(primary_path, radius)

    result = {
        "schema": "janus.cosmos.love_edem.gaia_mirror_cluster.result.v1",
        "experiment_id": prereg["experiment_id"],
        "mode": prereg["mode"],
        "anomaly_scoring_used": False,
        "gaia_mirror_catalog": prereg["gaia_mirror"]["catalog"],
        "targets": {},
        "firewall": prereg["firewall"],
    }

    for label, target in prereg["targets"].items():
        ra = float(target["ra_deg_icrs"])
        dec = float(target["dec_deg_icrs"])
        try:
            gaia = gaia_vizier(ra, dec, radius)
        except Exception as exc:
            gaia = {"status": "ERROR", "count": 0, "rows": [], "error": f"{type(exc).__name__}: {exc}"}
        p = primary["targets"][label]
        detections = _collect(p, gaia)
        clusters = _cluster(detections, prereg["source_clustering"]["pair_tolerance_arcsec"], ra, dec)
        result["targets"][label] = {
            "semantic_alias": target["semantic_alias"],
            "center_icrs": {"ra_deg": ra, "dec_deg": dec},
            "search_radius_arcsec": radius,
            "direct_gaia_service_status": p.get("gaia_dr3", {}).get("status"),
            "gaia_vizier": gaia,
            "total_catalog_detections": len(detections),
            "cluster_count": len(clusters),
            "multi_catalog_cluster_count": sum(c["catalog_count"] >= 2 for c in clusters),
            "nearest_cluster": clusters[0] if clusters else None,
            "clusters": clusters,
            "interpretation": {
                "question": "what catalogued sources are at or near the coordinate",
                "anomaly_gate_used": False,
                "cluster_is_not_semantic_identity": True,
            },
        }

    out = output_dir / "gaia-mirror-cluster-result.json"
    out.write_text(json.dumps(result, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps({
        label: {
            "gaia_vizier_status": t["gaia_vizier"]["status"],
            "gaia_vizier_count": t["gaia_vizier"]["count"],
            "total_catalog_detections": t["total_catalog_detections"],
            "cluster_count": t["cluster_count"],
            "multi_catalog_cluster_count": t["multi_catalog_cluster_count"],
            "nearest_cluster": t["nearest_cluster"],
        }
        for label, t in result["targets"].items()
    }, indent=2, default=str))
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prereg", default="data/love/LOVE_EDEM_GAIA_MIRROR_CLUSTER_PREREG.json")
    ap.add_argument("--output-dir", default="results/love_edem_gaia_mirror_cluster")
    args = ap.parse_args()
    run(Path(args.prereg), Path(args.output_dir))


if __name__ == "__main__":
    main()
