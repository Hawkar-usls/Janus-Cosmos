from __future__ import annotations

import argparse
import csv
import io
import json
import math
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from astropy import units as u
from astropy.coordinates import SkyCoord


def _jsonable(value: Any):
    if value is None:
        return None
    try:
        if getattr(value, "mask", False) is True:
            return None
    except Exception:
        pass
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    if isinstance(value, (str, int, float, bool)):
        return value
    try:
        if math.isnan(float(value)):
            return None
    except Exception:
        pass
    return str(value)


def _distance_pc_from_parallax(parallax_mas: float) -> float:
    return 1000.0 / parallax_mas


def _separation_3d_pc(ra1_deg: float, dec1_deg: float, d1_pc: float,
                      ra2_deg: float, dec2_deg: float, d2_pc: float) -> tuple[float, float]:
    c1 = SkyCoord(ra1_deg * u.deg, dec1_deg * u.deg, distance=d1_pc * u.pc, frame="icrs")
    c2 = SkyCoord(ra2_deg * u.deg, dec2_deg * u.deg, distance=d2_pc * u.pc, frame="icrs")
    return float(c1.separation_3d(c2).pc), float(c1.separation(c2).deg)


def _k_like(teff: float | None, bp_rp: float | None, prereg: dict) -> bool:
    t0, t1 = prereg["ranking"]["k_like_temperature_k"]
    c0, c1 = prereg["ranking"]["k_like_bp_rp_fallback"]
    if teff is not None and math.isfinite(teff):
        return t0 <= teff <= t1
    return bp_rp is not None and math.isfinite(bp_rp) and c0 <= bp_rp <= c1


def query_gaia(prereg: dict) -> tuple[list[dict], str | None]:
    try:
        from astroquery.gaia import Gaia

        o = prereg["origin"]["edsm_icrs"]
        s = prereg["search_volume"]
        query = f"""
        SELECT TOP 5000
          g.source_id, g.ra, g.dec, g.parallax, g.parallax_error,
          g.parallax_over_error, g.pmra, g.pmdec, g.ruwe,
          g.phot_g_mean_mag, g.bp_rp,
          ap.teff_gspphot, ap.logg_gspphot, ap.mh_gspphot
        FROM gaiadr3.gaia_source AS g
        LEFT OUTER JOIN gaiadr3.astrophysical_parameters AS ap
          ON g.source_id = ap.source_id
        WHERE 1 = CONTAINS(
          POINT('ICRS', g.ra, g.dec),
          CIRCLE('ICRS', {o['ra_deg']}, {o['dec_deg']}, {s['angular_radius_deg']})
        )
          AND g.parallax BETWEEN {s['parallax_min_mas']} AND {s['parallax_max_mas']}
          AND g.parallax_over_error >= {s['gaia_parallax_over_error_min']}
        """
        job = Gaia.launch_job_async(query, dump_to_file=False)
        table = job.get_results()
        target_d = float(prereg["origin"]["target_distance_pc"])
        rows: list[dict] = []
        for row in table:
            p = float(row["parallax"])
            d = _distance_pc_from_parallax(p)
            sep3d, sepang = _separation_3d_pc(
                float(row["ra"]), float(row["dec"]), d,
                float(o["ra_deg"]), float(o["dec_deg"]), target_d,
            )
            teff = _jsonable(row["teff_gspphot"])
            bp_rp = _jsonable(row["bp_rp"])
            teff_f = float(teff) if teff is not None else None
            bp_rp_f = float(bp_rp) if bp_rp is not None else None
            ruwe = _jsonable(row["ruwe"])
            ruwe_f = float(ruwe) if ruwe is not None else None
            rows.append({
                "source_id": str(row["source_id"]),
                "ra_deg": float(row["ra"]),
                "dec_deg": float(row["dec"]),
                "parallax_mas": p,
                "distance_pc": d,
                "distance_from_lave_point_pc": sep3d,
                "angular_separation_deg": sepang,
                "parallax_over_error": float(row["parallax_over_error"]),
                "ruwe": ruwe_f,
                "phot_g_mean_mag": _jsonable(row["phot_g_mean_mag"]),
                "bp_rp": bp_rp_f,
                "teff_gspphot_k": teff_f,
                "logg_gspphot": _jsonable(row["logg_gspphot"]),
                "mh_gspphot": _jsonable(row["mh_gspphot"]),
                "k_like": _k_like(teff_f, bp_rp_f, prereg),
                "preferred_astrometry": ruwe_f is None or ruwe_f <= float(s["gaia_ruwe_max_preferred"]),
            })
        rows.sort(key=lambda x: (x["distance_from_lave_point_pc"], not x["k_like"], not x["preferred_astrometry"]))
        return rows[: int(s["top_ranked_sources"])], None
    except Exception as exc:
        return [], f"{type(exc).__name__}: {exc}"


def query_nasa_exoplanets(prereg: dict) -> tuple[list[dict], str | None]:
    try:
        o = prereg["origin"]["edsm_icrs"]
        radius = float(prereg["search_volume"]["angular_radius_deg"])
        query = f"""
        select pl_name,hostname,ra,dec,sy_dist,st_spectype,disc_year,discoverymethod
        from pscomppars
        where contains(point('ICRS',ra,dec),circle('ICRS',{o['ra_deg']},{o['dec_deg']},{radius}))=1
          and sy_dist between 20 and 60
        """.strip()
        url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urllib.parse.urlencode({
            "query": query,
            "format": "csv",
        })
        req = urllib.request.Request(url, headers={"User-Agent": "Janus-Cosmos-LOVE/1.0"})
        with urllib.request.urlopen(req, timeout=60) as response:
            text = response.read().decode("utf-8")
        rows = []
        for r in csv.DictReader(io.StringIO(text)):
            rows.append({
                "pl_name": r.get("pl_name"),
                "hostname": r.get("hostname"),
                "ra_deg": float(r["ra"]) if r.get("ra") else None,
                "dec_deg": float(r["dec"]) if r.get("dec") else None,
                "distance_pc": float(r["sy_dist"]) if r.get("sy_dist") else None,
                "spectral_type": r.get("st_spectype") or None,
                "discovery_year": int(float(r["disc_year"])) if r.get("disc_year") else None,
                "discovery_method": r.get("discoverymethod") or None,
                "known_before_love_search": True,
            })
        return rows, None
    except Exception as exc:
        return [], f"{type(exc).__name__}: {exc}"


def probe_tic_tess(candidates: list[dict], prereg: dict) -> tuple[list[dict], list[str]]:
    errors: list[str] = []
    if not candidates:
        return candidates, errors
    try:
        from astroquery.mast import Catalogs, Observations
    except Exception as exc:
        return candidates, [f"MAST_IMPORT: {type(exc).__name__}: {exc}"]

    n = int(prereg["search_volume"]["tess_probe_top_sources"])
    for item in candidates[:n]:
        coord = SkyCoord(item["ra_deg"] * u.deg, item["dec_deg"] * u.deg)
        try:
            tic = Catalogs.query_region(coord, radius=5 * u.arcsec, catalog="TIC")
            if len(tic) == 0:
                item["tic"] = {"status": "NO_TIC_MATCH"}
                continue
            tic.sort("dstArcSec") if "dstArcSec" in tic.colnames else None
            t = tic[0]
            tic_id = str(t["ID"])
            item["tic"] = {
                "status": "MATCH",
                "tic_id": tic_id,
                "separation_arcsec": float(t["dstArcSec"]) if "dstArcSec" in tic.colnames else None,
                "tmag": _jsonable(t["Tmag"]) if "Tmag" in tic.colnames else None,
            }
            try:
                obs = Observations.query_criteria(
                    obs_collection="TESS",
                    target_name=tic_id,
                    dataproduct_type="timeseries",
                )
                item["tess"] = {
                    "status": "AVAILABLE" if len(obs) else "NO_TARGETED_TIMESERIES",
                    "observation_count": int(len(obs)),
                }
            except Exception as exc:
                item["tess"] = {"status": "QUERY_ERROR", "error": f"{type(exc).__name__}: {exc}"}
        except Exception as exc:
            item["tic"] = {"status": "QUERY_ERROR", "error": f"{type(exc).__name__}: {exc}"}
            errors.append(f"Gaia {item['source_id']}: {type(exc).__name__}: {exc}")
    return candidates, errors


def run(prereg_path: Path, output_path: Path) -> dict:
    prereg = json.loads(prereg_path.read_text(encoding="utf-8"))
    gaia, gaia_error = query_gaia(prereg)
    known_planets, nasa_error = query_nasa_exoplanets(prereg)
    gaia, mast_errors = probe_tic_tess(gaia, prereg)

    k_like = [x for x in gaia if x.get("k_like")]
    tess_ready = [x for x in gaia if x.get("tess", {}).get("status") == "AVAILABLE"]
    result = {
        "schema": "janus.cosmos.love_search.result.v1",
        "experiment_id": prereg["experiment_id"],
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "prereg": str(prereg_path),
        "target": prereg["origin"],
        "gaia": {
            "status": "OK" if gaia else ("ERROR" if gaia_error else "NO_MATCHES"),
            "error": gaia_error,
            "ranked_source_count": len(gaia),
            "k_like_count": len(k_like),
            "top_sources": gaia,
        },
        "nasa_exoplanet_archive": {
            "status": "OK" if nasa_error is None else "ERROR",
            "error": nasa_error,
            "known_planet_count_in_cone_and_distance_slice": len(known_planets),
            "known_planets": known_planets,
        },
        "mast_tess": {
            "probed_source_count": min(len(gaia), int(prereg["search_volume"]["tess_probe_top_sources"])),
            "targeted_timeseries_available_count": len(tess_ready),
            "errors": mast_errors,
        },
        "love_gate": {
            "candidate_codename": "Love",
            "candidate_activated": False,
            "reason": "This pass resolves real stellar counterparts and known-planet/TESS coverage only; it does not detect a new planetary signal.",
            "next_gate": "TESS_BLS_TRANSIT_SEARCH_ON_TOP_UNCATALOGUED_K_LIKE_COUNTERPARTS",
        },
        "claim_ceiling": prereg["claim_ceiling"],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description="Resolve the Elite/EDSM Lave point into real Gaia/NASA/MAST targets.")
    ap.add_argument("--prereg", default="data/love/LAVE_TO_LOVE_PREREG.json")
    ap.add_argument("--output", default="results/love/love-real-sky-result.json")
    args = ap.parse_args()
    result = run(Path(args.prereg), Path(args.output))
    print(json.dumps({
        "experiment_id": result["experiment_id"],
        "gaia_ranked": result["gaia"]["ranked_source_count"],
        "gaia_k_like": result["gaia"]["k_like_count"],
        "known_planets": result["nasa_exoplanet_archive"]["known_planet_count_in_cone_and_distance_slice"],
        "tess_ready": result["mast_tess"]["targeted_timeseries_available_count"],
        "love_candidate_activated": result["love_gate"]["candidate_activated"],
    }, indent=2))
    return 0 if result["gaia"]["status"] != "ERROR" else 2


if __name__ == "__main__":
    raise SystemExit(main())
