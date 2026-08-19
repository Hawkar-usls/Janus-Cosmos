from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from janus_cosmos.discovery import product_rank
from janus_cosmos.luci import read_luci_fits_image


ANALYSIS_LABEL = "TARGET_A"


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
        f = float(value)
        if math.isnan(f):
            return None
        return f
    except Exception:
        return str(value)


def _distance_pc_from_parallax(parallax_mas: float) -> float:
    return 1000.0 / parallax_mas


def _separation_3d_pc(ra1: float, dec1: float, d1: float, ra2: float, dec2: float, d2: float) -> tuple[float, float]:
    a = SkyCoord(ra1 * u.deg, dec1 * u.deg, distance=d1 * u.pc, frame="icrs")
    b = SkyCoord(ra2 * u.deg, dec2 * u.deg, distance=d2 * u.pc, frame="icrs")
    return float(a.separation_3d(b).pc), float(a.separation(b).deg)


def _k_like(teff: float | None, bp_rp: float | None, cfg: dict) -> bool:
    t0, t1 = cfg["k_like_temperature_k"]
    c0, c1 = cfg["k_like_bp_rp_fallback"]
    if teff is not None and math.isfinite(teff):
        return t0 <= teff <= t1
    return bp_rp is not None and math.isfinite(bp_rp) and c0 <= bp_rp <= c1


def query_gaia(prereg: dict) -> tuple[list[dict], str | None]:
    try:
        from astroquery.gaia import Gaia

        o = prereg["origin"]["edsm_icrs"]
        s = prereg["gaia"]
        adql = f"""
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
          AND g.parallax_over_error >= {s['parallax_over_error_min']}
        """
        table = Gaia.launch_job_async(adql, dump_to_file=False).get_results()
        target_d = float(prereg["origin"]["target_distance_pc"])
        rows: list[dict] = []
        for row in table:
            parallax = float(row["parallax"])
            distance = _distance_pc_from_parallax(parallax)
            sep3d, sepang = _separation_3d_pc(
                float(row["ra"]), float(row["dec"]), distance,
                float(o["ra_deg"]), float(o["dec_deg"]), target_d,
            )
            teff = _jsonable(row["teff_gspphot"])
            bp_rp = _jsonable(row["bp_rp"])
            ruwe = _jsonable(row["ruwe"])
            teff_f = float(teff) if teff is not None else None
            bp_rp_f = float(bp_rp) if bp_rp is not None else None
            ruwe_f = float(ruwe) if ruwe is not None else None
            rows.append({
                "blind_label": ANALYSIS_LABEL,
                "source_id": str(row["source_id"]),
                "ra_deg": float(row["ra"]),
                "dec_deg": float(row["dec"]),
                "parallax_mas": parallax,
                "distance_pc": distance,
                "distance_from_target_point_pc": sep3d,
                "angular_separation_from_target_deg": sepang,
                "parallax_over_error": float(row["parallax_over_error"]),
                "ruwe": ruwe_f,
                "phot_g_mean_mag": _jsonable(row["phot_g_mean_mag"]),
                "bp_rp": bp_rp_f,
                "teff_gspphot_k": teff_f,
                "logg_gspphot": _jsonable(row["logg_gspphot"]),
                "mh_gspphot": _jsonable(row["mh_gspphot"]),
                "k_like": _k_like(teff_f, bp_rp_f, s),
                "preferred_astrometry": ruwe_f is None or ruwe_f <= float(s["ruwe_max_preferred"]),
            })
        rows.sort(key=lambda x: (
            x["distance_from_target_point_pc"],
            not x["k_like"],
            not x["preferred_astrometry"],
        ))
        return rows[: int(s["top_sources"])], None
    except Exception as exc:
        return [], f"{type(exc).__name__}: {exc}"


def query_known_planets(prereg: dict) -> tuple[list[dict], str | None]:
    try:
        o = prereg["origin"]["edsm_icrs"]
        radius = float(prereg["gaia"]["angular_radius_deg"])
        q = f"""
        select pl_name,hostname,ra,dec,sy_dist,st_spectype,disc_year,discoverymethod
        from pscomppars
        where contains(point('ICRS',ra,dec),circle('ICRS',{o['ra_deg']},{o['dec_deg']},{radius}))=1
          and sy_dist between 20 and 60
        """.strip()
        url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urllib.parse.urlencode({"query": q, "format": "csv"})
        req = urllib.request.Request(url, headers={"User-Agent": "Janus-Cosmos-LOVE-MULTIEYE/2.0"})
        with urllib.request.urlopen(req, timeout=90) as response:
            text = response.read().decode("utf-8")
        out = []
        for r in csv.DictReader(io.StringIO(text)):
            out.append({
                "pl_name": r.get("pl_name"),
                "hostname": r.get("hostname"),
                "ra_deg": float(r["ra"]) if r.get("ra") else None,
                "dec_deg": float(r["dec"]) if r.get("dec") else None,
                "distance_pc": float(r["sy_dist"]) if r.get("sy_dist") else None,
                "spectral_type": r.get("st_spectype") or None,
                "discovery_year": int(float(r["disc_year"])) if r.get("disc_year") else None,
                "discovery_method": r.get("discoverymethod") or None,
                "known_before_search": True,
            })
        return out, None
    except Exception as exc:
        return [], f"{type(exc).__name__}: {exc}"


def _row_text(row, key: str, default: str = "") -> str:
    try:
        value = row[key]
    except Exception:
        return default
    if value is None:
        return default
    try:
        if getattr(value, "mask", False):
            return default
    except Exception:
        pass
    return str(value).strip()


def _fits_contains_coordinate(path: Path, ra_deg: float, dec_deg: float) -> tuple[bool, dict]:
    try:
        with fits.open(path, memmap=False, lazy_load_hdus=True) as hdul:
            for index, hdu in enumerate(hdul):
                data = getattr(hdu, "data", None)
                if data is None or np.ndim(data) < 2:
                    continue
                shape = np.shape(data)[-2:]
                try:
                    w = WCS(hdu.header).celestial
                    if w.pixel_n_dim != 2 or w.world_n_dim != 2:
                        continue
                    x, y = w.world_to_pixel(SkyCoord(ra_deg * u.deg, dec_deg * u.deg, frame="icrs"))
                    inside = bool(np.isfinite(x) and np.isfinite(y) and 0 <= x < shape[1] and 0 <= y < shape[0])
                    if inside:
                        return True, {
                            "hdu": int(index),
                            "pixel_x": float(x),
                            "pixel_y": float(y),
                            "shape_yx": [int(shape[0]), int(shape[1])],
                        }
                except Exception:
                    continue
        return False, {"reason": "NO_CELESTIAL_WCS_CONTAINMENT"}
    except Exception as exc:
        return False, {"error": f"{type(exc).__name__}: {exc}"}


def _safe_filename(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text.strip())
    return text[:180] or "archive.fits"


def probe_hst(candidates: list[dict], prereg: dict, raw_root: Path) -> tuple[list[dict], list[str]]:
    errors: list[str] = []
    try:
        from astroquery.mast import Observations
    except Exception as exc:
        return [], [f"MAST_IMPORT: {type(exc).__name__}: {exc}"]

    cfg = prereg["hst_mast"]
    n = min(len(candidates), int(prereg["gaia"]["archive_probe_sources"]))
    results = []
    total_downloads = 0
    global_cap = 12
    max_bytes = 250_000_000

    for rank, source in enumerate(candidates[:n], start=1):
        coord = SkyCoord(source["ra_deg"] * u.deg, source["dec_deg"] * u.deg, frame="icrs")
        entry = {
            "source_rank": rank,
            "gaia_source_id": source["source_id"],
            "ra_deg": source["ra_deg"],
            "dec_deg": source["dec_deg"],
            "query_radius_arcsec": cfg["query_radius_arcsec"],
            "status": "NO_HST_COVERAGE",
            "observation_count": 0,
            "downloaded": [],
        }
        try:
            obs = Observations.query_region(coord, radius=float(cfg["query_radius_arcsec"]) * u.arcsec)
            if len(obs):
                mask = np.ones(len(obs), dtype=bool)
                if "obs_collection" in obs.colnames:
                    mask &= np.asarray([_row_text(r, "obs_collection").upper() == "HST" for r in obs])
                if "dataproduct_type" in obs.colnames:
                    mask &= np.asarray([_row_text(r, "dataproduct_type").lower() == "image" for r in obs])
                if "dataRights" in obs.colnames:
                    mask &= np.asarray([_row_text(r, "dataRights", "PUBLIC").upper() == "PUBLIC" for r in obs])
                obs = obs[mask]
            entry["observation_count"] = int(len(obs))
            if len(obs) == 0:
                results.append(entry)
                continue

            products = Observations.get_unique_product_list(obs)
            admitted = []
            for row in products:
                filename = _row_text(row, "productFilename")
                uri = _row_text(row, "dataURI")
                ptype = _row_text(row, "productType").upper()
                size = int(float(_row_text(row, "size", "0") or 0))
                if not filename.lower().endswith(tuple(cfg["accepted_product_suffixes"])):
                    continue
                if ptype and ptype != "SCIENCE":
                    continue
                if not uri:
                    continue
                if size > max_bytes:
                    continue
                admitted.append(row)
            admitted.sort(key=product_rank, reverse=True)
            if not admitted:
                entry["status"] = "HST_OBSERVATIONS_WITHOUT_ADMITTED_FITS"
                results.append(entry)
                continue

            per_source = int(cfg["max_downloads_per_source"])
            for row in admitted[:per_source]:
                if total_downloads >= global_cap:
                    break
                filename = _safe_filename(_row_text(row, "productFilename"))
                target_dir = raw_root / "hst" / f"gaia_{source['source_id']}"
                target_dir.mkdir(parents=True, exist_ok=True)
                local = target_dir / filename
                status, message, remote_url = Observations.download_file(
                    _row_text(row, "dataURI"), local_path=str(local), cache=True, verbose=False
                )
                record = {
                    "product_filename": filename,
                    "data_uri": _row_text(row, "dataURI"),
                    "product_subgroup": _row_text(row, "productSubGroupDescription"),
                    "download_status": str(status),
                    "download_message": _jsonable(message),
                    "remote_url": _jsonable(remote_url),
                    "local_path": str(local),
                    "size_bytes": local.stat().st_size if local.exists() else None,
                }
                if str(status).upper() == "COMPLETE" and local.exists():
                    contained, wcs_meta = _fits_contains_coordinate(local, source["ra_deg"], source["dec_deg"])
                    record["wcs_contains_gaia_coordinate"] = contained
                    record["wcs"] = wcs_meta
                    total_downloads += 1
                entry["downloaded"].append(record)
            entry["status"] = "HST_FITS_DOWNLOADED" if entry["downloaded"] else "HST_COVERAGE_NO_DOWNLOAD"
        except Exception as exc:
            entry["status"] = "HST_QUERY_ERROR"
            entry["error"] = f"{type(exc).__name__}: {exc}"
            errors.append(f"Gaia {source['source_id']}: {entry['error']}")
        results.append(entry)
    return results, errors


def _norm_col(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def _pick_coord_col(columns: list[str], axis: str) -> str | None:
    priority_ra = ["ra", "radeg", "raj2000", "ra2000", "objra", "targetra", "alpha", "crval1"]
    priority_dec = ["dec", "decdeg", "dej2000", "decj2000", "dec2000", "objdec", "targetdec", "delta", "crval2"]
    wanted = priority_ra if axis == "ra" else priority_dec
    norm_map = {_norm_col(c): c for c in columns}
    for token in wanted:
        if token in norm_map:
            return norm_map[token]
    for c in columns:
        n = _norm_col(c)
        if axis == "ra" and (n.startswith("ra") or n.endswith("ra")) and "rate" not in n:
            return c
        if axis == "dec" and (n.startswith("dec") or n.endswith("dec")):
            return c
    return None


def _luci_schema_coordinate_binding():
    from experiments.luci.build_luci_archive_manifest import _s, tap_query

    bindings = []
    schema = {}
    for table in ("lbt.luci", "lbt.lbt"):
        t = tap_query(
            "SELECT column_name, datatype FROM TAP_SCHEMA.columns "
            f"WHERE table_name='{table}'"
        )
        cols = [_s(row["column_name"]) for row in t]
        schema[table] = cols
        ra = _pick_coord_col(cols, "ra")
        dec = _pick_coord_col(cols, "dec")
        if ra and dec:
            bindings.append({"table": table, "ra": ra, "dec": dec})
    return bindings, schema


def _download_http(url: str, local: Path, timeout: int = 180) -> None:
    local.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "Janus-Cosmos-LOVE-MULTIEYE/2.0"})
    with urllib.request.urlopen(req, timeout=timeout) as response, local.open("wb") as fh:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            fh.write(chunk)


def probe_luci(candidates: list[dict], prereg: dict, raw_root: Path) -> tuple[list[dict], dict, list[str]]:
    from experiments.luci.build_luci_archive_manifest import _s, tap_query

    errors: list[str] = []
    cfg = prereg["luci_lbt"]
    n = min(len(candidates), int(prereg["gaia"]["archive_probe_sources"]))
    results = []
    try:
        bindings, schema = _luci_schema_coordinate_binding()
    except Exception as exc:
        return [], {"status": "SCHEMA_ERROR", "error": f"{type(exc).__name__}: {exc}"}, []
    schema_report = {
        "status": "COORDINATE_BINDING_FOUND" if bindings else "NO_RA_DEC_BINDING_IN_TAP_SCHEMA",
        "bindings": bindings,
        "luci_columns": schema.get("lbt.luci", []),
        "lbt_columns": schema.get("lbt.lbt", []),
    }
    if not bindings:
        return results, schema_report, errors

    binding = bindings[0]
    table = binding["table"]
    qra = f"{table}.{binding['ra']}"
    qdec = f"{table}.{binding['dec']}"
    radius_deg = float(cfg["query_radius_arcsec"]) / 3600.0
    total_downloads = 0
    global_cap = 12

    for rank, source in enumerate(candidates[:n], start=1):
        entry = {
            "source_rank": rank,
            "gaia_source_id": source["source_id"],
            "ra_deg": source["ra_deg"],
            "dec_deg": source["dec_deg"],
            "query_radius_arcsec": cfg["query_radius_arcsec"],
            "status": "NO_LUCI_COVERAGE",
            "tap_rows": 0,
            "downloaded": [],
        }
        adql = (
            "SELECT TOP 100 "
            "lbt.luci.instrument, lbt.luci.telescope, lbt.luci.object, "
            "lbt.luci.filters, lbt.luci.filter1, lbt.luci.filter2, "
            "lbt.luci.gratname, lbt.luci.imagetyp, lbt.luci.file_name, "
            "lbt.luci.date_obs, lbt.lbt.file_url, lbt.lbt.policy, "
            f"{qra} AS archive_ra, {qdec} AS archive_dec "
            "FROM lbt.luci JOIN lbt.lbt ON lbt.luci.file_name=lbt.lbt.file_name "
            "WHERE lbt.luci.imagetyp='SCIENCE' "
            "AND (lbt.luci.gratname='Mirror' OR lbt.luci.gratname='mirror') "
            "AND lbt.lbt.policy='FREE' "
            f"AND 1=CONTAINS(POINT('ICRS',{qra},{qdec}),"
            f"CIRCLE('ICRS',{source['ra_deg']},{source['dec_deg']},{radius_deg}))"
        )
        try:
            rows = tap_query(adql, timeout=120)
            entry["tap_rows"] = int(len(rows))
            if len(rows) == 0:
                results.append(entry)
                continue
            selected = []
            seen_filters = set()
            for row in rows:
                url = _s(row["file_url"])
                if not url.startswith(("http://", "https://")):
                    continue
                f2 = _s(row["filter2"])
                f1 = _s(row["filter1"])
                filt = f2 if f2.lower() not in {"", "clear", "blind", "open", "none"} else f1
                key = filt.casefold()
                if key in seen_filters:
                    continue
                seen_filters.add(key)
                selected.append((row, filt))
                if len(selected) >= int(cfg["max_downloads_per_source"]):
                    break
            for row, filt in selected:
                if total_downloads >= global_cap:
                    break
                url = _s(row["file_url"])
                filename = _safe_filename(_s(row["file_name"]) or Path(urllib.parse.urlparse(url).path).name)
                target_dir = raw_root / "luci" / f"gaia_{source['source_id']}"
                local = target_dir / filename
                record = {
                    "archive_file_name": filename,
                    "archive_object": _s(row["object"]),
                    "archive_filter": filt,
                    "archive_date_obs": _s(row["date_obs"]),
                    "instrument": _s(row["instrument"]),
                    "url": url,
                    "local_path": str(local),
                }
                try:
                    _download_http(url, local)
                    image, provenance = read_luci_fits_image(local, require_imaging=True)
                    contained, wcs_meta = _fits_contains_coordinate(local, source["ra_deg"], source["dec_deg"])
                    record.update({
                        "download_status": "COMPLETE",
                        "size_bytes": local.stat().st_size,
                        "luci_provenance": provenance,
                        "image_shape": [int(x) for x in image.shape],
                        "wcs_contains_gaia_coordinate": contained,
                        "wcs": wcs_meta,
                    })
                    total_downloads += 1
                except Exception as exc:
                    record["download_status"] = "ERROR"
                    record["error"] = f"{type(exc).__name__}: {exc}"
                entry["downloaded"].append(record)
            entry["status"] = "LUCI_FITS_DOWNLOADED" if entry["downloaded"] else "LUCI_COVERAGE_NO_DOWNLOAD"
        except Exception as exc:
            entry["status"] = "LUCI_QUERY_ERROR"
            entry["error"] = f"{type(exc).__name__}: {exc}"
            errors.append(f"Gaia {source['source_id']}: {entry['error']}")
        results.append(entry)
    return results, schema_report, errors


def probe_tess(candidates: list[dict], prereg: dict) -> tuple[list[dict], list[str]]:
    errors: list[str] = []
    try:
        from astroquery.mast import Catalogs, Observations
    except Exception as exc:
        return [], [f"MAST_IMPORT: {type(exc).__name__}: {exc}"]
    n = min(len(candidates), int(prereg["tess"]["probe_top_sources"]))
    out = []
    for rank, source in enumerate(candidates[:n], start=1):
        coord = SkyCoord(source["ra_deg"] * u.deg, source["dec_deg"] * u.deg, frame="icrs")
        entry = {"source_rank": rank, "gaia_source_id": source["source_id"], "status": "NO_TIC_MATCH"}
        try:
            tic = Catalogs.query_region(coord, radius=5 * u.arcsec, catalog="TIC")
            if len(tic) == 0:
                out.append(entry)
                continue
            if "dstArcSec" in tic.colnames:
                tic.sort("dstArcSec")
            row = tic[0]
            tic_id = str(row["ID"])
            obs = Observations.query_criteria(obs_collection="TESS", target_name=tic_id, dataproduct_type="timeseries")
            entry.update({
                "status": "TESS_AVAILABLE" if len(obs) else "TIC_MATCH_NO_TARGETED_TIMESERIES",
                "tic_id": tic_id,
                "tmag": _jsonable(row["Tmag"]) if "Tmag" in tic.colnames else None,
                "tess_observation_count": int(len(obs)),
            })
        except Exception as exc:
            entry["status"] = "TESS_QUERY_ERROR"
            entry["error"] = f"{type(exc).__name__}: {exc}"
            errors.append(f"Gaia {source['source_id']}: {entry['error']}")
        out.append(entry)
    return out, errors


def run(prereg_path: Path, output_dir: Path) -> dict:
    prereg = json.loads(prereg_path.read_text(encoding="utf-8"))
    raw_root = output_dir / "raw"
    raw_root.mkdir(parents=True, exist_ok=True)

    gaia, gaia_error = query_gaia(prereg)
    known_planets, nasa_error = query_known_planets(prereg)
    hst, hst_errors = probe_hst(gaia, prereg, raw_root) if gaia else ([], [])
    luci, luci_schema, luci_errors = probe_luci(gaia, prereg, raw_root) if gaia else ([], {"status": "SKIPPED_NO_GAIA"}, [])
    tess, tess_errors = probe_tess(gaia, prereg) if gaia else ([], [])

    hst_contained = sum(
        1 for e in hst for d in e.get("downloaded", []) if d.get("wcs_contains_gaia_coordinate") is True
    )
    luci_contained = sum(
        1 for e in luci for d in e.get("downloaded", []) if d.get("wcs_contains_gaia_coordinate") is True
    )
    tess_available = sum(1 for e in tess if e.get("status") == "TESS_AVAILABLE")
    k_like = sum(1 for x in gaia if x.get("k_like"))

    result = {
        "schema": "janus.cosmos.love.multieye.result.v2",
        "experiment_id": prereg["experiment_id"],
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_label": ANALYSIS_LABEL,
        "target": prereg["origin"],
        "blind_identity_firewall": prereg["blind_identity_firewall"],
        "gaia": {
            "status": "OK" if gaia else ("ERROR" if gaia_error else "NO_MATCHES"),
            "error": gaia_error,
            "ranked_source_count": len(gaia),
            "k_like_count": k_like,
            "top_sources": gaia,
        },
        "known_planet_firewall": {
            "status": "OK" if nasa_error is None else "ERROR",
            "error": nasa_error,
            "known_planet_count": len(known_planets),
            "known_planets": known_planets,
        },
        "hst_mast": {
            "probed_source_count": len(hst),
            "wcs_contained_download_count": hst_contained,
            "sources": hst,
            "errors": hst_errors,
        },
        "luci_lbt": {
            "schema_coordinate_binding": luci_schema,
            "probed_source_count": len(luci),
            "wcs_contained_download_count": luci_contained,
            "sources": luci,
            "errors": luci_errors,
        },
        "tess": {
            "probed_source_count": len(tess),
            "timeseries_available_count": tess_available,
            "sources": tess,
            "errors": tess_errors,
        },
        "love_gate": {
            "reserved_codename": prereg["blind_identity_firewall"]["reserved_post_gate_codename"],
            "candidate_activated": False,
            "reason": "Archive discovery and real FITS acquisition do not by themselves establish a new planet. The codename remains sealed until a new repeatable planetary signal survives controls and independent replication.",
            "next_gate": "RUN_SIGNAL_SEARCH_ON_COORDINATE_BOUND_HST_LUCI_TESS_DATA_AND_REQUIRE_INDEPENDENT_REPLICATION",
        },
        "claim_ceiling": prereg["claim_ceiling"],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "love-multieye-result.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description="Blind coordinate-bound Gaia + HST/MAST + LUCI/LBT + TESS search around the mapped Lave point.")
    ap.add_argument("--prereg", default="data/love/LAVE_TO_LOVE_MULTIEYE_PREREG.json")
    ap.add_argument("--output-dir", default="results/love_multieye")
    args = ap.parse_args()
    result = run(Path(args.prereg), Path(args.output_dir))
    print(json.dumps({
        "experiment_id": result["experiment_id"],
        "analysis_label": result["analysis_label"],
        "gaia_ranked": result["gaia"]["ranked_source_count"],
        "gaia_k_like": result["gaia"]["k_like_count"],
        "known_planets": result["known_planet_firewall"]["known_planet_count"],
        "hst_wcs_downloads": result["hst_mast"]["wcs_contained_download_count"],
        "luci_wcs_downloads": result["luci_lbt"]["wcs_contained_download_count"],
        "tess_available": result["tess"]["timeseries_available_count"],
        "love_candidate_activated": result["love_gate"]["candidate_activated"],
    }, indent=2))
    return 0 if result["gaia"]["status"] != "ERROR" else 2


if __name__ == "__main__":
    raise SystemExit(main())
