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

import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares


def _safe_float(value, default=None):
    try:
        v = float(value)
        return v if math.isfinite(v) else default
    except Exception:
        return default


def _safe_name(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in str(text))[:180]


def _robust_sde(power: np.ndarray) -> tuple[float, float, float]:
    p = np.asarray(power, dtype=float)
    p = p[np.isfinite(p)]
    if p.size < 10:
        return float("nan"), float("nan"), float("nan")
    med = float(np.median(p))
    mad = float(np.median(np.abs(p - med)))
    scale = 1.4826 * mad
    if not math.isfinite(scale) or scale <= 0:
        scale = float(np.std(p))
    if not math.isfinite(scale) or scale <= 0:
        return float("nan"), med, scale
    return float((np.max(p) - med) / scale), med, scale


def _estimate_transits(time: np.ndarray, period: float, duration: float, t0: float) -> int:
    if time.size == 0 or period <= 0:
        return 0
    phase_index = np.round((time - t0) / period).astype(int)
    half = duration / 2.0
    centers = t0 + phase_index * period
    in_transit = np.abs(time - centers) <= half
    if not np.any(in_transit):
        return 0
    return int(len(np.unique(phase_index[in_transit])))


def _read_tess_lightcurve(path: Path) -> dict:
    with fits.open(path, memmap=False) as hdul:
        hdr0 = hdul[0].header
        if len(hdul) < 2 or hdul[1].data is None:
            raise RuntimeError("No TESS lightcurve table extension")
        data = hdul[1].data
        names = {str(n).upper(): n for n in data.names}
        if "TIME" not in names:
            raise RuntimeError("TIME column missing")
        flux_col = "PDCSAP_FLUX" if "PDCSAP_FLUX" in names else ("SAP_FLUX" if "SAP_FLUX" in names else None)
        if flux_col is None:
            raise RuntimeError("PDCSAP_FLUX/SAP_FLUX missing")
        err_col = "PDCSAP_FLUX_ERR" if flux_col == "PDCSAP_FLUX" and "PDCSAP_FLUX_ERR" in names else (
            "SAP_FLUX_ERR" if "SAP_FLUX_ERR" in names else None
        )
        time = np.asarray(data[names["TIME"]], dtype=float)
        flux = np.asarray(data[names[flux_col]], dtype=float)
        err = np.asarray(data[names[err_col]], dtype=float) if err_col else np.full_like(flux, np.nan)
        quality = np.asarray(data[names["QUALITY"]], dtype=int) if "QUALITY" in names else np.zeros_like(time, dtype=int)
        good = np.isfinite(time) & np.isfinite(flux) & (quality == 0) & (flux > 0)
        time = time[good]
        flux = flux[good]
        err = err[good]
        if time.size < 200:
            raise RuntimeError(f"Too few good cadences: {time.size}")
        med = float(np.nanmedian(flux))
        if not math.isfinite(med) or med <= 0:
            raise RuntimeError("Invalid flux median")
        flux = flux / med
        err = err / med
        finite_err = np.isfinite(err) & (err > 0)
        if np.any(finite_err):
            fallback = float(np.nanmedian(err[finite_err]))
            err[~finite_err] = fallback
        else:
            scatter = 1.4826 * float(np.nanmedian(np.abs(flux - np.nanmedian(flux))))
            err[:] = max(scatter, 1e-5)
        return {
            "sector": int(hdr0.get("SECTOR", -1)),
            "camera": int(hdr0.get("CAMERA", -1)),
            "ccd": int(hdr0.get("CCD", -1)),
            "tic_id_header": str(hdr0.get("TICID", hdr0.get("TIC_ID", ""))),
            "ra_obj": _safe_float(hdr0.get("RA_OBJ", hdr0.get("RA_TARG"))),
            "dec_obj": _safe_float(hdr0.get("DEC_OBJ", hdr0.get("DEC_TARG"))),
            "flux_column": flux_col,
            "cadence_count": int(time.size),
            "time": time,
            "flux": flux,
            "err": err,
        }


def _sector_bls(lightcurve: dict, cfg: dict) -> dict:
    time = lightcurve["time"]
    flux = lightcurve["flux"]
    err = lightcurve["err"]
    baseline = float(np.max(time) - np.min(time))
    min_period = float(cfg["minimum_period_days"])
    max_period = min(float(cfg["maximum_period_days"]), baseline / 2.0)
    if max_period <= min_period:
        return {"status": "BASELINE_TOO_SHORT", "baseline_days": baseline}
    durations = np.asarray(cfg["duration_hours"], dtype=float) / 24.0
    durations = durations[durations < min_period]
    model = BoxLeastSquares(time, flux, dy=err)
    result = model.autopower(
        durations,
        objective="snr",
        minimum_n_transit=2,
        minimum_period=min_period,
        maximum_period=max_period,
        frequency_factor=float(cfg["frequency_factor"]),
    )
    power = np.asarray(result.power, dtype=float)
    idx = int(np.nanargmax(power))
    period = float(result.period[idx])
    duration = float(result.duration[idx])
    transit_time = float(result.transit_time[idx])
    depth = float(result.depth[idx])
    depth_snr = float(result.depth_snr[idx]) if hasattr(result, "depth_snr") else None
    sde, median_power, power_scale = _robust_sde(power)
    n_transits = _estimate_transits(time, period, duration, transit_time)
    return {
        "status": "OK",
        "baseline_days": baseline,
        "best_period_days": period,
        "best_duration_hours": duration * 24.0,
        "best_transit_time_btjd": transit_time,
        "best_depth_fraction": depth,
        "best_depth_snr": depth_snr,
        "power_sde": sde,
        "power_median": median_power,
        "power_robust_scale": power_scale,
        "estimated_transits": n_transits,
        "screen_pass": bool(math.isfinite(sde) and sde >= float(cfg["sector_screen_sde_min"])),
    }


def _repeat_clusters(sectors: list[dict], cfg: dict) -> list[dict]:
    tol = float(cfg["repeat_period_relative_tolerance"])
    eligible = [s for s in sectors if s.get("bls", {}).get("status") == "OK" and s["bls"].get("screen_pass")]
    clusters = []
    used = set()
    for i, a in enumerate(eligible):
        if i in used:
            continue
        pa = float(a["bls"]["best_period_days"])
        members = [a]
        for j, b in enumerate(eligible):
            if j == i:
                continue
            pb = float(b["bls"]["best_period_days"])
            rel = abs(pa - pb) / ((pa + pb) / 2.0)
            if rel <= tol:
                members.append(b)
                used.add(j)
        if len(members) >= 2:
            periods = [float(m["bls"]["best_period_days"]) for m in members]
            clusters.append({
                "period_days_median": float(np.median(periods)),
                "matching_sector_count": len(members),
                "sectors": [int(m["sector"]) for m in members],
                "periods_days": periods,
                "max_relative_spread": float((max(periods) - min(periods)) / np.mean(periods)),
            })
        used.add(i)
    clusters.sort(key=lambda c: (-c["matching_sector_count"], c["max_relative_spread"]))
    return clusters


def _combined_local_bls(lightcurves: list[dict], center_period: float, cfg: dict) -> dict:
    if not lightcurves:
        return {"status": "NO_DATA"}
    times = []
    fluxes = []
    errs = []
    for lc in lightcurves:
        times.append(lc["time"])
        fluxes.append(lc["flux"])
        errs.append(lc["err"])
    time = np.concatenate(times)
    flux = np.concatenate(fluxes)
    err = np.concatenate(errs)
    order = np.argsort(time)
    time, flux, err = time[order], flux[order], err[order]
    periods = np.linspace(center_period * 0.985, center_period * 1.015, 4000)
    durations = np.asarray(cfg["duration_hours"], dtype=float) / 24.0
    durations = durations[durations < np.min(periods)]
    model = BoxLeastSquares(time, flux, dy=err)
    result = model.power(periods, durations, objective="snr")
    power = np.asarray(result.power, dtype=float)
    idx = int(np.nanargmax(power))
    period = float(result.period[idx])
    duration = float(result.duration[idx])
    transit_time = float(result.transit_time[idx])
    depth = float(result.depth[idx])
    depth_snr = float(result.depth_snr[idx]) if hasattr(result, "depth_snr") else None
    sde, median_power, power_scale = _robust_sde(power)
    n_transits = _estimate_transits(time, period, duration, transit_time)
    try:
        stats = model.compute_stats(period, duration, transit_time)
        odd = float(stats["depth_odd"][0]) if "depth_odd" in stats else None
        even = float(stats["depth_even"][0]) if "depth_even" in stats else None
    except Exception:
        odd, even = None, None
    return {
        "status": "OK",
        "best_period_days": period,
        "best_duration_hours": duration * 24.0,
        "best_transit_time_btjd": transit_time,
        "best_depth_fraction": depth,
        "best_depth_snr": depth_snr,
        "power_sde": sde,
        "power_median": median_power,
        "power_robust_scale": power_scale,
        "estimated_transits": n_transits,
        "odd_depth_fraction": odd,
        "even_depth_fraction": even,
        "cadence_count": int(time.size),
        "time_baseline_days": float(np.max(time) - np.min(time)),
    }


def _known_planets_at(ra_deg: float | None, dec_deg: float | None) -> tuple[list[str], str | None]:
    if ra_deg is None or dec_deg is None:
        return [], "NO_COORDINATES"
    try:
        radius_deg = 30.0 / 3600.0
        query = (
            "select pl_name,hostname,ra,dec from pscomppars where "
            f"contains(point('ICRS',ra,dec),circle('ICRS',{ra_deg},{dec_deg},{radius_deg}))=1"
        )
        url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urllib.parse.urlencode({"query": query, "format": "csv"})
        req = urllib.request.Request(url, headers={"User-Agent": "Janus-Cosmos-LOVE-TESS-BLS/1.0"})
        with urllib.request.urlopen(req, timeout=60) as response:
            text = response.read().decode("utf-8")
        names = [r.get("pl_name") for r in csv.DictReader(io.StringIO(text)) if r.get("pl_name")]
        return sorted(set(names)), None
    except Exception as exc:
        return [], f"{type(exc).__name__}: {exc}"


def _download_tess_products(tic_id: str, raw_dir: Path) -> tuple[list[Path], list[dict], list[str]]:
    from astroquery.mast import Observations

    errors = []
    records = []
    paths = []
    obs = Observations.query_criteria(obs_collection="TESS", target_name=str(tic_id), dataproduct_type="timeseries")
    if len(obs) == 0:
        return paths, records, errors
    products = Observations.get_product_list(obs)
    candidates = []
    for row in products:
        try:
            filename = str(row["productFilename"]).strip()
        except Exception:
            continue
        lower = filename.lower()
        try:
            ptype = str(row["productType"]).strip().upper()
        except Exception:
            ptype = ""
        try:
            subgroup = str(row["productSubGroupDescription"]).strip().upper()
        except Exception:
            subgroup = ""
        try:
            uri = str(row["dataURI"]).strip()
        except Exception:
            uri = ""
        is_lc = lower.endswith("lc.fits") or lower.endswith("lc.fits.gz") or subgroup == "LC"
        if not is_lc or (ptype and ptype != "SCIENCE") or not uri:
            continue
        candidates.append(("fast-lc" in lower, filename, uri))
    candidates.sort(key=lambda x: (x[0], x[1]))
    raw_dir.mkdir(parents=True, exist_ok=True)
    for is_fast, filename, uri in candidates:
        local = raw_dir / _safe_name(filename)
        try:
            status, msg, remote = Observations.download_file(uri, local_path=str(local), cache=True, verbose=False)
            rec = {
                "filename": filename,
                "data_uri": uri,
                "fast_cadence_product": bool(is_fast),
                "download_status": str(status),
                "message": str(msg) if msg else None,
                "remote_url": str(remote) if remote else None,
                "local_path": str(local),
                "size_bytes": local.stat().st_size if local.exists() else None,
            }
            records.append(rec)
            if str(status).upper() == "COMPLETE" and local.exists():
                paths.append(local)
        except Exception as exc:
            errors.append(f"{filename}: {type(exc).__name__}: {exc}")
    return paths, records, errors


def run(prereg_path: Path, output_dir: Path) -> dict:
    prereg = json.loads(prereg_path.read_text(encoding="utf-8"))
    cfg = prereg["bls_search"]
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_root = output_dir / "raw_tess"
    source_results = []

    for src in prereg["frozen_sources"]:
        tic_id = str(src["tic_id"])
        label = f"TARGET_A_R{src['rank']}"
        paths, downloads, download_errors = _download_tess_products(tic_id, raw_root / f"tic_{tic_id}")
        sector_map = {}
        parse_errors = []
        for path in paths:
            try:
                lc = _read_tess_lightcurve(path)
                sector = int(lc["sector"])
                current = sector_map.get(sector)
                prefer_new = current is None or ("fast-lc" in current["path"].name.lower() and "fast-lc" not in path.name.lower())
                if prefer_new:
                    sector_map[sector] = {"path": path, "lc": lc}
            except Exception as exc:
                parse_errors.append(f"{path.name}: {type(exc).__name__}: {exc}")

        sectors = []
        lightcurves = []
        ra_obj = dec_obj = None
        for sector in sorted(sector_map):
            item = sector_map[sector]
            lc = item["lc"]
            if ra_obj is None and lc.get("ra_obj") is not None:
                ra_obj, dec_obj = lc.get("ra_obj"), lc.get("dec_obj")
            bls = _sector_bls(lc, cfg)
            sectors.append({
                "sector": sector,
                "camera": lc["camera"],
                "ccd": lc["ccd"],
                "filename": item["path"].name,
                "flux_column": lc["flux_column"],
                "cadence_count": lc["cadence_count"],
                "bls": bls,
            })
            lightcurves.append(lc)

        known_planets, known_error = _known_planets_at(ra_obj, dec_obj)
        clusters = _repeat_clusters(sectors, cfg)
        combined = {"status": "NO_REPEAT_CLUSTER"}
        if clusters:
            combined = _combined_local_bls(lightcurves, clusters[0]["period_days_median"], cfg)

        repeat_count = clusters[0]["matching_sector_count"] if clusters else 0
        depth = combined.get("best_depth_fraction")
        sde = combined.get("power_sde")
        n_transits = combined.get("estimated_transits", 0)
        watchlist = bool(
            not known_planets
            and combined.get("status") == "OK"
            and sde is not None and math.isfinite(float(sde)) and float(sde) >= float(cfg["combined_sde_min"])
            and depth is not None and float(cfg["depth_fraction_min"]) <= float(depth) <= float(cfg["depth_fraction_max"])
            and int(n_transits) >= int(cfg["minimum_estimated_transits"])
            and int(repeat_count) >= int(cfg["minimum_matching_sectors"])
        )
        source_results.append({
            "blind_label": label,
            "rank": src["rank"],
            "gaia_source_id": src["gaia_source_id"],
            "tic_id": tic_id,
            "k_like": bool(src["k_like"]),
            "coordinates_from_tess_header": {"ra_deg": ra_obj, "dec_deg": dec_obj},
            "downloaded_product_count": len(paths),
            "download_records": downloads,
            "download_errors": download_errors,
            "parse_errors": parse_errors,
            "distinct_sector_count": len(sectors),
            "sectors": sectors,
            "repeat_clusters": clusters,
            "combined_repeat_fit": combined,
            "known_planet_control": {"matched_planets": known_planets, "query_error": known_error},
            "tess_repeat_watchlist": watchlist,
            "semantic_status": "UNCONFIRMED",
        })

    watch = [s for s in source_results if s["tess_repeat_watchlist"]]
    result = {
        "schema": "janus.cosmos.love.tess_bls.result.v1",
        "experiment_id": prereg["experiment_id"],
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "parent_receipt": prereg["parent_receipt"],
        "blind_identity_firewall": prereg["blind_identity_firewall"],
        "frozen_source_count": len(prereg["frozen_sources"]),
        "sources": source_results,
        "watchlist_count": len(watch),
        "watchlist_blind_labels": [s["blind_label"] for s in watch],
        "love_gate": {
            "reserved_codename": prereg["blind_identity_firewall"]["reserved_post_gate_codename"],
            "candidate_activated": False,
            "reason": prereg["admission"]["reason"],
            "next_gate": "INDEPENDENT_NON_TESS_REPLICATION_AND_FALSE_POSITIVE_VETTING" if watch else "NO_TESS_REPEAT_WATCHLIST__EXPAND_ARCHIVE_OR_OBSERVATIONAL_BASELINE",
        },
        "claim_ceiling": prereg["claim_ceiling"],
    }
    (output_dir / "love-tess-bls-result.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "experiment_id": result["experiment_id"],
        "sources": result["frozen_source_count"],
        "watchlist_count": result["watchlist_count"],
        "watchlist_labels": result["watchlist_blind_labels"],
        "love_candidate_activated": result["love_gate"]["candidate_activated"],
    }, indent=2))
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description="Blind TESS BLS repeatability screen for frozen LOVE target counterparts")
    ap.add_argument("--prereg", default="data/love/LAVE_TO_LOVE_TESS_BLS_PREREG.json")
    ap.add_argument("--output-dir", default="results/love_tess_bls")
    args = ap.parse_args()
    run(Path(args.prereg), Path(args.output_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
