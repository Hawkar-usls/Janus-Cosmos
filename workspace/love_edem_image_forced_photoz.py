from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astroquery.skyview import SkyView
from astroquery.vizier import Vizier

OUT = Path("data/love/IMAGE_LEVEL_FORCED_PHOTOMETRY_AND_EDEM_GALAXY_PHOTOZ-v1-LATEST-RECEIPT.json")

TARGETS = {
    "EDEM_NEAREST_WISE_ONLY": {"ra": 139.2249571, "dec": 30.262781, "anchor": "AllWISE J091653.98+301546.0"},
    "LOVE_NEAREST_WISE_ONLY": {"ra": 204.2960668, "dec": -36.7804774, "anchor": "AllWISE J133711.05-364649.7"},
}
GALAXY = SkyCoord(139.222681 * u.deg, 30.256918 * u.deg, frame="icrs")
SURVEYS = {
    "WISE 3.4": 9.0,
    "WISE 4.6": 9.0,
    "2MASS-J": 4.0,
    "DSS2 Red": 4.0,
}


def _finite_float(v: Any) -> float | None:
    try:
        if getattr(v, "mask", False):
            return None
    except Exception:
        pass
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except Exception:
        return None


def _aperture_flux(data: np.ndarray, x: float, y: float, r_pix: float) -> tuple[float | None, dict[str, Any]]:
    yy, xx = np.indices(data.shape)
    rr = np.hypot(xx - x, yy - y)
    ap = (rr <= r_pix) & np.isfinite(data)
    ann = (rr >= 2.0 * r_pix) & (rr <= 3.0 * r_pix) & np.isfinite(data)
    if ap.sum() < 4 or ann.sum() < 20:
        return None, {"aperture_pixels": int(ap.sum()), "annulus_pixels": int(ann.sum())}
    bg = float(np.nanmedian(data[ann]))
    flux = float(np.nansum(data[ap] - bg))
    return flux, {"aperture_pixels": int(ap.sum()), "annulus_pixels": int(ann.sum()), "background_median": bg}


def _forced_measurement(hdu, coord: SkyCoord, aperture_arcsec: float) -> dict[str, Any]:
    data = np.asarray(hdu.data, dtype=float)
    while data.ndim > 2:
        data = data[0]
    w = WCS(hdu.header).celestial
    x, y = w.world_to_pixel(coord)
    scales = proj_plane_pixel_scales(w) * 3600.0
    pixscale = float(np.sqrt(abs(scales[0] * scales[1])))
    r_pix = aperture_arcsec / pixscale
    source_flux, meta = _aperture_flux(data, x, y, r_pix)
    controls = []
    for k in range(24):
        pa = (360.0 * k / 24.0) * u.deg
        c = coord.directional_offset_by(pa, 45.0 * u.arcsec)
        cx, cy = w.world_to_pixel(c)
        if cx < 3*r_pix or cy < 3*r_pix or cx > data.shape[1]-3*r_pix or cy > data.shape[0]-3*r_pix:
            continue
        f, _ = _aperture_flux(data, cx, cy, r_pix)
        if f is not None and math.isfinite(f):
            controls.append(float(f))
    z = None
    robust_sigma = None
    median_control = None
    if source_flux is not None and len(controls) >= 8:
        median_control = float(np.median(controls))
        mad = float(np.median(np.abs(np.asarray(controls) - median_control)))
        robust_sigma = 1.4826 * mad
        if robust_sigma > 0:
            z = (source_flux - median_control) / robust_sigma
    digest = hashlib.sha256(np.nan_to_num(data, nan=0.0).astype(np.float32).tobytes()).hexdigest()
    return {
        "pixel_scale_arcsec": pixscale,
        "aperture_radius_arcsec": aperture_arcsec,
        "source_pixel": [float(x), float(y)],
        "source_flux_native_units_background_subtracted": source_flux,
        "control_count": len(controls),
        "control_median_flux": median_control,
        "control_robust_sigma": robust_sigma,
        "forced_control_z": z,
        "passes_frozen_5sigma_gate": bool(z is not None and z >= 5.0),
        "image_sha256_float32_nan0": digest,
        "bunit": hdu.header.get("BUNIT"),
        **meta,
    }


def query_images(name: str, info: dict[str, Any]) -> dict[str, Any]:
    coord = SkyCoord(info["ra"] * u.deg, info["dec"] * u.deg, frame="icrs")
    out = {"anchor": info["anchor"], "icrs": {"ra_deg": info["ra"], "dec_deg": info["dec"]}, "surveys": {}}
    for survey, aperture in SURVEYS.items():
        item: dict[str, Any] = {"survey": survey, "status": "UNKNOWN"}
        try:
            images = SkyView.get_images(
                position=coord,
                survey=[survey],
                radius=1.5 * u.arcmin,
                pixels=256,
                cache=False,
                show_progress=False,
            )
            if not images:
                item["status"] = "NO_IMAGE"
            else:
                item["status"] = "OK"
                item["measurement"] = _forced_measurement(images[0][0], coord, aperture)
        except Exception as exc:
            item["status"] = "QUERY_ERROR"
            item["error"] = f"{type(exc).__name__}: {exc}"
        out["surveys"][survey] = item
    passes = [s for s, item in out["surveys"].items() if item.get("measurement", {}).get("passes_frozen_5sigma_gate")]
    out["summary"] = {
        "surveys_passing_forced_5sigma": passes,
        "wise_pixel_persistence": any(s.startswith("WISE") for s in passes),
        "independent_nonwise_pixel_counterpart": any(not s.startswith("WISE") for s in passes),
    }
    return out


def _find_col(names: list[str], candidates: list[str]) -> str | None:
    low = {n.lower(): n for n in names}
    for c in candidates:
        if c.lower() in low:
            return low[c.lower()]
    return None


def query_photoz(catalog: str) -> dict[str, Any]:
    out: dict[str, Any] = {"catalog": catalog, "status": "UNKNOWN", "matches": []}
    try:
        viz = Vizier(columns=["*"], row_limit=50)
        tabs = viz.query_region(GALAXY, radius=2.0 * u.arcsec, catalog=catalog)
        rows_out = []
        for table in tabs:
            names = list(table.colnames)
            ra_col = _find_col(names, ["RAdeg", "RA_ICRS", "RAJ2000", "RA"])
            de_col = _find_col(names, ["DEdeg", "DE_ICRS", "DEJ2000", "DEC", "DE"])
            z_col = _find_col(names, ["zphot", "zph", "zph1", "photoz", "z_phot"])
            ez_col = _find_col(names, ["e_zphot", "e_zph", "e_zph1", "photoz_err", "zphot_err"])
            fclean_col = _find_col(names, ["fclean", "flag_clean"])
            fqual_col = _find_col(names, ["fqual", "flag_qual"])
            type_col = _find_col(names, ["type", "class"])
            pstar_col = _find_col(names, ["pstar"])
            if not (ra_col and de_col):
                continue
            for row in table:
                ra = _finite_float(row[ra_col]); de = _finite_float(row[de_col])
                if ra is None or de is None:
                    continue
                c = SkyCoord(ra*u.deg, de*u.deg, frame="icrs")
                sep = float(GALAXY.separation(c).to_value(u.arcsec))
                z = _finite_float(row[z_col]) if z_col else None
                ez = _finite_float(row[ez_col]) if ez_col else None
                fclean = _finite_float(row[fclean_col]) if fclean_col else None
                fqual = _finite_float(row[fqual_col]) if fqual_col else None
                q_ok = True
                if fclean is not None:
                    q_ok = q_ok and int(fclean) == 1
                if fqual is not None:
                    q_ok = q_ok and int(fqual) == 1
                accepted = bool(sep <= 1.0 and z is not None and 0 < z < 7 and q_ok)
                rows_out.append({
                    "table_id": str(getattr(table, "meta", {}).get("name", "")),
                    "ra_deg": ra, "dec_deg": de, "separation_arcsec": sep,
                    "zphot": z, "e_zphot": ez,
                    "fclean": fclean, "fqual": fqual,
                    "type": str(row[type_col]) if type_col else None,
                    "pstar": _finite_float(row[pstar_col]) if pstar_col else None,
                    "accepted_high_quality": accepted,
                    "available_columns": {"z": z_col, "e_z": ez_col, "fclean": fclean_col, "fqual": fqual_col},
                })
        out["matches"] = sorted(rows_out, key=lambda r: r["separation_arcsec"])
        out["status"] = "OK"
    except Exception as exc:
        out["status"] = "QUERY_ERROR"
        out["error"] = f"{type(exc).__name__}: {exc}"
    return out


def main() -> None:
    image_results = {name: query_images(name, info) for name, info in TARGETS.items()}
    photoz_queries = {cat: query_photoz(cat) for cat in ["VII/292", "V/147/sdss12"]}
    accepted = []
    for cat, q in photoz_queries.items():
        for row in q.get("matches", []):
            if row.get("accepted_high_quality"):
                accepted.append({"catalog": cat, **row})
    accepted.sort(key=lambda r: (r["separation_arcsec"], r.get("e_zphot") if r.get("e_zphot") is not None else 999))
    best = accepted[0] if accepted else None
    distance = None
    if best and best.get("zphot") is not None:
        z = best["zphot"]
        distance = {
            "zphot": z,
            "luminosity_distance_mpc": float(Planck18.luminosity_distance(z).to_value(u.Mpc)),
            "comoving_distance_mpc": float(Planck18.comoving_distance(z).to_value(u.Mpc)),
            "lookback_time_gyr": float(Planck18.lookback_time(z).to_value(u.Gyr)),
            "warning": "Distances are conditional on a photometric redshift and Planck18 cosmology; not a spectroscopic distance measurement."
        }

    payload = {
        "schema": "janus.cosmos.love_edem.image_level_forced_photometry_photoz.receipt.v1",
        "experiment_id": "IMAGE_LEVEL_FORCED_PHOTOMETRY_AND_EDEM_GALAXY_PHOTOZ",
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETE",
        "mode": "IMAGE_LEVEL_FORCED_PHOTOMETRY_AND_PHOTOZ_WITHOUT_ANOMALY_SCORING",
        "anomaly_scoring_used": False,
        "image_level": image_results,
        "edem_galaxy_photoz": {
            "icrs": {"ra_deg": float(GALAXY.ra.deg), "dec_deg": float(GALAXY.dec.deg)},
            "queries": photoz_queries,
            "accepted_high_quality_matches": accepted,
            "best_photoz": best,
            "conditional_cosmology": distance,
        },
        "frozen_rules": {
            "forced_detection_threshold_sigma": 5.0,
            "control_apertures": 24,
            "control_ring_arcsec": 45.0,
            "recenter_for_detection": False,
            "photoz_match_radius_arcsec": 2.0,
            "accepted_photoz_match_max_arcsec": 1.0,
            "no_post_hoc_radius_or_threshold_changes": True,
        },
        "firewall": {
            "WISE_image_recovery_is_independent_of_AllWISE": False,
            "pixel_excess_is_anomaly": False,
            "photoz_is_spectroscopic_redshift": False,
            "near_source_is_semantic_LOVE_or_EDEM_identity": False,
            "physical_LOVE_EDEM_ORION_association_established": False,
            "planet_claim": False,
            "claim_ceiling": "PIXEL_LEVEL_PERSISTENCE_AND_CATALOG_PHOTOZ_ONLY",
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "image_summary": {k: v["summary"] for k, v in image_results.items()},
        "accepted_photoz_count": len(accepted),
        "best_photoz": best,
        "output": str(OUT),
    }, ensure_ascii=False, indent=2))

if __name__ == "__main__":
    main()
