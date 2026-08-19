from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from janus_cosmos import love_edem_center_object_probe as primary_probe
from janus_cosmos import love_edem_open_sky_scan as bypass_scan


TARGET_MAP = {
    "TARGET_A": "LOVE",
    "TARGET_B": "EDEM_ZAPORIZHZHIA_DIRECTION_CANDIDATE",
}

CATALOG_MAP = {
    "gaia_dr3": ("gaia_dr3", "gaia_dr3"),
    "simbad": ("simbad", "simbad"),
    "allwise": ("allwise", "allwise"),
    "2mass_psc": ("2mass_psc", "2mass_psc"),
    "panstarrs": ("panstarrs_vizier", "panstarrs_dr2"),
}


def _sep_arcsec(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    a = SkyCoord(float(ra1) * u.deg, float(dec1) * u.deg, frame="icrs")
    b = SkyCoord(float(ra2) * u.deg, float(dec2) * u.deg, frame="icrs")
    return float(a.separation(b).arcsec)


def _position(row: dict | None) -> tuple[float, float] | None:
    if not row:
        return None
    for ra_key, dec_key in (
        ("ra", "dec"),
        ("raMean", "decMean"),
        ("RAJ2000", "DEJ2000"),
        ("RAdeg", "DEdeg"),
    ):
        ra, dec = row.get(ra_key), row.get(dec_key)
        try:
            if ra is not None and dec is not None and math.isfinite(float(ra)) and math.isfinite(float(dec)):
                return float(ra), float(dec)
        except (TypeError, ValueError):
            pass
    return None


def _nearest(block: dict | None) -> dict | None:
    if not isinstance(block, dict):
        return None
    rows = block.get("rows") or []
    return rows[0] if rows else None


def _catalog_calibration(primary: dict, bypass: dict, primary_name: str, bypass_name: str, tolerance_arcsec: float) -> dict:
    a = primary.get(primary_name, {})
    b = bypass.get(bypass_name, {})
    row_a = _nearest(a)
    row_b = _nearest(b)
    pos_a = _position(row_a)
    pos_b = _position(row_b)
    delta = None
    agrees = None
    if pos_a and pos_b:
        delta = _sep_arcsec(pos_a[0], pos_a[1], pos_b[0], pos_b[1])
        agrees = bool(delta <= tolerance_arcsec)
    return {
        "primary_status": a.get("status"),
        "secondary_status": b.get("status"),
        "primary_count": a.get("count"),
        "secondary_count": b.get("count"),
        "primary_nearest": row_a,
        "secondary_nearest": row_b,
        "nearest_position_delta_arcsec": delta,
        "nearest_position_agrees_within_tolerance": agrees,
        "tolerance_arcsec": tolerance_arcsec,
    }


def _wcs_center_checks(skyview_dir: Path, target_ra: float, target_dec: float) -> list[dict]:
    rows: list[dict] = []
    for path in sorted(skyview_dir.rglob("*.fits")):
        item = {"file": str(path)}
        try:
            with fits.open(path, memmap=False) as hdul:
                chosen = None
                for hdu in hdul:
                    data = getattr(hdu, "data", None)
                    if data is None or getattr(data, "ndim", 0) < 2:
                        continue
                    wcs = WCS(hdu.header).celestial
                    if wcs.pixel_n_dim != 2 or wcs.world_n_dim != 2:
                        continue
                    chosen = (data, wcs)
                    break
                if chosen is None:
                    raise RuntimeError("no 2D celestial WCS HDU")
                data, wcs = chosen
                ny, nx = data.shape[-2], data.shape[-1]
                cx, cy = (nx - 1) / 2.0, (ny - 1) / 2.0
                ra, dec = wcs.pixel_to_world_values(cx, cy)
                residual = _sep_arcsec(target_ra, target_dec, float(ra), float(dec))
                item.update(
                    {
                        "status": "OK",
                        "center_ra_deg": float(ra),
                        "center_dec_deg": float(dec),
                        "requested_center_residual_arcsec": residual,
                        "shape": [int(ny), int(nx)],
                    }
                )
        except Exception as exc:
            item.update({"status": "ERROR", "error": f"{type(exc).__name__}: {exc}"})
        rows.append(item)
    return rows


def _image_channel_summary(entry: dict, skyview_dir: Path, target_ra: float, target_dec: float) -> dict:
    surveys = entry.get("skyview", {})
    survey_rows = {}
    for name, row in surveys.items():
        downloads = row.get("downloaded", []) if isinstance(row, dict) else []
        successful = [x for x in downloads if isinstance(x, dict) and not x.get("error")]
        survey_rows[name] = {
            "status": row.get("status") if isinstance(row, dict) else None,
            "links_found": len(row.get("links_found", [])) if isinstance(row, dict) else 0,
            "downloads_successful": len(successful),
            "downloads_total": len(downloads),
        }
    wcs_rows = _wcs_center_checks(skyview_dir, target_ra, target_dec)
    residuals = [r["requested_center_residual_arcsec"] for r in wcs_rows if r.get("status") == "OK"]
    return {
        "surveys": survey_rows,
        "fits_wcs_center_checks": wcs_rows,
        "fits_wcs_check_count": len(residuals),
        "fits_wcs_center_residual_arcsec_min": min(residuals) if residuals else None,
        "fits_wcs_center_residual_arcsec_max": max(residuals) if residuals else None,
    }


def run(prereg_path: Path, output_dir: Path) -> dict:
    prereg = json.loads(prereg_path.read_text(encoding="utf-8"))
    radius = float(prereg["primary_channel"]["radius_arcsec"])
    tolerance = float(prereg["cross_check_rules"]["same_source_positional_tolerance_arcsec"])

    output_dir.mkdir(parents=True, exist_ok=True)
    primary_dir = output_dir / "primary_catalog_resolver"
    bypass_dir = output_dir / "skyview_bypass_calibration"
    primary_path = primary_dir / "result.json"

    primary_result = primary_probe.run(primary_path, radius)
    bypass_scan.run(bypass_dir, radius)
    bypass_result = json.loads((bypass_dir / "open-sky-scan.json").read_text(encoding="utf-8"))

    result = {
        "schema": "janus.cosmos.love_edem.coordinate_resolver.result.v1",
        "experiment_id": prereg["experiment_id"],
        "scientific_question": prereg["scientific_question"],
        "mode": "WHAT_IS_AT_THESE_COORDINATES",
        "anomaly_scoring_used": False,
        "legacy_janus_v21_morphology_role": "DIAGNOSTIC_ONLY",
        "primary_channel": prereg["primary_channel"]["name"],
        "secondary_channel": prereg["secondary_channel"]["name"],
        "targets": {},
        "firewall": {
            "catalog_proximity_is_not_identity": True,
            "query_error_is_not_zero_sources": True,
            "secondary_channel_can_calibrate_but_not_override_authoritative_catalog_identity": True,
            "anomaly_gate_required": False,
            "edem_identity_confirmed": False,
            "love_candidate_activated": False,
            "claim_ceiling": prereg["claim_ceiling"],
        },
    }

    for primary_label, bypass_label in TARGET_MAP.items():
        p = primary_result["targets"][primary_label]
        b = bypass_result["targets"][bypass_label]
        ra = float(p["ra_deg"])
        dec = float(p["dec_deg"])
        catalog_checks = {
            common: _catalog_calibration(p, b, primary_name, bypass_name, tolerance)
            for common, (primary_name, bypass_name) in CATALOG_MAP.items()
        }
        usable_agreements = [
            v["nearest_position_agrees_within_tolerance"]
            for v in catalog_checks.values()
            if v["nearest_position_agrees_within_tolerance"] is not None
        ]
        image_dir = bypass_dir / bypass_label.lower() / "skyview"
        result["targets"][primary_label] = {
            "semantic_alias": p["semantic_alias"],
            "center_icrs": {"ra_deg": ra, "dec_deg": dec},
            "radius_arcsec": radius,
            "primary_nearest_sources": p["summary"],
            "primary_allwise_crossmatches": p.get("nearest_allwise_crossmatches"),
            "cross_service_catalog_calibration": catalog_checks,
            "cross_service_agreement_count": sum(bool(x) for x in usable_agreements),
            "cross_service_comparable_catalog_count": len(usable_agreements),
            "skyview_bypass": _image_channel_summary(b, image_dir, ra, dec),
            "interpretation": {
                "purpose": "resolve and characterize sources at the coordinate",
                "anomaly_interpretation_forbidden": True,
                "secondary_channel_role": "calibration_and_independent_cross_check",
            },
        }

    output = output_dir / "coordinate-resolver-result.json"
    output.write_text(json.dumps(result, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, default=str))
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prereg", default="data/love/LOVE_EDEM_COORDINATE_RESOLVER_PREREG.json")
    ap.add_argument("--output-dir", default="results/love_edem_coordinate_resolver")
    args = ap.parse_args()
    run(Path(args.prereg), Path(args.output_dir))


if __name__ == "__main__":
    main()
