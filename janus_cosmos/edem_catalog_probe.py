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

from astroquery.gaia import Gaia


def _safe_float(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except Exception:
        return None


def _query_known_planets(ra_deg: float, dec_deg: float, radius_deg: float) -> tuple[list[dict], str | None]:
    try:
        query = (
            "select pl_name,hostname,ra,dec,discoverymethod,disc_year from pscomppars where "
            f"contains(point('ICRS',ra,dec),circle('ICRS',{ra_deg},{dec_deg},{radius_deg}))=1"
        )
        url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urllib.parse.urlencode({"query": query, "format": "csv"})
        req = urllib.request.Request(url, headers={"User-Agent": "Janus-Cosmos-EDEM-Catalog-Probe/1.0"})
        with urllib.request.urlopen(req, timeout=90) as response:
            text = response.read().decode("utf-8")
        rows = []
        for row in csv.DictReader(io.StringIO(text)):
            rows.append({
                "pl_name": row.get("pl_name"),
                "hostname": row.get("hostname"),
                "ra_deg": _safe_float(row.get("ra")),
                "dec_deg": _safe_float(row.get("dec")),
                "discoverymethod": row.get("discoverymethod"),
                "disc_year": row.get("disc_year"),
            })
        return rows, None
    except Exception as exc:
        return [], f"{type(exc).__name__}: {exc}"


def run(receipt_path: Path, output_path: Path, radius_deg: float = 0.25) -> dict:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    center = receipt["edem_geometry_direction_candidate_icrs"]
    ra = float(center["ra_deg"])
    dec = float(center["dec_deg"])
    radius = float(radius_deg)

    adql = f"""
SELECT TOP 50
  source_id, ra, dec, parallax, parallax_error, pmra, pmdec,
  phot_g_mean_mag, bp_rp, ruwe, teff_gspphot,
  DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra}, {dec})) AS sep_deg
FROM gaiadr3.gaia_source
WHERE 1 = CONTAINS(
  POINT('ICRS', ra, dec),
  CIRCLE('ICRS', {ra}, {dec}, {radius})
)
ORDER BY sep_deg ASC
""".strip()

    job = Gaia.launch_job_async(adql, dump_to_file=False)
    table = job.get_results()
    gaia_rows = []
    for row in table:
        gaia_rows.append({
            "source_id": str(row["source_id"]),
            "ra_deg": _safe_float(row["ra"]),
            "dec_deg": _safe_float(row["dec"]),
            "separation_arcsec": (_safe_float(row["sep_deg"]) or 0.0) * 3600.0,
            "parallax_mas": _safe_float(row["parallax"]),
            "parallax_error_mas": _safe_float(row["parallax_error"]),
            "pmra_masyr": _safe_float(row["pmra"]),
            "pmdec_masyr": _safe_float(row["pmdec"]),
            "phot_g_mean_mag": _safe_float(row["phot_g_mean_mag"]),
            "bp_rp": _safe_float(row["bp_rp"]),
            "ruwe": _safe_float(row["ruwe"]),
            "teff_gspphot_k": _safe_float(row["teff_gspphot"]),
        })

    planets, planet_error = _query_known_planets(ra, dec, radius)
    nearest = gaia_rows[0] if gaia_rows else None
    result = {
        "schema": "janus.cosmos.edem.catalog_probe.result.v1",
        "experiment_id": "EDEM-GIZA-CATALOG-PROBE-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "parent_receipt": str(receipt_path),
        "center_icrs": center,
        "search_radius_deg": radius,
        "gaia_dr3": {
            "returned_top_count": len(gaia_rows),
            "nearest_source": nearest,
            "sources": gaia_rows,
        },
        "nasa_exoplanet_archive": {
            "known_planets_in_cone_count": len(planets),
            "known_planets": planets,
            "query_error": planet_error,
        },
        "interpretation": "Catalog population around the geometry-derived direction only. Nearest catalog source is not automatically Edem.",
        "edem_identity_confirmed": False,
        "claim_ceiling": "CATALOG_PROBE_AROUND_GIZA_DERIVED_DIRECTION_ONLY__NO_EDEM_IDENTIFICATION",
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "center": center,
        "gaia_count": len(gaia_rows),
        "nearest_gaia": nearest,
        "known_planets_in_cone": len(planets),
        "planet_query_error": planet_error,
    }, indent=2))
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description="Gaia/NASA catalog probe around the Giza-derived Edem direction")
    ap.add_argument("--receipt", default="data/love/EDEM-GIZA-REVERSE-SPEAR-v1-RUN-001-RECEIPT.json")
    ap.add_argument("--output", default="results/edem_catalog_probe/edem-catalog-probe.json")
    ap.add_argument("--radius-deg", type=float, default=0.25)
    args = ap.parse_args()
    run(Path(args.receipt), Path(args.output), args.radius_deg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
