from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18
from astroquery.vizier import Vizier

# Frozen source-group positions from the prior Janus-Cosmos source-nature receipt.
TARGETS = {
    "EDEM_NEAREST_WISE_ONLY": {
        "coord": SkyCoord(139.2249571 * u.deg, 30.262781 * u.deg, frame="icrs"),
        "allwise_id": "J091653.98+301546.0",
        "optical_catalogs": ["II/349/ps1"],
    },
    "LOVE_NEAREST_WISE_ONLY": {
        "coord": SkyCoord(204.2960668 * u.deg, -36.7804774 * u.deg, frame="icrs"),
        "allwise_id": "J133711.05-364649.7",
        "optical_catalogs": ["II/379/smssdr4", "II/371/des_dr2"],
    },
}

EDEM_GALAXY = SkyCoord(139.222681 * u.deg, 30.256918 * u.deg, frame="icrs")

CATALOGS = {
    "CATWISE2020": "II/365/catwise",
    "UNWISE": "II/363/unwise",
    "PANSTARRS_DR1": "II/349/ps1",
    "SKYMAPPER_DR4": "II/379/smssdr4",
    "DES_DR2": "II/371/des_dr2",
    "SDSS_DR16": "V/154/sdss16",
    "DESI_DR1": "V/161/zcatdr1",
}

POSITION_PAIRS = [
    ("RA_ICRS", "DE_ICRS"),
    ("RAICRS", "DEICRS"),
    ("RAJ2000", "DEJ2000"),
    ("RAdeg", "DEdeg"),
    ("RA", "DEC"),
]

FEATURE_NAMES = [
    "Name", "objID", "ObjectId", "DES", "SMSS", "TargetID",
    "W1mproPM", "e_W1mproPM", "W2mproPM", "e_W2mproPM", "snrW1pm", "snrW2pm",
    "pmRA", "pmDE", "e_pmRA", "e_pmDE", "pmQual", "abf",
    "FW1", "FW2", "e_FW1", "e_FW2", "q_W1", "q_W2", "fFW1", "fFW2",
    "gmag", "rmag", "imag", "zmag", "ymag", "umag",
    "class", "zsp", "e_zsp", "f_zsp", "spCl", "subCl", "SpObjID",
    "z", "e_z", "ZWARN", "SPECTYPE", "SUBTYPE",
]


def jsonable(v):
    if v is None:
        return None
    try:
        if getattr(v, "mask", False):
            return None
    except Exception:
        pass
    try:
        x = v.item()
    except Exception:
        x = v
    if isinstance(x, bytes):
        return x.decode("utf-8", errors="replace")
    if isinstance(x, (str, int, float, bool)):
        if isinstance(x, float) and not math.isfinite(x):
            return None
        return x
    return str(x)


def position_columns(names):
    lower = {n.lower(): n for n in names}
    for ra, dec in POSITION_PAIRS:
        if ra.lower() in lower and dec.lower() in lower:
            return lower[ra.lower()], lower[dec.lower()]
    return None, None


def query_nearest(catalog_id: str, center: SkyCoord, radius_arcsec: float = 2.0):
    result = {"catalog_id": catalog_id, "radius_arcsec": radius_arcsec, "status": "UNKNOWN", "count": 0, "nearest": None}
    try:
        viz = Vizier(columns=["*"], row_limit=50)
        tables = viz.query_region(center, radius=radius_arcsec * u.arcsec, catalog=catalog_id)
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

        candidates = []
        for row, table in zip(rows, table_for_row):
            names = list(table.colnames)
            ra_col, dec_col = position_columns(names)
            if not ra_col or not dec_col:
                continue
            try:
                ra = float(row[ra_col])
                dec = float(row[dec_col])
            except Exception:
                continue
            c = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
            sep = float(center.separation(c).to_value(u.arcsec))
            features = {}
            for name in FEATURE_NAMES:
                if name in names:
                    features[name] = jsonable(row[name])
            candidates.append({"ra_deg": ra, "dec_deg": dec, "separation_arcsec": sep, "features": features})
        if candidates:
            result["nearest"] = min(candidates, key=lambda x: x["separation_arcsec"])
        else:
            result["status"] = "POSITION_COLUMNS_NOT_FOUND"
        return result
    except Exception as exc:
        result["status"] = "QUERY_ERROR"
        result["error"] = f"{type(exc).__name__}: {exc}"
        return result


def catwise_motion(nearest):
    if not nearest:
        return None
    f = nearest.get("features", {})
    vals = [f.get("pmRA"), f.get("pmDE"), f.get("e_pmRA"), f.get("e_pmDE")]
    if any(v is None for v in vals):
        return None
    pmra, pmde, epmra, epmde = map(float, vals)
    if epmra <= 0 or epmde <= 0:
        return None
    sig = math.sqrt((pmra / epmra) ** 2 + (pmde / epmde) ** 2)
    return {
        "pmra_arcsec_per_year": pmra,
        "pmdec_arcsec_per_year": pmde,
        "pmra_error_arcsec_per_year": epmra,
        "pmdec_error_arcsec_per_year": epmde,
        "vector_significance_sigma": sig,
        "pm_quality": f.get("pmQual"),
        "significant_at_frozen_5sigma_gate": sig >= 5.0,
        "interpretation": "motion candidate only; CatWISE motion alone is not promoted to a stellar identity in this gate",
    }


def deep_source_pass(name: str, spec: dict):
    center = spec["coord"]
    catwise = query_nearest(CATALOGS["CATWISE2020"], center, 2.0)
    unwise = query_nearest(CATALOGS["UNWISE"], center, 2.0)
    optical = {}
    for catalog_id in spec["optical_catalogs"]:
        key = next(k for k, v in CATALOGS.items() if v == catalog_id)
        optical[key] = query_nearest(catalog_id, center, 2.0)

    deep_ir_matches = [
        key for key, item in (("CATWISE2020", catwise), ("UNWISE", unwise))
        if item.get("nearest") is not None and item["nearest"]["separation_arcsec"] <= 2.0
    ]
    optical_matches = [
        key for key, item in optical.items()
        if item.get("nearest") is not None and item["nearest"]["separation_arcsec"] <= 2.0
    ]
    motion = catwise_motion(catwise.get("nearest"))

    if len(deep_ir_matches) == 2 and not optical_matches:
        nature = "PERSISTENT_DEEP_WISE_INFRARED_SOURCE__OPTICAL_COUNTERPART_NOT_RECOVERED"
        confidence = "MEDIUM_HIGH_FOR_IR_REALITY__LOW_FOR_ASTROPHYSICAL_CLASS"
    elif optical_matches:
        nature = "INFRARED_SOURCE_WITH_DEEP_OPTICAL_COUNTERPART"
        confidence = "MEDIUM_HIGH"
    else:
        nature = "DEEP_COUNTERPART_INCONCLUSIVE"
        confidence = "LOW_MEDIUM"

    return {
        "allwise_anchor": spec["allwise_id"],
        "representative_icrs": {"ra_deg": float(center.ra.deg), "dec_deg": float(center.dec.deg)},
        "queries": {"CATWISE2020": catwise, "UNWISE": unwise, **optical},
        "deep_ir_match_catalogs": deep_ir_matches,
        "optical_match_catalogs": optical_matches,
        "catwise_motion_test": motion,
        "result": {
            "nature_gate": nature,
            "confidence": confidence,
            "plain_ru": (
                "Источник повторно виден в более глубоких WISE-каталогах, но оптический counterpart в замороженном радиусе не найден; это укрепляет реальность ИК-источника, но не делает его автоматически звездой, галактикой или планетой."
                if len(deep_ir_matches) == 2 and not optical_matches else
                "У ИК-источника найден глубокий оптический counterpart; следующий шаг — объединённая фотометрическая классификация."
                if optical_matches else
                "Глубокий проход не дал достаточного нового counterpart для более узкой классификации."
            )
        }
    }


def extract_redshift(item, source_name):
    nearest = item.get("nearest")
    if not nearest:
        return None
    f = nearest.get("features", {})
    if source_name == "SDSS_DR16":
        z = f.get("zsp")
        zerr = f.get("e_zsp")
        zwarn = f.get("f_zsp")
        sclass = f.get("spCl")
        good = z is not None and (zwarn in (None, 0, "0"))
    else:
        z = f.get("z")
        zerr = f.get("e_z")
        zwarn = f.get("ZWARN")
        sclass = f.get("SPECTYPE")
        good = z is not None and (zwarn in (None, 0, "0"))
    if z is None:
        return None
    zf = float(z)
    out = {
        "source": source_name,
        "z": zf,
        "z_error": None if zerr is None else float(zerr),
        "warning": zwarn,
        "spectral_class": sclass,
        "accepted_by_frozen_quality_gate": bool(good),
        "separation_arcsec": nearest["separation_arcsec"],
    }
    if good and zf > 0:
        try:
            out["planck18_luminosity_distance_mpc"] = float(Planck18.luminosity_distance(zf).to_value(u.Mpc))
            out["planck18_comoving_distance_mpc"] = float(Planck18.comoving_distance(zf).to_value(u.Mpc))
        except Exception:
            pass
    return out


def galaxy_redshift_pass():
    sdss = query_nearest(CATALOGS["SDSS_DR16"], EDEM_GALAXY, 2.0)
    desi = query_nearest(CATALOGS["DESI_DR1"], EDEM_GALAXY, 2.0)
    redshifts = []
    for src, item in (("SDSS_DR16", sdss), ("DESI_DR1", desi)):
        r = extract_redshift(item, src)
        if r is not None:
            redshifts.append(r)
    accepted = [r for r in redshifts if r["accepted_by_frozen_quality_gate"]]
    if accepted:
        result = "SPECTROSCOPIC_REDSHIFT_RECOVERED"
        plain = "Для EDEM-галактики найден качественный спектроскопический redshift; расстояние можно оценивать космологически."
    else:
        result = "NO_ACCEPTED_SPECTROSCOPIC_REDSHIFT_IN_SDSS_DR16_OR_DESI_DR1"
        plain = "Морфологическая классификация галактики сохраняется, но качественный спектроскопический redshift в SDSS DR16/DESI DR1 в этом проходе не найден."
    return {
        "representative_icrs": {"ra_deg": float(EDEM_GALAXY.ra.deg), "dec_deg": float(EDEM_GALAXY.dec.deg)},
        "queries": {"SDSS_DR16": sdss, "DESI_DR1": desi},
        "redshift_candidates": redshifts,
        "accepted_redshifts": accepted,
        "result": result,
        "plain_ru": plain,
    }


def main():
    out = Path("data/love/LOVE-EDEM-DEEP-COUNTERPART-REDSHIFT-v1-LATEST-RECEIPT.json")
    payload = {
        "schema": "janus.cosmos.love_edem.deep_counterpart_redshift.receipt.v1",
        "experiment_id": "LOVE-EDEM-DEEP-COUNTERPART-REDSHIFT-v1",
        "status": "COMPLETE",
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "DEEP_COUNTERPART_AND_REDSHIFT_WITHOUT_ANOMALY_SCORING",
        "anomaly_scoring_used": False,
        "frozen_rules": {
            "counterpart_radius_arcsec": 2.0,
            "catwise_motion_significance_gate_sigma": 5.0,
            "catwise_motion_alone_does_not_establish_stellar_identity": True,
            "spectroscopic_redshift_acceptance": "non-null z and warning bitmask == 0 (or absent in mirror row)",
            "no_post_hoc_threshold_changes": True,
        },
        "deep_sources": {name: deep_source_pass(name, spec) for name, spec in TARGETS.items()},
        "edem_optical_galaxy_redshift": galaxy_redshift_pass(),
        "provenance": {
            "catalog_frontend": "CDS VizieR via astroquery.vizier",
            "catalogs": CATALOGS,
            "cosmology_if_redshift_available": "Astropy Planck18",
        },
        "firewall": {
            "deep_ir_persistence_is_not_anomaly": True,
            "deep_ir_persistence_is_not_planet": True,
            "catwise_motion_candidate_is_not_automatic_star": True,
            "near_source_is_not_semantic_LOVE_or_EDEM_identity": True,
            "physical_LOVE_EDEM_ORION_association_not_established": True,
            "claim_ceiling": "DEEP_CATALOG_COUNTERPART_AND_SPECTROSCOPIC_REDSHIFT_CHECK_ONLY",
        },
    }
    failures = []
    for source in payload["deep_sources"].values():
        for key, item in source["queries"].items():
            if item.get("status") == "QUERY_ERROR":
                failures.append({"catalog": key, "error": item.get("error")})
    for key, item in payload["edem_optical_galaxy_redshift"]["queries"].items():
        if item.get("status") == "QUERY_ERROR":
            failures.append({"catalog": key, "error": item.get("error")})
    payload["service_failures"] = failures
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "deep_source_results": {k: v["result"] for k, v in payload["deep_sources"].items()},
        "galaxy_redshift_result": payload["edem_optical_galaxy_redshift"]["result"],
        "service_failure_count": len(failures),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
