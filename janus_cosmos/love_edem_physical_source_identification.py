from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from astropy import units as u
from astropy.coordinates import SkyCoord

from janus_cosmos import love_edem_center_object_probe as primary_probe
from janus_cosmos import love_edem_gaia_mirror_cluster as cluster_probe


def _sep_arcsec(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    a = SkyCoord(float(ra1) * u.deg, float(dec1) * u.deg, frame="icrs")
    b = SkyCoord(float(ra2) * u.deg, float(dec2) * u.deg, frame="icrs")
    return float(a.separation(b).arcsec)


def _find(rows: list[dict], key: str, value: str) -> dict | None:
    for row in rows or []:
        if str(row.get(key)) == str(value):
            return row
    return None


def _finite(v) -> bool:
    try:
        return math.isfinite(float(v))
    except Exception:
        return False


def _member_cluster(clusters: list[dict], catalog: str, source_id: str) -> dict | None:
    for cluster in clusters:
        for member in cluster.get("members", []):
            if member.get("catalog") == catalog and str(member.get("source_id")) == str(source_id):
                return cluster
    return None


def _wise_assessment(row: dict | None, cluster: dict | None) -> dict:
    if not row:
        return {"status": "MISSING", "classification": "INSUFFICIENT_DATA"}
    cc = str(row.get("cc_flags") or "")
    ext = row.get("ext_flg")
    nb = row.get("nb")
    na = row.get("na")
    w1snr = row.get("w1snr")
    clean_flags = cc == "0000" and ext == 0 and nb == 1 and na == 0 and _finite(w1snr) and float(w1snr) > 0
    catalogs = set(cluster.get("catalogs", [])) if cluster else {"allwise"}
    counterparts = sorted(c for c in catalogs if c != "allwise")
    if clean_flags and not counterparts:
        label = "WISE_CLEAN_POINTLIKE_IR_ONLY_UNCLASSIFIED"
    elif not clean_flags:
        label = "WISE_FLAGGED_OR_BLEND_SUSPECT"
    else:
        label = "INSUFFICIENT_DATA"
    return {
        "status": "OK",
        "classification": label,
        "allwise_flags": {
            "ph_qual": row.get("ph_qual"), "cc_flags": row.get("cc_flags"),
            "ext_flg": row.get("ext_flg"), "nb": row.get("nb"), "na": row.get("na"),
            "w1_mag": row.get("w1mpro"), "w1_snr": row.get("w1snr"),
            "w2_snr": row.get("w2snr"), "w3_snr": row.get("w3snr"), "w4_snr": row.get("w4snr")
        },
        "catalog_level_clean_pointlike_flags": clean_flags,
        "cluster_catalogs": sorted(catalogs),
        "counterpart_catalogs_under_frozen_tolerances": counterparts,
        "interpretation": "AllWISE has no catalog-level artifact/extended/blend flag for this source, but no frozen-radius optical/NIR/Gaia counterpart means its physical class remains unresolved." if label == "WISE_CLEAN_POINTLIKE_IR_ONLY_UNCLASSIFIED" else "Current catalog metadata does not support a clean IR-only classification."
    }


def _optical_assessment(sdss: dict | None, ps1: dict | None) -> dict:
    if not sdss or not ps1:
        return {"status": "MISSING", "classification": "INSUFFICIENT_DATA"}
    sep = _sep_arcsec(sdss["ra"], sdss["dec"], ps1["RAJ2000"], ps1["DEJ2000"])
    label = "SDSS_GALAXY_MORPHOLOGY_SUPPORTED" if sdss.get("type") == 3 and sep <= 1.0 else "INSUFFICIENT_DATA"
    return {
        "status": "OK",
        "classification": label,
        "sdss_type": sdss.get("type"),
        "sdss_photometry": {k: sdss.get(k) for k in ("u", "g", "r", "i", "z")},
        "panstarrs_photometry": {k: ps1.get(k) for k in ("gmag", "rmag", "imag", "zmag", "ymag")},
        "sdss_panstarrs_separation_arcsec": sep,
        "interpretation": "SDSS type=3 supports galaxy morphology and Pan-STARRS independently matches the same optical position within 1 arcsec." if label == "SDSS_GALAXY_MORPHOLOGY_SUPPORTED" else "Optical evidence does not satisfy the frozen v1 galaxy-support rule."
    }


def _gaia_assessment(row: dict | None, cluster: dict | None, snr_cut: float) -> dict:
    if not row:
        return {"status": "MISSING", "classification": "INSUFFICIENT_DATA"}
    plx, eplx = row.get("Plx"), row.get("e_Plx")
    snr = float(plx) / float(eplx) if _finite(plx) and _finite(eplx) and float(eplx) != 0 else None
    distance_pc = 1000.0 / float(plx) if snr is not None and snr >= snr_cut and float(plx) > 0 else None
    return {
        "status": "OK",
        "classification": "GAIA_LINKED_STELLAR_SOURCE_GROUP" if cluster else "INSUFFICIENT_DATA",
        "source_id": str(row.get("Source")), "Gmag": row.get("Gmag"), "BP_RP": row.get("BP-RP"), "RUWE": row.get("RUWE"),
        "parallax_mas": plx, "parallax_error_mas": eplx, "parallax_snr": snr,
        "pmRA_masyr": row.get("pmRA"), "pmDE_masyr": row.get("pmDE"),
        "naive_inverse_parallax_distance_pc": distance_pc,
        "distance_rule": f"reported only when parallax SNR >= {snr_cut}",
        "cluster_catalogs": sorted(set(cluster.get("catalogs", []))) if cluster else ["gaia_dr3_vizier"]
    }


def run(prereg_path: Path, output_dir: Path) -> dict:
    prereg = json.loads(prereg_path.read_text(encoding="utf-8"))
    output_dir.mkdir(parents=True, exist_ok=True)
    primary = primary_probe.run(output_dir / "primary-catalog-result.json", 30.0)
    f = prereg["frozen_matching_rules"]
    pair_cfg = {"default": float(f["default_optical_arcsec"]), "with_2mass": float(f["with_2mass_arcsec"]), "with_allwise": float(f["with_allwise_arcsec"])}
    target_results = {}
    for label, target_key in (("LOVE", "TARGET_A"), ("EDEM", "TARGET_B")):
        center = prereg["targets"][label]["center_icrs"]
        ra, dec = float(center["ra_deg"]), float(center["dec_deg"])
        p = primary["targets"][target_key]
        try:
            gaia = cluster_probe.gaia_vizier(ra, dec, 30.0)
        except Exception as exc:
            gaia = {"status": "ERROR", "count": 0, "rows": [], "error": f"{type(exc).__name__}: {exc}"}
        detections = cluster_probe._collect(p, gaia)
        clusters = cluster_probe._cluster(detections, pair_cfg, ra, dec)
        target_results[label] = {"primary": p, "gaia": gaia, "clusters": clusters}

    cut = float(prereg["distance_rule"]["naive_inverse_parallax_allowed_only_if_snr_gte"])
    edem, ef = target_results["EDEM"], prereg["targets"]["EDEM"]["focus_sources"]
    ewise = _find(edem["primary"]["allwise"].get("rows", []), "designation", ef["WISE_NEAREST"])
    esdss = _find(edem["primary"]["sdss_dr17"].get("rows", []), "objid", ef["OPTICAL_SDSS"])
    eps1 = _find(edem["primary"]["panstarrs_vizier"].get("rows", []), "objID", ef["OPTICAL_PS1"])
    egaia = _find(edem["gaia"].get("rows", []), "Source", ef["GAIA_CONTROL"])
    ewise_c = _member_cluster(edem["clusters"], "allwise", ef["WISE_NEAREST"])
    esdss_c = _member_cluster(edem["clusters"], "sdss_dr17", ef["OPTICAL_SDSS"])
    egaia_c = _member_cluster(edem["clusters"], "gaia_dr3_vizier", ef["GAIA_CONTROL"])
    wise_opt_sep = _sep_arcsec(ewise["ra"], ewise["dec"], esdss["ra"], esdss["dec"]) if ewise and esdss else None

    love, lf = target_results["LOVE"], prereg["targets"]["LOVE"]["focus_sources"]
    lwise = _find(love["primary"]["allwise"].get("rows", []), "designation", lf["WISE_NEAREST"])
    lgaia = _find(love["gaia"].get("rows", []), "Source", lf["GAIA_MULTICATALOG_CONTROL"])
    lwise_c = _member_cluster(love["clusters"], "allwise", lf["WISE_NEAREST"])
    lgaia_c = _member_cluster(love["clusters"], "gaia_dr3_vizier", lf["GAIA_MULTICATALOG_CONTROL"])

    result = {
        "schema": "janus.cosmos.love_edem.physical_source_identification.result.v1",
        "experiment_id": prereg["experiment_id"], "mode": prereg["mode"],
        "anomaly_scoring_used": False, "anomaly_gate_required": False,
        "EDEM": {
            "center_icrs": prereg["targets"]["EDEM"]["center_icrs"],
            "wise_nearest": _wise_assessment(ewise, ewise_c),
            "optical_13arcsec_group": _optical_assessment(esdss, eps1),
            "gaia_control_group": _gaia_assessment(egaia, egaia_c, cut),
            "wise_to_optical_group_separation_arcsec": wise_opt_sep,
            "wise_and_optical_same_frozen_cluster": bool(ewise_c and esdss_c and ewise_c.get("cluster_id") == esdss_c.get("cluster_id")),
            "edem_identity_confirmed": False
        },
        "LOVE": {
            "center_icrs": prereg["targets"]["LOVE"]["center_icrs"],
            "wise_nearest": _wise_assessment(lwise, lwise_c),
            "gaia_multicatalog_control": _gaia_assessment(lgaia, lgaia_c, cut)
        },
        "firewall": prereg["firewall"],
        "next_gate": "DEEPER_WISE_COUNTERPART_CHECK_OR_SED_ONLY_IF_CURRENT_CLASSIFICATION_REMAINS_UNRESOLVED"
    }
    out = output_dir / "physical-source-identification-result.json"
    out.write_text(json.dumps(result, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, default=str))
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prereg", default="data/love/LOVE_EDEM_PHYSICAL_SOURCE_IDENTIFICATION_PREREG.json")
    ap.add_argument("--output-dir", default="results/love_edem_physical_source_identification")
    args = ap.parse_args()
    run(Path(args.prereg), Path(args.output_dir))


if __name__ == "__main__":
    main()
