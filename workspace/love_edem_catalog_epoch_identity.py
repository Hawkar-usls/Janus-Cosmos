#!/usr/bin/env python3
"""LIVE LOVE/EDEM catalog-epoch identity gate.

Resolves exact/reference epochs and catalog-specific positional uncertainty for
already-frozen source IDs, then asks whether those independent detections are
compatible with the corresponding Gaia DR3 stellar worldline.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable

import requests

from janus_catalog_epoch_identity import (
    CatalogMeasurement,
    conservative_unknown_corr_bound_mas2,
    covariance_from_cosigma_mas2,
    covariance_from_sigmas_corr_mas2,
    ellipse_covariance_mas2,
    evaluate_identity,
    jd_to_mjd,
)
from love_edem_epoch_worldline_test import query_gaia_live, state_from_features

SCHEMA = "janus.cosmos.love-edem.catalog-epoch-identity.v1"
SDSS_SQL_ENDPOINT = "https://skyserver.sdss.org/dr16/SkyServerWS/SearchTools/SqlSearch"
PS1_SYSTEMATIC_FLOOR_MAS = 20.0  # official MAST DR1 FAQ: likely ~20 mas 1-sigma 2-d systematic

GROUPS = {
    "LOVE_GAIA_WISE_2MASS": {
        "target": "LOVE",
        "gaia_source_id": "6163586620213012352",
        "center": [204.29573215827, -36.77911588007],
        "catalogs": {
            "ALLWISE": "J133710.98-364644.9",
            "2MASS_PSC": "13371100-3646446",
        },
    },
    "EDEM_GAIA_WISE_SDSS_PS1": {
        "target": "EDEM_SEARCH_CENTER_ZP",
        "gaia_source_id": "699051350998534656",
        "center": [139.21802643661, 30.26276387854],
        "catalogs": {
            "ALLWISE": "J091652.31+301545.8",
            "PANSTARRS_DR1": "144311392180445928",
            "SDSS_DR16": "1237664320537362761",
        },
    },
}

CATALOG_META = {
    "ALLWISE": {
        "vizier": "II/328/allwise", "id_columns": ["AllWISE"],
        "constraint_names": ["AllWISE"],
    },
    "2MASS_PSC": {
        "vizier": "II/246/out", "id_columns": ["2MASS", "_2MASS"],
        "constraint_names": ["2MASS", "_2MASS"],
    },
    "PANSTARRS_DR1": {
        "vizier": "II/349/ps1", "id_columns": ["objID"],
        "constraint_names": ["objID"],
    },
}


def _jsonable(v: Any) -> Any:
    if v is None:
        return None
    try:
        if getattr(v, "mask", False):
            return None
    except Exception:
        pass
    if isinstance(v, (str, int, float, bool)):
        return v
    try:
        return v.item()
    except Exception:
        return str(v)


def _get(d: Dict[str, Any], *names: str) -> Any:
    low = {str(k).lower(): k for k in d}
    for name in names:
        k = low.get(name.lower())
        if k is not None:
            v = d[k]
            if v is not None and str(v).strip() not in ("", "--", "nan"):
                return v
    return None


def _float(d: Dict[str, Any], *names: str) -> float | None:
    v = _get(d, *names)
    if v is None:
        return None
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except Exception:
        return None


def _id_matches(row: Dict[str, Any], expected: str, id_columns: Iterable[str]) -> bool:
    e = str(expected).strip().replace("WISEA ", "").replace("WISE ", "")
    for c in id_columns:
        v = _get(row, c)
        if v is None:
            continue
        s = str(v).strip().replace("WISEA ", "").replace("WISE ", "")
        if s == e:
            return True
    return False


def query_vizier_exact(catalog_key: str, source_id: str, center: list[float]) -> Dict[str, Any]:
    from astroquery.vizier import Vizier
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    cfg = CATALOG_META[catalog_key]
    viz = Vizier(columns=["**"], row_limit=100)
    errors = []
    for constraint in cfg["constraint_names"]:
        try:
            tables = viz.query_constraints(catalog=cfg["vizier"], **{constraint: str(source_id)})
            for table in tables:
                for row in table:
                    snap = {name: _jsonable(row[name]) for name in table.colnames}
                    if _id_matches(snap, source_id, cfg["id_columns"]):
                        snap["_janus_query_mode"] = f"EXACT_CONSTRAINT:{constraint}"
                        return snap
        except Exception as exc:
            errors.append(f"{constraint}:{type(exc).__name__}:{exc}")

    # Exact-ID matching after a tiny positional query is a recovery route, not a
    # nearest-neighbour substitution: only the frozen source ID is accepted.
    try:
        c = SkyCoord(center[0]*u.deg, center[1]*u.deg, frame="icrs")
        tables = viz.query_region(c, radius=2*u.arcsec, catalog=cfg["vizier"])
        for table in tables:
            for row in table:
                snap = {name: _jsonable(row[name]) for name in table.colnames}
                if _id_matches(snap, source_id, cfg["id_columns"]):
                    snap["_janus_query_mode"] = "POSITION_RECOVERY_WITH_EXACT_ID_ACCEPTANCE"
                    snap["_janus_constraint_errors"] = errors
                    return snap
    except Exception as exc:
        errors.append(f"region:{type(exc).__name__}:{exc}")
    raise RuntimeError(f"EXACT_SOURCE_NOT_RETURNED:{catalog_key}:{source_id}:{errors}")


def allwise_measurement(source_id: str, row: Dict[str, Any]) -> CatalogMeasurement:
    # Prefer the catalog's explicit PM-model reference position: AllWISE defines
    # RA_pm/DE_pm at MJD 55400.0. This avoids pretending a release date or a band
    # exposure mean is the position epoch.
    ra = _float(row, "RA_pm", "ra_pm")
    dec = _float(row, "DE_pm", "dec_pm")
    era = _float(row, "e_RA_pm", "sigra_pm")
    edec = _float(row, "e_DE_pm", "sigdec_pm")
    co = _float(row, "cosig_pm", "sigradec_pm")
    if None in (ra, dec, era, edec, co):
        raise ValueError("ALLWISE_PM_REFERENCE_SOLUTION_INCOMPLETE")
    cov = covariance_from_cosigma_mas2(era, edec, co)
    return CatalogMeasurement(
        catalog="ALLWISE", source_id=source_id, ra_deg=ra, dec_deg=dec,
        epoch_mjd=55400.0, covariance_mas2=cov.tolist(),
        covariance_status="MEASURED_FULL_2D_COVARIANCE_FROM_ALLWISE_PM_COSIGMA",
        epoch_status="EXACT_CATALOG_REFERENCE_EPOCH_MJD_55400",
        position_semantics="ALLWISE_MOTION_FIT_POSITION_AT_STANDARD_REFERENCE_EPOCH",
        provenance={
            "catalog_id":"II/328/allwise", "query_mode":row.get("_janus_query_mode"),
            "raw_epoch_evidence":{"RA_pm":ra,"DE_pm":dec,"reference_mjd":55400.0},
            "raw_uncertainty":{"e_RA_pm_arcsec":era,"e_DE_pm_arcsec":edec,"cosig_pm_arcsec":co},
            "stationary_position_not_used": True,
        },
    )


def twomass_measurement(source_id: str, row: Dict[str, Any]) -> CatalogMeasurement:
    ra, dec = _float(row,"RAJ2000"), _float(row,"DEJ2000")
    maj, minor, pa = _float(row,"errMaj"), _float(row,"errMin"), _float(row,"errPA")
    raw_jd = _float(row,"JD","jdate")
    if None in (ra,dec,maj,minor,pa,raw_jd):
        raise ValueError("2MASS_POSITION_EPOCH_OR_ELLIPSE_INCOMPLETE")
    mjd = jd_to_mjd(raw_jd)
    return CatalogMeasurement(
        catalog="2MASS_PSC", source_id=source_id, ra_deg=ra, dec_deg=dec,
        epoch_mjd=mjd, covariance_mas2=ellipse_covariance_mas2(maj,minor,pa).tolist(),
        covariance_status="MEASURED_FULL_2D_ERROR_ELLIPSE",
        epoch_status="EXACT_OBSERVATION_EPOCH_FROM_2MASS_JD",
        position_semantics="2MASS_PSC_POSITION_AT_SCAN_OBSERVATION_EPOCH",
        provenance={
            "catalog_id":"II/246/out", "query_mode":row.get("_janus_query_mode"),
            "raw_epoch_field":"JD", "raw_epoch_value":raw_jd, "interpreted_mjd":mjd,
            "date_precision_contract":"2MASS_JDATE_DOCUMENTED_TO_APPROX_PLUS_MINUS_30_SECONDS",
            "raw_error_ellipse":{"major_arcsec":maj,"minor_arcsec":minor,"pa_deg_east_of_north":pa},
        },
    )


def ps1_measurement(source_id: str, row: Dict[str, Any]) -> CatalogMeasurement:
    ra, dec = _float(row,"RAJ2000"), _float(row,"DEJ2000")
    era, edec, epoch = _float(row,"e_RAJ2000"), _float(row,"e_DEJ2000"), _float(row,"Epoch")
    if None in (ra,dec,era,edec,epoch):
        raise ValueError("PS1_MEAN_POSITION_OR_EPOCH_INCOMPLETE")
    cov = conservative_unknown_corr_bound_mas2(era, edec, PS1_SYSTEMATIC_FLOOR_MAS)
    return CatalogMeasurement(
        catalog="PANSTARRS_DR1", source_id=source_id, ra_deg=ra, dec_deg=dec,
        epoch_mjd=epoch, covariance_mas2=cov.tolist(),
        covariance_status="CONSERVATIVE_TRACE_ISOTROPIC_UPPER_BOUND__UNKNOWN_RA_DEC_CORRELATION__20MAS_SYSTEMATIC_FLOOR",
        epoch_status="EXACT_WEIGHTED_MEAN_EPOCH_MJD_FROM_PS1_EPOCH",
        position_semantics="PS1_WEIGHTED_MEAN_SINGLE_EPOCH_DETECTION_POSITION_AT_EPOCHMEAN",
        provenance={
            "catalog_id":"II/349/ps1", "query_mode":row.get("_janus_query_mode"),
            "raw_epoch_field":"Epoch", "raw_epoch_mjd":epoch,
            "raw_statistical_errors_arcsec":{"ra":era,"dec":edec},
            "systematic_floor_mas":PS1_SYSTEMATIC_FLOOR_MAS,
            "unknown_correlation_policy":"TRACE_TIMES_IDENTITY_UPPER_BOUND_NOT_ZERO_CORRELATION",
        },
    )


def query_sdss_exact(objid: str) -> Dict[str, Any]:
    sql = f"""SELECT TOP 1 p.objID,p.ra,p.dec,p.raErr,p.decErr,p.raDecCorr,p.fieldID,
 f.mjd_u,f.mjd_g,f.mjd_r,f.mjd_i,f.mjd_z
 FROM PhotoObjAll p JOIN Field f ON p.fieldID=f.fieldID
 WHERE p.objID={int(objid)}"""
    resp = requests.get(SDSS_SQL_ENDPOINT, params={"cmd":sql,"format":"csv"}, timeout=45)
    resp.raise_for_status()
    text = resp.text.strip()
    if not text:
        raise RuntimeError("SDSS_EMPTY_RESPONSE")
    # SkyServer CSV is normally a header followed by rows; tolerate comment lines.
    lines = [line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    rows = list(csv.DictReader(io.StringIO("\n".join(lines))))
    if not rows:
        raise RuntimeError(f"SDSS_NO_ROW:{text[:500]}")
    row = rows[0]
    if str(_get(row,"objID")).strip() != str(objid):
        raise RuntimeError("SDSS_OBJID_MISMATCH")
    row["_janus_query_mode"] = "EXACT_OBJID_SQL_JOIN_PHOTOOBJALL_FIELD"
    row["_janus_endpoint"] = SDSS_SQL_ENDPOINT
    return row


def sdss_measurement(source_id: str, row: Dict[str, Any]) -> CatalogMeasurement:
    ra,dec = _float(row,"ra"),_float(row,"dec")
    era,edec,rho = _float(row,"raErr"),_float(row,"decErr"),_float(row,"raDecCorr")
    mjdr = _float(row,"mjd_r")
    if None in (ra,dec,era,edec,rho,mjdr):
        raise ValueError("SDSS_RBAND_POSITION_EPOCH_OR_COVARIANCE_INCOMPLETE")
    cov = covariance_from_sigmas_corr_mas2(era,edec,rho)
    return CatalogMeasurement(
        catalog="SDSS_DR16", source_id=source_id, ra_deg=ra, dec_deg=dec,
        epoch_mjd=mjdr, covariance_mas2=cov.tolist(),
        covariance_status="MEASURED_FULL_2D_COVARIANCE_FROM_RAERR_DECERR_RADEC_CORR",
        epoch_status="EXACT_R_BAND_FIELD_EPOCH_MJD_R",
        position_semantics="SDSS_PHOTOOBJ_FINAL_RA_DEC_ARE_R_BAND_ASTROMETRY",
        provenance={
            "query_mode":row.get("_janus_query_mode"), "endpoint":row.get("_janus_endpoint"),
            "field_id":_get(row,"fieldID"),
            "band_mjd":{b:_float(row,f"mjd_{b}") for b in "ugriz"},
            "raw_uncertainty":{"raErr_arcsec":era,"decErr_arcsec":edec,"raDecCorr":rho},
        },
    )


def _resolve_measurement(catalog: str, source_id: str, center: list[float]) -> CatalogMeasurement:
    if catalog == "SDSS_DR16":
        return sdss_measurement(source_id, query_sdss_exact(source_id))
    row = query_vizier_exact(catalog, source_id, center)
    if catalog == "ALLWISE": return allwise_measurement(source_id,row)
    if catalog == "2MASS_PSC": return twomass_measurement(source_id,row)
    if catalog == "PANSTARRS_DR1": return ps1_measurement(source_id,row)
    raise KeyError(catalog)


def run() -> Dict[str, Any]:
    groups: Dict[str, Any] = {}
    resolved_epoch_count = 0
    required_count = sum(len(g["catalogs"]) for g in GROUPS.values())
    for group_name,cfg in GROUPS.items():
        gaia_error = None
        try:
            gf = query_gaia_live(cfg["gaia_source_id"])
            state = state_from_features(cfg["gaia_source_id"], gf)
        except Exception as exc:
            groups[group_name] = {"status":"I_DO_NOT_KNOW","reason":f"GAIA_QUERY_FAILED:{type(exc).__name__}:{exc}"}
            continue

        measurements=[]; audit={}
        for catalog,source_id in cfg["catalogs"].items():
            try:
                m=_resolve_measurement(catalog,source_id,cfg["center"])
                measurements.append(m); resolved_epoch_count += 1
                audit[catalog]={
                    "source_id":source_id,"status":"EPOCH_AND_POSITION_COVARIANCE_RESOLVED",
                    "epoch_mjd":m.epoch_mjd,"epoch_jyear":m.epoch_jyear,
                    "epoch_status":m.epoch_status,"covariance_status":m.covariance_status,
                    "position_semantics":m.position_semantics,
                }
            except Exception as exc:
                audit[catalog]={"source_id":source_id,"status":"BLOCKED","reason":f"{type(exc).__name__}:{exc}"}

        identity = evaluate_identity(state, measurements) if measurements else {
            "status":"I_DO_NOT_KNOW","reason":"NO_EXTERNAL_MEASUREMENTS_RESOLVED"
        }
        groups[group_name] = {
            "status":"EVALUATED" if measurements else "I_DO_NOT_KNOW",
            "target":cfg["target"],
            "gaia_source_id":cfg["gaia_source_id"],
            "catalog_epoch_audit":audit,
            "identity_test":identity,
        }

    return {
        "schema":SCHEMA,
        "experiment":"LOVE_EDEM_CATALOG_EPOCH_IDENTITY_GEN3",
        "formula":"RESPICIENS_ET_PROSPICIENS_GEN3",
        "run_time_utc":datetime.now(timezone.utc).isoformat(),
        "groups":groups,
        "summary":{
            "required_external_catalog_measurements":required_count,
            "resolved_external_catalog_measurements":resolved_epoch_count,
            "all_required_epochs_resolved":resolved_epoch_count==required_count,
            "release_year_substitutions":0,
            "simulation_count_increased":False,
        },
        "archive_search_bridge":{
            "lineage":"TOPA tools/topa_arxiv_gateway.py",
            "role":"DISCOVERY_AND_PROVENANCE_INDEX_ONLY",
            "catalog_values_must_come_from_primary_catalog_or_archive":True,
            "search_rank_is_truth":False,
        },
        "epistemic_firewall":{
            "catalog_agreement_is_love_edem_identity":False,
            "catalog_agreement_is_anomaly":False,
            "missing_epoch_may_use_release_year":False,
            "missing_covariance_may_be_zeroed":False,
            "mahalanobis_pass_is_identity_proof":False,
            "negative_result_is_valid":True,
        },
        "claim_ceiling":"CROSS_CATALOG_COMPATIBILITY_OF_ALREADY_GROUPED_STELLAR_DETECTIONS_ONLY__NOT_LOVE_EDEM_IDENTITY",
    }


def self_test() -> None:
    # No network: test the exact epoch/covariance constructors against synthetic rows.
    a=allwise_measurement("W",{"RA_pm":10.0,"DE_pm":20.0,"e_RA_pm":.1,"e_DE_pm":.2,"cosig_pm":.01})
    t=twomass_measurement("T",{"RAJ2000":10.0,"DEJ2000":20.0,"errMaj":.2,"errMin":.1,"errPA":30.0,"JD":2451545.0})
    p=ps1_measurement("P",{"RAJ2000":10.0,"DEJ2000":20.0,"e_RAJ2000":.01,"e_DEJ2000":.02,"Epoch":56000.0})
    assert a.epoch_mjd==55400.0 and abs(t.epoch_mjd-51544.5)<1e-9 and p.epoch_mjd==56000.0
    assert "CONSERVATIVE" in p.covariance_status
    print("LOVE_EDEM_CATALOG_EPOCH_IDENTITY_SELF_TEST=PASS")


def main() -> int:
    p=argparse.ArgumentParser(); p.add_argument("--output"); p.add_argument("--self-test",action="store_true"); a=p.parse_args()
    if a.self_test: self_test(); return 0
    out=run(); text=json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)
    if a.output: Path(a.output).write_text(text+"\n",encoding="utf-8")
    else: print(text)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
