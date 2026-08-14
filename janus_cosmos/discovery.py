from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

DEFAULT_FILTERS = ("F435W", "F555W", "F606W", "F814W")
DEFAULT_TARGETS = (
    "NGC1365", "NGC1425", "NGC1637", "NGC2841", "NGC3031", "NGC3627", "NGC4321",
    "M51", "M81", "M82", "NGC253", "NGC4258", "NGC4565", "NGC5194", "NGC5195",
    "NGC2403", "NGC6946", "NGC5457", "NGC3351", "NGC3621", "NGC6744", "NGC7793",
)
FOCUS_TARGETS = {"NGC1425", "NGC1637"}


@dataclass(frozen=True)
class DiscoveryConfig:
    radius_deg: float = 0.20
    retries: int = 3
    retry_sleep_seconds: float = 2.0
    filters: tuple[str, ...] = DEFAULT_FILTERS


def canonical_name(target: str) -> str:
    t = target.strip()
    if re.fullmatch(r"NGC\s*\d+", t, flags=re.I):
        number = re.search(r"\d+", t).group(0)
        return f"NGC {number}"
    return t


def _value(row, key: str, default=""):
    try:
        v = row[key]
    except Exception:
        try:
            v = row.get(key, default)
        except Exception:
            return default
    if v is None:
        return default
    return v


def _filter_match(value: str, wanted: str) -> bool:
    text = str(value or "").upper()
    parts = {p.strip() for p in re.split(r"[;,/ ]+", text) if p.strip()}
    return wanted.upper() in parts or wanted.upper() in text


def product_rank(row) -> tuple:
    filename = str(_value(row, "productFilename", "")).lower()
    subgroup = str(_value(row, "productSubGroupDescription", "")).upper()
    group = str(_value(row, "productGroupDescription", "")).lower()
    ptype = str(_value(row, "productType", "")).upper()
    data_uri = str(_value(row, "dataURI", "")).lower()
    try:
        calib = int(_value(row, "calib_level", 0) or 0)
    except Exception:
        calib = 0
    try:
        size = int(_value(row, "size", 0) or 0)
    except Exception:
        size = 0

    subgroup_score = {
        "DRZ": 7,
        "DRC": 7,
        "SCI": 6,
        "FLC": 5,
        "FLT": 4,
        "CAL": 3,
        "RAW": 1,
    }.get(subgroup, 0)
    calibrated_name = 1 if any(x in filename for x in ("_drz", "_drc", "_sci", "mosaic")) else 0
    science = 1 if ptype == "SCIENCE" else 0
    mrp = 1 if "minimum recommended" in group else 0
    fits = 1 if filename.endswith((".fits", ".fits.gz")) or data_uri.endswith((".fits", ".fits.gz")) else 0
    public = 1 if str(_value(row, "dataRights", "PUBLIC")).upper() == "PUBLIC" else 0
    return (public, science, subgroup_score, calibrated_name, calib, mrp, fits, size, filename)


def select_products(products: Iterable, filters: Iterable[str] = DEFAULT_FILTERS) -> list[dict]:
    rows = list(products)
    selected = []
    for filt in filters:
        candidates = []
        for row in rows:
            if not _filter_match(_value(row, "filters", ""), filt):
                continue
            filename = str(_value(row, "productFilename", ""))
            data_uri = str(_value(row, "dataURI", ""))
            if not (filename.lower().endswith((".fits", ".fits.gz")) or data_uri.lower().endswith((".fits", ".fits.gz"))):
                continue
            candidates.append(row)
        if not candidates:
            continue
        row = max(candidates, key=product_rank)
        selected.append({
            "filter": filt,
            "band": {"F435W": "B", "F555W": "V", "F606W": "R", "F814W": "I"}.get(filt, filt),
            "dataURI": str(_value(row, "dataURI", "")),
            "productFilename": str(_value(row, "productFilename", "")),
            "productType": str(_value(row, "productType", "")),
            "productSubGroupDescription": str(_value(row, "productSubGroupDescription", "")),
            "obs_collection": str(_value(row, "obs_collection", "")),
            "obs_id": str(_value(row, "obs_id", "")),
            "obsid": str(_value(row, "obsid", _value(row, "obsID", ""))),
            "calib_level": int(_value(row, "calib_level", 0) or 0),
            "size": int(_value(row, "size", 0) or 0),
        })
    return selected


def _query_observations(observations, query_name: str, cfg: DiscoveryConfig):
    return observations.query_criteria(
        objectname=query_name,
        radius=f"{cfg.radius_deg} deg",
        obs_collection="HST",
        dataproduct_type="image",
        intentType="science",
        dataRights="PUBLIC",
    )


def discover_target(target: str, cfg: DiscoveryConfig | None = None, observations=None) -> tuple[list[dict], dict]:
    cfg = cfg or DiscoveryConfig()
    if observations is None:
        from astroquery.mast import Observations
        observations = Observations

    query_name = canonical_name(target)
    last_error = None
    for attempt in range(1, cfg.retries + 1):
        try:
            obs = _query_observations(observations, query_name, cfg)
            if len(obs) == 0:
                return [], {
                    "target": target,
                    "query_name": query_name,
                    "status": "NO_OBSERVATIONS",
                    "attempts": attempt,
                    "observation_count": 0,
                }
            products = observations.get_unique_product_list(obs)
            selected = select_products(products, cfg.filters)
            return selected, {
                "target": target,
                "query_name": query_name,
                "status": "OK" if len(selected) >= 2 else "INSUFFICIENT_FILTERS",
                "attempts": attempt,
                "observation_count": int(len(obs)),
                "product_count": int(len(products)),
                "selected_filter_count": int(len(selected)),
                "selected_filters": [x["filter"] for x in selected],
            }
        except Exception as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            if attempt < cfg.retries:
                time.sleep(cfg.retry_sleep_seconds * attempt)
    return [], {
        "target": target,
        "query_name": query_name,
        "status": "MAST_QUERY_ERROR",
        "attempts": cfg.retries,
        "error": last_error,
    }


def build_manifest(targets: Iterable[str], *, cfg: DiscoveryConfig | None = None, observations=None) -> dict:
    cfg = cfg or DiscoveryConfig()
    unique_targets = list(dict.fromkeys(targets))
    target_rows = []
    statuses = []
    for target in unique_targets:
        selected, status = discover_target(target, cfg=cfg, observations=observations)
        statuses.append(status)
        if len(selected) >= 2:
            target_rows.append({
                "target": target,
                "query_name": canonical_name(target),
                "class": "galaxy",
                "focus": target in FOCUS_TARGETS,
                "filters": selected,
            })

    return {
        "schema": "janus.cosmos.hst.live_manifest.v1",
        "status": "LIVE_MAST_DISCOVERY",
        "source": "MAST / STScI",
        "selection": "Deterministic HST public science FITS selection; at least two requested filters required.",
        "filters": list(cfg.filters),
        "requested_targets": len(unique_targets),
        "selected_targets": len(target_rows),
        "selected_products": sum(len(x["filters"]) for x in target_rows),
        "targets": target_rows,
        "target_status": statuses,
        "blind_restrictions": {
            "ocr": False,
            "face_search": False,
            "semantic_analysis": False,
            "cipher_search": False,
            "post_hoc_tuning": False,
        },
    }


def write_manifest(path: str | Path, manifest: dict) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
