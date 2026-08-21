#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import io
import json
import re
import tarfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin, urlparse

import requests
from bs4 import BeautifulSoup

from cousteau_ea_tphase_blind_cluster import (
    EXPECTED_EVENT_COUNT,
    DBSCAN_GRID,
    NULL_DOMAIN,
    NULL_SAMPLES,
    NULL_SEED,
    blind_cluster,
    reveal_and_score,
    parse_catalog,
    sha256_bytes,
)

DATASET_UID = "30497"
FILE_UID = "2504732"
DOI = "10.26022/IEDA/330497"
LANDING = f"https://www.marine-geo.org/tools/files/{DATASET_UID}"
UID_URL = f"https://www.marine-geo.org/tools/search/file_uids.php?data_set_uid={DATASET_UID}"
MODAL_URL = "https://www.marine-geo.org/services/download/download_modal.php"
ACCEPT_URL = "https://api.marine-geo.org/services/download/download_accept.php"
ANCHOR_LAT = -3.865418
ANCHOR_LON = 3.854924
MIN_COORDS = 1000
MAX_FOLLOW = 40


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def response_meta(r: requests.Response) -> dict:
    return {
        "requested_url": r.request.url if r.request else None,
        "final_url": r.url,
        "status": r.status_code,
        "bytes": len(r.content),
        "content_type": r.headers.get("content-type"),
        "content_disposition": r.headers.get("content-disposition"),
        "sha256": sha256_bytes(r.content),
    }


def looks_text(raw: bytes) -> bool:
    if not raw:
        return False
    sample = raw[:8192]
    if sample.count(b"\x00") > max(2, len(sample) // 100):
        return False
    printable = sum((b in b"\t\n\r") or (32 <= b <= 126) or b >= 128 for b in sample)
    return printable / max(1, len(sample)) >= 0.72


def parse_candidate(raw: bytes, label: str, trace: list[dict]):
    if len(raw) < 1000 or not looks_text(raw):
        return None
    low = raw[:500].lower()
    if b"<html" in low or b"<!doctype" in low:
        return None
    try:
        coords, meta = parse_catalog(raw)
    except Exception as exc:
        trace.append({
            "stage": "parse_candidate",
            "label": label,
            "bytes": len(raw),
            "sha256": sha256_bytes(raw),
            "error": f"{type(exc).__name__}: {exc}",
        })
        return None
    row = {
        "stage": "parse_candidate",
        "label": label,
        "bytes": len(raw),
        "sha256": sha256_bytes(raw),
        "valid_coordinates": int(len(coords)),
        "expected_event_count": EXPECTED_EVENT_COUNT,
        "expected_count_exact_match": int(len(coords)) == EXPECTED_EVENT_COUNT,
    }
    trace.append(row)
    if len(coords) >= MIN_COORDS:
        return raw, coords, meta, label
    return None


def unpack_candidates(raw: bytes, label: str, trace: list[dict]):
    """Yield direct or archive-member payloads without changing scientific parsing rules."""
    yielded: set[str] = set()

    def emit(name: str, b: bytes):
        h = sha256_bytes(b)
        if h in yielded:
            return None
        yielded.add(h)
        return name, b

    direct = emit(label, raw)
    if direct:
        yield direct

    # ZIP
    if raw.startswith(b"PK\x03\x04"):
        try:
            with zipfile.ZipFile(io.BytesIO(raw)) as zf:
                members = [x for x in zf.infolist() if not x.is_dir()]
                trace.append({"stage": "archive_detected", "format": "zip", "label": label,
                              "members": [{"name": x.filename, "bytes": x.file_size} for x in members[:200]]})
                for info in members:
                    if info.file_size <= 0 or info.file_size > 20_000_000:
                        continue
                    b = zf.read(info)
                    item = emit(f"{label}::zip::{info.filename}", b)
                    if item:
                        yield item
        except Exception as exc:
            trace.append({"stage": "archive_error", "format": "zip", "label": label,
                          "error": f"{type(exc).__name__}: {exc}"})

    # TAR, including tar.gz/tgz when tarfile can identify it directly.
    try:
        with tarfile.open(fileobj=io.BytesIO(raw), mode="r:*") as tf:
            members = [m for m in tf.getmembers() if m.isfile()]
            trace.append({"stage": "archive_detected", "format": "tar", "label": label,
                          "members": [{"name": m.name, "bytes": m.size} for m in members[:200]]})
            for m in members:
                if m.size <= 0 or m.size > 20_000_000:
                    continue
                f = tf.extractfile(m)
                if f is None:
                    continue
                b = f.read()
                item = emit(f"{label}::tar::{m.name}", b)
                if item:
                    yield item
    except (tarfile.ReadError, EOFError):
        pass
    except Exception as exc:
        trace.append({"stage": "archive_error", "format": "tar", "label": label,
                      "error": f"{type(exc).__name__}: {exc}"})

    # Standalone gzip that is not a tar archive.
    if raw.startswith(b"\x1f\x8b"):
        try:
            b = gzip.decompress(raw)
            trace.append({"stage": "archive_detected", "format": "gzip", "label": label,
                          "uncompressed_bytes": len(b), "uncompressed_sha256": sha256_bytes(b)})
            item = emit(f"{label}::gzip", b)
            if item:
                yield item
            # gzip may contain a tar stream that tarfile could not see in the original wrapper.
            try:
                with tarfile.open(fileobj=io.BytesIO(b), mode="r:*") as tf:
                    for m in tf.getmembers():
                        if not m.isfile() or m.size <= 0 or m.size > 20_000_000:
                            continue
                        f = tf.extractfile(m)
                        if f is None:
                            continue
                        mb = f.read()
                        item = emit(f"{label}::gzip-tar::{m.name}", mb)
                        if item:
                            yield item
            except Exception:
                pass
        except Exception as exc:
            trace.append({"stage": "archive_error", "format": "gzip", "label": label,
                          "error": f"{type(exc).__name__}: {exc}"})


def allowed_url(url: str) -> bool:
    try:
        host = (urlparse(url).hostname or "").lower()
    except Exception:
        return False
    return host == "marine-geo.org" or host.endswith(".marine-geo.org")


def discover_links(r: requests.Response) -> list[str]:
    ctype = (r.headers.get("content-type") or "").lower()
    text = r.text if ("text" in ctype or "html" in ctype or "json" in ctype or len(r.content) < 1_000_000) else ""
    out: list[str] = []
    if text:
        try:
            soup = BeautifulSoup(text, "html.parser")
            for tag, attr in [("a", "href"), ("form", "action"), ("iframe", "src")]:
                for el in soup.find_all(tag):
                    u = el.get(attr)
                    if u:
                        out.append(urljoin(r.url, u.replace("&amp;", "&")))
        except Exception:
            pass
        for u in re.findall(r'https?://[^"\'<>\s]+', text, flags=re.I):
            out.append(u.replace("&amp;", "&"))
        # Also inspect JSON string values recursively.
        if "json" in ctype:
            try:
                obj = r.json()
                stack = [obj]
                while stack:
                    x = stack.pop()
                    if isinstance(x, dict):
                        stack.extend(x.values())
                    elif isinstance(x, list):
                        stack.extend(x)
                    elif isinstance(x, str) and (x.startswith("http://") or x.startswith("https://")):
                        out.append(x)
            except Exception:
                pass
    return list(dict.fromkeys(u for u in out if allowed_url(u)))


def try_response(r: requests.Response, label: str, trace: list[dict]):
    trace.append({"stage": "response", "label": label, **response_meta(r)})
    for cand_label, payload in unpack_candidates(r.content, label, trace):
        got = parse_candidate(payload, cand_label, trace)
        if got:
            inner_raw, coords, meta, source_label = got
            return {
                "archive_raw": r.content,
                "catalog_raw": inner_raw,
                "coords": coords,
                "parse_meta": meta,
                "source_label": source_label,
                "download_meta": response_meta(r),
            }
    return None


def acquire():
    s = requests.Session()
    s.headers.update({
        "User-Agent": "Janus-Echo-Cousteau/1.4 archive-aware scientific reproducibility audit",
        "Referer": LANDING,
    })
    trace: list[dict] = []

    # PREBIND: authoritative website dataset -> file UID, before catalog bytes are read.
    r = s.get(UID_URL, timeout=45)
    r.raise_for_status()
    uids = [str(x) for x in r.json()]
    trace.append({"stage": "authoritative_prebind", **response_meta(r), "dataset_uid": DATASET_UID, "file_uids": uids})
    if FILE_UID not in uids:
        raise RuntimeError(f"authoritative file UID drift: expected {FILE_UID}, got {uids}")

    # Verify the site's own download modal and its accept endpoint contract.
    modal = s.post(MODAL_URL, data={"FileDownload": FILE_UID, "data_set_uid": DATASET_UID}, timeout=60, allow_redirects=True)
    modal.raise_for_status()
    trace.append({"stage": "modal", **response_meta(modal), "body_prefix": modal.text[:5000]})
    soup = BeautifulSoup(modal.text, "html.parser")
    form = soup.find("form", id="data_link") or soup.find("form")
    if form is None:
        raise RuntimeError("MGDS download modal no longer exposes a download form")
    action = urljoin(modal.url, form.get("action") or ACCEPT_URL)
    if not allowed_url(action):
        raise RuntimeError(f"refusing non-MGDS download action: {action}")

    # Critical v5 repair: include the radio field required by the website UI.
    payload = {
        "purpose": "Research",
        "client": "DataLink",
        "force_download": "1",
        "data_uids": FILE_UID,
    }
    trace.append({"stage": "accept_contract", "action": action, "payload": payload,
                  "note": "purpose=Research reproduces the mandatory UI radio selection; this changes acquisition only, not scientific thresholds"})
    accepted = s.post(action, data=payload, timeout=180, allow_redirects=True)
    accepted.raise_for_status()
    got = try_response(accepted, "download_accept", trace)
    if got:
        return got, trace

    # If the accept endpoint returns a landing page/JSON pointer, follow only marine-geo links.
    queue = discover_links(accepted)
    seen = {accepted.url}
    followed = 0
    while queue and followed < MAX_FOLLOW:
        u = queue.pop(0)
        if u in seen:
            continue
        seen.add(u)
        low = u.lower()
        # Prefer file/archive/download routes; ignore decorative site links.
        if not any(k in low for k in ["download", "archive", "file", "data", "2504732"]):
            continue
        followed += 1
        try:
            rr = s.get(u, timeout=180, allow_redirects=True)
            trace.append({"stage": "follow", "source": accepted.url, "requested": u, **response_meta(rr)})
            if rr.status_code != 200:
                continue
            got = try_response(rr, f"follow:{u}", trace)
            if got:
                return got, trace
            queue.extend(v for v in discover_links(rr) if v not in seen)
        except Exception as exc:
            trace.append({"stage": "follow_error", "requested": u, "error": f"{type(exc).__name__}: {exc}"})

    raise RuntimeError("MGDS v5 archive-aware accept flow yielded no parseable >=1000-coordinate catalog")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--status-output", required=True)
    args = ap.parse_args()
    out = Path(args.output)
    status_path = Path(args.status_output)
    out.parent.mkdir(parents=True, exist_ok=True)
    artifact = "JANUS-ECHO-COUSTEAU-EA-TPHASE-BLIND-CLUSTER-RUN-005-2026-08-21-v1.4"
    started = utcnow()

    try:
        acq, trace = acquire()
        coords = acq["coords"]

        # BLIND PHASE: exact original frozen DBSCAN/grid/null configuration. Anchor absent here.
        blind = blind_cluster(coords)
        blind_hash = blind["freeze_sha256"]

        # REVEAL only after blind JSON state is frozen and hashed.
        reveal = reveal_and_score(coords, blind, ANCHOR_LAT, ANCHOR_LON)
        nearest = [x["nearest_cluster"]["anchor_to_center_km"] for x in reveal["configs"] if x.get("nearest_cluster")]
        p95 = any(x.get("nearest_cluster") and x["nearest_cluster"]["anchor_inside_cluster_p95_radius"] for x in reveal["configs"])
        maxr = any(x.get("nearest_cluster") and x["nearest_cluster"]["anchor_inside_cluster_max_radius"] for x in reveal["configs"])
        verdict = "ANCHOR_OVERLAPS_BLIND_CLUSTER_P95__TECTONIC_CONTROL_REQUIRED" if p95 else "NO_BLIND_CLUSTER_P95_OVERLAP_WITH_FROZEN_ANCHOR"

        result = {
            "artifact_id": artifact,
            "research_branch": "Janus-Echo-Кусто",
            "started_at_utc": started,
            "completed_at_utc": utcnow(),
            "source": {
                "doi": DOI,
                "dataset": "EA_Hydroacoustics",
                "data_set_uid": int(DATASET_UID),
                "file_uid": FILE_UID,
                "prebind_endpoint": UID_URL,
                "download_modal_endpoint": MODAL_URL,
                "download_accept_endpoint": action if False else ACCEPT_URL,
                "acquisition_repair": "PURPOSE_RESEARCH_PLUS_ARCHIVE_UNPACK",
                "archive_sha256": sha256_bytes(acq["archive_raw"]),
                "archive_bytes": len(acq["archive_raw"]),
                "catalog_member_sha256": sha256_bytes(acq["catalog_raw"]),
                "catalog_member_bytes": len(acq["catalog_raw"]),
                "catalog_member_label": acq["source_label"],
                "raw_committed": False,
                "expected_event_count": EXPECTED_EVENT_COUNT,
                "parsed_valid_coordinate_count": int(len(coords)),
                "expected_count_exact_match": int(len(coords)) == EXPECTED_EVENT_COUNT,
                "parse": acq["parse_meta"],
                "download": acq["download_meta"],
                "acquisition_trace": trace,
            },
            "preregistration": {
                "anchor_hidden_during_clustering": True,
                "clustering_parameters_frozen_before_anchor_reveal": True,
                "dbscan_grid": DBSCAN_GRID,
                "null_domain": NULL_DOMAIN,
                "null_samples": NULL_SAMPLES,
                "null_seed": NULL_SEED,
                "parameter_change_from_runs_001_004": False,
            },
            "blind_phase": blind,
            "post_reveal": reveal,
            "summary": {
                "blind_freeze_sha256": blind_hash,
                "parsed_count": int(len(coords)),
                "expected_count_exact_match": int(len(coords)) == EXPECTED_EVENT_COUNT,
                "nearest_catalog_event_to_anchor_km": reveal["nearest_event"]["distance_km"],
                "nearest_blind_cluster_center_across_grid_km": round(min(nearest), 3) if nearest else None,
                "anchor_inside_any_blind_cluster_p95_radius": p95,
                "anchor_inside_any_blind_cluster_max_radius": maxr,
                "verdict": verdict,
                "semantic_status": "UNCONFIRMED",
                "tectonic_control_required": True,
            },
            "hard_rules": [
                "AUTHORITATIVE_DATASET_AND_FILE_UID_PREBOUND_BEFORE_CATALOG_READ",
                "PURPOSE_FIELD_REPAIRS_ACQUISITION_ONLY",
                "ARCHIVE_UNPACK_REPAIRS_TRANSPORT_ONLY",
                "BLIND_CLUSTER_BEFORE_ANCHOR_REVEAL",
                "NO_PARAMETER_RETUNING_AFTER_REVEAL",
                "MID_ATLANTIC_RIDGE_SEISMICITY_IS_MANDATORY_TECTONIC_CONTROL",
                "RECTANGULAR_LOOK_ELSEWHERE_NULL_IS_DIAGNOSTIC_NOT_FORMAL",
                "DISTANCE_IS_NOT_CAUSATION",
                "NO_RECENTERING",
                "NO_UNDERWATER_PYRAMID_DETECTED_YET",
            ],
            "status": "BLIND_CLUSTER_RUN_COMPLETE",
        }
        out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        status_path.write_text(json.dumps({
            "artifact_id": artifact,
            "status": "SUCCESS",
            "completed_at_utc": utcnow(),
            "parsed_count": int(len(coords)),
            "expected_count_exact_match": int(len(coords)) == EXPECTED_EVENT_COUNT,
            "archive_sha256": sha256_bytes(acq["archive_raw"]),
            "catalog_member_sha256": sha256_bytes(acq["catalog_raw"]),
            "blind_freeze_sha256": blind_hash,
            "verdict": verdict,
            "result_path": str(out),
        }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(json.dumps(result["summary"], indent=2, ensure_ascii=False))
        return 0
    except Exception as exc:
        payload = {
            "artifact_id": artifact,
            "status": "BLOCKED_DATA_ACQUISITION_OR_PARSE",
            "started_at_utc": started,
            "completed_at_utc": utcnow(),
            "error_type": type(exc).__name__,
            "error": str(exc),
            "scientific_interpretation": "ACQUISITION_BLOCKER_ONLY__NOT_A_NEGATIVE_CLUSTER_RESULT",
        }
        if "trace" in locals():
            payload["acquisition_trace"] = trace
        status_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
