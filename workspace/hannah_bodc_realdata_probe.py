#!/usr/bin/env python3
"""Ephemeral real-data probe for Hannah/BODC Cousteau lane.

Downloads BODC request data to a temporary GitHub Actions runner, hashes original
bytes, derives edge-window measurement features, and passes those features to
the frozen Cousteau Synesthetic Memory Core. Raw BODC files are NOT committed
or uploaded as artifacts; only derived JSON receipts/passports are emitted.

Scientific boundary:
    sensory similarity -> retrieval priority only
    sensory similarity != scientific convergence
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import json
import math
import os
import posixpath
import statistics
import sys
import tempfile
from collections import deque
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Iterable

try:
    import paramiko
except Exception as exc:  # pragma: no cover
    raise SystemExit(f"paramiko is required for SFTP probe: {exc}")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cousteau_synesthetic_memory_core import build_passport, compare_passports  # noqa: E402
from cousteau_synesthetic_semantic_overlay import enrich_passport  # noqa: E402

HOST = "livftp.noc.ac.uk"
PORT = 22
USERNAME = "anonymous"
REMOTE_9408 = "/bodc/bodc/data/BODCREQ-9408"
REMOTE_9406 = "/bodc/bodc/data/BODCREQ-9406"
TARGET_NAMES = {
    "em122.aco",
    "em122.tpl",
    "seatex-gga.aco",
    "seatex-gga.tpl",
    "gyro.aco",
    "gyro.tpl",
    "ea600.aco",
    "ea600.tpl",
    "furuno-vtg.aco",
    "furuno-vtg.tpl",
}
EDGE_SECONDS = (60, 300, 1800, 7200)
SCALE_BY_SECONDS = {60: "60s", 300: "300s", 1800: "1800s", 7200: "7200s"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def blake2_file(path: Path) -> str:
    h = hashlib.blake2b(digest_size=32)
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sftp_walk(sftp: paramiko.SFTPClient, root: str, max_depth: int = 8) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    stack = [(root, 0)]
    seen: set[str] = set()
    while stack:
        current, depth = stack.pop()
        if current in seen or depth > max_depth:
            continue
        seen.add(current)
        try:
            entries = sftp.listdir_attr(current)
        except Exception as exc:
            out.append({"path": current, "type": "ERROR", "error": str(exc)})
            continue
        for ent in entries:
            p = posixpath.join(current, ent.filename)
            mode = ent.st_mode
            import stat
            if stat.S_ISDIR(mode):
                out.append({"path": p, "type": "dir"})
                stack.append((p, depth + 1))
            else:
                out.append({
                    "path": p,
                    "type": "file",
                    "size_bytes": int(ent.st_size),
                    "mtime_epoch": int(ent.st_mtime),
                })
    return sorted(out, key=lambda x: x.get("path", ""))


def download(sftp: paramiko.SFTPClient, remote: str, local: Path) -> None:
    local.parent.mkdir(parents=True, exist_ok=True)
    sftp.get(remote, str(local))


def split_fields(raw: bytes) -> list[str]:
    text = raw.decode("utf-8", errors="replace").strip()
    if not text:
        return []
    if "," in text:
        return [x.strip() for x in text.split(",")]
    return text.split()


def parse_timestamp(fields: list[str]) -> float | None:
    if len(fields) < 4:
        return None
    try:
        year = int(float(fields[0]))
        julian_with_fraction = float(fields[1])
        if not (1900 <= year <= 2200 and 1.0 <= julian_with_fraction < 367.0):
            return None
        dt = datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=julian_with_fraction - 1.0)
        return dt.timestamp()
    except Exception:
        return None


def parse_tpl(tpl_path: Path | None) -> dict[str, Any]:
    result: dict[str, Any] = {"present": False, "variables": [], "raw_line_count": 0}
    if tpl_path is None or not tpl_path.exists():
        return result
    variables = []
    lines = tpl_path.read_text(encoding="utf-8", errors="replace").splitlines()
    for line in lines:
        s = line.strip()
        if not s or s.startswith(("#", ";")):
            continue
        parts = [p.strip() for p in s.split(",")]
        if len(parts) >= 2:
            variables.append({
                "id": parts[0],
                "name": parts[1],
                "unit": parts[2] if len(parts) >= 3 else None,
            })
    return {"present": True, "variables": variables, "raw_line_count": len(lines)}


def infer_depth_index(sample_fields: list[str], tpl: dict[str, Any]) -> dict[str, Any]:
    ncols = len(sample_fields)
    vars_ = tpl.get("variables") or []
    depth_candidates = [i for i, v in enumerate(vars_) if "depth" in str(v.get("name", "")).lower()]
    if len(depth_candidates) == 1:
        i = depth_candidates[0]
        if len(vars_) == ncols:
            return {"resolved": True, "index": i, "basis": "TPL_DIRECT_COLUMN_ALIGNMENT", "tpl_variable": vars_[i]}
        if len(vars_) == max(0, ncols - 4):
            return {"resolved": True, "index": 4 + i, "basis": "TPL_AFTER_FOUR_COMMON_TIME_COLUMNS", "tpl_variable": vars_[i]}
    if ncols == 5:
        return {"resolved": True, "index": 4, "basis": "SOLE_VALUE_AFTER_FOUR_COMMON_TIME_COLUMNS", "warning": "TPL_NOT_USED_FOR_COLUMN_ID"}
    return {
        "resolved": False,
        "index": None,
        "basis": "UNRESOLVED",
        "ncols": ncols,
        "tpl_variable_count": len(vars_),
        "depth_candidate_count": len(depth_candidates),
    }


def to_float(x: str) -> float | None:
    s = x.strip().lower()
    if s in {"", "nan", "null", "none", "-999", "-9999", "-99999", "99999"}:
        return None
    try:
        v = float(s)
        return v if math.isfinite(v) else None
    except Exception:
        return None


def median_or_none(xs: list[float]) -> float | None:
    return float(statistics.median(xs)) if xs else None


def mad(xs: list[float]) -> float | None:
    if not xs:
        return None
    m = statistics.median(xs)
    return float(statistics.median(abs(x - m) for x in xs))


def max_run(values: Iterable[Any], predicate_equal) -> int:
    best = cur = 0
    prev = object()
    have_prev = False
    for v in values:
        if have_prev and predicate_equal(prev, v):
            cur += 1
        else:
            cur = 1
        best = max(best, cur)
        prev = v
        have_prev = True
    return best


def max_null_run(depths: list[float | None]) -> int:
    best = cur = 0
    for d in depths:
        if d is None:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


def max_identical_valid_run(depths: list[float | None]) -> int:
    best = cur = 0
    prev: float | None = None
    for d in depths:
        if d is None:
            cur = 0
            prev = None
            continue
        if prev is not None and d == prev:
            cur += 1
        else:
            cur = 1
        best = max(best, cur)
        prev = d
    return best


def records_for_window(records: list[tuple[float, float | None, bytes]], seconds: int, side: str) -> list[tuple[float, float | None, bytes]]:
    if not records:
        return []
    if side == "HEAD":
        t0 = records[0][0]
        return [r for r in records if r[0] <= t0 + seconds]
    t1 = records[-1][0]
    return [r for r in records if r[0] >= t1 - seconds]


def window_payload(records: list[tuple[float, float | None, bytes]]) -> tuple[dict[str, Any], bytes, dict[str, Any]]:
    if not records:
        return {}, b"", {"record_count": 0}
    ts = [r[0] for r in records]
    ds = [r[1] for r in records]
    valid = [(t, d) for t, d, _ in records if d is not None]
    vals = [d for _, d in valid]
    diffs_t = [b - a for a, b in zip(ts, ts[1:]) if b > a]
    slopes = []
    if len(valid) >= 2:
        for (t0, d0), (t1, d1) in zip(valid, valid[1:]):
            dt = t1 - t0
            if dt > 0:
                slopes.append((d1 - d0) / dt)
    accelerations = []
    if len(slopes) >= 2 and diffs_t:
        for a, b in zip(slopes, slopes[1:]):
            accelerations.append(b - a)
    med = median_or_none(vals)
    m_mad = mad(vals)
    local_range = (max(vals) - min(vals)) if vals else None
    slope_global = None
    if len(valid) >= 2 and valid[-1][0] > valid[0][0]:
        slope_global = (valid[-1][1] - valid[0][1]) / (valid[-1][0] - valid[0][0])
    cadence = median_or_none(diffs_t)
    cadence_mad = mad(diffs_t)
    outlier_score = None
    if vals and m_mad is not None:
        if m_mad > 0:
            outlier_score = max(abs(x - med) for x in vals) / (1.4826 * m_mad)
        else:
            outlier_score = 0.0 if local_range == 0 else 99.0
    payload = {
        "em122_depth": med,
        "rolling_depth_median": med,
        "rolling_depth_mad": m_mad,
        "depth_local_range": local_range,
        "depth_local_slope": slope_global,
        "delta_depth_dt": median_or_none(slopes),
        "delta2_depth_dt2": median_or_none(accelerations),
        "timestamp_cadence": cadence,
        "cadence_jitter": cadence_mad,
        "missing_fraction": (sum(d is None for d in ds) / len(ds)) if ds else None,
        "null_run_length": max_null_run(ds),
        "identical_value_run_length": max_identical_valid_run(ds),
        "outlier_score": outlier_score,
    }
    payload = {k: v for k, v in payload.items() if v is not None}
    raw_bytes = b"".join(r[2] for r in records)
    meta = {
        "record_count": len(records),
        "valid_depth_count": len(vals),
        "start_utc": datetime.fromtimestamp(ts[0], timezone.utc).isoformat(),
        "end_utc": datetime.fromtimestamp(ts[-1], timezone.utc).isoformat(),
        "duration_seconds": ts[-1] - ts[0],
    }
    return payload, raw_bytes, meta


def parse_em122_edges(path: Path, depth_index: int) -> dict[str, Any]:
    head: list[tuple[float, float | None, bytes]] = []
    tail: deque[tuple[float, float | None, bytes]] = deque()
    first_ts: float | None = None
    last_ts: float | None = None
    valid_rows = 0
    parsed_rows = 0
    malformed_rows = 0
    column_counts: dict[int, int] = {}

    with path.open("rb") as f:
        for raw in f:
            fields = split_fields(raw)
            column_counts[len(fields)] = column_counts.get(len(fields), 0) + 1
            ts = parse_timestamp(fields)
            if ts is None or depth_index >= len(fields):
                malformed_rows += 1
                continue
            depth = to_float(fields[depth_index])
            parsed_rows += 1
            valid_rows += int(depth is not None)
            if first_ts is None:
                first_ts = ts
            last_ts = ts
            rec = (ts, depth, raw)
            if first_ts is not None and ts <= first_ts + max(EDGE_SECONDS):
                head.append(rec)
            tail.append(rec)
            while tail and ts - tail[0][0] > max(EDGE_SECONDS) + 5:
                tail.popleft()

    tail_list = list(tail)
    return {
        "head_records": head,
        "tail_records": tail_list,
        "summary": {
            "parsed_rows": parsed_rows,
            "valid_depth_rows": valid_rows,
            "malformed_or_unparsed_rows": malformed_rows,
            "first_utc": datetime.fromtimestamp(first_ts, timezone.utc).isoformat() if first_ts is not None else None,
            "last_utc": datetime.fromtimestamp(last_ts, timezone.utc).isoformat() if last_ts is not None else None,
            "column_count_histogram": {str(k): v for k, v in sorted(column_counts.items())},
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--password-env", default="BODC_ANON_PASSWORD")
    args = ap.parse_args()
    password = os.environ.get(args.password_env)
    if not password:
        raise SystemExit(f"missing {args.password_env}")

    result: dict[str, Any] = {
        "schema": "janus.cosmos.cousteau.hannah_bodc.realdata_probe.v1",
        "status": "STARTED",
        "source": {
            "host": HOST,
            "port": PORT,
            "username": USERNAME,
            "requests": ["BODCREQ-9408", "BODCREQ-9406"],
            "password_stored_or_emitted": False,
        },
        "scientific_claim": False,
        "raw_files_redistributed": False,
        "hard_rules": [
            "RAW_BYTES_OUTRANK_MNEMONIC",
            "SENSORY_MATCH_IS_RETRIEVAL_ONLY",
            "NO_TIME_TO_SPACE_CLAIM_FROM_ACO_DEPTH_ALONE",
        ],
    }

    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        client.connect(HOST, port=PORT, username=USERNAME, password=password, timeout=20, auth_timeout=20, banner_timeout=20)
        transport = client.get_transport()
        if transport is None:
            raise RuntimeError("no SSH transport")
        key = transport.get_remote_server_key()
        result["source"]["server_host_key"] = {
            "type": key.get_name(),
            "sha256_b64": base64.b64encode(hashlib.sha256(key.asbytes()).digest()).decode("ascii").rstrip("="),
            "verification": "RECORDED_NOT_PREPINNED",
        }
        sftp = client.open_sftp()
        inv9408 = sftp_walk(sftp, REMOTE_9408)
        inv9406 = sftp_walk(sftp, REMOTE_9406)
        result["inventory"] = {"BODCREQ-9408": inv9408, "BODCREQ-9406": inv9406}

        with tempfile.TemporaryDirectory(prefix="janus_bodc_") as td:
            tmp = Path(td)
            files9408 = [x for x in inv9408 if x.get("type") == "file"]
            by_name = {posixpath.basename(x["path"]).lower(): x for x in files9408}
            downloaded: dict[str, dict[str, Any]] = {}
            local_paths: dict[str, Path] = {}
            for name in sorted(TARGET_NAMES):
                info = by_name.get(name)
                if not info:
                    continue
                local = tmp / posixpath.basename(info["path"])
                download(sftp, info["path"], local)
                local_paths[name] = local
                downloaded[name] = {
                    "remote_path": info["path"],
                    "size_bytes": local.stat().st_size,
                    "sha256": sha256_file(local),
                    "blake2b_256": blake2_file(local),
                }
            result["downloaded_originals"] = downloaded

            em = local_paths.get("em122.aco")
            tplp = local_paths.get("em122.tpl")
            if em is None:
                result["status"] = "CONNECTED_BUT_EM122_ACO_NOT_FOUND"
            else:
                # sample one valid-looking line solely for schema inference; raw sample is not emitted.
                sample_fields: list[str] = []
                with em.open("rb") as f:
                    for raw in f:
                        fields = split_fields(raw)
                        if parse_timestamp(fields) is not None:
                            sample_fields = fields
                            break
                tpl = parse_tpl(tplp)
                mapping = infer_depth_index(sample_fields, tpl)
                result["em122_schema"] = {
                    "sample_column_count": len(sample_fields),
                    "tpl": tpl,
                    "depth_mapping": mapping,
                    "raw_sample_emitted": False,
                }
                if not mapping.get("resolved"):
                    result["status"] = "RAW_HASHED_DEPTH_COLUMN_UNRESOLVED"
                else:
                    edge = parse_em122_edges(em, int(mapping["index"]))
                    result["em122_parse_summary"] = edge["summary"]
                    result["edge_windows"] = {}
                    result["head_tail_comparisons"] = {}
                    for sec in EDGE_SECONDS:
                        scale = SCALE_BY_SECONDS[sec]
                        hrec = records_for_window(edge["head_records"], sec, "HEAD")
                        trec = records_for_window(edge["tail_records"], sec, "TAIL")
                        hp, hraw, hmeta = window_payload(hrec)
                        tp, traw, tmeta = window_payload(trec)
                        hpass = enrich_passport(build_passport(
                            hp,
                            direction="HEAD_FORWARD",
                            scale=scale,
                            provenance={"source": "BODCREQ-9408/em122.ACO", "window": "HEAD", "seconds": sec},
                            raw_bytes=hraw,
                        ))
                        tpass = enrich_passport(build_passport(
                            tp,
                            direction="TAIL_REVERSE",
                            scale=scale,
                            provenance={"source": "BODCREQ-9408/em122.ACO", "window": "TAIL", "seconds": sec},
                            raw_bytes=traw,
                        ))
                        result["edge_windows"][scale] = {
                            "HEAD": {"meta": hmeta, "measurement_payload": hp, "passport": hpass},
                            "TAIL": {"meta": tmeta, "measurement_payload": tp, "passport": tpass},
                        }
                        result["head_tail_comparisons"][scale] = compare_passports(hpass, tpass)
                    result["status"] = "REAL_EM122_EDGE_MNEMONICS_READY"

        result["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    except Exception as exc:
        result["status"] = "SFTP_OR_PROBE_FAILED"
        result["error_type"] = type(exc).__name__
        result["error"] = str(exc)
        result["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    finally:
        try:
            client.close()
        except Exception:
            pass

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result.get("status"),
        "output": str(args.output),
        "downloaded": sorted((result.get("downloaded_originals") or {}).keys()),
        "em122_parse_summary": result.get("em122_parse_summary"),
        "scientific_claim": False,
    }, indent=2))
    return 0 if result.get("status") not in {"SFTP_OR_PROBE_FAILED"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
