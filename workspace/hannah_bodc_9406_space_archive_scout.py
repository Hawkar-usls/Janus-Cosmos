#!/usr/bin/env python3
"""Metadata-only scout for BODCREQ-9406 before sample-level space replay.

The scout inventories the request tree and reads short prefixes only from likely
navigation/readme/index files. It does not inspect sonar samples or rank any
measurement feature, preserving the preregistered SPACE REPLAY order.
"""
from __future__ import annotations

import argparse
import ftplib
import hashlib
import json
import posixpath
import re
from datetime import datetime, timezone
from pathlib import Path

HOST = "livftp.noc.ac.uk"
ROOT = "/bodc/bodc/data/BODCREQ-9406"
MAX_ENTRIES = 5000
MAX_DEPTH = 6
PREVIEW_BYTES = 65536

KEYWORDS = (
    "nav", "navi", "position", "track", "tow", "fish", "cable", "layback",
    "readme", "info", "index", "format", "deploy", "log", "bathy", "swath",
    "sonar", "tobi", "cdf", "processed", "mosaic",
)
TEXT_EXTS = {".txt", ".doc", ".readme", ".info", ".lst", ".csv", ".nav", ".log", ".asc"}
BINARY_NAV_EXTS = {".cdf", ".nc", ".mat", ".bin"}


def connect() -> ftplib.FTP:
    f = ftplib.FTP(timeout=60)
    f.connect(HOST, 21)
    f.login("anonymous", "janus-probe@example.invalid")
    f.voidcmd("TYPE I")
    return f


def basename(x: str) -> str:
    return posixpath.basename(x.rstrip("/"))


def ext(name: str) -> str:
    n = name.lower()
    i = n.rfind(".")
    return n[i:] if i >= 0 else ""


def likely_relevant(path: str) -> bool:
    s = path.lower()
    return any(k in s for k in KEYWORDS) or ext(s) in TEXT_EXTS | BINARY_NAV_EXTS


def may_preview(path: str) -> bool:
    s = path.lower()
    return ext(s) in TEXT_EXTS or any(k in basename(s) for k in ("readme", "info", "index", "format", "nav", "track", "tow", "cable", "layback"))


def list_names(f: ftplib.FTP, path: str) -> list[str]:
    names = f.nlst(path)
    out = []
    for n in names:
        if n in {".", "..", path, path + "/"}:
            continue
        if not n.startswith("/"):
            n = posixpath.join(path, n)
        out.append(posixpath.normpath(n))
    return sorted(set(out))


def is_dir(f: ftplib.FTP, path: str) -> bool:
    old = f.pwd()
    try:
        f.cwd(path)
        return True
    except Exception:
        return False
    finally:
        try:
            f.cwd(old)
        except Exception:
            pass


def size_of(f: ftplib.FTP, path: str) -> int | None:
    try:
        v = f.size(path)
        return int(v) if v is not None else None
    except Exception:
        return None


def read_prefix(path: str, n: int = PREVIEW_BYTES) -> bytes:
    f = connect()
    out = bytearray()
    try:
        sock = f.transfercmd("RETR " + path)
        try:
            while len(out) < n:
                b = sock.recv(min(65536, n - len(out)))
                if not b:
                    break
                out.extend(b)
        finally:
            try:
                sock.close()
            except Exception:
                pass
    finally:
        try:
            f.close()
        except Exception:
            pass
    return bytes(out)


def text_preview(raw: bytes) -> dict:
    text = raw.decode("utf-8", errors="replace")
    lines = text.splitlines()
    printable = sum(1 for c in text if c.isprintable() or c in "\r\n\t") / max(1, len(text))
    return {
        "prefix_sha256": hashlib.sha256(raw).hexdigest(),
        "prefix_bytes": len(raw),
        "printable_fraction": round(printable, 6),
        "text_lines": [line[:500] for line in lines[:120]],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    result = {
        "schema": "janus.cosmos.cousteau.hannah_bodc.space_archive_scout.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "root": ROOT,
        "status": "STARTED",
        "sample_data_inspected": False,
        "scientific_claim": False,
        "entries": [],
        "relevant_files": [],
    }
    f = connect()
    try:
        stack = [(ROOT, 0)]
        seen_dirs = set()
        while stack and len(result["entries"]) < MAX_ENTRIES:
            path, depth = stack.pop()
            if path in seen_dirs or depth > MAX_DEPTH:
                continue
            seen_dirs.add(path)
            try:
                children = list_names(f, path)
            except Exception as e:
                result["entries"].append({"path": path, "type": "dir", "list_error": f"{type(e).__name__}: {e}"})
                continue
            for child in children:
                if len(result["entries"]) >= MAX_ENTRIES:
                    break
                d = is_dir(f, child)
                if d:
                    row = {"path": child, "relative_path": posixpath.relpath(child, ROOT), "type": "dir", "depth": depth + 1}
                    result["entries"].append(row)
                    stack.append((child, depth + 1))
                else:
                    sz = size_of(f, child)
                    row = {
                        "path": child,
                        "relative_path": posixpath.relpath(child, ROOT),
                        "type": "file",
                        "size_bytes": sz,
                        "relevant_by_name": likely_relevant(child),
                    }
                    result["entries"].append(row)
                    if row["relevant_by_name"]:
                        rel = dict(row)
                        if may_preview(child) and (sz is None or sz <= 8_000_000):
                            try:
                                rel["preview"] = text_preview(read_prefix(child))
                            except Exception as e:
                                rel["preview_error"] = f"{type(e).__name__}: {e}"
                        result["relevant_files"].append(rel)
        result["entry_count"] = len(result["entries"])
        result["relevant_file_count"] = len(result["relevant_files"])
        result["directory_count"] = sum(1 for x in result["entries"] if x.get("type") == "dir")
        result["file_count"] = sum(1 for x in result["entries"] if x.get("type") == "file")
        result["truncated"] = len(result["entries"]) >= MAX_ENTRIES
        result["status"] = "ARCHIVE_METADATA_SCOUT_READY"
    except Exception as e:
        result["status"] = "ARCHIVE_METADATA_SCOUT_FAILED"
        result["error_type"] = type(e).__name__
        result["error"] = str(e)
    finally:
        try:
            f.quit()
        except Exception:
            try:
                f.close()
            except Exception:
                pass
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "entry_count": result.get("entry_count"),
        "relevant_file_count": result.get("relevant_file_count"),
        "truncated": result.get("truncated"),
    }, indent=2))
    return 0 if result["status"] == "ARCHIVE_METADATA_SCOUT_READY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
