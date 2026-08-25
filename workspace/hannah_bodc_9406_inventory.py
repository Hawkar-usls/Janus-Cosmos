#!/usr/bin/env python3
"""Read-only inventory of real BODCREQ-9406 / CD169_TOBI archive.

Uses conventional anonymous FTP with a generic invalid contact address; no
requester credential is required. Emits only file metadata, hashes for small
metadata files, and keyword-derived README hints. Raw scientific data are not
uploaded or committed by this script.
"""
from __future__ import annotations

import argparse
import ftplib
import hashlib
import io
import json
import posixpath
import re
from datetime import datetime, timezone
from pathlib import Path

HOST = "livftp.noc.ac.uk"
USER = "anonymous"
PASSWORD = "janus-probe@example.invalid"
ROOT = "/bodc/bodc/data/BODCREQ-9406/CD169_TOBI"
MAX_DEPTH = 12
MAX_FILES = 20000
META_MAX_BYTES = 2_000_000
META_NAME_RE = re.compile(r"(^|[-_.])(readme|info|format|index|manifest|notes?|log)([-_.]|$)|\.(txt|md|csv|lst)$", re.I)
KEYWORDS = re.compile(r"JD0?59|line\s*19|TOBI[-_ ]?02|bath|nav|magnet|chirp|sidescan|side.scan|PRISM|MAPR|altitude|layback", re.I)


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def ftp_bytes(ftp: ftplib.FTP, path: str) -> bytes:
    buf = io.BytesIO()
    ftp.retrbinary(f"RETR {path}", buf.write)
    return buf.getvalue()


def list_dir(ftp: ftplib.FTP, path: str):
    entries = []
    try:
        for name, facts in ftp.mlsd(path, facts=["type", "size", "modify"]):
            if name in {".", ".."}:
                continue
            entries.append((name, facts))
        return entries, "MLSD"
    except Exception:
        # Fallback: NLST plus CWD/SIZE tests.
        names = ftp.nlst(path)
        base = path.rstrip("/")
        for raw in names:
            name = posixpath.basename(raw.rstrip("/"))
            p = raw if raw.startswith("/") else posixpath.join(base, name)
            typ = "file"
            size = None
            current = ftp.pwd()
            try:
                ftp.cwd(p)
                typ = "dir"
            except Exception:
                try:
                    size = ftp.size(p)
                except Exception:
                    pass
            finally:
                try:
                    ftp.cwd(current)
                except Exception:
                    pass
            entries.append((name, {"type": typ, "size": str(size) if size is not None else None}))
        return entries, "NLST_FALLBACK"


def normalize_modify(s: str | None) -> str | None:
    if not s:
        return None
    try:
        dt = datetime.strptime(s, "%Y%m%d%H%M%S").replace(tzinfo=timezone.utc)
        return dt.isoformat()
    except Exception:
        return s


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    result = {
        "schema": "janus.cosmos.cousteau.hannah_bodc.cd169_tobi_inventory.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {"host": HOST, "root": ROOT, "auth": "ANONYMOUS_GENERIC_NON_REQUESTER"},
        "status": "STARTED",
        "raw_scientific_bytes_redistributed": False,
        "entries": [],
        "metadata_hints": [],
        "limits": {"max_depth": MAX_DEPTH, "max_files": MAX_FILES, "metadata_max_bytes": META_MAX_BYTES},
    }

    ftp = ftplib.FTP(timeout=30)
    try:
        ftp.connect(HOST, 21, timeout=30)
        ftp.login(USER, PASSWORD)
        ftp.cwd(ROOT)
        resolved = ftp.pwd()
        result["source"]["resolved_root"] = resolved
        stack = [(resolved, 0)]
        count = 0
        listing_methods = set()
        while stack and count < MAX_FILES:
            directory, depth = stack.pop()
            if depth > MAX_DEPTH:
                continue
            try:
                children, method = list_dir(ftp, directory)
                listing_methods.add(method)
            except Exception as exc:
                result["entries"].append({"path": directory, "type": "ERROR", "error": str(exc)})
                continue
            for name, facts in children:
                path = posixpath.join(directory.rstrip("/"), name)
                typ = facts.get("type") or "unknown"
                size = None
                try:
                    size = int(facts.get("size")) if facts.get("size") not in {None, "None"} else None
                except Exception:
                    pass
                entry = {
                    "path": path,
                    "relative_path": posixpath.relpath(path, resolved),
                    "type": typ,
                    "size_bytes": size,
                    "modified_utc": normalize_modify(facts.get("modify")),
                }
                result["entries"].append(entry)
                count += 1
                if typ in {"dir", "cdir", "pdir"}:
                    if typ == "dir":
                        stack.append((path, depth + 1))
                    continue
                if typ != "file":
                    continue
                if size is not None and size <= META_MAX_BYTES and META_NAME_RE.search(name):
                    try:
                        raw = ftp_bytes(ftp, path)
                        text = raw.decode("utf-8", errors="replace")
                        hits = []
                        for lineno, line in enumerate(text.splitlines(), start=1):
                            if KEYWORDS.search(line):
                                # Keep only compact derived line hints; cap length/count.
                                hits.append({"line": lineno, "hint": line.strip()[:300]})
                                if len(hits) >= 50:
                                    break
                        result["metadata_hints"].append({
                            "relative_path": entry["relative_path"],
                            "size_bytes": len(raw),
                            "sha256": sha256_bytes(raw),
                            "keyword_hit_count_capped": len(hits),
                            "keyword_hints": hits,
                        })
                    except Exception as exc:
                        result["metadata_hints"].append({"relative_path": entry["relative_path"], "error": str(exc)})
        result["listing_methods"] = sorted(listing_methods)
        result["entry_count"] = len(result["entries"])
        result["truncated"] = count >= MAX_FILES
        result["status"] = "REAL_ARCHIVE_INVENTORY_READY"
    except Exception as exc:
        result["status"] = "FTP_INVENTORY_FAILED"
        result["error_type"] = type(exc).__name__
        result["error"] = str(exc)
    finally:
        try:
            ftp.quit()
        except Exception:
            try:
                ftp.close()
            except Exception:
                pass

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result.get("status"),
        "entry_count": result.get("entry_count"),
        "metadata_hint_files": len(result.get("metadata_hints", [])),
        "raw_scientific_bytes_redistributed": False,
    }, indent=2))
    return 0 if result.get("status") == "REAL_ARCHIVE_INVENTORY_READY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
