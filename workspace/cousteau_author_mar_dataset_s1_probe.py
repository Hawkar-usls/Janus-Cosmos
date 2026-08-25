#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import io
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

DOI = "10.1029/2022JB024008"
FILENAME = "2022JB024008-sup-0003-Data Set SI-S01.zip"
URLS = [
    "https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2022JB024008&file=2022JB024008-sup-0003-Data+Set+SI-S01.zip",
    "https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2022JB024008&file=2022JB024008-sup-0003-Data+Set+SI-S01.zip",
    "https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2022JB024008&file=2022JB024008-sup-0003-Data%20Set%20SI-S01.zip",
]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def response_record(response: requests.Response) -> dict[str, Any]:
    return {
        "requested_url": response.request.url,
        "final_url": response.url,
        "status": int(response.status_code),
        "content_type": response.headers.get("content-type"),
        "content_disposition": response.headers.get("content-disposition"),
        "bytes": len(response.content),
        "sha256": sha256_bytes(response.content),
    }


def inspect_zip(payload: bytes) -> list[dict[str, Any]]:
    members: list[dict[str, Any]] = []
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        infos = archive.infolist()
        names = [info.filename for info in infos]
        if len(names) != len(set(names)):
            raise ValueError("DUPLICATE_ZIP_MEMBER_NAMES")
        for info in infos:
            if info.is_dir():
                continue
            raw = archive.read(info)
            member: dict[str, Any] = {
                "name": info.filename,
                "bytes": len(raw),
                "sha256": sha256_bytes(raw),
            }
            if len(raw) < 2_000_000:
                text = raw.decode("utf-8", errors="replace")
                lines = text.splitlines()
                member["line_count"] = len(lines)
                member["first_30_lines"] = lines[:30]
                member["last_5_lines"] = lines[-5:]
            members.append(member)
    return members


def run(output: Path) -> dict[str, Any]:
    report: dict[str, Any] = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-AUTHOR-MAR-DATASET-S1-ACQUISITION-PROBE-2026-08-22-v1.1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "paper_doi": DOI,
        "expected_filename": FILENAME,
        "responses": [],
        "status": "BLOCKED",
        "claim_ceiling": "Acquisition/provenance probe only. SUCCESS means one configured publisher URL returned a readable ZIP and the member inventory was recorded. It does not establish scientific replication, event-count reconciliation, waveform validity, or any target claim.",
    }

    session = requests.Session()
    session.headers["User-Agent"] = (
        "Mozilla/5.0 Janus-Echo-Cousteau/1.1 independent MAR Data Set S1 replication"
    )
    session.headers["Referer"] = (
        "https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022JB024008"
    )

    for url in URLS:
        try:
            response = session.get(url, timeout=90, allow_redirects=True)
            item = response_record(response)
            report["responses"].append(item)
            if response.status_code != 200 or len(response.content) < 1000:
                continue
            try:
                members = inspect_zip(response.content)
            except Exception as exc:
                item["zip_error"] = f"{type(exc).__name__}:{exc}"
                item["prefix_hex"] = response.content[:32].hex()
                continue
            item["zip_members"] = members
            report["status"] = "SUCCESS"
            report["selected_download"] = {
                "url": response.url,
                "bytes": len(response.content),
                "sha256": sha256_bytes(response.content),
                "members": [member["name"] for member in members],
            }
            break
        except Exception as exc:
            report["responses"].append(
                {
                    "requested_url": url,
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                }
            )

    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"NO_OVERWRITE:{output}")
    output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    report = run(args.output)
    print(
        json.dumps(
            {"status": report["status"], "selected": report.get("selected_download")},
            indent=2,
            ensure_ascii=False,
        )
    )
    return 0 if report["status"] == "SUCCESS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
