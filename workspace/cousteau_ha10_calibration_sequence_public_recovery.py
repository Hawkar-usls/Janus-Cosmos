#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path
from typing import Any
from urllib.parse import urljoin, urlparse

import requests
from bs4 import BeautifulSoup
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cousteau"
PROTOCOL = DATA / "JANUS-ECHO-COUSTEAU-HA10-CALIBRATION-SEQUENCE-PUBLIC-RECOVERY-PROTOCOL-2026-08-22-v1.0.json"
FROZEN_119 = DATA / "JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0.json"
EXPECTED_PROTOCOL_BLOB = "44d14209974c31c010b1f6a9ba502ea558013a57"
EXPECTED_FROZEN_119_BLOB = "eb8b48fb7f043160c057f9df6264a781412ed854"
EXPECTED_FROZEN_119_VERDICT = "NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE"
UA = "JANUS-Cousteau-Calibration-Recovery/1.0 (+https://github.com/Hawkar-usls/Janus-Cosmos)"

CAL_MARKERS = re.compile(r"calibrat|transfer function|RBC|random broadband code|calibration flag", re.I)
H10_MARKERS = re.compile(r"H10S|H10N|HA10|Ascension", re.I)
RAW_MARKERS = re.compile(r"raw|sequence|sample|waveform|RBC|random broadband code|CW sinus|zero-input|silence", re.I)
ACCESS_MARKERS = re.compile(r"vDEC|request access|request data|contract|authorized|confidentiality", re.I)
DERIVED_MARKERS = re.compile(r"0\.8396|4\.7832|2\.5340|H10S electronic|onset of electrical noise|cross-talk", re.I)
TABLE_MARKERS = re.compile(r"start/end times|start time|end time|calibration flag|flag bit", re.I)
RELEVANT_LINK = re.compile(r"calib|H10|hydro|attachment|material|contribution|presentation|\.pdf$|\.csv$|\.json$|\.txt$|\.zip$|\.dat$|\.seed$|\.mseed$|vdec", re.I)
NUMERIC_TOKEN = re.compile(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def git_blob(path: Path) -> str:
    b = path.read_bytes()
    return hashlib.sha1(f"blob {len(b)}\0".encode() + b).hexdigest()


def verify_contracts() -> tuple[dict[str, Any], dict[str, Any]]:
    if git_blob(PROTOCOL) != EXPECTED_PROTOCOL_BLOB:
        raise RuntimeError("PROTOCOL_BLOB_DRIFT")
    if git_blob(FROZEN_119) != EXPECTED_FROZEN_119_BLOB:
        raise RuntimeError("FROZEN_119_BLOB_DRIFT")
    p = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    f = json.loads(FROZEN_119.read_text(encoding="utf-8"))
    if p.get("status") != "PREREGISTERED_BEFORE_GATE_ACQUISITION":
        raise RuntimeError("PROTOCOL_NOT_PREREGISTERED")
    if p.get("authority", {}).get("authority_delta_for_119hz") != 0:
        raise RuntimeError("119_AUTHORITY_DRIFT")
    if f.get("summary", {}).get("verdict") != EXPECTED_FROZEN_119_VERDICT:
        raise RuntimeError("119_VERDICT_DRIFT")
    return p, f


def host_allowed(url: str, domains: list[str]) -> bool:
    host = (urlparse(url).hostname or "").lower().strip(".")
    if host.startswith("www."):
        host = host[4:]
    for d in domains:
        d = d.lower().strip(".")
        if d.startswith("www."):
            d = d[4:]
        if host == d or host.endswith("." + d):
            return True
    return False


def fetch_bounded(session: requests.Session, url: str, domains: list[str], max_bytes: int) -> dict[str, Any]:
    row: dict[str, Any] = {"requested_url": url, "status": "UNFETCHED"}
    if not host_allowed(url, domains):
        row["status"] = "REJECTED_DOMAIN"
        return row
    if re.search(r"(?:login|signin|oauth|token=|auth=)", url, re.I):
        row["status"] = "REJECTED_AUTH_SURFACE"
        return row
    try:
        with session.get(url, timeout=25, allow_redirects=True, stream=True) as r:
            row["http_status"] = r.status_code
            row["final_url"] = r.url
            if not host_allowed(r.url, domains):
                row["status"] = "REJECTED_REDIRECT_DOMAIN"
                return row
            if r.status_code >= 400:
                row["status"] = "HTTP_ERROR"
                return row
            chunks: list[bytes] = []
            total = 0
            for chunk in r.iter_content(65536):
                if not chunk:
                    continue
                total += len(chunk)
                if total > max_bytes:
                    row["status"] = "REJECTED_TOO_LARGE"
                    row["bytes_seen_before_stop"] = total
                    return row
                chunks.append(chunk)
            data = b"".join(chunks)
            row["status"] = "FETCHED"
            row["bytes"] = len(data)
            row["sha256"] = hashlib.sha256(data).hexdigest()
            row["content_type"] = (r.headers.get("content-type") or "").split(";", 1)[0].lower()
            row["content_disposition"] = r.headers.get("content-disposition") or ""
            row["data"] = data
            return row
    except Exception as exc:
        row["status"] = "FETCH_FAILED"
        row["error"] = f"{type(exc).__name__}:{exc}"[:1000]
        return row


def pdf_text(data: bytes) -> str:
    try:
        reader = PdfReader(io.BytesIO(data))
        return "\n".join((page.extract_text() or "") for page in reader.pages[:80])[:200000]
    except Exception:
        return ""


def decode_text(data: bytes) -> str:
    return data.decode("utf-8", errors="replace")[:200000]


def payload_text(row: dict[str, Any]) -> str:
    data = row.get("data") or b""
    ctype = row.get("content_type") or ""
    final = str(row.get("final_url") or row.get("requested_url") or "").lower()
    if "pdf" in ctype or final.endswith(".pdf"):
        return pdf_text(data)
    if "html" in ctype or final.endswith((".html", "/")):
        decoded = decode_text(data)
        soup = BeautifulSoup(decoded, "html.parser")
        for tag in soup(["script", "style", "noscript"]):
            tag.decompose()
        return soup.get_text("\n", strip=True)[:200000]
    if ctype.startswith("text/") or "json" in ctype or final.endswith((".csv", ".txt", ".json", ".dat")):
        return decode_text(data)
    return ""


def discover_links(row: dict[str, Any], domains: list[str], limit: int) -> list[str]:
    if row.get("status") != "FETCHED":
        return []
    ctype = row.get("content_type") or ""
    if "html" not in ctype:
        return []
    base = str(row.get("final_url") or row.get("requested_url"))
    soup = BeautifulSoup(decode_text(row.get("data") or b""), "html.parser")
    out: list[str] = []
    for a in soup.find_all("a", href=True):
        href = urljoin(base, a["href"])
        label = f"{href} {a.get_text(' ', strip=True)}"
        if not host_allowed(href, domains) or not RELEVANT_LINK.search(label):
            continue
        if re.search(r"(?:login|signin|oauth|token=|auth=)", href, re.I):
            continue
        if href not in out:
            out.append(href)
        if len(out) >= limit:
            break
    return out


def numeric_density(text: str) -> float:
    if not text:
        return 0.0
    nums = NUMERIC_TOKEN.findall(text)
    return len(nums) / max(1, len(text.split()))


def classify(row: dict[str, Any]) -> tuple[str, list[str]]:
    if row.get("status") != "FETCHED":
        return "FETCH_FAILED_OR_REJECTED", []
    text = payload_text(row)
    final = str(row.get("final_url") or row.get("requested_url") or "")
    ctype = row.get("content_type") or ""
    lower = final.lower()
    evidence: list[str] = []
    if CAL_MARKERS.search(text): evidence.append("CALIBRATION_MARKER")
    if H10_MARKERS.search(text): evidence.append("H10_MARKER")
    if DERIVED_MARKERS.search(text): evidence.append("H10_DERIVED_RESULT_MARKER")
    if ACCESS_MARKERS.search(text): evidence.append("CONTROLLED_ACCESS_MARKER")
    if TABLE_MARKERS.search(text): evidence.append("CALIBRATION_TABLE_DESCRIPTION_MARKER")

    archive_members: list[str] = []
    if lower.endswith(".zip") or "zip" in ctype:
        try:
            with zipfile.ZipFile(io.BytesIO(row["data"])) as zf:
                archive_members = zf.namelist()
            row["archive_members"] = archive_members
        except Exception:
            pass

    binary_raw_ext = lower.endswith((".mseed", ".seed", ".sac", ".wav", ".bin", ".raw"))
    raw_archive = any(re.search(r"(?:H10S|HA10).*(?:calib|RBC|sequence|\.mseed|\.seed|\.raw)", n, re.I) for n in archive_members)
    if (binary_raw_ext or raw_archive) and (H10_MARKERS.search(final) or H10_MARKERS.search(text)):
        return "RAW_CALIBRATION_SEQUENCE_BYTES", evidence + ["RAW_SAMPLE_PAYLOAD_SHAPE"]

    machine_ext = lower.endswith((".csv", ".json", ".txt", ".dat")) or "json" in ctype or "csv" in ctype
    if machine_ext and H10_MARKERS.search(text) and CAL_MARKERS.search(text) and numeric_density(text) >= 0.12:
        if TABLE_MARKERS.search(text):
            return "CALIBRATION_FLAG_OR_TIMESTAMP_TABLE", evidence + ["MACHINE_READABLE_NUMERIC_TABLE"]
        return "MACHINE_READABLE_NUMERIC_RESPONSE", evidence + ["MACHINE_READABLE_NUMERIC_TABLE"]

    is_pdf_or_slides = "pdf" in ctype or lower.endswith((".pdf", ".ppt", ".pptx"))
    if is_pdf_or_slides and H10_MARKERS.search(text) and (CAL_MARKERS.search(text) or DERIVED_MARKERS.search(text)):
        return "PUBLISHED_DERIVED_CALIBRATION_FIGURE_OR_PRESENTATION", evidence + ["DERIVED_PUBLICATION_CONTAINER"]

    if ACCESS_MARKERS.search(text) and re.search(r"vDEC|request", final, re.I):
        return "ACCESS_POINTER_ONLY", evidence
    if CAL_MARKERS.search(text) or H10_MARKERS.search(text):
        return "ABSTRACT_OR_DESCRIPTIVE_TEXT_ONLY", evidence
    return "FETCHED_OTHER", evidence


def historical_fact_flags(texts: list[str]) -> dict[str, bool]:
    joined = "\n".join(texts)
    patterns = {
        "sequence_rbc_cw": r"random broadband code|\bRBC\b.*CW|pseudo-random.*CW",
        "stored_at_idc": r"stored at the CTBTO.*International Data Centre|stored at.*IDC",
        "flag_start_end_table": r"Table from the CTBTO International Data Centre.*start/end times|calibration flag bit",
        "earliest_2013": r"Earliest calibration responses.*2013|calibration responses.*2013",
        "h10s1_08396_db": r"0\.8396\s*dB",
        "h10s23_large_differences": r"4\.7832\s*dB|2\.5340\s*dB",
        "fault_20130719": r"2013[/\-]07[/\-]19",
        "crosstalk": r"cross[- ]talk",
    }
    return {k: bool(re.search(p, joined, re.I | re.S)) for k, p in patterns.items()}


def run(output: Path) -> dict[str, Any]:
    try:
        output.resolve().relative_to(DATA.resolve())
        raise RuntimeError("CANONICAL_DATA_WRITE_FORBIDDEN")
    except ValueError:
        pass
    protocol, frozen = verify_contracts()
    domains = list(protocol["allowed_domains"])
    cc = protocol["crawl_contract"]
    max_bytes = int(cc["max_bytes_per_resource"])
    max_total = int(cc["max_total_fetches"])
    per_seed = int(cc["max_discovered_links_per_seed"])
    session = requests.Session()
    session.headers["User-Agent"] = UA

    queue = list(protocol["frozen_candidate_urls"])
    fetched_urls: set[str] = set()
    resources: list[dict[str, Any]] = []
    public_texts: list[str] = []
    while queue and len(resources) < max_total:
        url = queue.pop(0)
        if url in fetched_urls:
            continue
        fetched_urls.add(url)
        row = fetch_bounded(session, url, domains, max_bytes)
        cls, evidence = classify(row)
        row["evidence_class"] = cls
        row["classification_evidence"] = evidence
        text = payload_text(row) if row.get("status") == "FETCHED" else ""
        if text:
            row["extracted_text_sha256"] = hashlib.sha256(text.encode("utf-8")).hexdigest()
            row["extracted_text_chars"] = len(text)
            row["excerpt"] = text[:3000]
            public_texts.append(text)
        row.pop("data", None)
        resources.append(row)
        if row.get("status") == "FETCHED":
            raw_again = fetch_bounded(session, url, domains, max_bytes)
            links = discover_links(raw_again, domains, per_seed)
            for link in links:
                if link not in fetched_urls and link not in queue:
                    queue.append(link)

    counts: dict[str, int] = {}
    for r in resources:
        counts[r["evidence_class"]] = counts.get(r["evidence_class"], 0) + 1
    raw_count = counts.get("RAW_CALIBRATION_SEQUENCE_BYTES", 0)
    numeric_count = counts.get("MACHINE_READABLE_NUMERIC_RESPONSE", 0) + counts.get("CALIBRATION_FLAG_OR_TIMESTAMP_TABLE", 0)
    derived_count = counts.get("PUBLISHED_DERIVED_CALIBRATION_FIGURE_OR_PRESENTATION", 0)
    descriptive_count = counts.get("ABSTRACT_OR_DESCRIPTIVE_TEXT_ONLY", 0)
    access_count = counts.get("ACCESS_POINTER_ONLY", 0)

    if raw_count:
        verdict = "PUBLIC_H10S_CALIBRATION_SEQUENCE_BYTES_RECOVERED"
    elif numeric_count:
        verdict = "PUBLIC_NUMERIC_H10S_RESPONSE_RECOVERED__RAW_SEQUENCE_NOT_PUBLICLY_RECOVERED"
    elif derived_count or descriptive_count:
        verdict = "PUBLISHED_H10S_CALIBRATION_RESULTS_RECOVERED__RAW_SEQUENCE_PUBLIC_ACCESS_NOT_FOUND"
    else:
        verdict = "BLOCKED_CALIBRATION_PUBLIC_RECOVERY"

    receipt = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-HA10-CALIBRATION-SEQUENCE-PUBLIC-RECOVERY-RUN",
        "gate_id": protocol["gate_id"],
        "protocol_git_blob_sha1": EXPECTED_PROTOCOL_BLOB,
        "frozen_119hz_result_git_blob_sha1": EXPECTED_FROZEN_119_BLOB,
        "frozen_119hz_verdict": frozen["summary"]["verdict"],
        "selected_by_janus": protocol["selected_by_janus"],
        "verdict": verdict,
        "resource_count": len(resources),
        "evidence_class_counts": counts,
        "raw_sequence_bytes_recovered": raw_count > 0,
        "machine_readable_numeric_or_flag_table_recovered": numeric_count > 0,
        "derived_public_calibration_material_recovered": derived_count > 0,
        "controlled_access_pointer_recovered": access_count > 0,
        "historical_fact_flags": historical_fact_flags(public_texts),
        "resources": resources,
        "claim_ceiling": protocol["claim_ceiling"],
        "authority_delta_for_119hz": 0,
        "target_identity": "UNCONFIRMED",
        "target_evidence_delta": "NONE_FROM_CALIBRATION_PUBLIC_RECOVERY",
        "source_writeback": False,
        "janus_next_step_required": True,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    r = run(Path(args.output))
    print(json.dumps({
        "gate_id": r["gate_id"],
        "verdict": r["verdict"],
        "resource_count": r["resource_count"],
        "evidence_class_counts": r["evidence_class_counts"],
        "historical_fact_flags": r["historical_fact_flags"],
    }, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
