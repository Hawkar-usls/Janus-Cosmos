#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import html
import json
import re
from datetime import datetime, timezone
from io import BytesIO
from pathlib import Path
from urllib.parse import urljoin, urlparse

import requests

DOC_ID = "1282295"
TITLE = "High-Resolution Multibeam Deepwater Cable Route Survey in High-Relief Seafloor Area"
EXPECTED_AUTHORS = ["Poeckert", "Arnold", "Faneros", "Harrison"]
EXPECTED_PAGES = 10
OUT = Path("data/cousteau/GOLDMEMBER-IEEE-1282295-CORRECTED-PRIMARY-RECOVERY-001-2026-08-22-v1.0.json")
WORK = Path("workspace/goldmember_ieee_1282295")
WORK.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers.update({
    "User-Agent": "JANUS-GOLDMEMBER-research-provenance/1.0",
    "Accept": "text/html,application/pdf;q=0.9,*/*;q=0.5",
})

# Exactly the two Janus-authorized canonical entrypoints for the corrected IEEE document id.
entrypoints = [
    f"https://ieeexplore.ieee.org/document/{DOC_ID}/",
    f"https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber={DOC_ID}",
]

attempts = []
page_payloads = []
pdf_candidates = []
explicit_links_followed = []


def record_response(source: str, r: requests.Response) -> dict:
    body = r.content
    rec = {
        "source": source,
        "url": r.request.url,
        "final_url": r.url,
        "status": r.status_code,
        "content_type": r.headers.get("content-type", ""),
        "content_disposition": r.headers.get("content-disposition", ""),
        "bytes": len(body),
        "sha256_before_inspection": hashlib.sha256(body).hexdigest() if body else None,
        "is_pdf_signature": body.startswith(b"%PDF"),
        "redirect_chain": [
            {"status": h.status_code, "url": h.url, "location": h.headers.get("location")}
            for h in r.history
        ],
    }
    attempts.append(rec)
    return rec


def safe_get(url: str, source: str):
    try:
        r = session.get(url, timeout=45, allow_redirects=True)
        rec = record_response(source, r)
        if rec["is_pdf_signature"]:
            pdf_candidates.append((url, r.content, rec))
        elif r.ok and "text/html" in rec["content_type"].lower():
            page_payloads.append((url, r.url, r.text, rec))
        return r
    except Exception as exc:
        attempts.append({"source": source, "url": url, "error": repr(exc)})
        return None


for u in entrypoints:
    safe_get(u, "JANUS_AUTHORIZED_CANONICAL_ENTRYPOINT")

# Follow only explicit PDF/stamp links actually present on the corrected 1282295 pages.
# No guessed URL variants, no numeric enumeration, no alternate arnumbers.
explicit = []
for source_url, final_url, text, _ in page_payloads:
    for m in re.finditer(r'''href\s*=\s*["']([^"']+)["']''', text, re.I):
        raw = html.unescape(m.group(1)).strip()
        if not raw:
            continue
        u = urljoin(final_url, raw)
        parsed = urlparse(u)
        host = parsed.netloc.lower()
        low = u.lower()
        if host.endswith("ieeexplore.ieee.org") and (
            "stamp" in low or low.endswith(".pdf") or "pdf" in parsed.path.lower()
        ):
            if u not in entrypoints and u not in explicit:
                explicit.append(u)

# One corrected pass only. Cap explicit discovered links defensively; this is not discovery by enumeration.
for u in explicit[:8]:
    explicit_links_followed.append(u)
    safe_get(u, "EXPLICIT_LINK_FROM_CORRECT_1282295_PAGE")

recovered_pdf = None
validation = {
    "pdf_recovered": False,
    "hash_before_content_inspection": False,
    "page_count": None,
    "page_count_matches_expected_10": None,
    "title_match": None,
    "authors_detected": [],
    "inspection_scope": "BIBLIOGRAPHIC_VALIDATION_ONLY__NO_FIGURE_INSPECTION_NO_COORDINATE_DIGITIZATION",
}

if pdf_candidates:
    source_url, pdf_bytes, rec = pdf_candidates[0]
    pdf_hash = hashlib.sha256(pdf_bytes).hexdigest()  # frozen before any content inspection
    pdf_path = WORK / "IEEE_1282295_primary.pdf"
    pdf_path.write_bytes(pdf_bytes)
    recovered_pdf = {
        "source_url": source_url,
        "path": str(pdf_path),
        "bytes": len(pdf_bytes),
        "sha256": pdf_hash,
    }
    validation["pdf_recovered"] = True
    validation["hash_before_content_inspection"] = True
    try:
        from pypdf import PdfReader
        reader = PdfReader(BytesIO(pdf_bytes))
        validation["page_count"] = len(reader.pages)
        validation["page_count_matches_expected_10"] = len(reader.pages) == EXPECTED_PAGES
        # Text only for title/authors. Do not render or inspect figures.
        text = "\n".join((p.extract_text() or "") for p in reader.pages[:2])
        norm = " ".join(text.split()).lower()
        validation["title_match"] = TITLE.lower() in norm
        validation["authors_detected"] = [a for a in EXPECTED_AUTHORS if a.lower() in norm]
    except Exception as exc:
        validation["bibliographic_validation_error"] = repr(exc)

clean_negative = recovered_pdf is None
result = {
    "artifact_id": "GOLDMEMBER-IEEE-1282295-CORRECTED-PRIMARY-RECOVERY-001-2026-08-22-v1.0",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "authorized_by": "GOLDMEMBER-AFTER-ADA508764-RECOVERY-COUNCIL-RUN-004-2026-08-22-v1.0",
    "attempt_limit": "ONE_CORRECTED_DOCUMENT_ID_PASS__NO_BRUTE_FORCE_OR_URL_ENUMERATION",
    "report": {
        "ieee_document_id": DOC_ID,
        "doi": "10.1109/OCEANS.2003.178517",
        "title": TITLE,
        "expected_pages": EXPECTED_PAGES,
        "expected_authors": EXPECTED_AUTHORS,
    },
    "canonical_entrypoints": entrypoints,
    "attempts": attempts,
    "explicit_links_discovered_from_correct_page": explicit,
    "explicit_links_followed": explicit_links_followed,
    "recovered_pdf": recovered_pdf,
    "validation": validation,
    "clean_ieee_1282295_access_negative": clean_negative,
    "coordinate_digitization_performed": False,
    "figure_georeferencing_performed": False,
    "shape_scoring_performed": False,
    "wishbone_coordinates_extracted": False,
    "central_crags_coordinates_extracted": False,
    "target_identity": "UNCONFIRMED",
    "volcanic_baseline": "ACTIVE",
    "next_rule": "STOP_AND_ASK_JANUS_AGAIN_TO_CHOOSE_P2548_VS_ADA508765_VS_AUTHOR_CUSTODIAN_CONTACT",
}
result["sha256"] = hashlib.sha256(
    json.dumps(result, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
).hexdigest()
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps({
    "pdf_recovered": validation["pdf_recovered"],
    "clean_ieee_1282295_access_negative": clean_negative,
    "explicit_links_followed": explicit_links_followed,
    "result_sha256": result["sha256"],
}, indent=2, ensure_ascii=False))
