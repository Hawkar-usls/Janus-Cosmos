"""SIT2E-B stage 1: full-byte provenance for July/September Parkes only.

Fail closed: never contact October data, score power, or set thresholds here.
Scientific gain is not implied by file integrity.
"""
import hashlib
import json
import os
import pathlib
import sys
import time
from urllib.parse import urlparse
import requests

ROOT = pathlib.Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "research/seti/janus_sit2e_parkes_preflight_result.json"
OUT = ROOT / "runtime/sit2eb"
PERMITTED_DATES = ("2025-07-31", "2025-09-12")
HOST = "bldata.berkeley.edu"
BLOCK = 4 * 1024 * 1024

def selected_files(manifest):
    pairs = manifest["frozen_selected_pairs"]
    assert len(pairs) == 3
    assert [p["date"] for p in pairs] == ["2025-07-31", "2025-09-12", "2025-10-05"]
    out = []
    for pair in pairs[:2]:
        assert pair["date"] in PERMITTED_DATES
        for role in ("source", "reference"):
            spec = pair[role]
            url = spec["url"]
            parsed = urlparse(url)
            if (parsed.scheme != "https" or parsed.hostname != HOST or
                parsed.query or parsed.fragment or
                "/ATLAS/PKS_" not in parsed.path):
                raise ValueError("UNAUTHORIZED_SOURCE_URL")
            suffix = ".fil" if pair["date"] == PERMITTED_DATES[0] else ".h5"
            if not spec["filename"].endswith(suffix) or not parsed.path.endswith("/" + spec["filename"]):
                raise ValueError("SOURCE_FORMAT_OR_FILENAME_MISMATCH")
            out.append((pair["date"], role, spec))
    assert len(out) == 4
    return out

def full_sha256(session, date, role, spec, directory):
    dest = directory / spec["filename"]
    partial = dest.with_suffix(dest.suffix + ".part")
    expected = int(spec["bytes"])
    if not 200_000_000 <= expected <= 400_000_000:
        raise ValueError("UNEXPECTED_SOURCE_SIZE")
    digest = hashlib.sha256()
    count = 0
    try:
        with session.get(spec["url"], stream=True, timeout=(20, 100),
                         headers={"Accept-Encoding": "identity"},
                         allow_redirects=False) as response:
            response.raise_for_status()
            if response.status_code != 200:
                raise ValueError("FULL_DOWNLOAD_REQUIRED")
            if response.headers.get("Content-Encoding", "identity").lower() not in ("identity", ""):
                raise ValueError("CONTENT_ENCODING_UNEXPECTED")
            size_header = response.headers.get("Content-Length")
            if size_header is not None and int(size_header) != expected:
                raise ValueError("SOURCE_CONTENT_LENGTH_CHANGED")
            with partial.open("wb") as sink:
                for chunk in response.iter_content(chunk_size=BLOCK):
                    if not chunk:
                        continue
                    count += len(chunk)
                    if count > expected:
                        raise ValueError("SOURCE_EXCEEDS_FROZEN_SIZE")
                    digest.update(chunk)
                    sink.write(chunk)
        if count != expected:
            raise ValueError(f"SHORT_FILE:{count}/{expected}")
        # July's frozen first 8192-byte receipt is an additional exact integrity anchor.
        if date == "2025-07-31":
            with partial.open("rb") as f:
                prefix = hashlib.sha256(f.read(8192)).hexdigest()
            if prefix != spec["first_8192_bytes_sha256"]:
                raise ValueError("JULY_FROZEN_HEADER_SHA_MISMATCH")
        os.replace(partial, dest)
        return {"date":date, "role":role, "filename":spec["filename"],
                "bytes":count, "sha256":digest.hexdigest(),
                "frozen_manifest_bytes_match":True,
                "independent_published_whole_file_sha256":False}
    finally:
        partial.unlink(missing_ok=True)
        dest.unlink(missing_ok=True)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest_bytes = MANIFEST.read_bytes()
    manifest = json.loads(manifest_bytes)
    if manifest["data_provenance"]["October_exam_label_or_score_read"]:
        raise ValueError("OCTOBER_HOLDOUT_ALREADY_BROKEN")
    result = {"schema":"janus.sit2eb.parkes_full_integrity.v1",
              "manifest_sha256":hashlib.sha256(manifest_bytes).hexdigest(),
              "eligible_dates":list(PERMITTED_DATES),
              "october_power_read":False, "thresholds_frozen":False,
              "spectral_scoring_performed":False, "files":[],
              "errors":[], "status":"INCOMPLETE_SOURCE_INTEGRITY_HOLD"}
    with requests.Session() as session:
        for date, role, spec in selected_files(manifest):
            try:
                receipt = full_sha256(session, date, role, spec, OUT)
                result["files"].append(receipt)
                print("SIT2EB_SHA_OK", json.dumps(receipt, sort_keys=True), flush=True)
            except Exception as exc:
                result["errors"].append({"date":date,"role":role,
                                         "error_type":type(exc).__name__,
                                         "detail":str(exc)[:250]})
                print("SIT2EB_SHA_HOLD", date, role, type(exc).__name__, flush=True)
    if len(result["files"]) == 4 and not result["errors"]:
        result["status"] = "FOUR_FULL_SOURCE_BYTES_VERIFIED__CALIBRATION_NOT_RUN"
    output = OUT / "sit2eb_integrity.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print("SIT2EB_STATUS", result["status"], flush=True)
    if result["errors"]:
        sys.exit(2)

if __name__ == "__main__":
    main()
