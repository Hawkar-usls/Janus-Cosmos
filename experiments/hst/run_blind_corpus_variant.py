from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import random
import tempfile
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "run_blind_corpus.py"
spec = importlib.util.spec_from_file_location("janus_cosmos_base", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError("Unable to load Janus Cosmos base runner")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


def download_variant(url: str, path: Path, backend: str) -> dict:
    if backend not in {"direct", "mast_api"}:
        raise ValueError(f"Unsupported backend: {backend}")
    source_url = url
    if backend == "mast_api":
        source_url = "https://mast.stsci.edu/api/v0.1/Download/file?" + urllib.parse.urlencode({"uri": url})
    req = urllib.request.Request(
        source_url,
        headers={
            "User-Agent": f"Janus-Cosmos/0.5-variant-{backend} (+https://github.com/Hawkar-usls/Janus-Cosmos)",
            "Accept": "application/fits,application/octet-stream,*/*",
            "Connection": "close",
        },
    )
    last_error = None
    for attempt in range(1, 4):
        try:
            h = hashlib.sha256(); total = 0; content_type = ""
            with urllib.request.urlopen(req, timeout=180) as r, open(path, "wb") as f:
                status = getattr(r, "status", 200)
                content_type = r.headers.get("Content-Type", "")
                if status != 200:
                    raise RuntimeError(f"HTTP status {status}")
                while True:
                    chunk = r.read(1024 * 1024)
                    if not chunk:
                        break
                    total += len(chunk)
                    h.update(chunk)
                    f.write(chunk)
            if total < 2880:
                raise RuntimeError(f"payload too small: {total}")
            with open(path, "rb") as f:
                magic = f.read(len(base.FITS_MAGIC))
            if magic != base.FITS_MAGIC:
                raise RuntimeError(f"not FITS: magic={magic!r}, content_type={content_type!r}")
            meta = {
                "backend": backend,
                "source_url": source_url,
                "original_url": url,
                "bytes": total,
                "content_type": content_type,
                "sha256": h.hexdigest(),
                "attempt": attempt,
            }
            base.emit("download_ok", **meta)
            return meta
        except Exception as exc:
            last_error = repr(exc)
            base.emit("download_retry", backend=backend, original_url=url, source_url=source_url, attempt=attempt, error=last_error)
            if attempt < 3:
                time.sleep(3 * attempt)
    base.emit("download_failed", backend=backend, original_url=url, source_url=source_url, error=last_error)
    raise RuntimeError(f"{backend} download failed after 3 attempts: {source_url}: {last_error}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=["direct", "mast_api"], required=True)
    ap.add_argument("--output-dir", default=".")
    args = ap.parse_args()

    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    base.EVENT_LOG = out_dir / "janus-cosmos-blind-events.jsonl"
    base.RECEIPT = out_dir / "janus-cosmos-blind-corpus-receipt.json"
    base.EVENT_LOG.unlink(missing_ok=True)
    seed = base.SEED
    rng = random.Random(seed)
    manifest = json.loads(Path("data/hst_blind_corpus.json").read_text(encoding="utf-8"))
    base.emit("run_started", backend=args.backend, seed=seed, nulls=base.NULLS, scales=base.SCALES, orientations_deg=base.ORIENTATIONS, image_size=base.SIZE, semantic_analysis=False, ocr=False, face_search=False, cipher_search=False, post_hoc_tuning=False)
    targets = []
    with tempfile.TemporaryDirectory() as td:
        for target in manifest["targets"]:
            t = {"target": target["target"], "target_class": target["class"], "filters": {}, "source_products": []}
            for item in target["filters"]:
                fitspath = Path(td) / f"{target['target']}_{item['filter']}.fits"
                base.emit("download_started", backend=args.backend, target=target["target"], filter=item["filter"], url=item["url"])
                metadata = download_variant(item["url"], fitspath, args.backend)
                image = base.read_image(fitspath)
                base.emit("fits_parsed", backend=args.backend, target=target["target"], filter=item["filter"], shape=list(image.shape))
                result = base.analyze(image, rng, target["target"], item["filter"])
                t["filters"][item["filter"]] = result
                t["source_products"].append({"filter": item["filter"], "band": item["band"], **metadata})
            passing = [f for f, r in t["filters"].items() if r["candidate_by_filter"]]
            t["passing_filters"] = passing
            t["cross_band_candidate"] = len(passing) >= 2
            targets.append(t)
            base.emit("target_completed", backend=args.backend, target=target["target"], passing_filters=passing, cross_band_candidate=t["cross_band_candidate"])
    receipt = {
        "schema": "janus.cosmos.hst.blind_corpus_variant_receipt.v0.1",
        "status": "BLIND_GEOMETRIC_CORPUS_VARIANT",
        "backend": args.backend,
        "source": manifest["source_archive"],
        "selection": manifest["selection"],
        "targets": targets,
        "candidate_count": sum(1 for t in targets if t["cross_band_candidate"]),
        "blind_gate": {"seed": seed, "nulls": base.NULLS, "scales": list(base.SCALES), "orientations_deg": list(base.ORIENTATIONS), "semantic_analysis": False, "ocr": False, "face_search": False, "cipher_search": False, "post_hoc_tuning": False, "cross_filter_required": True},
        "claim_ceiling": "Image-level geometric candidates only. No astronomical discovery claim; independent blind replication and scientific review are required.",
    }
    base.RECEIPT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    receipt_sha = hashlib.sha256(base.RECEIPT.read_bytes()).hexdigest()
    base.emit("run_completed", backend=args.backend, candidate_count=receipt["candidate_count"], receipt_sha256=receipt_sha)
    receipt["event_log"] = {"path": str(base.EVENT_LOG), "sha256": hashlib.sha256(base.EVENT_LOG.read_bytes()).hexdigest()}
    receipt["receipt_sha256"] = receipt_sha
    base.RECEIPT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
