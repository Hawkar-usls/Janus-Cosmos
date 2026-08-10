from __future__ import annotations

import hashlib
import json
import math
import random
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter, zoom

SEED = 20260810
NULLS = 256
SCALES = (1.0, 2.0, 4.0)
ORIENTATIONS = tuple(range(0, 180, 30))
SIZE = 256
USER_AGENT = "Janus-Cosmos/0.5-blind-corpus (+https://github.com/Hawkar-usls/Janus-Cosmos)"
FITS_MAGIC = b"SIMPLE  ="
EVENT_LOG = Path("janus-cosmos-blind-events.jsonl")
RECEIPT = Path("janus-cosmos-blind-corpus-receipt.json")


def emit(event: str, **fields: object) -> None:
    record = {"schema": "janus.cosmos.blind_event.v0.1", "event": event, "ts_unix": time.time(), **fields}
    with EVENT_LOG.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record, sort_keys=True) + "\n")


def download(url: str, path: Path) -> dict:
    last_error = None
    for attempt in range(1, 4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/fits,application/octet-stream,*/*", "Connection": "close"})
            h = hashlib.sha256(); total = 0; content_type = ""
            with urllib.request.urlopen(req, timeout=180) as r, open(path, "wb") as f:
                status = getattr(r, "status", 200); content_type = r.headers.get("Content-Type", "")
                if status != 200: raise RuntimeError(f"HTTP status {status}")
                while True:
                    chunk = r.read(1024 * 1024)
                    if not chunk: break
                    total += len(chunk); h.update(chunk); f.write(chunk)
            if total < 2880: raise RuntimeError(f"payload too small: {total}")
            with open(path, "rb") as f: magic = f.read(len(FITS_MAGIC))
            if magic != FITS_MAGIC: raise RuntimeError(f"not FITS: magic={magic!r}, content_type={content_type!r}")
            meta = {"bytes": total, "content_type": content_type, "sha256": h.hexdigest(), "attempt": attempt}
            emit("download_ok", url=url, **meta); return meta
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError, RuntimeError) as exc:
            last_error = repr(exc); emit("download_retry", url=url, attempt=attempt, error=last_error)
            if attempt < 3: time.sleep(3 * attempt)
    emit("download_failed", url=url, error=last_error)
    raise RuntimeError(f"download failed after 3 attempts: {url}: {last_error}")


def read_image(path: Path) -> np.ndarray:
    with fits.open(path, memmap=True, lazy_load_hdus=True) as hdul:
        planes = [np.asarray(h.data, dtype=np.float32) for h in hdul if getattr(h, "data", None) is not None and np.ndim(h.data) == 2]
    if not planes: raise RuntimeError(f"No 2-D image plane in {path}")
    image = max(planes, key=lambda a: a.size)
    image = np.nan_to_num(image, nan=0.0, posinf=0.0, neginf=0.0)
    h, w = image.shape; side = min(h, w); y0, x0 = (h-side)//2, (w-side)//2
    image = zoom(image[y0:y0+side, x0:x0+side], SIZE/side, order=1)
    lo, hi = np.percentile(image, [1.0, 99.5]); image = np.clip((image-lo)/max(hi-lo, 1e-9), 0, 1)
    return (np.arcsinh(6*image)/np.arcsinh(6)).astype(np.float32)


def directional_correlation(image: np.ndarray, sigma: float, degrees: int) -> float:
    sm = gaussian_filter(image, sigma=sigma, mode="reflect"); a = math.radians(degrees)
    dx, dy = int(round(2*sigma*math.cos(a))), int(round(2*sigma*math.sin(a)))
    ax0, ax1, bx0, bx1 = (dx, SIZE, 0, SIZE-dx) if dx >= 0 else (0, SIZE+dx, -dx, SIZE)
    ay0, ay1, by0, by1 = (dy, SIZE, 0, SIZE-dy) if dy >= 0 else (0, SIZE+dy, -dy, SIZE)
    if ax1-ax0 < 8 or ay1-ay0 < 8: return 0.0
    x, y = sm[ay0:ay1, ax0:ax1], sm[by0:by1, bx0:bx1]
    x -= float(x.mean()); y -= float(y.mean()); denom = float(np.sqrt(np.sum(x*x)*np.sum(y*y)))
    return abs(float(np.sum(x*y))/denom) if denom else 0.0


def score(image: np.ndarray) -> float:
    return float(np.mean([directional_correlation(image, s, o) for s in SCALES for o in ORIENTATIONS]))


def analyze(image: np.ndarray, rng: random.Random, target: str, filter_name: str) -> dict:
    observed = score(image); null = []
    for i in range(NULLS):
        shuffled = image.ravel().copy(); rng.shuffle(shuffled); value = score(shuffled.reshape(image.shape)); null.append(value)
        if i in (0, NULLS//2, NULLS-1): emit("null_checkpoint", target=target, filter=filter_name, null_index=i+1, null_score=value)
    ge = sum(v >= observed for v in null); p = (ge+1)/(NULLS+1)
    result = {"observed_score": observed, "null_median": float(np.median(null)), "null_min": float(np.min(null)), "null_max": float(np.max(null)), "p_empirical": p, "candidate_by_filter": p < 0.05}
    emit("filter_analyzed", target=target, filter=filter_name, **result); return result


def main() -> None:
    EVENT_LOG.unlink(missing_ok=True)
    emit("run_started", seed=SEED, nulls=NULLS, scales=SCALES, orientations_deg=ORIENTATIONS, image_size=SIZE, semantic_analysis=False, ocr=False, face_search=False, cipher_search=False, post_hoc_tuning=False)
    manifest = json.loads(Path("data/hst_blind_corpus.json").read_text())
    rng = random.Random(SEED); targets_out = []
    with tempfile.TemporaryDirectory() as td:
        for target in manifest["targets"]:
            target_out = {"target": target["target"], "target_class": target["class"], "filters": {}, "source_products": []}
            for item in target["filters"]:
                path = Path(td) / f"{target['target']}_{item['filter']}.fits"
                emit("download_started", target=target["target"], filter=item["filter"], url=item["url"])
                metadata = download(item["url"], path); image = read_image(path)
                emit("fits_parsed", target=target["target"], filter=item["filter"], shape=list(image.shape))
                target_out["filters"][item["filter"]] = analyze(image, rng, target["target"], item["filter"])
                target_out["source_products"].append({"filter": item["filter"], "band": item["band"], "url": item["url"], **metadata})
            passing = [f for f, r in target_out["filters"].items() if r["candidate_by_filter"]]
            target_out["passing_filters"] = passing; target_out["cross_band_candidate"] = len(passing) >= 2
            targets_out.append(target_out); emit("target_completed", target=target["target"], passing_filters=passing, cross_band_candidate=target_out["cross_band_candidate"])
    out = {"schema":"janus.cosmos.hst.blind_corpus_receipt.v0.1", "status":"BLIND_GEOMETRIC_CORPUS_PILOT", "source":manifest["source_archive"], "selection":manifest["selection"], "blind_gate":{"seed":SEED,"nulls":NULLS,"scales":SCALES,"orientations_deg":ORIENTATIONS,"semantic_analysis":False,"ocr":False,"face_search":False,"cipher_search":False,"post_hoc_tuning":False,"cross_filter_required":True}, "targets":targets_out, "candidate_count":sum(1 for t in targets_out if t["cross_band_candidate"]), "claim_ceiling":"Image-level geometric candidates only. No astronomical discovery claim; blind replication and independent scientific review are required.", "event_log":{"path":str(EVENT_LOG),"sha256":hashlib.sha256(EVENT_LOG.read_bytes()).hexdigest()}}
    RECEIPT.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    emit("run_completed", candidate_count=out["candidate_count"], receipt_sha256=hashlib.sha256(RECEIPT.read_bytes()).hexdigest())
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__": main()
