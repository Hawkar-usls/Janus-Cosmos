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
USER_AGENT = "Janus-Cosmos/0.4 (+https://github.com/Hawkar-usls/Janus-Cosmos)"
FITS_MAGIC = b"SIMPLE  ="
EVENT_LOG = Path("janus-cosmos-events.jsonl")


def emit(event: str, **fields: object) -> None:
    record = {"schema": "janus.cosmos.event.v0.1", "event": event, "ts_unix": time.time(), **fields}
    with EVENT_LOG.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record, sort_keys=True) + "\n")


def download(url: str, path: Path) -> dict:
    last_error = None
    for attempt in range(1, 4):
        try:
            request = urllib.request.Request(url, headers={
                "User-Agent": USER_AGENT,
                "Accept": "application/fits,application/octet-stream,*/*",
                "Connection": "close",
            })
            h = hashlib.sha256()
            total = 0
            content_type = ""
            with urllib.request.urlopen(request, timeout=180) as r, open(path, "wb") as f:
                status = getattr(r, "status", 200)
                content_type = r.headers.get("Content-Type", "")
                if status != 200:
                    raise RuntimeError(f"HTTP status {status} from {url}")
                while True:
                    chunk = r.read(1024 * 1024)
                    if not chunk:
                        break
                    total += len(chunk)
                    h.update(chunk)
                    f.write(chunk)
            if total < 2880:
                raise RuntimeError(f"Downloaded payload is too small to be FITS: {total} bytes")
            with open(path, "rb") as f:
                magic = f.read(len(FITS_MAGIC))
            if magic != FITS_MAGIC:
                raise RuntimeError(f"Payload is not a FITS primary HDU (magic={magic!r}, content_type={content_type!r})")
            metadata = {"bytes": total, "content_type": content_type, "sha256": h.hexdigest(), "attempt": attempt}
            emit("download_ok", url=url, **metadata)
            return metadata
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError, RuntimeError) as exc:
            last_error = exc
            emit("download_retry", url=url, attempt=attempt, error=repr(exc))
            print(f"Download attempt {attempt}/3 failed: {exc}")
            if attempt < 3:
                time.sleep(3 * attempt)
    emit("download_failed", url=url, error=repr(last_error))
    raise RuntimeError(f"MAST download failed after 3 attempts: {url}: {last_error}") from last_error


def read_image(path: Path) -> np.ndarray:
    with fits.open(path, memmap=True, lazy_load_hdus=True) as hdul:
        arrays = []
        for hdu in hdul:
            data = getattr(hdu, "data", None)
            if data is not None and np.ndim(data) == 2:
                arrays.append(np.asarray(data, dtype=np.float32))
        if not arrays:
            raise RuntimeError(f"No 2-D image plane in {path}")
        image = max(arrays, key=lambda a: a.size)
    image = np.nan_to_num(image, nan=0.0, posinf=0.0, neginf=0.0)
    h, w = image.shape
    side = min(h, w)
    y0, x0 = (h - side) // 2, (w - side) // 2
    image = image[y0:y0 + side, x0:x0 + side]
    image = zoom(image, SIZE / side, order=1)
    lo, hi = np.percentile(image, [1.0, 99.5])
    image = np.clip((image - lo) / max(hi - lo, 1e-9), 0, 1)
    image = np.arcsinh(6 * image) / np.arcsinh(6)
    return image.astype(np.float32)


def directional_correlation(image: np.ndarray, sigma: float, degrees: int) -> float:
    sm = gaussian_filter(image, sigma=sigma, mode="reflect")
    a = math.radians(degrees)
    dx = int(round(2.0 * sigma * math.cos(a)))
    dy = int(round(2.0 * sigma * math.sin(a)))
    if dx >= 0:
        ax0, ax1, bx0, bx1 = dx, SIZE, 0, SIZE - dx
    else:
        ax0, ax1, bx0, bx1 = 0, SIZE + dx, -dx, SIZE
    if dy >= 0:
        ay0, ay1, by0, by1 = dy, SIZE, 0, SIZE - dy
    else:
        ay0, ay1, by0, by1 = 0, SIZE + dy, -dy, SIZE
    if ax1 - ax0 < 8 or ay1 - ay0 < 8:
        return 0.0
    a1 = sm[ay0:ay1, ax0:ax1]
    b1 = sm[by0:by1, bx0:bx1]
    a1 = a1 - float(a1.mean())
    b1 = b1 - float(b1.mean())
    denom = float(np.sqrt(np.sum(a1 * a1) * np.sum(b1 * b1)))
    return abs(float(np.sum(a1 * b1)) / denom) if denom else 0.0


def score(image: np.ndarray) -> float:
    vals = [directional_correlation(image, s, o) for s in SCALES for o in ORIENTATIONS]
    return float(np.mean(vals))


def analyze(image: np.ndarray, rng: random.Random, filter_name: str) -> dict:
    observed = score(image)
    null = []
    for i in range(NULLS):
        shuffled = image.ravel().copy()
        rng.shuffle(shuffled)
        null_score = score(shuffled.reshape(image.shape))
        null.append(null_score)
        if i in (0, NULLS // 2, NULLS - 1):
            emit("null_checkpoint", filter=filter_name, null_index=i + 1, null_score=null_score)
    ge = sum(v >= observed for v in null)
    p = (ge + 1) / (NULLS + 1)
    result = {
        "observed_score": observed,
        "null_median": float(np.median(null)),
        "null_min": float(np.min(null)),
        "null_max": float(np.max(null)),
        "p_empirical": p,
        "candidate_by_filter": p < 0.05,
    }
    emit("filter_analyzed", filter=filter_name, **result)
    return result


def main():
    EVENT_LOG.unlink(missing_ok=True)
    emit("run_started", seed=SEED, nulls=NULLS, scales=SCALES, orientations_deg=ORIENTATIONS, image_size=SIZE)
    manifest = json.loads(Path("data/hst_real_pilot.json").read_text())
    out = {
        "schema": "janus.cosmos.hst.real_receipt.v0.3",
        "status": "REAL_HST_IMAGE_PILOT",
        "source": manifest["source_archive"],
        "target": manifest["target"],
        "target_class": manifest["target_class"],
        "blind_gate": manifest["blind_gate"],
        "source_products": [],
        "filters": {},
    }
    rng = random.Random(SEED)
    with tempfile.TemporaryDirectory() as td:
        for item in manifest["filters"]:
            path = Path(td) / (item["filter"] + ".fits")
            print(f"Downloading {item['filter']} from MAST: {item['url']}")
            emit("download_started", filter=item["filter"], band=item["band"], url=item["url"])
            metadata = download(item["url"], path)
            print(f"Downloaded {item['filter']}: {metadata['bytes']} bytes ({metadata['content_type']})")
            image = read_image(path)
            emit("fits_parsed", filter=item["filter"], shape=list(image.shape))
            result = analyze(image, rng, item["filter"])
            out["source_products"].append({"filter": item["filter"], "band": item["band"], "url": item["url"], **metadata})
            out["filters"][item["filter"]] = result
            print(f"Analyzed {item['filter']}: p_empirical={result['p_empirical']:.6f}")

    passing = [f for f, r in out["filters"].items() if r["candidate_by_filter"]]
    out["cross_band_candidate"] = len(passing) >= 2
    out["passing_filters"] = passing
    out["astronomical_discovery"] = False
    out["claim_ceiling"] = "Known M51 positive-control pilot only. A cross-band candidate is not a new discovery; it demonstrates that the image-level gate responds to real HST morphology."
    out["event_log"] = {"path": str(EVENT_LOG), "sha256": hashlib.sha256(EVENT_LOG.read_bytes()).hexdigest()}
    Path("janus-cosmos-real-hst-receipt.json").write_text(json.dumps(out, indent=2, sort_keys=True))
    emit("run_completed", cross_band_candidate=out["cross_band_candidate"], passing_filters=passing, receipt_sha256=hashlib.sha256(Path("janus-cosmos-real-hst-receipt.json").read_bytes()).hexdigest())
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
