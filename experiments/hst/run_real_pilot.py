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
USER_AGENT = "Janus-Cosmos/0.2 (+https://github.com/Hawkar-usls/Janus-Cosmos)"


def download(url: str, path: Path) -> str:
    """Download a public MAST product with retries and an explicit User-Agent."""
    last_error = None
    for attempt in range(1, 4):
        try:
            request = urllib.request.Request(
                url,
                headers={
                    "User-Agent": USER_AGENT,
                    "Accept": "application/fits,application/octet-stream,*/*",
                },
            )
            h = hashlib.sha256()
            with urllib.request.urlopen(request, timeout=180) as r, open(path, "wb") as f:
                while True:
                    chunk = r.read(1024 * 1024)
                    if not chunk:
                        break
                    h.update(chunk)
                    f.write(chunk)
            if path.stat().st_size < 1024:
                raise RuntimeError(f"Downloaded file is unexpectedly small: {path.stat().st_size} bytes")
            return h.hexdigest()
        except (urllib.error.URLError, TimeoutError, OSError, RuntimeError) as exc:
            last_error = exc
            if attempt < 3:
                time.sleep(3 * attempt)
    raise RuntimeError(f"MAST download failed after 3 attempts: {url}: {last_error}") from last_error


def read_image(path: Path) -> np.ndarray:
    with fits.open(path, memmap=True) as hdul:
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
    dx = max(1, int(round(2.0 * sigma * math.cos(a))))
    dy = max(1, int(round(2.0 * sigma * math.sin(a))))
    y1, y2 = max(0, dy), min(SIZE, SIZE + dy)
    x1, x2 = max(0, dx), min(SIZE, SIZE + dx)
    if y2 - y1 < 8 or x2 - x1 < 8:
        return 0.0
    a1 = sm[y1:y2, x1:x2]
    b1 = sm[y1-dy:y2-dy, x1-dx:x2-dx]
    a1 = a1 - float(a1.mean())
    b1 = b1 - float(b1.mean())
    denom = float(np.sqrt(np.sum(a1 * a1) * np.sum(b1 * b1)))
    return abs(float(np.sum(a1 * b1)) / denom) if denom else 0.0


def score(image: np.ndarray) -> float:
    vals = [directional_correlation(image, s, o) for s in SCALES for o in ORIENTATIONS]
    return float(np.mean(vals))


def analyze(image: np.ndarray, rng: random.Random) -> dict:
    observed = score(image)
    flat = image.ravel().copy()
    null = []
    for _ in range(NULLS):
        rng.shuffle(flat)
        null.append(score(flat.reshape(image.shape)))
    ge = sum(v >= observed for v in null)
    p = (ge + 1) / (NULLS + 1)
    return {
        "observed_score": observed,
        "null_median": float(np.median(null)),
        "null_min": float(np.min(null)),
        "null_max": float(np.max(null)),
        "p_empirical": p,
        "candidate_by_filter": p < 0.05,
    }


def main():
    manifest = json.loads(Path("data/hst_real_pilot.json").read_text())
    out = {
        "schema": "janus.cosmos.hst.real_receipt.v0.2",
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
            sha = download(item["url"], path)
            print(f"Downloaded {item['filter']}: {path.stat().st_size} bytes")
            image = read_image(path)
            result = analyze(image, rng)
            out["source_products"].append({
                "filter": item["filter"],
                "band": item["band"],
                "url": item["url"],
                "sha256": sha,
            })
            out["filters"][item["filter"]] = result
            print(f"Analyzed {item['filter']}: p_empirical={result['p_empirical']:.6f}")

    passing = [f for f, r in out["filters"].items() if r["candidate_by_filter"]]
    out["cross_band_candidate"] = len(passing) >= 2
    out["passing_filters"] = passing
    out["astronomical_discovery"] = False
    out["claim_ceiling"] = (
        "Known M51 positive-control pilot only. A cross-band candidate is not a new discovery; "
        "it demonstrates that the image-level gate responds to real HST morphology."
    )
    Path("janus-cosmos-real-hst-receipt.json").write_text(json.dumps(out, indent=2, sort_keys=True))
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
