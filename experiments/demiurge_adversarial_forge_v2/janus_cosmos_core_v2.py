#!/usr/bin/env python3
from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from scipy import ndimage

VERSION = "2.0.0"
IMAGE_SIZE = 128
IAAFT_ITER = 3

FEATURE_NAMES = [
    "directional_s1", "directional_s2", "directional_s4",
    "rot90_corr", "rot180_corr",
    "gradient_anisotropy", "high_frequency_energy", "fourier_angular_anisotropy",
    "component_count_q80", "component_count_q90", "component_count_q95",
    "largest_component_q80", "largest_component_q90", "largest_component_q95",
]
FEATURE_GROUPS = {
    "directional": [0, 1, 2],
    "rotation": [3, 4],
    "anisotropy": [5, 7],
    "high_frequency": [6],
    "component_count": [8, 9, 10],
    "largest_component": [11, 12, 13],
}


def canonical_json(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def stable_seed(*parts) -> int:
    return int.from_bytes(hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:8], "big")


def _parse_fits_value(raw: str):
    raw = raw.strip()
    if raw.startswith("'"):
        end = raw.find("'", 1)
        return raw[1:end].strip() if end >= 1 else raw.strip("' ")
    token = raw.split("/", 1)[0].strip()
    if token == "T":
        return True
    if token == "F":
        return False
    if not token:
        return None
    try:
        return float(token.replace("D", "E")) if any(c in token for c in ".EeDd") else int(token)
    except Exception:
        return token


def _read_header(f):
    cards = []
    while True:
        block = f.read(2880)
        if len(block) != 2880:
            raise RuntimeError("truncated FITS header")
        for i in range(0, 2880, 80):
            card = block[i:i+80].decode("ascii", "replace")
            cards.append(card)
            if card.startswith("END"):
                h = {}
                for c in cards:
                    k = c[:8].strip()
                    if k and len(c) > 8 and c[8] == "=":
                        h[k] = _parse_fits_value(c[10:])
                return h


def read_primary_fits(path: Path):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rb") as f:
        h = _read_header(f)
        bitpix = int(h.get("BITPIX", 0))
        naxis = int(h.get("NAXIS", 0))
        if naxis != 2:
            raise RuntimeError(f"expected NAXIS=2, got {naxis}")
        n1, n2 = int(h["NAXIS1"]), int(h["NAXIS2"])
        dtypes = {-32: ">f4", -64: ">f8", 16: ">i2", 32: ">i4", 8: ">u1"}
        if bitpix not in dtypes:
            raise RuntimeError(f"unsupported BITPIX={bitpix}")
        count = n1 * n2
        nbytes = count * (abs(bitpix) // 8)
        raw = f.read(nbytes)
        if len(raw) != nbytes:
            raise RuntimeError("truncated FITS data")
        arr = np.frombuffer(raw, dtype=dtypes[bitpix], count=count).reshape((n2, n1)).astype(np.float32)
        bscale = float(h.get("BSCALE", 1) or 1)
        bzero = float(h.get("BZERO", 0) or 0)
        if bscale != 1 or bzero != 0:
            arr = arr * bscale + bzero
        meta = {
            "fits_bitpix": bitpix,
            "fits_native_shape": [n2, n1],
            "fits_object": h.get("OBJECT", ""),
            "fits_header_filter": h.get("FILTER", ""),
            "fits_instrument": h.get("INSTRUME", ""),
            "fits_exptime": h.get("EXPTIME", None),
        }
        return arr, h, meta


def normalize(a: np.ndarray, genome: dict, image_size: int = IMAGE_SIZE) -> np.ndarray:
    x = np.nan_to_num(np.asarray(a, dtype=np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    h, w = x.shape
    side = max(h, w)
    py, px = side - h, side - w
    fill = float(np.median(x[np.isfinite(x)])) if np.any(np.isfinite(x)) else 0.0
    x = np.pad(x, ((py//2, py-py//2), (px//2, px-px//2)), mode="constant", constant_values=fill)
    if side != image_size:
        x = ndimage.zoom(x, (image_size/side, image_size/side), order=1)
    lo = float(np.percentile(x, 1.0))
    hi = float(np.percentile(x, float(genome["clip_high"])))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.zeros((image_size, image_size), np.float32)
    x = np.clip((x - lo) / (hi - lo), 0, 1)
    scale = float(genome["asinh_scale"])
    return (np.arcsinh(scale * x) / np.arcsinh(scale)).astype(np.float32)


def corr(a, b) -> float:
    aa = np.asarray(a, np.float64).ravel()
    bb = np.asarray(b, np.float64).ravel()
    aa -= aa.mean(); bb -= bb.mean()
    den = np.linalg.norm(aa) * np.linalg.norm(bb)
    return float(np.dot(aa, bb) / den) if den else 0.0


def shift_corr(sm, sigma, deg) -> float:
    h, w = sm.shape
    a = math.radians(deg)
    dx = int(round(2 * sigma * math.cos(a)))
    dy = int(round(2 * sigma * math.sin(a)))
    ax0, ax1, bx0, bx1 = (dx, w, 0, w-dx) if dx >= 0 else (0, w+dx, -dx, w)
    ay0, ay1, by0, by1 = (dy, h, 0, h-dy) if dy >= 0 else (0, h+dy, -dy, h)
    if ax1-ax0 < 8 or ay1-ay0 < 8:
        return 0.0
    return abs(corr(sm[ay0:ay1, ax0:ax1], sm[by0:by1, bx0:bx1]))


def geometry(x: np.ndarray) -> np.ndarray:
    ds = []
    for sigma in (1.0, 2.0, 4.0):
        sm = ndimage.gaussian_filter(x, sigma, mode="reflect")
        ds.append(float(np.mean([shift_corr(sm, sigma, d) for d in range(0, 180, 30)])))
    rot90 = corr(x, np.rot90(x, 1))
    rot180 = corr(x, np.rot90(x, 2))
    gy, gx = np.gradient(ndimage.gaussian_filter(x, 1.0, mode="reflect"))
    jxx, jyy, jxy = float(np.mean(gx*gx)), float(np.mean(gy*gy)), float(np.mean(gx*gy))
    tr = jxx + jyy
    disc = math.sqrt(max((jxx-jyy)**2 + 4*jxy*jxy, 0))
    l1, l2 = 0.5*(tr+disc), 0.5*(tr-disc)
    grad = float((l1-l2)/max(l1+l2, 1e-12))
    smooth = ndimage.gaussian_filter(x, 4.0, mode="reflect")
    hf = float(np.mean((x-smooth)**2)/max(np.mean(x**2), 1e-12))
    f = np.fft.fftshift(np.fft.fft2(x-x.mean()))
    p = np.abs(f)**2
    h, w = p.shape
    yy, xx = np.indices(p.shape)
    cy, cx = (h-1)/2, (w-1)/2
    dx, dy = xx-cx, yy-cy
    r = np.sqrt(dx*dx + dy*dy)
    th = np.mod(np.arctan2(dy, dx), np.pi)
    mask = r > max(h,w)*.03
    bins = np.linspace(0, np.pi, 13)
    vals = []
    for blo, bhi in zip(bins[:-1], bins[1:]):
        m = mask & (th >= blo) & (th < bhi)
        vals.append(float(p[m].mean()) if np.any(m) else 0.0)
    vals = np.asarray(vals)
    fa = float(vals.std()/max(vals.mean(),1e-12))
    counts, largest = [], []
    n = x.size
    for q in (.80, .90, .95):
        thr = float(np.quantile(x, q))
        labels, count = ndimage.label(x >= thr)
        sizes = np.bincount(labels.ravel())[1:] if count else np.array([], dtype=int)
        counts.append(float(np.log1p(count)))
        largest.append(float(sizes.max()/n) if sizes.size else 0.0)
    return np.asarray([*ds, rot90, rot180, grad, hf, fa, *counts, *largest], np.float64)


def feature_weights(genome: dict) -> np.ndarray:
    w = np.ones(len(FEATURE_NAMES), dtype=np.float64)
    for group, idxs in FEATURE_GROUPS.items():
        val = float(genome["group_weights"][group])
        for i in idxs:
            w[i] = val
    # Normalize so global weight scale cannot inflate the statistic.
    w /= max(float(np.mean(w)), 1e-12)
    return w


def weighted_std_dist(v, center, scale, genome: dict) -> float:
    z = (np.asarray(v) - np.asarray(center)) / np.maximum(np.asarray(scale), 1e-9)
    w = feature_weights(genome)
    return float(np.sqrt(np.sum(w * z * z) / np.sum(w)))


def quantile_remap(src, ref):
    s = np.asarray(src, np.float64).ravel()
    r = np.sort(np.asarray(ref, np.float64).ravel())
    order = np.argsort(s, kind="mergesort")
    out = np.empty_like(s)
    out[order] = r
    return out.reshape(src.shape).astype(np.float32)


def phase_iaaft(x, rng, genome: dict):
    h, w = x.shape
    orig = np.fft.rfft2(x)
    amp = np.abs(orig)
    op = np.angle(orig)
    phase_noise = np.angle(np.fft.rfft2(rng.normal(size=x.shape)))
    fy = np.fft.fftfreq(h)[:,None]
    fx = np.fft.rfftfreq(w)[None,:]
    rr = np.sqrt(fy*fy + fx*fx)
    cut = float(np.quantile(rr, float(genome["low_freq_frac"])))
    low = rr <= cut
    phase = np.where(low, op, phase_noise)
    s = np.fft.irfft2(amp*np.exp(1j*phase), s=x.shape).astype(np.float32)
    for _ in range(IAAFT_ITER):
        s = quantile_remap(s, x)
        sf = np.fft.rfft2(s)
        phase = np.where(low, op, np.angle(sf))
        s = np.fft.irfft2(amp*np.exp(1j*phase), s=x.shape).astype(np.float32)
    return quantile_remap(s, x)


def block_shuffle(x, rng, genome: dict):
    block = int(genome["block_size"])
    h, w = x.shape
    by, bx = math.ceil(h/block), math.ceil(w/block)
    p = np.pad(x, ((0,by*block-h),(0,bx*block-w)), mode="reflect")
    blocks = [p[i*block:(i+1)*block, j*block:(j+1)*block].copy() for i in range(by) for j in range(bx)]
    rng.shuffle(blocks)
    out = np.empty_like(p)
    k = 0
    for i in range(by):
        for j in range(bx):
            out[i*block:(i+1)*block, j*block:(j+1)*block] = blocks[k]
            k += 1
    return out[:h,:w]


def surrogate(x, rng, model: str, genome: dict):
    if model == "phase_iaaft":
        return phase_iaaft(x, rng, genome)
    if model == "block_shuffle":
        return block_shuffle(x, rng, genome)
    raise ValueError(model)


def calibrate(x, genome: dict, model: str, cal_nulls: int, seed_parts: Sequence) -> tuple[np.ndarray,np.ndarray]:
    rng = np.random.default_rng(stable_seed("cal", *seed_parts, model))
    arr = np.asarray([geometry(surrogate(x, rng, model, genome)) for _ in range(cal_nulls)])
    center = arr.mean(0)
    scale = arr.std(0, ddof=1) if cal_nulls > 1 else np.ones(arr.shape[1])
    return center, scale


def empirical_test(x, genome: dict, model: str, test_nulls: int, cal_nulls: int, seeds: Sequence[int], seed_parts: Sequence, progress=None) -> dict:
    center, scale = calibrate(x, genome, model, cal_nulls, seed_parts)
    obs = weighted_std_dist(geometry(x), center, scale, genome)
    base = test_nulls // len(seeds)
    extra = test_nulls % len(seeds)
    vals, chunks = [], []
    for si, seed in enumerate(seeds):
        n = base + (1 if si < extra else 0)
        rng = np.random.default_rng(stable_seed("test", *seed_parts, model, seed))
        vv = []
        for i in range(n):
            st = weighted_std_dist(geometry(surrogate(x, rng, model, genome)), center, scale, genome)
            vals.append(st); vv.append(st)
            if progress and ((i+1) % max(1,n//4) == 0 or i+1 == n):
                progress(seed, i+1, n)
        chunks.append({"seed": int(seed), "null_count": n, "median_stat": float(np.median(vv))})
    a = np.asarray(vals)
    ge = int(np.count_nonzero(a >= obs))
    return {
        "observed_stat": float(obs),
        "null_count": int(len(vals)),
        "ge_count": ge,
        "p_empirical": float((ge+1)/(len(vals)+1)),
        "null_min": float(np.min(a)),
        "null_median": float(np.median(a)),
        "null_max": float(np.max(a)),
        "calibration_nulls": int(cal_nulls),
        "seed_chunks": chunks,
    }


def _cd_matrix(h):
    if all(k in h for k in ("CD1_1","CD1_2","CD2_1","CD2_2")):
        return np.array([[float(h["CD1_1"]), float(h["CD1_2"])], [float(h["CD2_1"]), float(h["CD2_2"])]])
    c1, c2 = float(h.get("CDELT1",1.0)), float(h.get("CDELT2",1.0))
    if all(k in h for k in ("PC1_1","PC1_2","PC2_1","PC2_2")):
        pc = np.array([[float(h["PC1_1"]),float(h["PC1_2"])],[float(h["PC2_1"]),float(h["PC2_2"])]])
        return np.diag([c1,c2]) @ pc
    a = math.radians(float(h.get("CROTA2",0.0)))
    return np.array([[c1*math.cos(a),-c2*math.sin(a)],[c1*math.sin(a),c2*math.cos(a)]])


def world_to_pixel(h, ra_deg: float, dec_deg: float):
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    ra0, dec0 = math.radians(float(h["CRVAL1"])), math.radians(float(h["CRVAL2"]))
    dra = (ra-ra0+math.pi)%(2*math.pi)-math.pi
    cosc = math.sin(dec0)*math.sin(dec)+math.cos(dec0)*math.cos(dec)*math.cos(dra)
    if cosc <= 0:
        return float("nan"), float("nan")
    xi = math.cos(dec)*math.sin(dra)/cosc
    eta = (math.cos(dec0)*math.sin(dec)-math.sin(dec0)*math.cos(dec)*math.cos(dra))/cosc
    dx, dy = np.linalg.solve(_cd_matrix(h), np.degrees([xi,eta]))
    return float(h["CRPIX1"])+dx, float(h["CRPIX2"])+dy


def raw_to_norm(x1: float, y1: float, shape, n=IMAGE_SIZE):
    h,w = shape
    return float((x1-1)/(w-1)*(n-1)), float((y1-1)/(h-1)*(n-1))


def belt_corridor(x: np.ndarray, header: dict, native_shape, stars: list[dict], half_width: float = 10, margin: float = 8):
    pts=[]; rows=[]
    for s in stars:
        px,py=world_to_pixel(header,float(s["ra_deg"]),float(s["dec_deg"]))
        nx,ny=raw_to_norm(px,py,native_shape,x.shape[0])
        pts.append([nx,ny]); rows.append({"name":s["name"],"native_pixel_1based":[px,py],"norm_pixel_0based":[nx,ny]})
    pts=np.asarray(pts,float)
    c=pts.mean(0); q=pts-c
    _,_,vh=np.linalg.svd(q,full_matrices=False); axis=vh[0]
    if axis[0]<0: axis=-axis
    perp=np.array([-axis[1],axis[0]])
    along=q@axis
    amin=float(along.min()-margin); amax=float(along.max()+margin)
    width=max(16,int(math.ceil(amax-amin))+1); height=max(8,int(math.ceil(2*half_width))+1)
    av=np.linspace(amin,amax,width); pv=np.linspace(-half_width,half_width,height)
    A,P=np.meshgrid(av,pv)
    X=c[0]+A*axis[0]+P*perp[0]; Y=c[1]+A*axis[1]+P*perp[1]
    crop=ndimage.map_coordinates(x,[Y,X],order=1,mode="reflect")
    # Crop is already normalized; re-square it without re-reading original values.
    g={"clip_high":99.5,"asinh_scale":6.0}
    # Use geometric resize only to avoid a second photometric transform.
    ch,cw=crop.shape; side=max(ch,cw); py=side-ch; px=side-cw
    crop=np.pad(crop,((py//2,py-py//2),(px//2,px-px//2)),mode="reflect")
    if side!=IMAGE_SIZE:
        crop=ndimage.zoom(crop,(IMAGE_SIZE/side,IMAGE_SIZE/side),order=1)
    crop=np.clip(crop,0,1).astype(np.float32)
    return crop,{"stars":rows,"center_xy":c.tolist(),"axis_xy":axis.tolist(),"sample_shape":[height,width]}


def smooth_correlation(a: np.ndarray, b: np.ndarray, sigma=2.0) -> float:
    aa=ndimage.gaussian_filter(a,sigma,mode="reflect")
    bb=ndimage.gaussian_filter(b,sigma,mode="reflect")
    return corr(aa,bb)
