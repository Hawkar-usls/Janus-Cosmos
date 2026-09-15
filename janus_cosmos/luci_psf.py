from __future__ import annotations

import math
from dataclasses import dataclass, asdict

import numpy as np
from scipy import ndimage


@dataclass(frozen=True)
class PsfSource:
    y: float
    x: float
    peak_snr: float
    area_px: int
    fwhm_minor_px: float
    fwhm_major_px: float
    fwhm_geom_px: float
    elongation: float

    def to_dict(self) -> dict:
        return asdict(self)


def robust_background(image: np.ndarray) -> tuple[float, float]:
    a = np.asarray(image, dtype=float)
    x = a[np.isfinite(a)]
    if x.size < 128:
        raise ValueError("too few finite pixels")
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    sigma = 1.4826 * mad
    if not math.isfinite(sigma) or sigma <= 0:
        sigma = float(np.std(x))
    if not math.isfinite(sigma) or sigma <= 0:
        raise ValueError("non-positive background sigma")
    return med, sigma


def _measure_candidate(signal: np.ndarray, y: int, x: int, sigma: float, radius: int = 6) -> PsfSource | None:
    y0, y1 = y - radius, y + radius + 1
    x0, x1 = x - radius, x + radius + 1
    if y0 < 0 or x0 < 0 or y1 > signal.shape[0] or x1 > signal.shape[1]:
        return None
    sub = np.asarray(signal[y0:y1, x0:x1], dtype=float)
    peak = float(sub[radius, radius])
    if not math.isfinite(peak) or peak < 4.0 * sigma:
        return None
    yy, xx = np.indices(sub.shape, dtype=float)
    rr = np.hypot(xx - radius, yy - radius)
    core = (sub > 2.0 * sigma) & (rr <= radius)
    labels, _ = ndimage.label(core)
    lab = int(labels[radius, radius])
    if lab <= 0:
        return None
    comp = labels == lab
    area = int(np.count_nonzero(comp))
    if area < 5:
        return None
    weights = np.where(comp, np.clip(sub, 0.0, None), 0.0)
    sw = float(weights.sum())
    if sw <= 0:
        return None
    cx = float((weights * xx).sum() / sw)
    cy = float((weights * yy).sum() / sw)
    dx, dy = xx - cx, yy - cy
    cxx = float((weights * dx * dx).sum() / sw)
    cyy = float((weights * dy * dy).sum() / sw)
    cxy = float((weights * dx * dy).sum() / sw)
    vals = np.linalg.eigvalsh(np.array([[cxx, cxy], [cxy, cyy]], dtype=float))
    if not np.all(np.isfinite(vals)) or float(vals[0]) <= 0:
        return None
    fmin = 2.354820045 * math.sqrt(float(vals[0]))
    fmax = 2.354820045 * math.sqrt(float(vals[1]))
    elong = fmax / fmin
    if fmin < 0.8 or fmax > 12.0 or elong > 4.0:
        return None
    return PsfSource(
        y=float(y0 + cy), x=float(x0 + cx), peak_snr=peak / sigma, area_px=area,
        fwhm_minor_px=fmin, fwhm_major_px=fmax,
        fwhm_geom_px=math.sqrt(fmin * fmax), elongation=elong,
    )


def detect_psf_sources(
    image: np.ndarray,
    *,
    threshold_sigma: float = 5.0,
    smooth_sigma_px: float = 1.0,
    edge_px: int = 12,
    max_sources: int = 256,
) -> list[PsfSource]:
    a = np.asarray(image, dtype=float)
    med, sigma = robust_background(a)
    signal = np.where(np.isfinite(a), a - med, 0.0)
    sm = ndimage.gaussian_filter(signal, smooth_sigma_px, mode="nearest")
    _, sm_sigma = robust_background(sm)
    local_max = sm == ndimage.maximum_filter(sm, size=5, mode="nearest")
    mask = local_max & (sm >= threshold_sigma * sm_sigma)
    if edge_px > 0:
        mask[:edge_px, :] = False; mask[-edge_px:, :] = False
        mask[:, :edge_px] = False; mask[:, -edge_px:] = False
    ys, xs = np.nonzero(mask)
    order = np.argsort(sm[ys, xs])[::-1]
    out: list[PsfSource] = []
    for oi in order:
        src = _measure_candidate(signal, int(ys[oi]), int(xs[oi]), sigma)
        if src is None:
            continue
        if any((src.x - q.x) ** 2 + (src.y - q.y) ** 2 < 16.0 for q in out):
            continue
        out.append(src)
        if len(out) >= max_sources:
            break
    return out


def _inject_gaussian(a: np.ndarray, y: float, x: float, fwhm_px: float, peak: float) -> None:
    sig = fwhm_px / 2.354820045
    rad = max(4, int(math.ceil(4.0 * sig)))
    y0, y1 = max(0, int(y) - rad), min(a.shape[0], int(y) + rad + 1)
    x0, x1 = max(0, int(x) - rad), min(a.shape[1], int(x) + rad + 1)
    yy, xx = np.indices((y1 - y0, x1 - x0), dtype=float)
    yy += y0; xx += x0
    a[y0:y1, x0:x1] += peak * np.exp(-0.5 * (((xx - x) / sig) ** 2 + ((yy - y) / sig) ** 2))


def _safe_positions(shape: tuple[int, int], n: int, rng: np.random.Generator, occupied: list[tuple[float, float]], edge: int = 20) -> list[tuple[float, float]]:
    h, w = shape
    if h <= 2 * edge or w <= 2 * edge:
        raise ValueError("image too small for injection gate")
    out: list[tuple[float, float]] = []
    for _ in range(10000):
        if len(out) >= n:
            break
        y = float(rng.uniform(edge, h - edge)); x = float(rng.uniform(edge, w - edge))
        if all((x - ox) ** 2 + (y - oy) ** 2 >= 12.0 ** 2 for oy, ox in occupied + out):
            out.append((y, x))
    if len(out) != n:
        raise RuntimeError(f"could place only {len(out)}/{n} injections")
    return out


def _recovery_fraction(injected: list[dict], detected: list[PsfSource], radius_px: float = 2.5) -> tuple[float, list[bool]]:
    flags = []
    r2 = radius_px * radius_px
    for inj in injected:
        ok = any((q.x - inj["x"]) ** 2 + (q.y - inj["y"]) ** 2 <= r2 for q in detected)
        flags.append(bool(ok))
    return (float(sum(flags) / len(flags)) if flags else 0.0), flags


def injection_recovery_gate(
    image: np.ndarray,
    *,
    seed: int = 20260815,
    fwhm_grid_px: tuple[float, ...] = (1.5, 2.0, 3.0),
    snr_grid: tuple[float, ...] = (6.0, 8.0, 12.0),
    replicates: int = 4,
    min_all_recovery: float = 0.80,
    min_high_snr_recovery: float = 0.90,
    max_hot_pixel_acceptance: float = 0.05,
) -> dict:
    a = np.asarray(image, dtype=float)
    med, sigma = robust_background(a)
    native = detect_psf_sources(a)
    occupied = [(q.y, q.x) for q in native]
    rng = np.random.default_rng(int(seed))
    combos = [(f, s) for f in fwhm_grid_px for s in snr_grid for _ in range(int(replicates))]
    pos = _safe_positions(a.shape, len(combos), rng, occupied)
    injected: list[dict] = []
    star_frame = np.array(a, copy=True)
    for (fwhm, snr), (y, x) in zip(combos, pos):
        _inject_gaussian(star_frame, y, x, fwhm, snr * sigma)
        injected.append({"y": y, "x": x, "fwhm_px": fwhm, "peak_snr": snr})
    star_det = detect_psf_sources(star_frame)
    all_rec, flags = _recovery_fraction(injected, star_det)
    high_idx = [i for i, z in enumerate(injected) if z["peak_snr"] >= 8.0]
    high_rec = float(sum(flags[i] for i in high_idx) / len(high_idx)) if high_idx else 0.0

    hot_pos = _safe_positions(a.shape, len(combos), rng, occupied + pos)
    hot = np.array(a, copy=True)
    hot_inj = []
    for (_, snr), (y, x) in zip(combos, hot_pos):
        iy, ix = int(round(y)), int(round(x))
        hot[iy, ix] += snr * sigma
        hot_inj.append({"y": float(iy), "x": float(ix), "peak_snr": snr})
    hot_det = detect_psf_sources(hot)
    hot_accept, _ = _recovery_fraction(hot_inj, hot_det)

    passed = bool(all_rec >= min_all_recovery and high_rec >= min_high_snr_recovery and hot_accept <= max_hot_pixel_acceptance)
    return {
        "schema": "janus.cosmos.luci_psf_injection_gate.v1",
        "passed": passed,
        "seed": int(seed),
        "background_median": med,
        "background_sigma": sigma,
        "native_psf_source_count": len(native),
        "injected_star_count": len(injected),
        "star_recovery_fraction_all": all_rec,
        "star_recovery_fraction_snr_ge_8": high_rec,
        "hot_pixel_injected_count": len(hot_inj),
        "hot_pixel_acceptance_fraction": hot_accept,
        "thresholds": {
            "min_all_recovery": min_all_recovery,
            "min_high_snr_recovery": min_high_snr_recovery,
            "max_hot_pixel_acceptance": max_hot_pixel_acceptance,
            "min_component_area_px": 5,
            "min_minor_fwhm_px": 0.8,
            "max_elongation": 4.0,
        },
        "fwhm_grid_px": list(fwhm_grid_px),
        "snr_grid": list(snr_grid),
        "replicates_per_grid_cell": int(replicates),
    }
