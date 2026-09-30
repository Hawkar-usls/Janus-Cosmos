from __future__ import annotations

import math

import numpy as np

from .luci_psf import PsfSource, _inject_gaussian, _measure_candidate, detect_psf_sources, robust_background


def measure_psf_at(image: np.ndarray, y: float, x: float, *, search_radius_px: int = 2) -> PsfSource | None:
    a = np.asarray(image, dtype=float)
    med, sigma = robust_background(a)
    signal = np.where(np.isfinite(a), a - med, 0.0)
    iy, ix = int(round(y)), int(round(x))
    y0, y1 = max(0, iy-search_radius_px), min(a.shape[0], iy+search_radius_px+1)
    x0, x1 = max(0, ix-search_radius_px), min(a.shape[1], ix+search_radius_px+1)
    if y0 >= y1 or x0 >= x1:
        return None
    sub = signal[y0:y1, x0:x1]
    k = int(np.argmax(sub))
    dy, dx = np.unravel_index(k, sub.shape)
    py, px = y0 + int(dy), x0 + int(dx)
    src = _measure_candidate(signal, py, px, sigma)
    if src is None:
        return None
    if (src.x-x)**2 + (src.y-y)**2 > 2.5**2:
        return None
    return src


def _clean_positions(
    image: np.ndarray,
    n: int,
    rng: np.random.Generator,
    *,
    edge_px: int = 20,
    min_separation_px: float = 16.0,
) -> list[tuple[float, float]]:
    a = np.asarray(image, dtype=float)
    med, sigma = robust_background(a)
    signal = np.where(np.isfinite(a), a-med, 0.0)
    h, w = a.shape
    out: list[tuple[float, float]] = []
    for _ in range(50000):
        if len(out) >= n:
            break
        y = int(rng.integers(edge_px, h-edge_px)); x = int(rng.integers(edge_px, w-edge_px))
        patch = signal[y-5:y+6, x-5:x+6]
        if patch.shape != (11,11) or not np.all(np.isfinite(patch)):
            continue
        if float(np.max(patch)) >= 3.0*sigma:
            continue
        if any((x-ox)**2 + (y-oy)**2 < min_separation_px**2 for oy, ox in out):
            continue
        out.append((float(y), float(x)))
    if len(out) != n:
        raise RuntimeError(f"clean-background placement failed: {len(out)}/{n}")
    return out


def estimate_native_psf_fwhm(image: np.ndarray) -> tuple[float | None, int]:
    src = detect_psf_sources(image, max_sources=2048)
    good = [q.fwhm_geom_px for q in src if q.peak_snr >= 8.0 and q.elongation <= 2.5 and 0.8 <= q.fwhm_geom_px <= 10.0]
    if len(good) < 2:
        good = [q.fwhm_geom_px for q in src if 0.8 <= q.fwhm_geom_px <= 10.0]
    if len(good) < 2:
        return None, len(src)
    return float(np.median(good)), len(src)


def psf_relative_injection_recovery_gate(
    image: np.ndarray,
    *,
    seed: int = 20260815,
    fwhm_ratios: tuple[float, ...] = (0.8, 1.0, 1.2),
    snr_grid: tuple[float, ...] = (6.0, 8.0, 12.0),
    replicates: int = 4,
    min_all_recovery: float = 0.80,
    min_high_snr_recovery: float = 0.90,
    max_hot_pixel_acceptance: float = 0.05,
) -> dict:
    a = np.asarray(image, dtype=float)
    med, sigma = robust_background(a)
    psf_median, native_n = estimate_native_psf_fwhm(a)
    if psf_median is None:
        return {
            "schema": "janus.cosmos.luci_psf_relative_injection_gate.v1",
            "passed": False,
            "reason": "INSUFFICIENT_NATIVE_PSF_REFERENCE",
            "native_psf_source_count": native_n,
            "background_median": med,
            "background_sigma": sigma,
        }
    base_fwhm = max(2.0, min(8.0, float(psf_median)))
    combos = [(r, s) for r in fwhm_ratios for s in snr_grid for _ in range(int(replicates))]
    rng = np.random.default_rng(int(seed))
    star_pos = _clean_positions(a, len(combos), rng)
    star = np.array(a, copy=True)
    injected = []
    for (ratio, snr), (y, x) in zip(combos, star_pos):
        fwhm = max(2.0, min(10.0, base_fwhm*ratio))
        _inject_gaussian(star, y, x, fwhm, snr*sigma)
        injected.append({"y":y,"x":x,"ratio":ratio,"fwhm_px":fwhm,"peak_snr":snr})
    flags = [measure_psf_at(star, z["y"], z["x"]) is not None for z in injected]
    all_rec = float(sum(flags)/len(flags)) if flags else 0.0
    hi = [i for i,z in enumerate(injected) if z["peak_snr"] >= 8.0]
    high_rec = float(sum(flags[i] for i in hi)/len(hi)) if hi else 0.0

    hot_pos = _clean_positions(a, len(combos), rng)
    hot = np.array(a, copy=True)
    hot_flags = []
    for (_, snr), (y, x) in zip(combos, hot_pos):
        iy, ix = int(round(y)), int(round(x))
        hot[iy, ix] += snr*sigma
        hot_flags.append((float(iy), float(ix)))
    accepted_hot = sum(measure_psf_at(hot, y, x) is not None for y,x in hot_flags)
    hot_accept = float(accepted_hot/len(hot_flags)) if hot_flags else 0.0

    passed = bool(all_rec >= min_all_recovery and high_rec >= min_high_snr_recovery and hot_accept <= max_hot_pixel_acceptance)
    return {
        "schema": "janus.cosmos.luci_psf_relative_injection_gate.v1",
        "passed": passed,
        "seed": int(seed),
        "background_median": med,
        "background_sigma": sigma,
        "native_psf_source_count": native_n,
        "native_psf_median_fwhm_px": float(psf_median),
        "injection_base_fwhm_px": base_fwhm,
        "fwhm_ratios": list(fwhm_ratios),
        "snr_grid": list(snr_grid),
        "replicates_per_grid_cell": int(replicates),
        "injected_star_count": len(injected),
        "star_recovery_fraction_all": all_rec,
        "star_recovery_fraction_snr_ge_8": high_rec,
        "hot_pixel_injected_count": len(hot_flags),
        "hot_pixel_acceptance_fraction": hot_accept,
        "thresholds": {
            "min_all_recovery": min_all_recovery,
            "min_high_snr_recovery": min_high_snr_recovery,
            "max_hot_pixel_acceptance": max_hot_pixel_acceptance,
        },
    }
