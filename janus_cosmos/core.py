from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Iterable

import numpy as np
from scipy import ndimage

IMAGE_SIZE = 128
SCALES = (1.0, 2.0, 4.0)
ORIENTATIONS = tuple(range(0, 180, 30))
BLOCK_SIZE = 16
LOW_FREQ_FRACTION = 0.10
IAAFT_ITERATIONS = 3
CALIBRATION_NULLS = 128
NULL_MODELS = ("phase_iaaft", "block_shuffle")

FEATURE_NAMES = (
    "directional_s1",
    "directional_s2",
    "directional_s4",
    "rot90_corr",
    "rot180_corr",
    "gradient_anisotropy",
    "high_frequency_energy",
    "fourier_angular_anisotropy",
    "component_count_q80",
    "component_count_q90",
    "component_count_q95",
    "largest_component_q80",
    "largest_component_q90",
    "largest_component_q95",
)


@dataclass(frozen=True)
class GateConfig:
    alpha_family: float = 0.05
    image_size: int = IMAGE_SIZE
    calibration_nulls: int = CALIBRATION_NULLS
    block_size: int = BLOCK_SIZE
    low_freq_fraction: float = LOW_FREQ_FRACTION
    iaaft_iterations: int = IAAFT_ITERATIONS


def stable_seed(*parts: object) -> int:
    blob = "|".join(str(x) for x in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(blob).digest()[:8], "big", signed=False)


def normalize_image(image: np.ndarray, size: int = IMAGE_SIZE) -> np.ndarray:
    x = np.asarray(image, dtype=np.float32)
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    if x.ndim != 2:
        raise ValueError("Janus Cosmos expects a 2-D image")

    h, w = x.shape
    side = max(h, w)
    py = side - h
    px = side - w
    x = np.pad(
        x,
        ((py // 2, py - py // 2), (px // 2, px - px // 2)),
        mode="constant",
        constant_values=float(np.median(x)),
    )
    if side != size:
        x = ndimage.zoom(x, (size / side, size / side), order=1)

    lo, hi = np.percentile(x, [1.0, 99.5])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.zeros((size, size), dtype=np.float32)
    x = np.clip((x - lo) / (hi - lo), 0.0, 1.0)
    x = np.arcsinh(6.0 * x) / np.arcsinh(6.0)
    return x.astype(np.float32)


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    aa = np.asarray(a, dtype=np.float64).ravel()
    bb = np.asarray(b, dtype=np.float64).ravel()
    aa = aa - aa.mean()
    bb = bb - bb.mean()
    den = float(np.linalg.norm(aa) * np.linalg.norm(bb))
    return float(np.dot(aa, bb) / den) if den else 0.0


def _shift_corr(smoothed: np.ndarray, sigma: float, degrees: int) -> float:
    size_y, size_x = smoothed.shape
    angle = math.radians(degrees)
    dx = int(round(2 * sigma * math.cos(angle)))
    dy = int(round(2 * sigma * math.sin(angle)))

    ax0, ax1, bx0, bx1 = (dx, size_x, 0, size_x - dx) if dx >= 0 else (0, size_x + dx, -dx, size_x)
    ay0, ay1, by0, by1 = (dy, size_y, 0, size_y - dy) if dy >= 0 else (0, size_y + dy, -dy, size_y)
    if ax1 - ax0 < 8 or ay1 - ay0 < 8:
        return 0.0
    a = smoothed[ay0:ay1, ax0:ax1]
    b = smoothed[by0:by1, bx0:bx1]
    return abs(_corr(a, b))


def directional_by_scale(image: np.ndarray) -> tuple[float, ...]:
    out = []
    for sigma in SCALES:
        smoothed = ndimage.gaussian_filter(image, sigma=sigma, mode="reflect")
        vals = [_shift_corr(smoothed, sigma, deg) for deg in ORIENTATIONS]
        out.append(float(np.mean(vals)))
    return tuple(out)


def _gradient_anisotropy(image: np.ndarray) -> float:
    gy, gx = np.gradient(ndimage.gaussian_filter(image, 1.0, mode="reflect"))
    jxx = float(np.mean(gx * gx))
    jyy = float(np.mean(gy * gy))
    jxy = float(np.mean(gx * gy))
    trace = jxx + jyy
    disc = math.sqrt(max((jxx - jyy) ** 2 + 4 * jxy * jxy, 0.0))
    l1 = 0.5 * (trace + disc)
    l2 = 0.5 * (trace - disc)
    return float((l1 - l2) / max(l1 + l2, 1e-12))


def _fourier_angular_anisotropy(image: np.ndarray) -> float:
    f = np.fft.fftshift(np.fft.fft2(image - float(image.mean())))
    power = np.abs(f) ** 2
    h, w = power.shape
    yy, xx = np.indices(power.shape)
    cy, cx = (h - 1) / 2.0, (w - 1) / 2.0
    dx = xx - cx
    dy = yy - cy
    radius = np.sqrt(dx * dx + dy * dy)
    theta = np.mod(np.arctan2(dy, dx), np.pi)
    mask = radius > max(h, w) * 0.03
    bins = np.linspace(0.0, np.pi, 13)
    vals = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = mask & (theta >= lo) & (theta < hi)
        vals.append(float(power[m].mean()) if np.any(m) else 0.0)
    vals = np.asarray(vals, dtype=np.float64)
    mean = float(vals.mean())
    return float(vals.std() / max(mean, 1e-12))


def _component_features(image: np.ndarray) -> list[float]:
    counts = []
    largest = []
    npx = image.size
    for q in (0.80, 0.90, 0.95):
        threshold = float(np.quantile(image, q))
        labels, count = ndimage.label(image >= threshold)
        sizes = np.bincount(labels.ravel())[1:] if count else np.array([], dtype=np.int64)
        counts.append(float(np.log1p(count)))
        largest.append(float(sizes.max() / npx) if sizes.size else 0.0)
    return counts + largest


def geometry_vector(image: np.ndarray) -> np.ndarray:
    x = np.asarray(image, dtype=np.float32)
    ds = directional_by_scale(x)
    rot90 = _corr(x, np.rot90(x, 1))
    rot180 = _corr(x, np.rot90(x, 2))
    grad_aniso = _gradient_anisotropy(x)
    smooth = ndimage.gaussian_filter(x, sigma=4.0, mode="reflect")
    hf = float(np.mean((x - smooth) ** 2) / max(np.mean(x ** 2), 1e-12))
    fang = _fourier_angular_anisotropy(x)
    comps = _component_features(x)
    out = np.asarray([*ds, rot90, rot180, grad_aniso, hf, fang, *comps], dtype=np.float64)
    if out.shape != (len(FEATURE_NAMES),):
        raise AssertionError("geometry feature shape mismatch")
    return out


def quantile_remap(source: np.ndarray, reference: np.ndarray) -> np.ndarray:
    src = np.asarray(source, dtype=np.float64).ravel()
    ref = np.sort(np.asarray(reference, dtype=np.float64).ravel())
    order = np.argsort(src, kind="mergesort")
    out = np.empty_like(src)
    out[order] = ref
    return out.reshape(source.shape).astype(np.float32)


def phase_iaaft_surrogate(
    image: np.ndarray,
    rng: np.random.Generator,
    low_freq_fraction: float = LOW_FREQ_FRACTION,
    iterations: int = IAAFT_ITERATIONS,
) -> np.ndarray:
    """IAAFT-like surrogate preserving marginal values and approximate spectrum.

    Random phase is inherited from an FFT of a real noise field, so the rFFT
    representation remains compatible with a real-valued inverse transform.
    Low spatial-frequency phase is frozen to the observation to retain broad
    morphology while higher-frequency phase is randomized.
    """
    x = np.asarray(image, dtype=np.float32)
    h, w = x.shape
    original = np.fft.rfft2(x)
    amplitude = np.abs(original)
    original_phase = np.angle(original)

    noise_phase = np.angle(np.fft.rfft2(rng.normal(size=x.shape)))
    fy = np.fft.fftfreq(h)[:, None]
    fx = np.fft.rfftfreq(w)[None, :]
    radius = np.sqrt(fy * fy + fx * fx)
    cutoff = float(np.quantile(radius, low_freq_fraction))
    low = radius <= cutoff
    phase = np.where(low, original_phase, noise_phase)
    surrogate = np.fft.irfft2(amplitude * np.exp(1j * phase), s=x.shape).astype(np.float32)

    for _ in range(max(0, int(iterations))):
        surrogate = quantile_remap(surrogate, x)
        sf = np.fft.rfft2(surrogate)
        phase = np.angle(sf)
        phase = np.where(low, original_phase, phase)
        surrogate = np.fft.irfft2(amplitude * np.exp(1j * phase), s=x.shape).astype(np.float32)

    return quantile_remap(surrogate, x)


def block_shuffle_surrogate(
    image: np.ndarray,
    rng: np.random.Generator,
    block_size: int = BLOCK_SIZE,
) -> np.ndarray:
    x = np.asarray(image, dtype=np.float32)
    h, w = x.shape
    by = math.ceil(h / block_size)
    bx = math.ceil(w / block_size)
    padded = np.pad(x, ((0, by * block_size - h), (0, bx * block_size - w)), mode="reflect")
    blocks = [
        padded[iy * block_size:(iy + 1) * block_size, ix * block_size:(ix + 1) * block_size].copy()
        for iy in range(by)
        for ix in range(bx)
    ]
    rng.shuffle(blocks)
    out = np.empty_like(padded)
    k = 0
    for iy in range(by):
        for ix in range(bx):
            out[iy * block_size:(iy + 1) * block_size, ix * block_size:(ix + 1) * block_size] = blocks[k]
            k += 1
    return out[:h, :w]


def pixel_permutation_surrogate(image: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    flat = np.asarray(image, dtype=np.float32).ravel().copy()
    rng.shuffle(flat)
    return flat.reshape(image.shape)


def standardized_distance(vec: np.ndarray, center: np.ndarray, scale: np.ndarray) -> float:
    z = (vec - center) / np.maximum(scale, 1e-9)
    return float(np.sqrt(np.mean(z * z)))


def empirical_p(observed_stat: float, null_stats: Iterable[float]) -> dict:
    vals = np.asarray(list(null_stats), dtype=np.float64)
    if vals.size == 0:
        raise ValueError("null distribution is empty")
    ge = int(np.count_nonzero(vals >= observed_stat))
    p = float((ge + 1) / (vals.size + 1))
    return {
        "observed_stat": float(observed_stat),
        "null_count": int(vals.size),
        "ge_count": ge,
        "p_empirical": p,
        "null_min": float(vals.min()),
        "null_median": float(np.median(vals)),
        "null_max": float(vals.max()),
    }


def bonferroni_alpha(alpha_family: float, filter_count: int, null_model_count: int = len(NULL_MODELS)) -> float:
    tests = max(1, int(filter_count) * int(null_model_count))
    return float(alpha_family / tests)


def minimum_test_nulls(alpha: float) -> int:
    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must be between 0 and 1")
    return max(1, int(math.floor(1.0 / alpha)))


def distribute_nulls(total: int, seeds: list[int]) -> list[int]:
    if total <= 0 or not seeds:
        raise ValueError("positive total and at least one seed are required")
    base = total // len(seeds)
    extra = total % len(seeds)
    return [base + (1 if i < extra else 0) for i in range(len(seeds))]


def analyze_against_null_model(
    image: np.ndarray,
    *,
    target: str,
    filter_name: str,
    model: str,
    test_nulls: int,
    seeds: list[int],
    calibration_nulls: int = CALIBRATION_NULLS,
    config: GateConfig | None = None,
) -> dict:
    cfg = config or GateConfig()
    obs_vec = geometry_vector(image)

    if model == "phase_iaaft":
        surrogate_fn = lambda rng: phase_iaaft_surrogate(image, rng, cfg.low_freq_fraction, cfg.iaaft_iterations)
    elif model == "block_shuffle":
        surrogate_fn = lambda rng: block_shuffle_surrogate(image, rng, cfg.block_size)
    elif model == "pixel_permutation":
        surrogate_fn = lambda rng: pixel_permutation_surrogate(image, rng)
    else:
        raise ValueError(f"unknown null model: {model}")

    cal_rng = np.random.default_rng(stable_seed("cal", target, filter_name, model, seeds[0]))
    calibration = np.asarray(
        [geometry_vector(surrogate_fn(cal_rng)) for _ in range(calibration_nulls)],
        dtype=np.float64,
    )
    center = calibration.mean(axis=0)
    scale = calibration.std(axis=0, ddof=1)
    observed_stat = standardized_distance(obs_vec, center, scale)

    test_stats: list[float] = []
    chunks = []
    for seed, count in zip(seeds, distribute_nulls(test_nulls, seeds)):
        rng = np.random.default_rng(stable_seed("test", target, filter_name, model, seed))
        chunk_stats = []
        for _ in range(count):
            stat = standardized_distance(geometry_vector(surrogate_fn(rng)), center, scale)
            test_stats.append(stat)
            chunk_stats.append(stat)
        chunks.append({
            "seed": int(seed),
            "null_count": int(count),
            "median_stat": float(np.median(chunk_stats)) if chunk_stats else None,
        })

    result = empirical_p(observed_stat, test_stats)
    result.update({
        "model": model,
        "feature_names": list(FEATURE_NAMES),
        "observed_features": [float(x) for x in obs_vec],
        "calibration_nulls": int(calibration_nulls),
        "seed_chunks": chunks,
    })
    return result


def analyze_image(
    image: np.ndarray,
    *,
    target: str,
    filter_name: str,
    test_nulls: int,
    seeds: list[int],
    alpha: float,
    include_legacy: bool = True,
    config: GateConfig | None = None,
) -> dict:
    cfg = config or GateConfig()
    x = normalize_image(image, cfg.image_size)
    phase = analyze_against_null_model(
        x,
        target=target,
        filter_name=filter_name,
        model="phase_iaaft",
        test_nulls=test_nulls,
        seeds=seeds,
        calibration_nulls=cfg.calibration_nulls,
        config=cfg,
    )
    block = analyze_against_null_model(
        x,
        target=target,
        filter_name=filter_name,
        model="block_shuffle",
        test_nulls=test_nulls,
        seeds=seeds,
        calibration_nulls=cfg.calibration_nulls,
        config=cfg,
    )

    out = {
        "phase_iaaft": phase,
        "block_shuffle": block,
        "alpha_corrected": float(alpha),
        "robust_candidate": bool(phase["p_empirical"] < alpha and block["p_empirical"] < alpha),
        "robust_rule": "phase_iaaft p<alpha AND block_shuffle p<alpha",
    }
    if include_legacy:
        legacy = analyze_against_null_model(
            x,
            target=target,
            filter_name=filter_name,
            model="pixel_permutation",
            test_nulls=min(test_nulls, 512),
            seeds=seeds,
            calibration_nulls=min(cfg.calibration_nulls, 64),
            config=cfg,
        )
        out["legacy_pixel_permutation"] = legacy
    return out
