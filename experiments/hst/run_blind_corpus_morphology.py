from __future__ import annotations

import argparse
import hashlib
import json
import math
import tempfile
from pathlib import Path

import numpy as np

from run_blind_corpus_variant import download_variant, base

NULLS = 256
BLOCK = 16
LOW_FREQ_FRACTION = 0.10
PHASE_STRENGTH = 1.0


def quantile_remap(source: np.ndarray, reference: np.ndarray) -> np.ndarray:
    src = source.ravel()
    ref = np.sort(reference.ravel())
    order = np.argsort(src, kind="mergesort")
    ranks = np.empty_like(order, dtype=np.float64)
    ranks[order] = np.linspace(0.0, 1.0, len(src), endpoint=True)
    pos = ranks * (len(ref) - 1)
    lo = np.floor(pos).astype(np.int64)
    hi = np.ceil(pos).astype(np.int64)
    frac = pos - lo
    mapped = ref[lo] * (1.0 - frac) + ref[hi] * frac
    return mapped.reshape(source.shape).astype(np.float32)


def morphology_preserving_surrogate(image: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Preserve marginal intensity and the Fourier power envelope.

    The low spatial frequencies retain their original phase so the broad morphology
    remains represented; higher frequencies receive randomized phase. rfft2/irfft2
    preserves the real-valued conjugate symmetry exactly. Rank remapping restores the
    empirical intensity distribution. This is a statistical null, not a physical model.
    """
    h, w = image.shape
    fy = np.fft.fftfreq(h)[:, None]
    fx = np.fft.rfftfreq(w)[None, :]
    radius = np.sqrt(fy * fy + fx * fx)
    cutoff = np.quantile(radius, LOW_FREQ_FRACTION)

    F = np.fft.rfft2(image)
    amplitude = np.abs(F)
    phase = np.angle(F)
    random_phase = rng.uniform(-np.pi, np.pi, size=phase.shape)

    alpha = np.clip((radius - cutoff) / max(cutoff, 1e-9), 0.0, 1.0)
    alpha = alpha ** PHASE_STRENGTH
    new_phase = phase * (1.0 - alpha) + random_phase * alpha
    new_phase[0, 0] = phase[0, 0]
    out = np.fft.irfft2(amplitude * np.exp(1j * new_phase), s=image.shape).astype(np.float32)
    return quantile_remap(out, image)


def block_shuffle_surrogate(image: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Preserve local spatial correlation inside blocks while destroying global layout."""
    h, w = image.shape
    by = math.ceil(h / BLOCK)
    bx = math.ceil(w / BLOCK)
    padded = np.pad(image, ((0, by * BLOCK - h), (0, bx * BLOCK - w)), mode="reflect")
    blocks = []
    for iy in range(by):
        for ix in range(bx):
            blocks.append(padded[iy * BLOCK:(iy + 1) * BLOCK, ix * BLOCK:(ix + 1) * BLOCK].copy())
    rng.shuffle(blocks)
    out = np.zeros_like(padded)
    k = 0
    for iy in range(by):
        for ix in range(bx):
            out[iy * BLOCK:(iy + 1) * BLOCK, ix * BLOCK:(ix + 1) * BLOCK] = blocks[k]
            k += 1
    return out[:h, :w]


def power_spectrum_rmse(a: np.ndarray, b: np.ndarray) -> float:
    pa = np.abs(np.fft.fftshift(np.fft.fft2(a))) ** 2
    pb = np.abs(np.fft.fftshift(np.fft.fft2(b))) ** 2
    scale = max(float(np.mean(pa)), 1e-12)
    return float(np.sqrt(np.mean((pa - pb) ** 2)) / scale)


def local_corr_delta(a: np.ndarray, b: np.ndarray) -> float:
    def corr(x: np.ndarray, dx: int, dy: int) -> float:
        x1 = x[max(0, dy): min(x.shape[0], x.shape[0] + dy), max(0, dx): min(x.shape[1], x.shape[1] + dx)]
        x2 = x[max(0, -dy): min(x.shape[0], x.shape[0] - dy), max(0, -dx): min(x.shape[1], x.shape[1] - dx)]
        x1 = x1.ravel() - float(x1.mean()); x2 = x2.ravel() - float(x2.mean())
        den = float(np.sqrt(np.sum(x1 * x1) * np.sum(x2 * x2)))
        return float(np.sum(x1 * x2) / den) if den else 0.0
    vals_a = [corr(a, 1, 0), corr(a, 0, 1), corr(a, 1, 1)]
    vals_b = [corr(b, 1, 0), corr(b, 0, 1), corr(b, 1, 1)]
    return float(abs(np.mean(vals_a) - np.mean(vals_b)))


def smooth_residual_delta(a: np.ndarray, b: np.ndarray) -> float:
    sa = base.gaussian_filter(a, sigma=8.0, mode="reflect")
    sb = base.gaussian_filter(b, sigma=8.0, mode="reflect")
    scale = max(float(np.mean(np.abs(sa))), 1e-9)
    return float(np.mean(np.abs(sa - sb)) / scale)


def empirical_result(observed: float, null_scores: list[float]) -> dict:
    ge = sum(v >= observed for v in null_scores)
    p = (ge + 1) / (len(null_scores) + 1)
    return {
        "observed_score": float(observed),
        "null_min": float(np.min(null_scores)),
        "null_median": float(np.median(null_scores)),
        "null_max": float(np.max(null_scores)),
        "ge_count": int(ge),
        "p_empirical": float(p),
        "candidate": bool(p < 0.05),
    }


def analyze_filter(image: np.ndarray, rng: np.random.Generator, target: str, filter_name: str) -> dict:
    observed = base.score(image)
    permutation: list[float] = []
    morphology: list[float] = []
    blocks: list[float] = []
    diagnostics = []
    for i in range(NULLS):
        p = image.ravel().copy()
        rng.shuffle(p)
        p = p.reshape(image.shape)
        m = morphology_preserving_surrogate(image, rng)
        b = block_shuffle_surrogate(image, rng)
        permutation.append(base.score(p))
        morphology.append(base.score(m))
        blocks.append(base.score(b))
        if i in (0, NULLS // 2, NULLS - 1):
            diagnostics.append({
                "null_index": i + 1,
                "morphology_power_spectrum_rmse": power_spectrum_rmse(image, m),
                "morphology_local_corr_delta": local_corr_delta(image, m),
                "morphology_smooth_residual_delta": smooth_residual_delta(image, m),
            })

    legacy = empirical_result(observed, permutation)
    morph = empirical_result(observed, morphology)
    block = empirical_result(observed, blocks)
    robust = bool(morph["candidate"] and block["candidate"])
    result = {
        "legacy_pixel_permutation": legacy,
        "morphology_preserving_phase": morph,
        "local_block_shuffle": block,
        "morphology_preserving_diagnostics": diagnostics,
        "robust_candidate": robust,
        "robust_rule": "phase_preserving_null AND local_block_shuffle both p<0.05",
    }
    base.emit("morphology_filter_analyzed", target=target, filter=filter_name, **result)
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=["direct", "mast_api"], default="direct")
    ap.add_argument("--output-dir", default="results/morphology")
    args = ap.parse_args()

    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    base.EVENT_LOG = out_dir / "janus-cosmos-morphology-events.jsonl"
    base.RECEIPT = out_dir / "janus-cosmos-morphology-receipt.json"
    base.EVENT_LOG.unlink(missing_ok=True)

    manifest = json.loads(Path("data/hst_blind_corpus.json").read_text(encoding="utf-8"))
    rng = np.random.default_rng(base.SEED)
    base.emit(
        "run_started",
        schema="janus.cosmos.morphology_event.v0.1",
        backend=args.backend,
        seed=base.SEED,
        nulls=NULLS,
        block_size=BLOCK,
        low_freq_fraction=LOW_FREQ_FRACTION,
        phase_strength=PHASE_STRENGTH,
        semantic_analysis=False,
        ocr=False,
        face_search=False,
        cipher_search=False,
        post_hoc_tuning=False,
    )

    targets = []
    with tempfile.TemporaryDirectory() as td:
        for target in manifest["targets"]:
            t = {"target": target["target"], "target_class": target["class"], "filters": {}, "source_products": []}
            for item in target["filters"]:
                path = Path(td) / f"{target['target']}_{item['filter']}.fits"
                metadata = download_variant(item["url"], path, args.backend)
                image = base.read_image(path)
                t["filters"][item["filter"]] = analyze_filter(image, rng, target["target"], item["filter"])
                t["source_products"].append({"filter": item["filter"], "band": item["band"], **metadata})
            robust = [f for f, r in t["filters"].items() if r["robust_candidate"]]
            t["robust_passing_filters"] = robust
            t["robust_cross_band_candidate"] = len(robust) >= 2
            targets.append(t)
            base.emit("target_completed", target=target["target"], robust_passing_filters=robust, robust_cross_band_candidate=t["robust_cross_band_candidate"])

    receipt = {
        "schema": "janus.cosmos.hst.morphology_null_receipt.v0.2",
        "status": "MORPHOLOGY_PRESERVING_NULL_PILOT",
        "backend": args.backend,
        "source": manifest["source_archive"],
        "selection": manifest["selection"],
        "null_models": {
            "legacy_pixel_permutation": {
                "preserves": ["marginal_intensity"],
                "destroys": ["local_correlation", "power_spectrum", "large_scale_morphology"],
            },
            "morphology_preserving_phase": {
                "preserves": ["marginal_intensity", "Fourier_power_envelope", "low_frequency_phase"],
                "controls": ["local_correlation", "large_scale_smooth_morphology"],
                "low_freq_fraction": LOW_FREQ_FRACTION,
                "phase_strength": PHASE_STRENGTH,
            },
            "local_block_shuffle": {
                "block_size": BLOCK,
                "preserves": ["within_block_local_correlation", "marginal_intensity"],
                "destroys": ["global_layout"],
            },
        },
        "targets": targets,
        "robust_candidate_count": sum(1 for t in targets if t["robust_cross_band_candidate"]),
        "claim_ceiling": "Image-level geometric robustness test only. No astronomical discovery claim without independent replication and scientific review.",
        "blind_gate": {
            "seed": base.SEED,
            "nulls": NULLS,
            "semantic_analysis": False,
            "ocr": False,
            "face_search": False,
            "cipher_search": False,
            "post_hoc_tuning": False,
            "cross_filter_required": True,
        },
    }
    base.RECEIPT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    receipt_sha = hashlib.sha256(base.RECEIPT.read_bytes()).hexdigest()
    base.emit("run_completed", robust_candidate_count=receipt["robust_candidate_count"], receipt_sha256=receipt_sha)
    receipt["event_log"] = {"path": str(base.EVENT_LOG), "sha256": hashlib.sha256(base.EVENT_LOG.read_bytes()).hexdigest()}
    receipt["receipt_sha256"] = receipt_sha
    base.RECEIPT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
