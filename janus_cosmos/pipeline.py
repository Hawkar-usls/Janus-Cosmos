from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path
from typing import Iterable

import numpy as np

from . import __version__
from .core import GateConfig, NULL_MODELS, analyze_image, bonferroni_alpha, minimum_test_nulls
from .exploratory import run_exploratory


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


class EventWriter:
    def __init__(self, path: Path):
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text("", encoding="utf-8")

    def emit(self, event: str, **fields):
        record = {"schema": "janus.cosmos.event.v1", "ts_unix": time.time(), "event": event, **fields}
        line = json.dumps(record, sort_keys=True, ensure_ascii=False)
        print(line, flush=True)
        with self.path.open("a", encoding="utf-8") as f:
            f.write(line + "\n")


def _source_for_item(item: dict) -> str:
    for key in ("dataURI", "url", "source_url"):
        value = str(item.get(key, "") or "")
        if value:
            return value
    raise ValueError(f"filter item lacks a source URI/URL: {item}")


def download_source(source: str, cache_dir: Path, events: EventWriter, *, target: str, filter_name: str) -> tuple[Path, dict]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    suffix = ".fits.gz" if source.lower().endswith(".fits.gz") else ".fits"
    path = cache_dir / (hashlib.sha256(source.encode("utf-8")).hexdigest() + suffix)
    if path.exists() and path.stat().st_size >= 2880:
        meta = {"source": source, "cache_hit": True, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
        events.emit("download_cache_hit", target=target, filter=filter_name, **meta)
        return path, meta

    last_error = None
    for attempt in range(1, 4):
        try:
            events.emit("download_started", target=target, filter=filter_name, source=source, attempt=attempt)
            if source.startswith("mast:"):
                from astroquery.mast import Observations
                status, msg, url = Observations.download_file(source, local_path=str(path), cache=False)
                if status != "COMPLETE":
                    raise RuntimeError(f"MAST download status={status!r}, message={msg!r}, url={url!r}")
            else:
                req = urllib.request.Request(
                    source,
                    headers={
                        "User-Agent": "Janus-Cosmos/1.0 (+https://github.com/Hawkar-usls/Janus-Cosmos)",
                        "Accept": "application/fits,application/octet-stream,*/*",
                    },
                )
                with urllib.request.urlopen(req, timeout=240) as response, path.open("wb") as f:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        f.write(chunk)
            if not path.exists() or path.stat().st_size < 2880:
                raise RuntimeError("downloaded payload is too small to be FITS")
            meta = {
                "source": source,
                "cache_hit": False,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "attempt": attempt,
            }
            events.emit("download_ok", target=target, filter=filter_name, **meta)
            return path, meta
        except Exception as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            events.emit("download_retry", target=target, filter=filter_name, source=source, attempt=attempt, error=last_error)
            path.unlink(missing_ok=True)
            if attempt < 3:
                time.sleep(2 * attempt)
    raise RuntimeError(f"download failed for {source}: {last_error}")


def read_fits_image(path: Path) -> tuple[np.ndarray, dict]:
    from astropy.io import fits

    with fits.open(path, memmap=False, lazy_load_hdus=True) as hdul:
        candidates = []
        for idx, hdu in enumerate(hdul):
            data = getattr(hdu, "data", None)
            if data is None or np.ndim(data) != 2:
                continue
            arr = np.asarray(data, dtype=np.float32)
            if arr.size < 1024:
                continue
            extname = str(getattr(hdu, "header", {}).get("EXTNAME", ""))
            candidates.append((arr.size, 1 if extname.upper() == "SCI" else 0, idx, extname, arr))
        if not candidates:
            raise RuntimeError(f"No usable 2-D FITS image plane in {path}")
        _, _, idx, extname, image = max(candidates, key=lambda x: (x[1], x[0]))
    meta = {
        "selected_hdu": int(idx),
        "selected_extname": extname,
        "native_shape": [int(x) for x in image.shape],
        "nan_fraction": float(np.mean(~np.isfinite(image))),
    }
    return image, meta


def load_manifest(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    targets = data.get("targets")
    if not isinstance(targets, list) or not targets:
        raise RuntimeError("manifest has no targets")
    for target in targets:
        if not target.get("target") or not isinstance(target.get("filters"), list):
            raise RuntimeError(f"invalid target entry: {target}")
    return data


def _select_targets(targets: list[dict], only: set[str] | None) -> list[dict]:
    if not only:
        return targets
    selected = [t for t in targets if t["target"] in only]
    missing = sorted(only - {t["target"] for t in selected})
    if missing:
        raise RuntimeError(f"requested targets missing from manifest: {missing}")
    return selected


def run_pipeline(
    *,
    manifest_path: Path,
    output_dir: Path,
    cache_dir: Path,
    test_nulls: int,
    seeds: list[int],
    only_targets: set[str] | None = None,
    exploratory: bool = False,
    allow_underpowered: bool = False,
    config: GateConfig | None = None,
) -> dict:
    cfg = config or GateConfig()
    output_dir.mkdir(parents=True, exist_ok=True)
    events = EventWriter(output_dir / "janus-cosmos-events.jsonl")
    manifest = load_manifest(manifest_path)
    targets = _select_targets(manifest["targets"], only_targets)
    filter_count = sum(len(t["filters"]) for t in targets)
    alpha = bonferroni_alpha(cfg.alpha_family, filter_count, len(NULL_MODELS))
    required = minimum_test_nulls(alpha)
    powered = test_nulls >= required
    if not powered and not allow_underpowered:
        raise RuntimeError(
            f"test_nulls={test_nulls} cannot resolve corrected alpha={alpha:.8g}; "
            f"minimum is {required}. Use --allow-underpowered only for smoke tests."
        )

    events.emit(
        "run_started",
        version=__version__,
        manifest=str(manifest_path),
        target_count=len(targets),
        filter_count=filter_count,
        test_nulls=test_nulls,
        seeds=seeds,
        alpha_family=cfg.alpha_family,
        alpha_corrected=alpha,
        minimum_test_nulls=required,
        powered=powered,
        exploratory=exploratory,
    )

    targets_out = []
    errors = []
    for ti, target in enumerate(targets, start=1):
        name = target["target"]
        print(f"\n[{ti}/{len(targets)}] {name}", flush=True)
        target_out = {
            "target": name,
            "target_class": target.get("class", "unknown"),
            "focus": bool(target.get("focus", name in {"NGC1425", "NGC1637"})),
            "filters": {},
            "source_products": [],
        }
        for fi, item in enumerate(target["filters"], start=1):
            filt = str(item.get("filter", f"FILTER_{fi}"))
            try:
                source = _source_for_item(item)
                print(f"  [{fi}/{len(target['filters'])}] {filt}: download/read/analyze", flush=True)
                path, source_meta = download_source(source, cache_dir, events, target=name, filter_name=filt)
                image, image_meta = read_fits_image(path)
                events.emit("fits_parsed", target=name, filter=filt, **image_meta)
                result = analyze_image(
                    image,
                    target=name,
                    filter_name=filt,
                    test_nulls=test_nulls,
                    seeds=seeds,
                    alpha=alpha,
                    include_legacy=True,
                    config=cfg,
                )
                result["powered"] = powered
                if exploratory:
                    result["exploratory"] = run_exploratory(image)
                target_out["filters"][filt] = result
                target_out["source_products"].append({
                    "filter": filt,
                    "band": item.get("band"),
                    "productFilename": item.get("productFilename"),
                    **source_meta,
                    **image_meta,
                })
                events.emit(
                    "filter_completed",
                    target=name,
                    filter=filt,
                    robust_candidate=result["robust_candidate"],
                    phase_p=result["phase_iaaft"]["p_empirical"],
                    block_p=result["block_shuffle"]["p_empirical"],
                )
            except Exception as exc:
                err = {"target": name, "filter": filt, "error": f"{type(exc).__name__}: {exc}"}
                errors.append(err)
                target_out["filters"][filt] = {"status": "ERROR", **err}
                events.emit("filter_error", **err)

        passing = [
            filt for filt, result in target_out["filters"].items()
            if result.get("robust_candidate") is True and result.get("powered") is True
        ]
        target_out["robust_passing_filters"] = passing
        target_out["robust_cross_filter_candidate"] = bool(powered and len(passing) >= 2)
        target_out["cross_filter_semantics"] = (
            "At least two independently tested filters pass the same geometric gate; "
            "this is not a WCS-localized proof that the same pixels/structure persist across bands."
        )
        targets_out.append(target_out)
        events.emit(
            "target_completed",
            target=name,
            passing_filters=passing,
            robust_cross_filter_candidate=target_out["robust_cross_filter_candidate"],
        )

    status = "PASS"
    if errors:
        status = "PARTIAL_ANALYSIS" if any(t["source_products"] for t in targets_out) else "DATA_FAILURE"
    if not powered:
        status = "SMOKE_ONLY" if not errors else "SMOKE_WITH_ERRORS"

    receipt = {
        "schema": "janus.cosmos.hst.canonical_receipt.v1",
        "status": status,
        "version": __version__,
        "source": manifest.get("source", manifest.get("source_archive", "unknown")),
        "selection": manifest.get("selection", "unknown"),
        "manifest_path": str(manifest_path),
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "target_count": len(targets),
        "filter_count": filter_count,
        "test_nulls_per_null_model": int(test_nulls),
        "calibration_nulls_per_null_model": int(cfg.calibration_nulls),
        "seeds": [int(x) for x in seeds],
        "alpha_family": float(cfg.alpha_family),
        "alpha_corrected": float(alpha),
        "minimum_test_nulls": int(required),
        "powered_for_corrected_gate": bool(powered),
        "null_models": {
            "phase_iaaft": {
                "role": "primary morphology-preserving null",
                "low_freq_fraction": cfg.low_freq_fraction,
                "iaaft_iterations": cfg.iaaft_iterations,
            },
            "block_shuffle": {
                "role": "primary local-correlation-preserving null",
                "block_size": cfg.block_size,
            },
            "pixel_permutation": {"role": "legacy diagnostic only; never sufficient for candidate status"},
        },
        "blind_gate": {
            "ocr": False,
            "face_search": False,
            "semantic_analysis": False,
            "cipher_search": False,
            "post_hoc_tuning": False,
            "human_label_inference": False,
            "cross_filter_required": True,
        },
        "exploratory_enrichment": {
            "requested": bool(exploratory),
            "affects_blind_gate": False,
            "post_hoc_tuning": "FORBIDDEN_ON_EVALUATION_DATA",
        },
        "targets": targets_out,
        "candidate_count": sum(1 for t in targets_out if t["robust_cross_filter_candidate"]),
        "candidate_targets": [t["target"] for t in targets_out if t["robust_cross_filter_candidate"]],
        "errors": errors,
        "claim_ceiling": (
            "Image-level geometric candidate only. No astronomical discovery, hidden-message, intelligence, "
            "or unknown-physics claim without independent data replication and scientific review."
        ),
    }
    events.emit("run_completed", status=status, candidate_count=receipt["candidate_count"], errors=len(errors))
    receipt["event_log"] = {"path": str(events.path), "sha256": sha256_file(events.path)}
    out_path = output_dir / "janus-cosmos-receipt.json"
    out_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    receipt["receipt_sha256_pre_field"] = sha256_file(out_path)
    out_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def parse_seeds(text: str) -> list[int]:
    seeds = [int(x.strip()) for x in text.split(",") if x.strip()]
    if not seeds:
        raise argparse.ArgumentTypeError("at least one seed is required")
    return seeds


def main(argv: Iterable[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Janus Cosmos canonical v1 blind geometric pipeline")
    ap.add_argument("--manifest", default="data/hst_blind_corpus.json")
    ap.add_argument("--output-dir", default="results/canonical_v1")
    ap.add_argument("--cache-dir", default=".cache/janus_cosmos")
    ap.add_argument("--nulls", type=int, default=1024, help="test nulls per primary null model")
    ap.add_argument("--seeds", default="20260810,20260811,20260812")
    ap.add_argument("--targets", default="", help="comma-separated target names from manifest")
    ap.add_argument("--exploratory", action="store_true")
    ap.add_argument("--allow-underpowered", action="store_true")
    args = ap.parse_args(list(argv) if argv is not None else None)

    seeds = parse_seeds(args.seeds)
    only = {x.strip() for x in args.targets.split(",") if x.strip()} or None
    try:
        receipt = run_pipeline(
            manifest_path=Path(args.manifest),
            output_dir=Path(args.output_dir),
            cache_dir=Path(args.cache_dir),
            test_nulls=args.nulls,
            seeds=seeds,
            only_targets=only,
            exploratory=args.exploratory,
            allow_underpowered=args.allow_underpowered,
        )
    except Exception as exc:
        print(f"JANUS COSMOS FATAL: {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        return 2

    print(json.dumps({
        "status": receipt["status"],
        "candidate_count": receipt["candidate_count"],
        "candidate_targets": receipt["candidate_targets"],
        "powered_for_corrected_gate": receipt["powered_for_corrected_gate"],
        "alpha_corrected": receipt["alpha_corrected"],
        "minimum_test_nulls": receipt["minimum_test_nulls"],
    }, indent=2))
    return 0 if receipt["status"] in {"PASS", "SMOKE_ONLY", "PARTIAL_ANALYSIS", "SMOKE_WITH_ERRORS"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
