from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Iterable

from . import __version__
from .core import GateConfig, NULL_MODELS, analyze_image, bonferroni_alpha, minimum_test_nulls
from .luci import LuciProvenanceError, read_luci_fits_image
from .pipeline import EventWriter, _select_targets, _source_for_item, download_source, load_manifest, parse_seeds, sha256_file


def run_luci_pipeline(
    *,
    manifest_path: Path,
    output_dir: Path,
    cache_dir: Path,
    test_nulls: int,
    seeds: list[int],
    only_targets: set[str] | None = None,
    allow_underpowered: bool = False,
    config: GateConfig | None = None,
) -> dict:
    """Run the JANUS geometry gate with a fail-closed LUCI/LUCIFER provenance boundary."""
    cfg = config or GateConfig()
    output_dir.mkdir(parents=True, exist_ok=True)
    events = EventWriter(output_dir / "janus-cosmos-luci-events.jsonl")
    manifest = load_manifest(manifest_path)
    targets = _select_targets(manifest["targets"], only_targets)

    filter_count = sum(len(t["filters"]) for t in targets)
    alpha = bonferroni_alpha(cfg.alpha_family, filter_count, len(NULL_MODELS))
    required = minimum_test_nulls(alpha)
    powered = test_nulls >= required
    if not powered and not allow_underpowered:
        raise RuntimeError(
            f"test_nulls={test_nulls} cannot resolve corrected alpha={alpha:.8g}; minimum is {required}. "
            "Use --allow-underpowered only for smoke tests."
        )

    events.emit(
        "luci_run_started",
        version=__version__,
        manifest=str(manifest_path),
        target_count=len(targets),
        filter_count=filter_count,
        test_nulls=test_nulls,
        seeds=seeds,
        alpha_corrected=alpha,
        minimum_test_nulls=required,
        powered=powered,
        instrument_boundary="LUCI/LUCIFER_ONLY",
        observation_mode="IMAGING_ONLY",
    )

    targets_out = []
    errors = []
    provenance_rejections = []
    for target in targets:
        name = target["target"]
        target_out = {
            "target": name,
            "target_class": target.get("class", "unknown"),
            "filters": {},
            "source_products": [],
        }
        for fi, item in enumerate(target["filters"], start=1):
            filt = str(item.get("filter", f"FILTER_{fi}"))
            try:
                source = _source_for_item(item)
                path, source_meta = download_source(source, cache_dir, events, target=name, filter_name=filt)
                image, image_meta = read_luci_fits_image(
                    path,
                    require_imaging=True,
                    expected_instrument=item.get("instrument"),
                )
                events.emit("luci_provenance_pass", target=name, filter=filt, **image_meta)
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
                result["instrument_provenance"] = image_meta
                target_out["filters"][filt] = result
                target_out["source_products"].append({
                    "filter": filt,
                    "band": item.get("band"),
                    "instrument_expected": item.get("instrument"),
                    **source_meta,
                    **image_meta,
                })
                events.emit(
                    "luci_filter_completed",
                    target=name,
                    filter=filt,
                    instrument=image_meta["instrument"],
                    robust_candidate=result["robust_candidate"],
                    phase_p=result["phase_iaaft"]["p_empirical"],
                    block_p=result["block_shuffle"]["p_empirical"],
                )
            except LuciProvenanceError as exc:
                err = {"target": name, "filter": filt, "error": f"LuciProvenanceError: {exc}"}
                errors.append(err)
                provenance_rejections.append(err)
                target_out["filters"][filt] = {"status": "PROVENANCE_REJECTED", **err}
                events.emit("luci_provenance_rejected", **err)
            except Exception as exc:
                err = {"target": name, "filter": filt, "error": f"{type(exc).__name__}: {exc}"}
                errors.append(err)
                target_out["filters"][filt] = {"status": "ERROR", **err}
                events.emit("luci_filter_error", **err)

        passing = [
            filt for filt, result in target_out["filters"].items()
            if result.get("robust_candidate") is True and result.get("powered") is True
        ]
        target_out["robust_passing_filters"] = passing
        target_out["robust_cross_filter_candidate"] = bool(powered and len(passing) >= 2)
        targets_out.append(target_out)

    status = "PASS"
    if errors:
        status = "PARTIAL_ANALYSIS" if any(t["source_products"] for t in targets_out) else "DATA_FAILURE"
    if not powered:
        status = "SMOKE_ONLY" if not errors else "SMOKE_WITH_ERRORS"

    receipt = {
        "schema": "janus.cosmos.luci_only.receipt.v1",
        "status": status,
        "version": __version__,
        "instrument_scope": ["LUCI1", "LUCI2", "legacy LUCIFER naming"],
        "telescope_scope": "Large Binocular Telescope (ground-based)",
        "spectral_scope": "near infrared; nominal LUCI coverage approximately 0.9-2.5 micrometres",
        "observation_mode_scope": "imaging only",
        "source": manifest.get("source", manifest.get("source_archive", "unknown")),
        "selection": manifest.get("selection", "unknown"),
        "manifest_path": str(manifest_path),
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "target_count": len(targets),
        "filter_count": filter_count,
        "test_nulls_per_null_model": int(test_nulls),
        "seeds": [int(x) for x in seeds],
        "alpha_corrected": float(alpha),
        "minimum_test_nulls": int(required),
        "powered_for_corrected_gate": bool(powered),
        "provenance_gate": {
            "required_instrument_header": True,
            "allowed_instruments": ["LUCI1", "LUCI2", "LUCIFER", "LUCIFER1", "LUCIFER2"],
            "reject_non_lbt_telescope_header": True,
            "reject_spectroscopy": True,
            "fail_closed": True,
        },
        "targets": targets_out,
        "candidate_count": sum(1 for t in targets_out if t["robust_cross_filter_candidate"]),
        "candidate_targets": [t["target"] for t in targets_out if t["robust_cross_filter_candidate"]],
        "provenance_rejections": provenance_rejections,
        "errors": errors,
        "claim_ceiling": (
            "LUCI/LUCIFER near-infrared image-level geometric candidate only. This does not establish an "
            "astronomical anomaly, unknown physics, artificial structure, hidden communication, or extraterrestrial intelligence."
        ),
    }
    events.emit("luci_run_completed", status=status, candidate_count=receipt["candidate_count"], errors=len(errors))
    receipt["event_log"] = {"path": str(events.path), "sha256": sha256_file(events.path)}
    out_path = output_dir / "janus-cosmos-luci-receipt.json"
    out_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    receipt["receipt_sha256_pre_field"] = sha256_file(out_path)
    out_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def main(argv: Iterable[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="JANUS COSMOS LUCI/LUCIFER-only near-IR imaging runner")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--output-dir", default="results/luci_only_v1")
    ap.add_argument("--cache-dir", default=".cache/janus_cosmos_luci")
    ap.add_argument("--nulls", type=int, default=1024)
    ap.add_argument("--seeds", default="20260815,20260816,20260817")
    ap.add_argument("--targets", default="")
    ap.add_argument("--allow-underpowered", action="store_true")
    args = ap.parse_args(list(argv) if argv is not None else None)

    seeds = parse_seeds(args.seeds)
    only = {x.strip() for x in args.targets.split(",") if x.strip()} or None
    try:
        receipt = run_luci_pipeline(
            manifest_path=Path(args.manifest),
            output_dir=Path(args.output_dir),
            cache_dir=Path(args.cache_dir),
            test_nulls=args.nulls,
            seeds=seeds,
            only_targets=only,
            allow_underpowered=args.allow_underpowered,
        )
    except Exception as exc:
        print(f"JANUS COSMOS LUCI FAILED: {type(exc).__name__}: {exc}")
        return 2
    print(json.dumps({"status": receipt["status"], "candidates": receipt["candidate_targets"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
