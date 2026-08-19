from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import importlib
import json
import sys
from pathlib import Path

import numpy as np


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_legacy(legacy_dir: Path):
    legacy_dir = legacy_dir.resolve()
    if not legacy_dir.is_dir():
        raise RuntimeError(f"legacy directory not found: {legacy_dir}")
    sys.path.insert(0, str(legacy_dir))
    dl = importlib.import_module("download_sky_v2_1")
    core = importlib.import_module("janus_cosmos_core_v2")
    specificity = importlib.import_module("janus_cosmos_specificity_v2_1")
    return dl, core, specificity


def _target_item(dl, survey: dict, target: dict, fov: float, pixels: int, destination: Path) -> dict:
    urls = [
        dl.hips_url(
            endpoint,
            survey["hips"],
            float(target["ra_deg_icrs"]),
            float(target["dec_deg_icrs"]),
            fov,
            pixels,
        )
        for endpoint in dl.HIPS_ENDPOINTS
    ]
    return {
        "kind": "TARGET",
        "id": f"{target['analysis_label']}_{survey['family']}_{survey['band']}",
        "dst": destination,
        "urls": urls,
        "family": survey["family"],
        "band": survey["band"],
        "analysis_label": target["analysis_label"],
    }


def _download_many(dl, items: list[dict], workers: int) -> list[dict]:
    rows: list[dict | None] = [None] * len(items)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers, thread_name_prefix="janus-v21-replay") as pool:
        future_map = {
            pool.submit(dl.download, item["urls"], item["dst"]): (i, item)
            for i, item in enumerate(items)
        }
        for future in concurrent.futures.as_completed(future_map):
            i, item = future_map[future]
            info = future.result()
            rows[i] = {
                "id": item["id"],
                "file": str(item["dst"]),
                "sha256": _sha256(item["dst"]),
                "bytes": item["dst"].stat().st_size,
                "source_status": info.get("status"),
                "source_url": info.get("url"),
            }
    return [r for r in rows if r is not None]


def _read_image(core, path: Path) -> np.ndarray:
    raw, _header, _meta = core.read_primary_fits(path)
    arr = np.asarray(raw, dtype=np.float32)
    while arr.ndim > 2:
        arr = arr[0]
    if arr.ndim != 2:
        raise RuntimeError(f"expected 2-D FITS image, got shape {arr.shape} for {path}")
    if not np.any(np.isfinite(arr)):
        raise RuntimeError(f"no finite pixels in {path}")
    return arr


def _score_four_band(specificity, protocol: dict, images: dict[tuple[str, str], np.ndarray]) -> dict:
    required = {("DSS2", "RED"), ("DSS2", "BLUE"), ("2MASS", "J"), ("2MASS", "K")}
    missing = sorted(required.difference(images))
    if missing:
        raise RuntimeError(f"missing frozen survey bands: {missing}")
    return specificity.orion_cross_survey_agreement(images, protocol["morphology_agreement"])


def run(prereg_path: Path, legacy_dir: Path, output_dir: Path, workers: int) -> dict:
    prereg = json.loads(prereg_path.read_text(encoding="utf-8"))
    dl, core, specificity = _load_legacy(legacy_dir)
    protocol = dl.PROTOCOL

    verification = specificity.verify_protocol_sources(protocol)
    if len(protocol["real_sky_controls"]["centers"]) != int(prereg["replay"]["frozen_real_sky_control_count"]):
        raise RuntimeError("legacy v2.1 control-count drift")

    output_dir.mkdir(parents=True, exist_ok=True)
    fov = float(prereg["replay"]["fov_deg"])
    pixels = int(prereg["replay"]["pixels"])
    surveys = list(protocol["real_sky_controls"]["surveys"])

    # Reuse the original frozen v2.1 download plan for its 20 deterministic controls.
    control_items = [item for item in dl.plan() if item["kind"] == "REAL_SKY_CONTROL"]
    expected_controls = len(protocol["real_sky_controls"]["centers"]) * len(surveys)
    if len(control_items) != expected_controls:
        raise RuntimeError(f"legacy control download-plan drift: {len(control_items)} != {expected_controls}")

    target_items: list[dict] = []
    target_paths: dict[str, dict[tuple[str, str], Path]] = {}
    for target in prereg["targets"]:
        label = target["analysis_label"]
        target_paths[label] = {}
        for survey in surveys:
            family = survey["family"]
            band = survey["band"]
            safe_family = family.lower().replace("2mass", "tmass")
            path = output_dir / "target_fits" / label.lower() / f"{safe_family}_{band.lower()}.fits"
            target_items.append(_target_item(dl, survey, target, fov, pixels, path))
            target_paths[label][(family, band)] = path

    download_rows = _download_many(dl, control_items + target_items, workers)

    # Calculate the exact historical v2.1 cross-survey morphology score for each frozen control.
    control_item_by_id = {item["id"]: item for item in control_items}
    control_scores: list[dict] = []
    for field in protocol["real_sky_controls"]["centers"]:
        images: dict[tuple[str, str], np.ndarray] = {}
        files: list[dict] = []
        for survey in surveys:
            item_id = f"{field['id']}_{survey['family']}_{survey['band']}"
            item = control_item_by_id[item_id]
            p = Path(item["dst"])
            images[(survey["family"], survey["band"])] = _read_image(core, p)
            files.append({"id": item_id, "sha256": _sha256(p), "bytes": p.stat().st_size})
        score = _score_four_band(specificity, protocol, images)
        control_scores.append({
            "field_id": field["id"],
            "ra_deg": float(field["ra_deg"]),
            "dec_deg": float(field["dec_deg"]),
            "score": float(score["score"]),
            "cross_family_correlations": [float(x) for x in score["cross_family_correlations"]],
            "files": files,
        })

    controls_vector = [row["score"] for row in control_scores]
    target_results: list[dict] = []
    for target in prereg["targets"]:
        label = target["analysis_label"]
        images = {
            key: _read_image(core, path)
            for key, path in target_paths[label].items()
        }
        score = _score_four_band(specificity, protocol, images)
        rank = specificity.real_field_rank(float(score["score"]), controls_vector)
        pass_gate = bool(rank["outperforms_all_controls"])
        target_results.append({
            "analysis_label": label,
            "semantic_alias_unblinded_after_scoring": target["semantic_alias"],
            "ra_deg_icrs": float(target["ra_deg_icrs"]),
            "dec_deg_icrs": float(target["dec_deg_icrs"]),
            "morphology_score": float(score["score"]),
            "cross_family_correlations": [float(x) for x in score["cross_family_correlations"]],
            "real_field_rank": rank,
            "gate_status": "LEGACY_V21_CROSS_SURVEY_MORPHOLOGY_PASS" if pass_gate else "LEGACY_V21_CROSS_SURVEY_MORPHOLOGY_FAIL",
            "target_files": [
                {
                    "family": family,
                    "band": band,
                    "file": str(path.relative_to(output_dir)),
                    "sha256": _sha256(path),
                    "bytes": path.stat().st_size,
                }
                for (family, band), path in target_paths[label].items()
            ],
        })

    result = {
        "schema": "janus.cosmos.love_edem.legacy_v21_replay.result.v1",
        "experiment_id": prereg["experiment_id"],
        "legacy_authority": prereg["legacy_authority"],
        "legacy_protocol_verification": verification,
        "replay": prereg["replay"],
        "downloaded_product_count": len(download_rows),
        "controls": {
            "count": len(control_scores),
            "score_min": float(min(controls_vector)),
            "score_median": float(np.median(controls_vector)),
            "score_max": float(max(controls_vector)),
            "rows": control_scores,
        },
        "targets": target_results,
        "firewall": {
            **prereg["firewall"],
            "full_legacy_v21_admission_claimed": False,
            "why": "Only the independent cross-survey morphology rank gate is replayed here. The original v2.1 admission stack contains additional whole-detector, corridor-local and family-specific gates."
        },
    }
    out = output_dir / "love-edem-legacy-v21-replay-result.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(json.dumps({
        "experiment_id": result["experiment_id"],
        "legacy_protocol_verification": verification,
        "control_score_range": [result["controls"]["score_min"], result["controls"]["score_max"]],
        "targets": [
            {
                "analysis_label": row["analysis_label"],
                "score": row["morphology_score"],
                "control_exceedances": row["real_field_rank"]["control_exceedances"],
                "p_empirical": row["real_field_rank"]["p_empirical"],
                "gate": row["gate_status"],
            }
            for row in target_results
        ],
    }, indent=2))
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description="Replay the original Janus-Cosmos v2.1 DSS2/2MASS cross-survey morphology gate on frozen LOVE/EDEM directions")
    ap.add_argument("--prereg", default="data/love/LOVE_EDEM_LEGACY_V21_REPLAY_PREREG.json")
    ap.add_argument("--legacy-dir", required=True)
    ap.add_argument("--output-dir", default="results/love_edem_legacy_v21_replay")
    ap.add_argument("--workers", type=int, default=10)
    args = ap.parse_args()
    run(Path(args.prereg), Path(args.legacy_dir), Path(args.output_dir), max(1, min(10, args.workers)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
