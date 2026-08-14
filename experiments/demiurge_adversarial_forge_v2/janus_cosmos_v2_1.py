#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import math
import multiprocessing
import os
import sys
import threading
import time
import traceback
from pathlib import Path

import numpy as np

import janus_cosmos_core_v2 as core
import janus_cosmos_specificity_v2_1 as specificity
import janus_cosmos_v2_0 as parent_runtime


VERSION = "2.1.1"
ROOT = Path(__file__).resolve().parent
PROTOCOL = json.loads((ROOT / "SPECIFICITY_PROTOCOL_v2_1.json").read_text(encoding="utf-8"))
SCIENTIFIC_VERSION = str(PROTOCOL["version"])
MAX_WORKERS = 10
DATA = ROOT / "external_data"
OUT = ROOT / "results_v2_1"
CHECKPOINTS = OUT / "checkpoints"
HISTORY = OUT / "history"
EVENTS = OUT / "janus-cosmos-v2.1-events.jsonl"
REPORT = OUT / "janus-cosmos-v2.1-report.json"
SUMMARY = OUT / "SUMMARY_v2.1.txt"
TERMINAL = OUT / "terminal_v2.1.log"
PROGRESS_STATE = OUT / "progress_v2.1.json"
LOG_HANDLE = None
LOG_LOCK = threading.RLock()
PROGRESS = None
WORKER_MESSAGE_QUEUE = None


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{threading.get_ident()}.tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def atomic_write_json(path: Path, value: dict) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def rotate_previous_output(path: Path, run_stamp: str) -> None:
    if not path.exists() or path.stat().st_size == 0:
        return
    HISTORY.mkdir(parents=True, exist_ok=True)
    destination = HISTORY / f"{path.stem}__{run_stamp}{path.suffix}"
    counter = 1
    while destination.exists():
        destination = HISTORY / f"{path.stem}__{run_stamp}__{counter}{path.suffix}"
        counter += 1
    os.replace(path, destination)


def format_duration(seconds: float | None) -> str:
    if seconds is None or not math.isfinite(seconds):
        return "--:--"
    seconds = max(0, int(seconds))
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}" if hours else f"{minutes:02d}:{seconds:02d}"


class LiveProgress:
    SPINNER = "|/-\\"

    def __init__(self, total_units: int, workers: int, state_path: Path) -> None:
        self.total_units = max(1, int(total_units))
        self.workers = int(workers)
        self.state_path = state_path
        self.started = time.monotonic()
        self.fractions: dict[str, float] = {}
        self.active: set[str] = set()
        self.update_count = 0
        self.last_state_write = 0.0
        self.last_non_tty_percent = -5
        self.rendered_width = 0
        self.closed = False

    def start_task(self, task_id: str) -> None:
        with LOG_LOCK:
            self.fractions.setdefault(task_id, 0.0)
            self.active.add(task_id)
            self._render_locked()

    def update(self, task_id: str, fraction: float) -> None:
        with LOG_LOCK:
            self.fractions[task_id] = max(self.fractions.get(task_id, 0.0), min(1.0, float(fraction)))
            self.active.add(task_id)
            self.update_count += 1
            self._render_locked()

    def complete(self, task_id: str) -> None:
        with LOG_LOCK:
            self.fractions[task_id] = 1.0
            self.active.discard(task_id)
            self.update_count += 1
            self._render_locked(force_state=True)

    def snapshot(self) -> dict:
        completed_units = min(float(self.total_units), sum(self.fractions.values()))
        elapsed = max(time.monotonic() - self.started, 1e-9)
        rate = completed_units / elapsed
        eta = (self.total_units - completed_units) / rate if rate > 0 else None
        return {
            "schema": "janus.cosmos.live_progress.v1",
            "runtime_version": VERSION,
            "scientific_protocol_version": SCIENTIFIC_VERSION,
            "workers": self.workers,
            "total_model_units": self.total_units,
            "completed_model_units": completed_units,
            "completed_percent": 100.0 * completed_units / self.total_units,
            "active_tasks": sorted(self.active),
            "elapsed_seconds": elapsed,
            "eta_seconds": eta,
        }

    def _line(self, state: dict) -> str:
        fraction = state["completed_model_units"] / state["total_model_units"]
        width = 36
        filled = min(width, int(round(width * fraction)))
        bar = "█" * filled + "░" * (width - filled)
        spinner = self.SPINNER[self.update_count % len(self.SPINNER)]
        done = int(state["completed_model_units"])
        return (
            f"{spinner} [{bar}] {state['completed_percent']:6.2f}%  "
            f"models {done}/{state['total_model_units']}  active {len(self.active)}/{self.workers}  "
            f"elapsed {format_duration(state['elapsed_seconds'])}  ETA {format_duration(state['eta_seconds'])}"
        )

    def clear_locked(self) -> None:
        if sys.stdout.isatty() and self.rendered_width:
            sys.stdout.write("\r" + " " * self.rendered_width + "\r")
            sys.stdout.flush()

    def redraw_locked(self) -> None:
        if not self.closed:
            self._render_locked(write_state=False)

    def _render_locked(self, force_state: bool = False, write_state: bool = True) -> None:
        if self.closed:
            return
        state = self.snapshot()
        line = self._line(state)
        if sys.stdout.isatty():
            self.clear_locked()
            sys.stdout.write("\r" + line)
            sys.stdout.flush()
            self.rendered_width = len(line)
        else:
            milestone = int(state["completed_percent"] // 5) * 5
            if milestone >= self.last_non_tty_percent + 5:
                print("[PROGRESS] " + line, flush=True)
                if LOG_HANDLE is not None:
                    LOG_HANDLE.write("[PROGRESS] " + line + "\n")
                    LOG_HANDLE.flush()
                self.last_non_tty_percent = milestone
        now = time.monotonic()
        if write_state and (force_state or now - self.last_state_write >= 0.5):
            atomic_write_json(self.state_path, state)
            self.last_state_write = now

    def close(self) -> None:
        with LOG_LOCK:
            if self.closed:
                return
            self.clear_locked()
            state = self.snapshot()
            atomic_write_json(self.state_path, state)
            if sys.stdout.isatty():
                print(self._line(state), flush=True)
            self.closed = True


class QueueProgress:
    """Small child-process proxy; rendering and file writes stay in the parent."""

    def __init__(self, message_queue) -> None:
        self.message_queue = message_queue

    def start_task(self, task_id: str) -> None:
        self.message_queue.put({"kind": "progress", "action": "start", "task_id": task_id})

    def update(self, task_id: str, fraction: float) -> None:
        self.message_queue.put(
            {"kind": "progress", "action": "update", "task_id": task_id, "fraction": float(fraction)}
        )

    def complete(self, task_id: str) -> None:
        self.message_queue.put({"kind": "progress", "action": "complete", "task_id": task_id})


def process_worker_init(message_queue) -> None:
    global LOG_HANDLE, PROGRESS, WORKER_MESSAGE_QUEUE
    LOG_HANDLE = None
    WORKER_MESSAGE_QUEUE = message_queue
    PROGRESS = QueueProgress(message_queue)


def log(message: str) -> None:
    with LOG_LOCK:
        if PROGRESS is not None:
            PROGRESS.clear_locked()
        print(message, flush=True)
        if LOG_HANDLE is not None:
            LOG_HANDLE.write(message + "\n")
            LOG_HANDLE.flush()
        if PROGRESS is not None:
            PROGRESS.redraw_locked()


def write_event_row(row: dict) -> None:
    encoded = json.dumps(row, sort_keys=True, ensure_ascii=False)
    with LOG_LOCK:
        log(encoded)
        with EVENTS.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write(encoded + "\n")


def emit(event: str, **fields) -> None:
    row = {"ts_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "event": event, **fields}
    if WORKER_MESSAGE_QUEUE is not None:
        WORKER_MESSAGE_QUEUE.put({"kind": "event", "row": row})
        return
    write_event_row(row)


def start_message_pump(message_queue) -> threading.Thread:
    def pump() -> None:
        while True:
            message = message_queue.get()
            if message is None:
                return
            if message["kind"] == "event":
                write_event_row(message["row"])
                continue
            if message["kind"] != "progress" or PROGRESS is None:
                continue
            action = message["action"]
            if action == "start":
                PROGRESS.start_task(message["task_id"])
            elif action == "update":
                PROGRESS.update(message["task_id"], message["fraction"])
            elif action == "complete":
                PROGRESS.complete(message["task_id"])

    thread = threading.Thread(target=pump, name="janus-v21-progress-pump", daemon=True)
    thread.start()
    return thread


def combined_sha256(paths: list[Path]) -> str:
    rows = [{"path": str(path.relative_to(ROOT)), "sha256": core.sha256_file(path)} for path in paths]
    return hashlib.sha256(specificity.canonical_json(rows).encode("utf-8")).hexdigest()


def checkpoint_key(file_sha: str, genome_sha: str, label: str, model: str, nulls: int, cal: int, seeds: list[int]) -> tuple[str, dict]:
    settings = {
        # Keep the scientific checkpoint identity at v2.1.0 so an interrupted
        # v2.1.0 run can resume under the v2.1.1 execution-only acceleration.
        "version": SCIENTIFIC_VERSION,
        "file_sha": file_sha,
        "genome_sha": genome_sha,
        "protocol_sha": specificity.canonical_sha256(PROTOCOL),
        "label": label,
        "variant": "WHOLE",
        "model": model,
        "nulls": nulls,
        "cal": cal,
        "seeds": seeds,
    }
    return hashlib.sha256(specificity.canonical_json(settings).encode("utf-8")).hexdigest(), settings


def run_tail_model(
    x: np.ndarray,
    genome: dict,
    file_sha: str,
    label: str,
    model: str,
    nulls: int,
    cal: int,
    seeds: list[int],
) -> dict:
    key, settings = checkpoint_key(file_sha, core.sha256_bytes(core.canonical_json(genome).encode("utf-8")), label, model, nulls, cal, seeds)
    checkpoint = CHECKPOINTS / f"{label}__WHOLE__{model}.json"
    task_id = f"{label}:WHOLE:{model}"
    if checkpoint.exists():
        try:
            old = json.loads(checkpoint.read_text(encoding="utf-8"))
            if old.get("settings_hash") == key:
                if PROGRESS is not None:
                    PROGRESS.start_task(task_id)
                    PROGRESS.complete(task_id)
                emit("model_cache_hit", label=label, variant="WHOLE", model=model)
                return old["result"]
        except (OSError, ValueError, TypeError, KeyError) as error:
            emit("checkpoint_invalid_recompute", checkpoint=checkpoint.name, error=f"{type(error).__name__}: {error}")

    if PROGRESS is not None:
        PROGRESS.start_task(task_id)
    seed_counts: dict[int, int] = {}
    base, extra = divmod(nulls, len(seeds))
    for index, seed in enumerate(seeds):
        seed_counts[int(seed)] = base + (1 if index < extra else 0)
    seed_completed = {int(seed): 0 for seed in seeds}

    def progress(seed: int, index: int, count: int) -> None:
        seed_completed[int(seed)] = int(index)
        completed = sum(min(seed_completed[s], seed_counts[s]) for s in seed_completed)
        if PROGRESS is not None:
            PROGRESS.update(task_id, completed / max(nulls, 1))

    emit("model_start", label=label, variant="WHOLE", model=model)
    result = specificity.empirical_test_with_tail(
        x,
        genome,
        model,
        nulls,
        cal,
        seeds,
        (label, "WHOLE"),
        progress=progress,
    )
    atomic_write_json(checkpoint, {"settings_hash": key, "settings": settings, "result": result})
    if PROGRESS is not None:
        PROGRESS.complete(task_id)
    emit(
        "model_complete",
        label=label,
        variant="WHOLE",
        model=model,
        p=result["p_empirical"],
        tail_ratio_q99=result["tail_ratio_q99"],
    )
    return result


def analyse_normalized(
    x: np.ndarray,
    genome: dict,
    file_sha: str,
    label: str,
    nulls: int,
    cal: int,
    seeds: list[int],
) -> dict:
    models = {
        model: run_tail_model(x, genome, file_sha, label, model, nulls, cal, seeds)
        for model in ("phase_iaaft", "block_shuffle")
    }
    return {"models": models, "band_tail_effect": specificity.band_tail_effect(models)}


def analyse_fits(path: Path, label: str, genome: dict, nulls: int, cal: int, seeds: list[int]) -> tuple[dict, np.ndarray, dict, tuple]:
    if not path.is_file():
        raise RuntimeError(f"missing frozen source product: {path.relative_to(ROOT)}; run download_sky_v2_1.py")
    raw, header, meta = core.read_primary_fits(path)
    x = core.normalize(raw, genome)
    file_sha = core.sha256_file(path)
    analysis = {
        "file": str(path.relative_to(ROOT)),
        "sha256": file_sha,
        "meta": meta,
        **analyse_normalized(x, genome, file_sha, label, nulls, cal, seeds),
    }
    return analysis, x, header, raw.shape


def control_file(field: dict, survey: dict) -> Path:
    filename = f"{field['id'].lower()}_{survey['family'].lower()}_{survey['band'].lower()}.fits".replace("2mass", "tmass")
    return DATA / "controls_v2_1" / filename


def real_sky_source_paths(field: dict, is_target: bool) -> list[Path]:
    surveys = PROTOCOL["orion_target"]["surveys"] if is_target else PROTOCOL["real_sky_controls"]["surveys"]
    if is_target:
        return [DATA / "orion" / survey["filename"] for survey in surveys]
    return [control_file(field, survey) for survey in surveys]


def field_checkpoint_key(
    analysis_kind: str,
    field_id: str,
    paths: list[Path],
    genome: dict,
    nulls: int,
    cal: int,
    seeds: list[int],
) -> tuple[str, dict]:
    missing = [str(path.relative_to(ROOT)) for path in paths if not path.is_file()]
    if missing:
        raise RuntimeError("missing frozen source product(s): " + ", ".join(missing) + "; run download_sky_v2_1.py")
    settings = {
        "scientific_version": SCIENTIFIC_VERSION,
        "protocol_sha256": specificity.canonical_sha256(PROTOCOL),
        "analysis_kind": analysis_kind,
        "field_id": field_id,
        "source_sha256": [
            {"path": str(path.relative_to(ROOT)), "sha256": core.sha256_file(path)} for path in paths
        ],
        "genome_sha256": core.sha256_bytes(core.canonical_json(genome).encode("utf-8")),
        "nulls": int(nulls),
        "cal": int(cal),
        "seeds": [int(seed) for seed in seeds],
    }
    return hashlib.sha256(specificity.canonical_json(settings).encode("utf-8")).hexdigest(), settings


def cached_field_analysis(
    analysis_kind: str,
    field_id: str,
    paths: list[Path],
    genome: dict,
    nulls: int,
    cal: int,
    seeds: list[int],
    compute,
) -> dict:
    key, settings = field_checkpoint_key(analysis_kind, field_id, paths, genome, nulls, cal, seeds)
    checkpoint = CHECKPOINTS / f"FIELD__{analysis_kind}__{field_id}.json"
    if checkpoint.exists():
        try:
            old = json.loads(checkpoint.read_text(encoding="utf-8"))
            if old.get("settings_hash") == key:
                if PROGRESS is not None:
                    unit_count = 8 if analysis_kind in {"REAL_SKY_CONTROL", "ORION_TARGET"} else 4
                    for index in range(unit_count):
                        task_id = f"FIELD_CACHE:{analysis_kind}:{field_id}:{index}"
                        PROGRESS.start_task(task_id)
                        PROGRESS.complete(task_id)
                emit("field_cache_hit", analysis_kind=analysis_kind, field_id=field_id)
                return old["result"]
        except (OSError, ValueError, TypeError, KeyError) as error:
            emit("checkpoint_invalid_recompute", checkpoint=checkpoint.name, error=f"{type(error).__name__}: {error}")
    result = compute()
    atomic_write_json(checkpoint, {"settings_hash": key, "settings": settings, "result": result})
    return result


def corridor_field_score(bands: dict) -> float:
    # The local p-value remains a gate, but the cross-field rank uses a
    # continuous q99 tail ratio so p-floor saturation cannot create ties by
    # construction.
    values = [float(item["corridor"]["local_rank"]["tail_ratio_q99"]) for item in bands.values()]
    return float(min(values))


def local_corridor_family_gate(bands: dict) -> dict:
    required = int(PROTOCOL["corridor_null"]["required_passing_bands_per_family"])
    result = {}
    for family in ("DSS2", "2MASS"):
        passing = sum(
            bool(item["corridor"]["local_rank"]["passes_local_alpha"])
            for item in bands.values()
            if item["family"] == family
        )
        result[family] = {"passing_bands": passing, "required": required, "pass": passing >= required}
    result["pass"] = all(result[family]["pass"] for family in ("DSS2", "2MASS"))
    return result


def analyse_real_sky_field(
    field: dict,
    genome: dict,
    nulls: int,
    cal: int,
    seeds: list[int],
    is_target: bool,
) -> dict:
    surveys = PROTOCOL["orion_target"]["surveys"] if is_target else PROTOCOL["real_sky_controls"]["surveys"]
    bands = {}
    images: dict[tuple[str, str], np.ndarray] = {}
    corridors: dict[tuple[str, str], np.ndarray] = {}
    candidate_cfg = PROTOCOL["corridor_null"]
    heldout_spec = None
    if not is_target:
        heldout_spec = specificity.corridor_spec(
            (candidate_cfg["candidate_seed_domain"], field["id"]),
            float(candidate_cfg["length_pixels_normalized"]),
            float(candidate_cfg["half_width_pixels_normalized"]),
        )
    for survey in surveys:
        family, band = survey["family"], survey["band"]
        if is_target:
            path = DATA / "orion" / survey["filename"]
            label = f"ORION_{family}_{band}"
        else:
            path = control_file(field, survey)
            label = f"{field['id']}_{family}_{band}"
        analysis, x, header, native_shape = analyse_fits(path, label, genome, nulls, cal, seeds)
        if is_target:
            candidate, diagnostics = core.belt_corridor(
                x,
                header,
                native_shape,
                PROTOCOL["orion_target"]["belt_stars_j2000"],
                half_width=float(candidate_cfg["half_width_pixels_normalized"]),
                margin=8,
            )
        else:
            candidate = specificity.extract_corridor(x, heldout_spec)
            diagnostics = {"heldout_spec": heldout_spec, "designation": "FROZEN_RANDOM_CONTROL_CANDIDATE"}
        local_rank = specificity.corridor_local_rank(x, candidate, genome, analysis["sha256"], label, candidate_cfg)
        analysis["corridor"] = {"diagnostics": diagnostics, "local_rank": local_rank}
        bands[label] = {"family": family, "band": band, "analysis": analysis, "corridor": analysis["corridor"]}
        images[(family, band)] = x
        corridors[(family, band)] = candidate
    whole_analyses = [item["analysis"] for item in bands.values()]
    morphology_cfg = PROTOCOL["morphology_agreement"]
    return {
        "center": field,
        "bands": bands,
        "scores": {
            "whole_detector": specificity.field_tail_effect(whole_analyses),
            "corridor_detector": corridor_field_score(bands),
            "whole_cross_survey_morphology": specificity.orion_cross_survey_agreement(images, morphology_cfg),
            "corridor_cross_survey_morphology": specificity.orion_cross_survey_agreement(corridors, morphology_cfg),
        },
        "corridor_local_family_gate": local_corridor_family_gate(bands),
    }


def analyse_real_sky_field_cached(
    field: dict,
    genome: dict,
    nulls: int,
    cal: int,
    seeds: list[int],
    is_target: bool,
) -> dict:
    field_id = "ORION" if is_target else str(field["id"])
    kind = "ORION_TARGET" if is_target else "REAL_SKY_CONTROL"
    paths = real_sky_source_paths(field, is_target)
    return cached_field_analysis(
        kind,
        field_id,
        paths,
        genome,
        nulls,
        cal,
        seeds,
        lambda: analyse_real_sky_field(field, genome, nulls, cal, seeds, is_target),
    )


def hst_paths(field_id: str, control: bool) -> dict[str, Path]:
    root = DATA / "hst" / ("controls_v2_1" if control else "")
    if control:
        root = root / field_id
    else:
        root = DATA / "hst" / field_id
    chip = PROTOCOL["hst_real_controls"]["canonical_chip"]
    paths = {}
    for band in PROTOCOL["hst_target"]["bands"]:
        paths[f"{band}_science"] = root / f"h_{field_id}_{band}_{chip}.fits"
        paths[f"{band}_weight"] = root / f"h_{field_id}_{band}_{chip}_wgt.fits"
    return paths


def analyse_hst_field(
    field_id: str,
    genome: dict,
    nulls: int,
    cal: int,
    seeds: list[int],
    control: bool,
) -> dict:
    paths = hst_paths(field_id, control)
    for path in paths.values():
        if not path.is_file():
            raise RuntimeError(f"missing frozen HST source product: {path.relative_to(ROOT)}; run download_sky_v2_1.py")
    f555_raw, _, f555_meta = core.read_primary_fits(paths["f555_science"])
    f555_weight, _, f555_weight_meta = core.read_primary_fits(paths["f555_weight"])
    f814_raw, _, f814_meta = core.read_primary_fits(paths["f814_science"])
    f814_weight, _, f814_weight_meta = core.read_primary_fits(paths["f814_weight"])
    x555, x814, mask_diagnostics = specificity.common_valid_support_pair(
        f555_raw,
        f814_raw,
        f555_weight,
        f814_weight,
        genome,
        PROTOCOL["hst_target"]["common_valid_support"],
    )
    all_paths = list(paths.values())
    field_sha = combined_sha256(all_paths)
    label_prefix = f"HSTCTRL_{field_id}" if control else "NGC1425"
    analyses = {
        "F555": analyse_normalized(x555, genome, field_sha + ":f555", label_prefix + "_F555_WF3_MASKED", nulls, cal, seeds),
        "F814": analyse_normalized(x814, genome, field_sha + ":f814", label_prefix + "_F814_WF3_MASKED", nulls, cal, seeds),
    }
    analyses["F555"]["source"] = {
        "science": str(paths["f555_science"].relative_to(ROOT)),
        "weight": str(paths["f555_weight"].relative_to(ROOT)),
        "science_meta": f555_meta,
        "weight_meta": f555_weight_meta,
    }
    analyses["F814"]["source"] = {
        "science": str(paths["f814_science"].relative_to(ROOT)),
        "weight": str(paths["f814_weight"].relative_to(ROOT)),
        "science_meta": f814_meta,
        "weight_meta": f814_weight_meta,
    }
    morphology = specificity.morphology_correlation(x555, x814, PROTOCOL["morphology_agreement"])
    return {
        "field_id": field_id,
        "combined_source_sha256": field_sha,
        "mask": mask_diagnostics,
        "bands": analyses,
        "scores": {
            "whole_detector": specificity.field_tail_effect(analyses.values()),
            "cross_filter_morphology": float(morphology),
        },
    }


def analyse_hst_field_cached(
    field_id: str,
    genome: dict,
    nulls: int,
    cal: int,
    seeds: list[int],
    control: bool,
) -> dict:
    kind = "HST_REAL_CONTROL" if control else "HST_TARGET"
    paths = list(hst_paths(field_id, control).values())
    return cached_field_analysis(
        kind,
        field_id,
        paths,
        genome,
        nulls,
        cal,
        seeds,
        lambda: analyse_hst_field(field_id, genome, nulls, cal, seeds, control),
    )


def process_field_task(
    analysis_kind: str,
    item,
    genome: dict,
    nulls: int,
    cal: int,
    seeds: list[int],
) -> dict:
    if analysis_kind == "REAL_SKY_CONTROL":
        return analyse_real_sky_field_cached(item, genome, nulls, cal, seeds, False)
    if analysis_kind == "HST_REAL_CONTROL":
        return analyse_hst_field_cached(str(item), genome, nulls, cal, seeds, True)
    if analysis_kind == "ORION_TARGET":
        return analyse_real_sky_field_cached(item, genome, nulls, cal, seeds, True)
    if analysis_kind == "HST_TARGET":
        return analyse_hst_field_cached(str(item), genome, nulls, cal, seeds, False)
    raise RuntimeError(f"unknown process analysis kind: {analysis_kind}")


def persist_running_report(report: dict) -> None:
    report["global_status"]["last_checkpoint_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    atomic_write_json(REPORT, report)


def parallel_ordered_fields(
    items: list,
    key_of,
    analyse,
    workers: int,
    start_event: str,
    complete_event: str,
    report: dict,
    report_bucket: str,
    executor=None,
    process_kind: str | None = None,
    process_args: tuple = (),
) -> dict:
    frozen_order = [str(key_of(item)) for item in items]
    rows_by_key: dict[str, dict] = {}
    worker_count = max(1, min(int(workers), len(items) or 1))
    emit(
        "parallel_batch_start",
        batch=report_bucket,
        workers=worker_count,
        item_count=len(items),
        result_order="FROZEN_INPUT_ORDER",
    )
    owned_executor = None
    if executor is None:
        owned_executor = concurrent.futures.ThreadPoolExecutor(
            max_workers=worker_count, thread_name_prefix="janus-v21-selftest"
        )
        active_executor = owned_executor
    else:
        if not process_kind:
            raise RuntimeError("process_kind is required with an external process executor")
        active_executor = executor
    try:
        futures = {}
        for item in items:
            field_id = str(key_of(item))
            emit(start_event, field_id=field_id)
            if executor is None:
                future = active_executor.submit(analyse, item)
            else:
                future = active_executor.submit(process_field_task, process_kind, item, *process_args)
            futures[future] = field_id
        for future in concurrent.futures.as_completed(futures):
            field_id = futures[future]
            try:
                rows_by_key[field_id] = future.result()
            except Exception as error:
                for pending in futures:
                    pending.cancel()
                emit(
                    "parallel_field_failed",
                    batch=report_bucket,
                    field_id=field_id,
                    error=f"{type(error).__name__}: {error}",
                )
                raise
            report[report_bucket] = {key: rows_by_key[key] for key in frozen_order if key in rows_by_key}
            report["global_status"][f"{report_bucket}_completed"] = len(rows_by_key)
            persist_running_report(report)
            emit(
                complete_event,
                field_id=field_id,
                completed=len(rows_by_key),
                total=len(items),
            )
    finally:
        if owned_executor is not None:
            owned_executor.shutdown(wait=True, cancel_futures=True)
    ordered = {key: rows_by_key[key] for key in frozen_order}
    emit("parallel_batch_complete", batch=report_bucket, completed=len(ordered), workers=worker_count)
    return ordered


def rank_gate(target_score: float, controls: list[float]) -> dict:
    rank = specificity.real_field_rank(target_score, controls)
    maximum = int(PROTOCOL["real_field_admission"]["maximum_control_exceedances"])
    rank["maximum_control_exceedances"] = maximum
    rank["pass"] = bool(rank["control_exceedances"] <= maximum)
    return rank


def build_summary(report: dict) -> str:
    lines = [
        f"JANUS COSMOS v{VERSION} — DETECTOR SPECIFICITY REPAIR",
        "",
        f"status: {report['status']}",
        f"smoke_only: {report['smoke_only']}",
        f"protocol_sha256: {report['protocol']['protocol_sha256']}",
        f"frozen_genome: {report['frozen_detector']['genome_sha256']}",
        f"workers: {report.get('execution', {}).get('workers')}",
    ]
    if report.get("global_status"):
        lines.extend(
            [
                "",
                f"REAL_SKY_CONTROL_FIELDS: {report['global_status'].get('real_sky_control_fields')}",
                f"HST_REAL_CONTROL_FIELDS: {report['global_status'].get('hst_real_control_fields')}",
                f"ORION: {report['orion'].get('status') if report.get('orion') else 'NOT_RUN'}",
                f"NGC1425: {report['ngc1425'].get('status') if report.get('ngc1425') else 'NOT_RUN'}",
            ]
        )
    if report["errors"]:
        lines.extend(["", "ERRORS:", *[f"- {item['error']}" for item in report["errors"]]])
    lines.extend(["", "Claim ceiling: " + PROTOCOL["claim_ceiling"]])
    return "\n".join(lines) + "\n"


def resolve_workers(requested: int | None) -> int:
    raw = requested if requested is not None else os.environ.get("JANUS_COSMOS_WORKERS", str(MAX_WORKERS))
    try:
        workers = int(raw)
    except (TypeError, ValueError) as error:
        raise RuntimeError(f"invalid worker count: {raw!r}") from error
    return max(1, min(MAX_WORKERS, workers))


def main() -> int:
    global LOG_HANDLE, PROGRESS
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--nulls", type=int)
    parser.add_argument("--cal-nulls", type=int)
    parser.add_argument("--workers", type=int, help=f"parallel field workers, 1..{MAX_WORKERS}; default {MAX_WORKERS}")
    args = parser.parse_args()
    run_started_monotonic = time.monotonic()
    run_started_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    frozen = parent_runtime.verify_forge()
    protocol_receipt = specificity.verify_protocol_sources(PROTOCOL)
    if frozen["genome_sha256"] != PROTOCOL["parent_detector"]["genome_sha256"] or frozen["freeze_sha256"] != PROTOCOL["parent_detector"]["freeze_sha256"]:
        raise RuntimeError("v2.1 protocol is not bound to the packaged frozen detector")
    if args.self_test:
        print(f"SELF-TEST PASS: v2.0.2 frozen detector + v2.1 protocol verified; runtime={VERSION}", flush=True)
        return 0
    cfg = PROTOCOL["synthetic_null_diagnostics"]
    if not args.smoke and (args.nulls is not None or args.cal_nulls is not None):
        raise RuntimeError("full v2.1 run forbids Monte Carlo overrides; use the frozen protocol")
    nulls = int(args.nulls or (24 if args.smoke else cfg["test_nulls_per_model"]))
    cal = int(args.cal_nulls or (12 if args.smoke else cfg["calibration_nulls_per_model"]))
    seeds = list(cfg["seeds"])
    workers = resolve_workers(args.workers)
    sky_fields = list(PROTOCOL["real_sky_controls"]["centers"][: (2 if args.smoke else None)])
    hst_fields = list(PROTOCOL["hst_real_controls"]["field_ids"][: (2 if args.smoke else None)])
    total_model_jobs = len(sky_fields) * 8 + len(hst_fields) * 4 + (0 if args.smoke else 12)

    OUT.mkdir(exist_ok=True)
    CHECKPOINTS.mkdir(exist_ok=True)
    run_stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    for path in (EVENTS, REPORT, SUMMARY, TERMINAL, PROGRESS_STATE):
        rotate_previous_output(path, run_stamp)
    atomic_write_text(EVENTS, "")
    LOG_HANDLE = TERMINAL.open("w", encoding="utf-8", buffering=1, newline="\n")
    PROGRESS = LiveProgress(total_model_jobs, workers, PROGRESS_STATE)
    report = {
        "schema": "janus.cosmos.v2.1.report",
        "version": VERSION,
        "scientific_protocol_version": SCIENTIFIC_VERSION,
        "status": "RUNNING",
        "smoke_only": bool(args.smoke),
        "execution": {
            "scheduler": "WINDOWS_SAFE_SPAWN_PROCESS_POOL_FROZEN_RESULT_ORDER",
            "workers": workers,
            "maximum_workers": MAX_WORKERS,
            "total_model_jobs": total_model_jobs,
            "started_utc": run_started_utc,
            "scientific_statistics_changed_by_scheduler": False,
            "legacy_v2_1_0_model_checkpoints_accepted_when_settings_hash_matches": True,
        },
        "protocol": protocol_receipt,
        "frozen_detector": {"genome_sha256": frozen["genome_sha256"], "freeze_sha256": frozen["freeze_sha256"]},
        "negative_parent_certificate": PROTOCOL["negative_parent_certificate"],
        "synthetic_null_diagnostics": {"test_nulls_per_model": nulls, "calibration_nulls": cal, "seeds": seeds},
        "real_sky_controls": {},
        "hst_real_controls": {},
        "orion": {},
        "ngc1425": {},
        "global_status": {},
        "errors": [],
        "claim_ceiling": PROTOCOL["claim_ceiling"],
    }
    genome = frozen["genome"]
    persist_running_report(report)
    emit(
        "run_start",
        runtime_version=VERSION,
        scientific_protocol_version=SCIENTIFIC_VERSION,
        workers=workers,
        total_model_jobs=total_model_jobs,
    )
    process_context = multiprocessing.get_context("spawn")
    message_queue = process_context.Queue()
    message_pump = start_message_pump(message_queue)
    process_pool = concurrent.futures.ProcessPoolExecutor(
        max_workers=workers,
        mp_context=process_context,
        initializer=process_worker_init,
        initargs=(message_queue,),
    )
    try:
        report["real_sky_controls"] = parallel_ordered_fields(
            sky_fields,
            lambda field: field["id"],
            None,
            workers,
            "real_sky_control_start",
            "real_sky_control_complete",
            report,
            "real_sky_controls",
            executor=process_pool,
            process_kind="REAL_SKY_CONTROL",
            process_args=(genome, nulls, cal, seeds),
        )
        report["hst_real_controls"] = parallel_ordered_fields(
            hst_fields,
            str,
            None,
            workers,
            "hst_real_control_start",
            "hst_real_control_complete",
            report,
            "hst_real_controls",
            executor=process_pool,
            process_kind="HST_REAL_CONTROL",
            process_args=(genome, nulls, cal, seeds),
        )
        report["global_status"]["real_sky_control_fields"] = len(report["real_sky_controls"])
        report["global_status"]["hst_real_control_fields"] = len(report["hst_real_controls"])
        persist_running_report(report)
        if args.smoke:
            report["status"] = "PASS"
            report["global_status"]["admission_disabled"] = True
            report["global_status"]["admission_disabled_reason"] = "SMOKE_ONLY_INCOMPLETE_REAL_CONTROL_COHORT"
        else:
            emit("orion_target_start")
            emit("ngc1425_target_start")
            target_futures = {
                process_pool.submit(
                    process_field_task,
                    "ORION_TARGET",
                    PROTOCOL["orion_target"]["center_j2000"],
                    genome,
                    nulls,
                    cal,
                    seeds,
                ): "ORION",
                process_pool.submit(
                    process_field_task,
                    "HST_TARGET",
                    PROTOCOL["hst_target"]["id"],
                    genome,
                    nulls,
                    cal,
                    seeds,
                ): "NGC1425",
            }
            target_results = {}
            for future in concurrent.futures.as_completed(target_futures):
                target_id = target_futures[future]
                target_results[target_id] = future.result()
                emit("orion_target_complete" if target_id == "ORION" else "ngc1425_target_complete")
                if target_id == "ORION":
                    report["orion"] = target_results[target_id]
                else:
                    report["ngc1425"] = target_results[target_id]
                persist_running_report(report)
            orion = target_results["ORION"]
            ngc = target_results["NGC1425"]
            sky_controls = list(report["real_sky_controls"].values())
            orion["real_field_gates"] = {
                "whole_detector": rank_gate(orion["scores"]["whole_detector"], [item["scores"]["whole_detector"] for item in sky_controls]),
                "corridor_detector": rank_gate(orion["scores"]["corridor_detector"], [item["scores"]["corridor_detector"] for item in sky_controls]),
                "whole_cross_survey_morphology": rank_gate(
                    orion["scores"]["whole_cross_survey_morphology"]["score"],
                    [item["scores"]["whole_cross_survey_morphology"]["score"] for item in sky_controls],
                ),
                "corridor_cross_survey_morphology": rank_gate(
                    orion["scores"]["corridor_cross_survey_morphology"]["score"],
                    [item["scores"]["corridor_cross_survey_morphology"]["score"] for item in sky_controls],
                ),
            }
            orion_pass = orion["corridor_local_family_gate"]["pass"] and all(
                item["pass"] for item in orion["real_field_gates"].values()
            )
            orion["admitted"] = bool(orion_pass)
            orion["status"] = "SKY_FIXED_MORPHOLOGY_CANDIDATE" if orion_pass else "DETECTOR_SPECIFICITY_BLOCKED"
            report["orion"] = orion

            hst_controls = list(report["hst_real_controls"].values())
            ngc["real_field_gates"] = {
                "whole_detector": rank_gate(ngc["scores"]["whole_detector"], [item["scores"]["whole_detector"] for item in hst_controls]),
                "cross_filter_morphology": rank_gate(
                    ngc["scores"]["cross_filter_morphology"],
                    [item["scores"]["cross_filter_morphology"] for item in hst_controls],
                ),
            }
            ngc_pass = ngc["mask"]["mask_gate_pass"] and all(item["pass"] for item in ngc["real_field_gates"].values())
            ngc["admitted"] = bool(ngc_pass)
            ngc["status"] = "HST_CROSS_FILTER_MORPHOLOGY_CANDIDATE" if ngc_pass else "DETECTOR_SPECIFICITY_BLOCKED"
            report["ngc1425"] = ngc
            report["global_status"].update(
                {
                    "admission_disabled": False,
                    "orion_candidate_admitted": bool(orion_pass),
                    "ngc1425_candidate_admitted": bool(ngc_pass),
                    "claim_ceiling_enforced": True,
                }
            )
            report["status"] = "PASS"
    except Exception as error:
        report["status"] = "FAIL"
        report["errors"].append({"error": f"{type(error).__name__}: {error}", "traceback": traceback.format_exc(limit=16)})
        emit("fatal_error", error=str(error))
    finally:
        process_pool.shutdown(wait=True, cancel_futures=True)
        message_queue.put(None)
        message_pump.join(timeout=30)
        message_queue.close()
        message_queue.join_thread()
    if PROGRESS is not None:
        PROGRESS.close()
    report["execution"]["elapsed_seconds"] = float(time.monotonic() - run_started_monotonic)
    report["execution"]["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    emit("run_complete", status=report["status"])
    atomic_write_json(REPORT, report)
    summary = build_summary(report)
    atomic_write_text(SUMMARY, summary)
    log(summary.rstrip())
    if LOG_HANDLE:
        LOG_HANDLE.close()
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    multiprocessing.freeze_support()
    raise SystemExit(main())
