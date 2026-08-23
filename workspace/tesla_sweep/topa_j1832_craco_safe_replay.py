#!/usr/bin/env python3
"""Safe acquisition gate for ASKAP J1832-0911 CRACO filterbank.

This script deliberately avoids all external pickle execution. It downloads the
public author-exported NumPy .npy product from Zenodo, verifies the published MD5,
checks the NPY header for object dtype, then loads with allow_pickle=False.

It freezes only array provenance and data-shape diagnostics. It does NOT infer
message content, source origin, or even time/frequency axis semantics unless those
are supplied independently by metadata.
"""
from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path

import numpy as np
import requests

URL = "https://zenodo.org/records/15228816/files/SB55237_CRACO_filterbank.npy?download=1"
EXPECTED_MD5 = "41fca6b1dccc464e65e948ed6c82695e"
OUT = Path("data/tesla-sweep/results/TOPA-HUNT-005C2-J1832-CRACO-SAFE-NPY-ACQUISITION-RUN-001.json")


def md5_file(path: Path, chunk: int = 1024 * 1024) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def npy_header(path: Path) -> dict:
    with path.open("rb") as f:
        version = np.lib.format.read_magic(f)
        if version == (1, 0):
            shape, fortran, dtype = np.lib.format.read_array_header_1_0(f)
        elif version == (2, 0):
            shape, fortran, dtype = np.lib.format.read_array_header_2_0(f)
        else:
            # NumPy uses the v2 parser for the current larger-header formats.
            shape, fortran, dtype = np.lib.format.read_array_header_2_0(f)
    return {
        "version": list(version),
        "shape": list(shape),
        "fortran_order": bool(fortran),
        "dtype": str(dtype),
        "dtype_has_object": bool(dtype.hasobject),
    }


def sampled_stats(a: np.ndarray) -> dict:
    # Bound work for large arrays. This is a provenance/health diagnostic only.
    if a.size == 0:
        return {"size": 0}
    flat = a.reshape(-1)
    step = max(1, flat.size // 1_000_000)
    sample = np.asarray(flat[::step][:1_000_000], dtype=np.float64)
    finite = np.isfinite(sample)
    if not finite.any():
        return {
            "size": int(a.size),
            "sample_size": int(sample.size),
            "finite_fraction_sample": 0.0,
        }
    s = sample[finite]
    return {
        "size": int(a.size),
        "sample_size": int(sample.size),
        "finite_fraction_sample": float(finite.mean()),
        "sample_min": float(np.min(s)),
        "sample_max": float(np.max(s)),
        "sample_median": float(np.median(s)),
        "sample_p01": float(np.percentile(s, 1)),
        "sample_p99": float(np.percentile(s, 99)),
    }


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    receipt = {
        "schema": "JANUS_TOPA_J1832_CRACO_SAFE_NPY_ACQUISITION",
        "version": "1.0",
        "hunt_id": "TOPA-TESLA-HUNT-005C2",
        "source": {
            "zenodo_record": "10.5281/zenodo.15228816",
            "filename": "SB55237_CRACO_filterbank.npy",
            "url": URL,
            "published_md5": EXPECTED_MD5,
            "author_readme_says_numpy_load_directly": True,
        },
        "security_rules": {
            "external_pickle_execution": False,
            "numpy_allow_pickle": False,
            "reject_object_dtype": True,
            "classification_ceiling": "SAFE_ARRAY_PROVENANCE_ONLY",
        },
        "status": "STARTED",
    }
    try:
        with tempfile.TemporaryDirectory(prefix="topa_j1832_") as td:
            path = Path(td) / "SB55237_CRACO_filterbank.npy"
            with requests.get(URL, stream=True, timeout=(20, 120)) as r:
                r.raise_for_status()
                with path.open("wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)
            receipt["downloaded_bytes"] = path.stat().st_size
            got = md5_file(path)
            receipt["computed_md5"] = got
            receipt["md5_match"] = got.lower() == EXPECTED_MD5.lower()
            if not receipt["md5_match"]:
                raise RuntimeError("MD5 mismatch; refuse to inspect array")

            header = npy_header(path)
            receipt["npy_header"] = header
            if header["dtype_has_object"]:
                raise RuntimeError("object dtype detected; refuse allow_pickle path")

            arr = np.load(path, mmap_mode="r", allow_pickle=False)
            receipt["array_loaded_safely"] = True
            receipt["sampled_health"] = sampled_stats(arr)
            receipt["axis_semantics"] = "NOT_INFERRED_FROM_SHAPE_ALONE"
            receipt["raw_telescope_voltage_processing_completed"] = False
            receipt["author_exported_filterbank_array_processing_completed"] = True
            receipt["morphology_scoring_completed"] = False
            receipt["status"] = "PASS_SAFE_NONPICKLE_ARRAY_ACQUISITION__READY_FOR_FROZEN_MORPHOLOGY_GATE"
            receipt["next_gate"] = "TOPA-HUNT-005C3-J1832-CRACO-MORPHOLOGY-WITH-METADATA"
    except Exception as exc:
        receipt["status"] = "BLOCKED_OR_FAILED_SAFE_ACQUISITION"
        receipt["error_type"] = type(exc).__name__
        receipt["error"] = str(exc)
        receipt["array_loaded_safely"] = False
        receipt["morphology_scoring_completed"] = False
        receipt["rule"] = "FAIL_CLOSED_AND_KEEP_RECEIPT"

    OUT.write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print(json.dumps(receipt, indent=2))
    return 0 if receipt["status"].startswith("PASS_") else 2


if __name__ == "__main__":
    raise SystemExit(main())
