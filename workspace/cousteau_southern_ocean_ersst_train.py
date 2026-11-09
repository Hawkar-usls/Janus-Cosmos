#!/usr/bin/env python3
"""JANUS Cousteau Southern Ocean SST: source-pinned retrospective one-step benchmark."""
import argparse
import concurrent.futures
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
import numpy as np
import requests
from netCDF4 import Dataset

BASE = "https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/netcdf"
FIRST = (2013, 1)
LAST = (2025, 12)
HOLDOUT = (2024, 1)
OUT = Path("data/cousteau/southern_ocean")
RAW = Path(os.environ.get("COUSTEAU_RAW_DIR", "/tmp/janus-cousteau-southern-raw"))
RUN_ID = os.environ.get("GITHUB_RUN_ID", "local")
USER_AGENT = "JANUS-Cousteau-SouthernOcean-Research/1.0"

def months(start, end):
    y, m = start
    while (y, m) <= end:
        yield (y, m)
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)

def file_for(ym):
    y, m = ym
    return "ersst.v5.%04d%02d.nc" % (y, m)

def southern_mean(blob):
    if not (blob.startswith(b"CDF") or blob.startswith(b"\\x89HDF")):
        # HDF5 magic is checked below with a true byte literal.
        if not blob.startswith(bytes.fromhex("89484446")):
            raise ValueError("Not a NetCDF3 or HDF5 file")
    with Dataset("inmemory.nc", memory=blob) as ds:
        lat = np.asarray(ds.variables["lat"][:], dtype=float)
        sea = np.ma.filled(np.ma.array(ds.variables["sst"][:]).squeeze(), np.nan)
        sea = np.asarray(sea, dtype=float)
        if sea.ndim != 2 or sea.shape[0] != len(lat):
            raise ValueError("Unexpected NOAA ERSST grid shape")
        where = lat <= -50.0
        if not np.any(where):
            raise ValueError("Southern Ocean latitude absent")
        selected = sea[where, :]
        weights = np.cos(np.deg2rad(lat[where]))[:, None]
        valid = np.isfinite(selected) & (selected >= -3.0) & (selected <= 45.0)
        if np.count_nonzero(valid) < 100:
            raise ValueError("Insufficient valid ocean SST cells")
        den = float(np.sum(weights * valid))
        return (float(np.sum(np.where(valid, selected, 0.0) * weights) / den),
                int(np.count_nonzero(valid)), den)

def acquire(ym):
    name = file_for(ym)
    path = RAW / name
    path.parent.mkdir(parents=True, exist_ok=True)
    url = BASE + "/" + name
    if path.exists():
        blob = path.read_bytes()
    else:
        error = None
        for attempt in range(3):
            try:
                response = requests.get(url, headers={"User-Agent": USER_AGENT},
                                        timeout=(20, 90))
                response.raise_for_status()
                blob = response.content
                if len(blob) < 50000 or len(blob) > 4000000:
                    raise ValueError("Unexpected file size: %s" % len(blob))
                path.write_bytes(blob)
                break
            except (requests.RequestException, ValueError) as exc:
                error = exc
                if attempt < 2:
                    time.sleep(2 * (attempt + 1))
        else:
            raise RuntimeError("%s: %s" % (name, error))
    value, cells, area_weight = southern_mean(blob)
    return {"year_month": "%04d-%02d" % ym,
            "source_url": url, "source_sha256": hashlib.sha256(blob).hexdigest(),
            "source_bytes": len(blob), "south_50_sst_c": round(value, 7),
            "valid_ocean_grid_cells": cells,
            "grid_cell_cosine_weight_sum": round(area_weight, 4)}

def features(series, i):
    m = i % 12
    return ([float(m == j) for j in range(12)]
            + [i / 12.0, series[i-1], series[i-12]])

def train_and_evaluate(rows):
    if len(rows) != len(list(months(FIRST, LAST))):
        raise ValueError("Incomplete time series")
    x = np.array([r["south_50_sst_c"] for r in rows], dtype=float)
    dates = [r["year_month"] for r in rows]
    split = dates.index("%04d-%02d" % HOLDOUT)
    idx_train = list(range(12, split))
    idx_test = list(range(split, len(x)))
    X = np.array([features(x, i) for i in idx_train])
    Y = x[idx_train]
    mean = X.mean(axis=0)
    scale = X.std(axis=0)
    scale[scale < 1e-10] = 1.0
    Z = (X-mean)/scale
    # Ridge regularization fixed before inspection of the holdout.
    alpha = 10.0
    beta = np.linalg.solve(Z.T @ Z + alpha * np.eye(Z.shape[1]), Z.T @ (Y-Y.mean()))
    bias = float(Y.mean())
    T = np.array([features(x, i) for i in idx_test])
    prediction = bias + ((T-mean)/scale) @ beta
    naive = np.array([x[i-12] for i in idx_test])
    truth = x[idx_test]
    def metrics(pred):
        error = pred-truth
        return {"rmse_c": round(float(np.sqrt(np.mean(error*error))), 7),
                "mae_c": round(float(np.mean(np.abs(error))), 7)}
    ridge_metrics, seasonal_metrics = metrics(prediction), metrics(naive)
    return {
        "model": "ridge_month_1hot_trend_lag1_lag12",
        "evaluation": "retrospective rolling one-step-ahead; lag1 may use prior holdout observation",
        "model_parameters": {"alpha": alpha, "mean": mean.tolist(),
                             "scale": scale.tolist(), "coefficients": beta.tolist(),
                             "bias": bias},
        "train_start": dates[12], "train_end": dates[split-1],
        "test_start": dates[split], "test_end": dates[-1],
        "n_train": len(idx_train), "n_test": len(idx_test),
        "metrics": {"ridge": ridge_metrics, "seasonal_naive_12m": seasonal_metrics},
        "holdout_predictions": [
            {"month": dates[i], "observed_sst_c": round(float(x[i]), 7),
             "ridge_sst_c": round(float(prediction[j]), 7),
             "seasonal_naive_sst_c": round(float(naive[j]), 7)}
            for j, i in enumerate(idx_test)],
        "beats_seasonal_baseline": ridge_metrics["rmse_c"] < seasonal_metrics["rmse_c"],
        "promotion_authorized": False,
        "claim_ceiling": "ONE_STEP_REGIONAL_SST_BENCHMARK_ONLY__NO_SONAR_OR_SEAFLOOR_ANOMALY_CLAIM"
    }

def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        sample = [{"year_month": "%04d-%02d" % ym,
                   "south_50_sst_c": round(2.0+math.sin(i*math.pi/6)+i*0.002, 7)}
                  for i, ym in enumerate(months(FIRST, LAST))]
        result = train_and_evaluate(sample)
        assert result["n_test"] == 24 and result["n_train"] == 120
        assert math.isfinite(result["metrics"]["ridge"]["rmse_c"])
        print("SELF_TEST_PASS (SYNTHETIC_ONLY)")
        return
    receipt_path = OUT / "training_runs" / ("run-" + RUN_ID + ".json")
    observations_path = OUT / "observations" / ("NOAA-ERSSTv5-2013-2025-run-" + RUN_ID + ".json")
    start = dt.datetime.now(dt.timezone.utc).isoformat()
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
            rows = list(pool.map(acquire, months(FIRST, LAST)))
        benchmark = train_and_evaluate(rows)
        write_json(observations_path, {
            "schema": "janus.cousteau.southern_ocean.ersst_observations.v1",
            "created_utc": start, "latitude_selection": "grid lat <= -50 degrees",
            "aggregation": "cos(latitude) weighted mean of valid ocean SST cells; NOT a depth/sonar product",
            "dataset": "NOAA ERSST v5 monthly NetCDF", "records": rows
        })
        obs_sha = hashlib.sha256(observations_path.read_bytes()).hexdigest()
        receipt = {"schema": "janus.cousteau.southern_ocean.train_receipt.v1",
                   "run_id": RUN_ID, "started_utc": start,
                   "finished_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                   "status": "TRAINED_AND_EVALUATED", "source": "NOAA ERSSTv5",
                   "source_monthly_files": len(rows), "source_bytes_total": sum(r["source_bytes"] for r in rows),
                   "observations_path": str(observations_path),
                   "observations_sha256": obs_sha,
                   "benchmark": benchmark, "source_access_snapshot": "mutable upstream; per-file SHA256 pinned",
                   "scientific_note": "Retrospective data may have been revised after observation dates; not prospective performance."}
        write_json(receipt_path, receipt)
        print(json.dumps({"status": receipt["status"], "run_id": RUN_ID,
                          "metrics": benchmark["metrics"], "beating_baseline": benchmark["beats_seasonal_baseline"],
                          "n_months": len(rows)}, sort_keys=True))
    except Exception as exc:
        write_json(receipt_path, {
            "schema": "janus.cousteau.southern_ocean.train_receipt.v1",
            "run_id": RUN_ID, "started_utc": start,
            "finished_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "status": "FAILED_NO_TRAINING_CLAIM",
            "error_type": type(exc).__name__, "error": str(exc)[:2000]})
        raise

if __name__ == "__main__":
    main()
