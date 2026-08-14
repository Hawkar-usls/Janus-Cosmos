<div align="center">

# JANUS COSMOS
### Proof-carrying blind geometry search in public astronomical imagery

![Status](https://img.shields.io/badge/status-canonical%20v1-1f6feb)
![Gate](https://img.shields.io/badge/gate-morphology--preserving-8957e5)
![Semantics](https://img.shields.io/badge/exploratory-post--gate%20only-6e7681)

</div>

## Mission

Janus Cosmos is an experimental, reproducible pipeline for finding **image-level geometric outliers** in public astronomical data. HST/MAST is the first corpus. A passing candidate is a reason to replicate and inspect, not an astronomical discovery claim.

The canonical v1 question is:

> Does an observation contain geometric structure that remains unusual under morphology-preserving null models, survives a family-wise statistical gate, and persists in at least two independently tested filters?

## Why v1 exists

The earlier pilots taught us where the first implementation was weak. A pixel-permutation null produced `7/7` candidates because it destroyed normal galaxy morphology. A later morphology-preserving pilot reduced that historical set to `2/7` (`NGC1425`, `NGC1637`), but it used 256 nulls. For the fixed seven-galaxy / 14-filter family with two primary null models, canonical v1 uses `alpha = 0.05 / 28 ≈ 0.0017857`; 256 nulls have a best possible empirical p-value of `1/257 ≈ 0.003891`, so those two objects are **focus targets, not confirmed anomalies**.

Canonical v1 fixes the null design, empirical-p resolution guard, multiple-testing denominator, seed semantics, source caching, FITS handling, live MAST discovery, runner/result schema mismatch, Windows launcher, and machine-readable logging.

## Canonical gate

```text
source manifest
   ↓
provenance-preserving download + SHA-256 cache
   ↓
FITS SCI-plane selection
   ↓
14-feature non-semantic geometry vector
   ↓
phase/IAAFT-like null  AND  local block-shuffle null
   ↓
empirical p-values
   ↓
Bonferroni family-wise gate
   ↓
≥2 filters pass both primary nulls
   ↓
CROSS-FILTER IMAGE CANDIDATE
```

Pixel permutation remains a **legacy diagnostic only** and can never admit a canonical candidate.

## Geometry vector

The fixed feature set contains multiscale directional correlation (σ=1/2/4), 90°/180° rotational agreement, gradient anisotropy, high-frequency energy, Fourier angular anisotropy, connected-component counts at q80/q90/q95, and largest-component fractions at q80/q90/q95.

## Statistical power gate

```text
alpha_corrected = 0.05 / (total_filter_tests × 2 primary null models)
```

The runner checks whether the requested empirical null count can actually resolve that threshold. Underpowered runs fail closed unless `--allow-underpowered` is explicitly supplied; then the receipt is `SMOKE_ONLY` and cannot contain a powered candidate.

Multiple seeds are deterministic Monte-Carlo chunks pooled into one null distribution. They are **not** treated as independent astrophysical replications.

## Fixed source-confirmed baseline

`data/hst_blind_corpus.json` is the first end-to-end baseline and contains public STScI Spiral Galaxies HLSP FITS URLs for F555W/F814W:

- NGC1365
- NGC1425 — historical focus target
- NGC1637 — historical focus target
- NGC2841
- NGC3031
- NGC3627
- NGC4321

Because these URLs are fixed, the baseline does not depend on live archive discovery.

## Live MAST expansion

`janus_cosmos.discovery` uses the high-level astroquery workflow:

```text
Observations.query_criteria(...)
  → Observations.get_unique_product_list(...)
  → deterministic public science FITS selection
  → mast: dataURI
```

Targets enter the expanded scoring manifest only when at least two requested filters are available.

## Exploratory layer

User-requested exploratory searches are isolated **after** the blind decision and have `affects_blind_gate=false`:

- OCR: optional local Tesseract/pytesseract;
- face-like regions: optional **non-identifying** OpenCV detector;
- cipher-like scan: descriptive entropy/repetition over OCR output only;
- semantic analysis: unavailable unless a separate local model is explicitly configured;
- post-hoc tuning on evaluation data: forbidden.

This prevents interesting-looking post-hoc results from changing the primary p-values.

## Run on Python

```bash
python -m pip install -r requirements.txt
python -m janus_cosmos.pipeline \
  --manifest data/hst_blind_corpus.json \
  --nulls 1024 \
  --seeds 20260810,20260811,20260812
```

Smoke test only:

```bash
python -m janus_cosmos.pipeline --nulls 64 --seeds 20260810 --allow-underpowered
```

## Run on Windows

From the repository folder, double-click:

```text
run_janus_cosmos_windows.bat
```

The BAT always changes to its own directory (`cd /d "%~dp0"`), streams output through PowerShell `Tee-Object`, writes a terminal log, and keeps the window open.

## Demiurge adversarial forge v2.0.2

The independent [v2.0.2 forge line](experiments/demiurge_adversarial_forge_v2/) adds a synthetic adversarial detector forge, a hard freeze-before-target wall, unrelated real-sky specificity controls, Orion cross-survey validation, and the historical NGC1425 HST gate.

Its portability contract deliberately separates:

- `freeze_sha256`: the portable detector identity (genome, validation decision, normalized source/manifests, blind-wall declaration and numerical contract);
- `metrics_sha256`: the exact raw synthetic metrics emitted by the current Python/NumPy/SciPy platform;
- `FORGE_METRIC_REFERENCE_v2_0.json`: a source-bound `2e-6` cross-platform conformance envelope.

The scientific gate always uses raw metrics. Metric drift outside the envelope, source drift, genome drift, manifest drift or validation failure blocks the real-sky run. CI covers Linux plus the reported Windows stack: Python 3.14, NumPy 2.5.2 and SciPy 1.18.0.

Windows entry point:

```text
experiments\demiurge_adversarial_forge_v2\run_janus_cosmos_v2_0.bat
```

## Tests and CI

```bash
python -m pip install -r requirements.txt pytest
python -m pytest -q
python -m py_compile janus_cosmos/*.py
```

`.github/workflows/canonical-v1.yml` runs unit/synthetic tests and an end-to-end fixed-HLSP focus smoke test. Manual dispatch also supports the powered fixed-seven baseline and live expanded-manifest discovery.

## Receipt semantics

Canonical runs write `janus-cosmos-receipt.json` and `janus-cosmos-events.jsonl`, including source SHA-256 hashes, FITS-plane metadata, corrected alpha, minimum null resolution, primary null results, errors, and claim ceiling.

`robust_cross_filter_candidate=true` means at least two filters pass the whole-image gate. It does **not** yet prove WCS-localized persistence of the same exact structure at the same sky coordinates.

## Claim ceiling

**CANONICAL v1 — RESEARCH PIPELINE. NO ASTRONOMICAL DISCOVERY CLAIMED.**

The pipeline cannot by itself establish unknown physics, artificial structures, hidden communication, extraterrestrial intelligence, or semantic meaning. The next scientific gate after a powered corpus result is WCS-aware localization and independent-data replication.
