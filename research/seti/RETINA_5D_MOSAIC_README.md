# JANUS COSMOS — RETINA × 5D Spiral × Mosaic Kaleidoscope (research only)

**Status: SYNTHETIC ENGINEERING PASS / SCIENCE HOLD. No production adoption.**

This is a **new, derived experimental mosaic operator** inspired by the canonical [JANUS 5D Spiral](https://github.com/Hawkar-usls/janus-meta-registry/blob/main/data/JANUS-DEEP-ANALYSIS-5D-SPIRAL-ORIGIN-PRIME-ORACLE-PROTOCOL-2026-08-21-v2.1.json) and [R41S preregistration](https://github.com/Hawkar-usls/Janus-Cosmos/blob/research/seti-r41s-5d-spiral-20261120/research/seti/r41s_5d_spiral_protocol.json). It does **not** claim the mosaic transform already existed in Meta Registry.

Native radio measurements and video Retina **are not the same input modality**.

- Radio: read-only numerical (time, frequency) arrays; fixed native-frequency paths; four temporal mosaic tiles; forward/reverse and frequency-mirror **parity checks only**; fixed path ON/OFF comparison; source digests; keep ALL baseline alerts.
- TOPA Retina: visual QC on source-bound decoded frames from the [shared module](https://github.com/Hawkar-usls/TOPA/blob/main/tools/topa_retina_video.py). Require the original TOPA frame manifest and validate frame SHA-256 + PTS before analyzing tiles. Sampled visual frames must never be treated as native radio power or independent astrophysical evidence.
- Neither mirrored copies nor nested mosaic pieces become independent observations. Shared ON/OFF structures are `advisory_only`, and an unresolved reference can **never** veto an ON candidate.
- `origin_prime_evidence_sha256` is an evidence-content state hash; repeated state returns `NO_INFORMATION_GAIN__STOP`. No time-travel or physical 5D claim.

## Run

```bash
python -m pip install numpy pillow
python -m unittest discover -s research/seti -p 'test_retina_5d_mosaic*.py' -v
python research/seti/benchmark_retina_5d_mosaic.py
python research/seti/retina_5d_mosaic_cli.py --native-on /path/native_on.npy \
  --native-off /path/native_off.npy --source-id SRC --reference-id REF \
  --output /tmp/retina_native_advisory.json
python research/seti/retina_5d_mosaic_cli.py \
  --retina-manifest /path/to/TOPA-unit/frame_manifest.json \
  --output /tmp/retina_visual_qc.json
```

For native radio arrays, matching dimensions alone do **not** authenticate instrument frequency semantics or true telescope pointing. `--native-grid-verified` is a **caller attestation**, not independent metadata proof.

## Frozen synthetic diagnostic (2026-11-23)

Calibration: 80 synthetic Gaussian nulls, independent seed from 200 evaluation cases (40 each). The two methods use the same three drift paths and rowwise robust normalization; calibration p95 is **synthetic-only**.

| Input group (n=40) | Direct exceedances | Mosaic exceedances | Shared pattern advisory |
|---|---:|---:|---:|
| Gaussian noise | 1 | 3 | 0 |
| Two impulsive outliers | 27 | 3 | 0 |
| Weak Source-only line | 35 | **31** | 0 |
| Strong Source-only line | 40 | 40 | 0 |
| Shared ON/OFF line | 40 | 40 | 39 |

**Critical negative result:** the mosaic loses **4/40 weak synthetic detections** versus the direct matched-mean comparator, and noise-only exceedances rise from 1 to 3 in this small test. The impulse reduction is exploratory, **not** a general or astrophysical false-alert reduction. The shared-line label does not remove candidates.

The real July/September Parkes spectra, the sealed October 2025 holdout, production models, thresholds, calibration files, and original TOPA source videos were **not accessed or modified** in this synthetic run.

Before any adoption: use authenticated native radio spectra and independent acquisition groups, ensure correct IF/pointing semantics, compare against frozen R41S/R42 baselines and full matched filter at equal measured false-alert burden, preserve every protected weak line, and report runtime and all counterexamples. External blind validation is still required.

See `retina_5d_mosaic_protocol.json` and `retina_5d_mosaic_synthetic_result.json` for immutable source pointers, method settings and the negative result.
