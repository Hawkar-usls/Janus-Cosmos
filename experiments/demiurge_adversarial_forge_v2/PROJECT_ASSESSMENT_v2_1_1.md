# Janus Cosmos v2.1.1 project assessment

## Current evidence

The supplied interrupted v2.1 run is structurally healthy: 12 of 20 real-sky
control fields completed and the thirteenth started. Of 103 completed synthetic
model tests, 77 reached the minimum empirical p-value (74.8%). This reinforces
the v2.0.2 negative result: synthetic-null significance is common in ordinary
sky imagery and must remain diagnostic rather than an admission rule. No target
result exists in the interrupted run, so it supports neither target admission
nor target rejection.

Completed whole-field tail effects in frozen field order were:

| Field | Effect | Field | Effect |
|---|---:|---|---:|
| SKYCTRL_01 | 0.866165 | SKYCTRL_07 | 0.342966 |
| SKYCTRL_02 | 0.983358 | SKYCTRL_08 | 0.501789 |
| SKYCTRL_03 | 1.346650 | SKYCTRL_09 | 1.924269 |
| SKYCTRL_04 | 1.146230 | SKYCTRL_10 | 0.669589 |
| SKYCTRL_05 | 1.362736 | SKYCTRL_11 | 0.635220 |
| SKYCTRL_06 | 0.843926 | SKYCTRL_12 | 0.644483 |

The median is 0.855046 and the observed range is 0.342966–1.924269. Unlike the
saturated synthetic p-values, this continuous real-field effect has usable
dispersion for the frozen cross-field rank.

The subsequently supplied checkpoint archive closes the resume uncertainty:
all 103 JSON files parse, all 103 settings hashes and frozen scientific
identities verify, and the v2.1.1 cache-hit path returns every stored result
exactly. Its key set and p/tail summaries also match all 103 `model_complete`
events. Fields 01–12 contain all eight expected checkpoints; field 13 contains
seven, with only `SKYCTRL_13_2MASS_K / block_shuffle` missing. The run therefore
preserves 103/160 sky-control model units and 103/252 model units in the entire
frozen run.

## Runtime validation

On a ten-field synthetic FITS benchmark using the unchanged numerical core,
24 test nulls and 12 calibration nulls per model, the Windows-safe spawned
process scheduler reduced elapsed time from 70.89 seconds at one worker to
12.25 seconds at ten workers (5.79×). All ten scientific field result objects
were exactly equal and report order remained frozen. A separate end-to-end
smoke comparison also produced exactly equal scientific report sections. These
are validation measurements, not a promise of identical speedup on every CPU.

A completed four-field smoke rerun from whole-field checkpoints finished in
0.60 seconds, demonstrating the new resume path.

## Strengths to preserve

- Freeze-before-target detector identity and fail-closed source/protocol hashes.
- Negative results are retained as first-class, hash-bound certificates.
- Admission uses complete independent real-field cohorts; ties block.
- Within-image random corridor positions/orientations control local selection.
- HST common-valid-pixel support removes the mosaic-footprint shortcut.
- Cross-survey/filter morphology is compared only after PSF/resolution and rank
  photometric alignment.
- Deterministic seeds, exact source provenance, resumable checkpoints and
  cross-platform numerical conformance make the experiment auditable.

## Weaknesses closed by v2.1.1

- Serial field execution and serial downloads; CPU work now uses spawned
  processes while network transfers use threads.
- Non-atomic checkpoints and final-report-only persistence.
- Completion-order nondeterminism under concurrency.
- Nested BLAS oversubscription on a many-thread desktop.
- Crashes on malformed checkpoint JSON.
- Accidental overwriting of prior logs/reports.
- User delivery bloated with historical raw evidence and obsolete launchers.

## Remaining scientific weaknesses

1. The full frozen v2.1 result has not completed, so there is no target verdict.
2. Twenty controls give a minimum real-field empirical p-value of 1/21,
   approximately 0.0476. That is enough for the frozen strict “beat all” rule,
   but gives very coarse tail resolution.
3. Sky controls are deterministic and blind but not yet matched/stratified by
   Galactic latitude, extinction, stellar density, survey plate or coverage.
4. The current morphology score is a global rank correlation after common
   smoothing/resizing. It does not yet quantify localized multiscale feature
   correspondence or registration uncertainty.
5. There is no pre-registered injection/recovery study measuring detector and
   morphology-gate power at controlled signal strengths.
6. The HST control archive is instrument-comparable but heterogeneous in object
   morphology, exposure and background.

## Recommended next frozen gate (only after v2.1 completes)

- Preserve v2.1 as-is; do not tune against its target outcome.
- Pre-register 50–100 real controls per family, or a sequentially valid cohort
  expansion rule, before downloading/evaluating the expanded blind set.
- Stratify controls on observable nuisance covariates and report within-stratum
  ranks as well as the aggregate rank.
- Add WCS-aware localized, multiscale cross-survey matching with uncertainty.
- Freeze negative and positive injection families and publish power curves.
- Replicate admitted image-level candidates in independent surveys/epochs.

## Claim ceiling

This is a promising reproducible anomaly-screening research pipeline, not an
astronomical discovery instrument yet. Its strongest contribution is currently
methodological: it makes detector specificity failures visible and preserves
them instead of converting saturated p-values into claims.
