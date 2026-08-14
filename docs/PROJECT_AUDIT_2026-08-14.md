# Janus Cosmos — canonical-v1 project audit

Date: 2026-08-14

## Executive finding

The project contained a useful research trajectory, but several experimental generations had accumulated incompatible assumptions. Canonical v1 keeps the useful evidence and removes the failure modes that could turn infrastructure errors or ordinary galaxy morphology into apparent discoveries.

## Historical results and corrected interpretation

### M51 positive control

A real HST pilot showed that the original image-level detector reacts strongly to structured HST morphology. This remains a positive-control result, not a discovery.

### Legacy fixed seven-galaxy gate: 7/7

The first seven-galaxy run compared real F555W/F814W imagery with pixel-permuted images and returned seven cross-filter candidates. The audit concludes that this primarily demonstrates that real galaxies contain spatial organization while pixel permutation destroys it. Pixel permutation is retained only as a diagnostic null.

### Morphology-preserving pilot: 2/7

A later phase-surrogate + block-shuffle pilot left NGC1425 and NGC1637 as historical survivors. It used 256 nulls. Canonical v1 treats 14 target/filter tests and two primary null models as one family, so alpha = 0.05/28 ≈ 0.0017857. With 256 nulls the smallest empirical p is 1/257 ≈ 0.003891; therefore that pilot could not satisfy the canonical family-wise gate. NGC1425 and NGC1637 remain focus targets, not confirmed anomalies.

## Defects found and repaired

1. **Destructive null overinterpretation.** Pixel permutation destroyed normal morphology and made ordinary structure easy to classify as unusual. Canonical admission now requires both a phase/IAAFT-like morphology-preserving null and a local block-shuffle null.

2. **Expanded result-schema mismatch.** Old expanded/focus runners accessed a top-level `observed_score` that the morphology analyzer did not return, creating a latent `KeyError` after expensive analysis. Canonical v1 owns one result schema in `janus_cosmos.core`.

3. **Wrong multiple-testing denominator.** The old expanded alpha multiplied by target count even though filter count already represented all target/filter tests. Canonical v1 uses `0.05 / (total_filters × primary_null_models)`.

4. **No empirical-p resolution gate.** Canonical v1 computes the minimum number of nulls needed to obtain p below the corrected alpha. Underpowered scientific runs fail closed. Explicit smoke runs are labelled `SMOKE_ONLY` and cannot create powered candidates.

5. **Seeds mistaken for replication.** Seeds now split and pool one Monte-Carlo null distribution. They are numerical chunks, not independent astrophysical replications.

6. **Repeated downloads.** Old expanded/focus code downloaded the same FITS product again for every seed. Canonical v1 downloads once, caches by source SHA-256 key, and records content SHA-256 provenance.

7. **Brittle live MAST discovery.** Several generations mixed low-level API calls and object-name resolution. Live expansion now uses the documented high-level astroquery path: `query_criteria` → `get_unique_product_list` → deterministic FITS product selection. The fixed seven-galaxy baseline remains independent of live discovery.

8. **Context-destroying preprocessing.** The legacy reader center-cropped to a square. Canonical v1 pads with the median before resizing so the complete 2-D field remains represented.

9. **Mutable overlapping correlation views.** Legacy shift correlation subtracted means in-place from overlapping source views. Canonical correlation copies/centers independent arrays and never mutates the input.

10. **Cross-band language exceeded evidence.** Two whole-image filters passing does not prove that the same localized sky structure persists. Canonical output is called `robust_cross_filter_candidate` and explicitly says WCS-localized persistence is not yet proven.

11. **Exploratory search contamination risk.** OCR, face-like detection, semantic analysis and cipher-like statistics are post-gate only, have `affects_blind_gate=false`, and cannot alter candidate p-values. Face processing is non-identifying. Post-hoc tuning on evaluation data remains forbidden.

12. **Windows launcher failures.** Historical local runners suffered from BOM corruption, execution from `C:\Windows\System32`, Unix `tee`, and missing working-directory normalization. The canonical BAT starts with `cd /d "%~dp0"`, uses PowerShell `Tee-Object`, logs to disk, and keeps the window open.

## Canonical geometry vector

The primary feature vector is frozen before evaluation and contains fourteen non-semantic measurements: directional correlation at scales 1/2/4, 90°/180° rotational agreement, gradient anisotropy, high-frequency energy, Fourier angular anisotropy, connected-component counts at q80/q90/q95, and largest-component fractions at q80/q90/q95.

The observed vector is standardized against a separate calibration-null ensemble. Held-out null vectors are then scored against the same calibration center/scale to produce empirical p-values.

## Candidate admission

A filter passes only when both primary null p-values are below the corrected alpha. A target passes only when at least two filters pass and the run is statistically powered. Missing downloads or FITS errors are recorded as errors and are never converted into a negative scientific result.

## Focus targets

NGC1425 and NGC1637 receive no special threshold, no feature weights and no admission shortcut. Their only special status is historical: they survived the earlier 256-null pilot and therefore make useful stress-test targets for the stronger gate.

## Operational acceptance gates

Canonical v1 is considered operational when offline unit/synthetic tests and module compilation pass; fixed public HLSP FITS can be downloaded and parsed; smoke runs are marked underpowered; powered runs refuse insufficient null counts; receipt/event logs are emitted; source SHA-256 values are recorded; and exploratory enrichment cannot affect the blind candidate flag.

## Claim ceiling

Canonical v1 can establish an image-level geometric candidate relative to its fixed feature/null family. It cannot by itself establish a new astrophysical object, unknown physics, artificial structure, hidden message or extraterrestrial intelligence.

The next scientific gate after a stable powered corpus result is WCS-aware localization: find where the excess originates, map candidate regions into common sky coordinates, and require localized persistence across independently calibrated filters/instruments and an independent dataset.
