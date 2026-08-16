# JANUS TACHYON STAR T1 — Look-back A→B→C Holdout Protocol

## Question
Can an isolated one-exposure point-source motif be reproduced at frozen Palomar coordinates when both the sky position and the temporal center are selected from metadata before any target pixels are read?

## Why this is a Tahyon Traveler descendant
The older Tahyon Traveler experiments separated preregistered target time, execution, Git write/materialization, Git timestamp, and later observation. T1 transfers only that epistemic architecture: **seal first, observe second**. Distributed scheduler/Git ordering effects are not treated as physical evidence.

## Training motif
TACHYON-STAR-001 in LUCI1 Br_gam: qualified absence at -232.5827 s, a marginal R1-compatible candidate at t0, and qualified absence at +432.2595 s. This motif defines timing weights only; it does not define a tachyon.

## Pixel-blind pool
Start from the 2F-E frozen exact-WCS table (443 pairs / 403 files / 42 sources). Remove all 92 files used by the 2F-E representative pixel stage and all 100 files frozen for 2F-F targeted recovery. The remaining pool is 239 rows / 211 files / 24 sources.

## Metadata-only selection
Within each `(src_id, instrument, filters)` group, sort by DATE-OBS. An eligible A-B-C triple is three consecutive untouched rows with both adjacent gaps in [120, 480] seconds. There are 24 eligible triples across 9 sources. For the primary analysis choose exactly one per source by minimizing `|Δpre-232.5827| + |Δpost-432.2595|`; ties are broken by B timestamp, instrument, filter, filename. The resulting 9 trials are frozen in the manifest SHA-256 `334fe0cacd0b056214d991be5256e05ab678fe5edb072d46a26ad07d4ae803b7`.

## Observation contract
Use unchanged LUCI reader, exact-WCS HDU binding, R1 frame gate, `measure_psf_at`, and the already-frozen local SNR 8/12 injection sensitivity contract. A and C must be sensitivity-qualified absences. B must contain an R1-compatible point-source candidate at the exact frozen coordinate. Failed sensitivity is UNRESOLVED, not absence.

## Independent detector witness
For six primary B epochs, metadata already provide an untouched same-filter exposure from the other LUCI instrument within 10 seconds. If a primary B candidate occurs, this paired frame is adjudicating only if its actual integration interval overlaps at least half of the shorter exposure. Same-coordinate detection in both instruments is a high-priority sky-origin candidate; single-instrument detection with a sensitive simultaneous absence is evidence against a sky-origin interpretation. Neither outcome identifies a tachyon.

## No post-pixel tuning
No retiming, coordinate reselection, source removal, threshold relaxation, or replacement of a blocked trial after pixel access. Exploratory analysis of all 24 eligible triples is reported separately and cannot upgrade the primary result.

## Claim ceiling
`PREREGISTERED_LUCI_HOLDOUT_ISOLATED_TRANSIENT_MOTIF_TEST_ONLY__NO_FTL__NO_RETROCAUSALITY__NO_TACHYON_IDENTITY__NO_NUCLEAR_CAUSALITY__NO_UAP_ORIGIN_CLAIM`
