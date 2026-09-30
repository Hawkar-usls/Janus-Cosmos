# LUCI-PALOMAR-JPFM-2F-D Exhaustive Coverage Protocol

## Purpose

2F-D asks a coverage/opportunity question before any astronomical interpretation:

> Does any member of the complete frozen POSS-I S0 candidate population fall inside an exact public LUCI/LUCIFER imaging FITS footprint?

The test is archive-first and fail-closed. A positive overlap is not an anomaly result. It only opens an independent decades-later near-IR counterpart test.

## Frozen parent

The source-extraction method is inherited from `LUCI-PALOMAR-JPFM-2F-C-R1`, whose six real LUCI validation frames passed preregistered injection-recovery and hot-pixel specificity thresholds.

## Full Palomar population

The complete S0 table is bound to:

- POSS repository commit: `4005e200541b321ead3d6608f0162a14430ef1c2`
- gzip SHA-256: `f19cf987756c62a68f55a472992d860e73ae63b3a4664189092b0e1fda77f7bb`
- decompressed CSV SHA-256: `2ff92f2210acb387ef9ef4b88d561595d3883e9aab27065042627272b96590f0`
- expected rows: `122820`

Only `src_id`, RA and Dec are used for coverage matching.

## Archive-first freeze order

1. Verify the fixed S0 bytes.
2. Write and SHA-256 the full 122,820-coordinate corpus.
3. Query the public LBT archive for all `FREE`, `SCIENCE`, `Mirror/mirror` LUCI rows **without using Palomar coordinates**.
4. Sort, write and hash that LUCI inventory.
5. Crossmatch frame centers against the entire Palomar corpus using unit-vector `cKDTree` geometry and frame-specific half-diagonals.
6. Freeze and hash all coarse pairs.
7. Apply archive-metadata WCS as a non-admitted download preflight.
8. Freeze and hash that metadata-WCS pair set.
9. Only then may any candidate FITS be downloaded.
10. Verify exact containment from the FITS celestial WCS.
11. Every exact-overlap frame must independently pass the same R1 injection-recovery/hot-pixel gate.
12. Only a passing frame may be inspected for an IR source at the Palomar coordinate.
13. If a source exists, compare its PSF morphology with same-frame local sources matched in peak SNR.

## Preregistered fail-closed limits

If the frozen metadata-WCS set requires more than 250 unique FITS files, the workflow stops before pixel inspection and preserves the overlap set for staged replay.

For an exact-overlap frame:

- all-star injection recovery >= 0.80;
- injection recovery for SNR >= 8 >= 0.90;
- hot-pixel acceptance <= 0.05.

Matched local morphology controls must be within 300 px, have peak-SNR ratio 0.5–2.0 relative to the target, and provide at least 8 controls (up to 20).

## Claim ceiling

No temporal/UAP labels enter archive enumeration, coordinate matching, extraction or morphology. Outcomes cannot establish anomaly, artificial origin, UAP origin, ETI, or causality.

The maximum admitted result is an independent near-IR counterpart/no-counterpart observation for exact Palomar coordinates with validated per-frame extraction performance and same-frame morphology controls.
