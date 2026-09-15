# LUCI × Palomar JPFM Transfer Protocol v1

## Purpose

This protocol asks whether the frozen Palomar `JPFM-2F-B` star-morphology measurement contract can be transported to an independent near-infrared instrument, LUCI/LUCIFER on the Large Binocular Telescope, without importing Palomar temporal/UAP labels into the LUCI analysis.

It also performs a separate sky-overlap preflight for the exact deterministic 64-source Palomar pilot sample. The two questions are intentionally not pooled.

## Frozen Palomar parent

The runner binds to the immutable `janus-meta-registry` commit `7890cd5c8f4650c02dd439dbf96f09bc45638654` and to the public POSS-I release commit `4005e200541b321ead3d6608f0162a14430ef1c2`.

Required parent hashes:

- structural manifest gzip: `166f5e6621ed2b065b7981b3c8208670f3c989b1394bd559c9005ab1fa6d07d9`
- structural manifest CSV: `34b0ccde7c3683d07626774e52dac0a197451f729242204e59aae81397bdbc2e`
- POSS-I S0 gzip: `f19cf987756c62a68f55a472992d860e73ae63b3a4664189092b0e1fda77f7bb`
- POSS-I S0 CSV: `2ff92f2210acb387ef9ef4b88d561595d3883e9aab27065042627272b96590f0`

The 64-source sample is rebuilt deterministically as in JPFM-2F-B: for each structural cluster 0..15, two sources closest to the cluster median anomaly score plus two with the largest anomaly score, with `src_id` tie-breaks.

## Arm A — morphology-contract transfer

Real public LUCI imaging FITS are admitted only through the existing LUCI fail-closed provenance gate. Spectroscopic frames or non-LUCI instruments remain excluded.

The following Palomar metric names are preserved:

- moment FWHM, major and minor axes
- elongation and orientation
- core concentration R2/R6
- threshold area
- circularity
- convex-defect fraction
- circle deviation
- local reference count and median FWHM
- FWHM ratio
- radial-profile difference

The LUCI detector polarity is not assumed blindly. The direct-positive image convention is used unless the robust negative high tail is >1.5× the positive high tail, in which case the frame is inverted and that choice is recorded.

A LUCI frame passes the transfer-feasibility gate only when at least 12 separated sources are measured, at least 95% of retained sources have finite FWHM and elongation, and at least 50% have three or more same-frame amplitude-matched references. The overall pilot is `TRANSFER_FEASIBILITY_PASS` only if at least 75% of tested frames pass. This is a measurement-feasibility gate, not an astrophysical anomaly gate.

## Arm B — direct Palomar-to-LUCI sky overlap preflight

The exact frozen 64 Palomar coordinates are queried against the official LBT TAP service. Only public (`FREE`), `SCIENCE`, `GRATNAME=Mirror` LUCI rows are considered.

A first-stage archive-centre search uses a deliberately generous 0.20 degree radius. Returned rows are then screened using archived `CRVAL1/2`, `NAXIS1/2`, and pixel scale to form a conservative half-diagonal image radius.

Even a surviving geometry overlap is **not** a counterpart detection. It only authorizes a later gate that must download the named FITS, prove that the frozen Palomar coordinate falls inside the actual LUCI WCS footprint, and apply a preregistered local-source/counterpart criterion.

## Firewalls and claim ceiling

The LUCI morphology arm does not receive nuclear dates, Blue Book/NUFORC labels, UAP labels, geomagnetic indices, lunar state, or the JPFM-2D/2E temporal outcome.

The Palomar 2F-B pilot itself previously failed closed on source recovery; this protocol does not reinterpret that failure as a Palomar morphology PASS.

Maximum claim after this protocol:

`METHOD_TRANSFER_AND_ARCHIVE_OVERLAP_PREFLIGHT_ONLY__NO_PALOMAR_ANOMALY_CONFIRMATION__NO_UAP_ORIGIN_IDENTIFICATION__NO_CAUSALITY`
