# JANUS COSMOS | JWST MOCSS HUDF 16004 Å — independent sky-registration replay

**Research date:** 2026-10-16  
**Status:** Independent WCS geometry validated; raw detector contrasts exploratory; no calibrated line flux, source identification, or physical trail length established.  
**Archive:** MAST JWST GO-3290, DOI `10.17909/tct8-jr92`.  
**Published case:** Huberty et al. (2026), arXiv `2610.12553v1`.

## 1. What was independently done

1. Retrieved **all 32 original NIRSpec detector RATE files (16 exposures × NRS1/NRS2)**, 83,966,400 bytes each, and two MSA metadata files containing the relevant eight configurations per orientation. The compact receipt records individual SHA-256 fingerprints for all 34 original FITS files.
2. Replayed each exposure with **official JWST pipeline 1.19.1, CRDS `jwst_1413.pmap`**, and exposure-specific `SCI` WCS pointing, using **reversible TEST-only virtual MSA shutters** close to the published sky position RA 53.13127°, Dec −27.78908°. The original FITS data were never modified.
3. Reconstructed source detector pixels in all 16 visits. August A/B map to **NRS2** and November C/D to **NRS1**, consistent with the paper's two-angle and 21-arcsec chip-gap-covering dithers. The source is near shutter bars in **7 of 16 dither pointings**, which forbids naive equal-throughput comparisons.
4. Independently projected three requested wavelengths into each detector: **16004 Å**, **15850.7933 Å** for [O III] 4960 if the primary is [O III] 5008, and **13656.7483 Å** for Paγ if the primary is Paβ.
5. Built four *WCS-registered 2D RATE* coadds A/B/C/D, first with constant per-row pedestal subtraction and then with a **41-pixel running median subtraction along dispersion** to control prominent smooth detector gradients. Repeated the comparison with 31- and 61-pixel windows. These tests are **not** standard STScI Stage 2 spectrophotometric reduction.
6. Reprojected the raw 16004 and 15850.8 Å cutouts to the same **3.0″ × 3.0″ local sky grid**, using independent WCS Jacobians for A/B and C/D. Sampled the published source and all three JADES neighbor locations as preliminary, *uncalibrated* spatial controls.

## 2. Independently reconstructed detector positions

For selected dither anchors; detector x/y are 0-based RATE pixels, not source sky offsets. All 16 exact coordinate sets are in the SHA-bearing machine receipt.

| Stack | Example exposure | Detector | x at 16004 Å | y at 16004 Å | Coverage |
|---|---:|---|---:|---:|---|
| A | 1 | NRS2 | 967.05 | 386.76 | inside the virtual open shutter |
| B | 8 | NRS2 | 754.02 | 366.35 | inside the virtual open shutter |
| C | 10 | NRS1 | 1749.66 | 289.06 | inside the virtual open shutter |
| D | 13 | NRS1 | 1556.74 | 295.22 | inside the virtual open shutter |

**Important metadata handling:** `MSAMETID` 37338, 66129 and other >32767 original header integers wrap in signed-16-bit MSA shutter tables. The JWST 1.19.1 MSA row selection compares the header ID directly to that table. The wrapper script normalizes the ID **only in a temporary TEST-rate copy** and records both original/normalized values; source FITS are immutable.

## 3. Raw detector-space consistency, *not significance*

The following values are **uncalibrated matched-filter contrasts divided by local empirical MAD** in the 4-frame coadd. They are **not Gaussian sigmas, detection significances, calibrated fluxes, or upper limits**. Correlations, unmodelled flat fields, background residuals, shutter throughput and source selection prohibit interpreting them that way.

| Stack | 16004 Å: row constant | 16004 Å: median 31 px | 41 px | 61 px | [O III] 15850.8 Å: 41 px | Paγ 13656.75 Å: 41 px |
|---|---:|---:|---:|---:|---:|---:|
| A | 0.10 | 3.28 | 2.93 | 2.64 | −0.44 | −0.62 |
| B | 0.11 | 3.82 | 3.78 | 2.91 | −0.23 | +2.04 |
| C | 0.56 | 0.38 | 0.51 | 1.02 | +0.02 | +1.59 |
| D | 2.56 | 2.23 | 2.49 | 2.29 | −0.19 | +1.62 |

**Interpretation:** A narrow detector feature centered at the independently predicted 16004 Å pixel becomes visible in **A, B and D** after removing a strongly wavelength-dependent smooth background. Its behavior is qualitatively stable across three high-pass scales. Stack C remains weak in this approximate, not flux-calibrated analysis. The dramatic changes from constant-pedestal to high-pass processing demonstrate that **a science-grade background/flat model is essential** before reporting a new statistical detection or deciding whether the source varied. The Paγ probe is mixed across groups and cannot yet confirm or rule out a Paβ model.

**Correct [O III] window:** If 16004 Å is [O III] 5008.239 Å, the [O III] 4960.295 Å partner must be at **16004 × 4960.295/5008.239 = 15850.793 Å**. The source paper describes a check at **15864 Å**, about **13.2 Å** away (roughly one reported FWHM). Its quoted 3σ limit at 15864 Å is **not** a demonstrated limit at 15850.8 Å. Our raw-rate check found no *stable strong* candidate at 15850.8 Å, but cannot set a meaningful flux upper limit without a calibrated spectrum and injection sensitivity.

## 4. Local sky-direction and neighboring-galaxy controls

This is a **two-position-angle WCS-projected RATE quick-look**, not an IFU data cube. A 0.1″ grid over ±1.5″ is centered on the published coordinates. Detector Jacobians are independently calculated from the official calibration for each group:

- A/B: approximately `d(x,y)/d(E,N) = [[−6.65,+7.26],[+6.95,+6.69]]` detector pixels per arcsecond.
- C/D: approximately `[[+7.40,+6.35],[+6.62,−7.03]]`, providing an almost perpendicular independent geometry.

When each group map is only *descriptively* normalized by local off-source MAD, the preselected source position has a combined relative contrast of **+2.28** at 16004 Å versus **−0.52** at 15850.8 Å. The 16004 source position is higher than approximately 95% of sky-grid cells; **this 5.1% rank is not a p-value or false-alarm probability**. The actual highest pixel lies off-target, showing why artifact rejection and a trials correction are indispensable.

| Neighbor sky position | Offset east | Offset north | Distance | Bearing east of north | Combined 16004-A relative raw contrast |
|---|---:|---:|---:|---:|---:|
| NW JADES galaxy (z_phot≈0.22) | −0.334″ | +0.540″ | 0.635″ | 328.2° | −0.58 |
| South JADES galaxy (z_phot≈0.19) | +0.099″ | −1.033″ | 1.038″ | 174.5° | +0.66 |
| East JADES galaxy (z_phot≈6.03) | +1.000″ | −0.090″ | 1.004″ | 95.1° | +0.34 |

These comparisons show **no repeatable same-wavelength extension specifically linking the candidate and any of the three galaxy centers** in the exploratory maps. They **do not** impose a valid physical-length or projected-flux bound. Astronomical sky coordinates are a photon **arrival direction**, not a measured beam or propagation axis. Distance requires a confirmed line identification/redshift.

## 5. Next falsifying experiment

- Reconstruct all four dispersed stacks with the MOCSS-specific **full-row detector background**, detector-level flats, Stage 2 spectroscopy calibration (paper skips `master_background`, `pathloss`, `barshadow`), error propagation, and sky-position-dependent shutter throughput. Perform sky/null controls and synthetic line injection/recovery with look-elsewhere accounting.
- Fit and publish **separate** line fluxes and credible limits at 16004±10 Å, **15850.8±10 Å**, **13656.75±10 Å**, plus unblinded off-line wavelengths, both epochs, and offset sky positions. Compare a pointlike line-spread function with extended morphology in A/B versus C/D.
- A pointed **NIRSpec IFU 3″×3″** field centered on the candidate nominally encompasses all three cataloged neighbors (maximum separation ≈1.04″). Use a **mosaic** for a hypothesized trail beyond that field. The IFU should check whether any photoionized cloud spectrum, velocity field and morphology point back to a galaxy/AGN; a transient needs contemporaneous imaging and spectroscopy.

**Bottom line:** This replay **actually did restore full 16×2 detector sky geometry and approximate independent 2D line/control maps**, but **did not identify the source, show a beam, or produce science-grade spectrophotometric detection statistics**. Preserve both the positive consistency and the unresolved limitations.
