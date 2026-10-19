# MOCSS HUDF 16004 Å — PSF injection into four actual October 2023 NIRCam images

**Date:** 2026-10-19. **Status:** original-MAST-pixel measurement; **instrument-anchored, analytic PSF sensitivity experiment**, *not* a STPSF-calibrated claim of intrinsic variability.

## New independent evidence

Use four archived, SHA-256-verified `JW01180` NIRCam F150W `i2d` SCI subsets acquired on **2023-10-01, 02, 03, 09**, in between the August and November NIRSpec detections of the 16004 Å feature. Measure directly at **RA 53.13127°, Dec −27.78908°** with the existing circular **0.12″** aperture and **0.24–0.38″** median-background annulus.

The official STScI NIRCam *resampled `i2d`* encircled-energy table reports **EE50 at 0.051″** and **EE80 at 0.134″** for F150W. We construct two mathematically distinct radial profiles, each reproducing both anchors to machine precision:

1. Single Moffat.
2. Two-Gaussian core/halo mixture, 72% core / 28% extended halo.

For every frame the source's sky coordinates are transformed with that frame's *actual* FITS WCS; each pixel receives the integral of the normalized PSF over a **7×7 subpixel grid**. We then redo the unchanged aperture/annulus estimator over the *actual original pixel values*. The original archive pixels are not edited. The injected source's total photometric F150W fluxes are **1.1, 2.5 and 5.5 nJy**; only the last is the paper's predicted broad-band flux from its tentative NIRSpec line flux.

We additionally inject the same 5.5 nJy point source into **95 matched-sky locations** in each of all four frames. All control positions, rejection mask, fixed epoch weights and empirical covariance are reproduced from the independent October covariance analysis (weighted covariance SD **0.662981 nJy**). Control positions share sky structure and thus are **not** 95 independent noise trials.

## Pointlike source hypothesis, under reported NIRSpec line calibration

| Measurement | Moffat | Core/halo Gaussian |
|---|---:|---:|
| Analytic EE inside 0.12″, before sky subtraction | 77.6% | 78.5% |
| **Actual FITS-pixel recovery after the local sky annulus** | **76.35%** | **76.76%** |
| Predicted source contribution if still 5.5 nJy | **4.199 nJy** | **4.222 nJy** |
| Matched-sky control injections with measured value at/below actual target | **0 / 95** | **0 / 95** |
| Difference / empirical matched-sky covariance SD | 6.33 (diagnostic only) | 6.37 (diagnostic only) |
| Difference / covariance plus 0.66 nJy model uncertainty | 5.04 (diagnostic only) | 5.06 (diagnostic only) |

The **real un-injected target value is 0.37098 nJy**, under the same 95-control, matched-grid median convention, and *not* 0.27850 nJy under the earlier, different per-image reference-median convention. Neither number implies a detection. Do **not** mix conventions.

The empirical 95-control *injected* weighted-aperture distribution spans approximately **2.586–6.588 nJy** (Moffat) or **2.615–6.626 nJy** (Gaussian), hence all injected controls exceed the un-injected target (0.371 nJy). This is a useful conditional sanity check, **not a calibrated p-value**, because matched sky cells are spatially correlated and the spectral-model flux has not been independently recalibrated.

## Astrometric miscentering changes the conclusion

We extract at the original NIRSpec sky coordinate while injecting 5.5 nJy at each of eight position angles for each specified sky offset (not observed physical motion):

| Deliberate injected–measured separation | Moffat median recovered fraction | Core/halo median recovered fraction |
|---:|---:|---:|
| 0.025″ | 75.5% | 76.2% |
| 0.050″ | 72.5% | 73.4% |
| 0.075″ | 66.0% | 65.8% |
| 0.100″ | 52.2% | 50.4% |
| 0.150″ | 12.6% | 13.6% |
| 0.200″ | 0.3% | 0.2% |

**Meaning:** if positional registration is good to well within 0.1″, the October non-detection strongly disfavors a persistent, compact source under the paper's F150W line-flux conversion. An otherwise constant source offset by ≈0.15″ might escape *this specific* forced aperture test. We have **not** determined that the actual NIRSpec–NIRCam registration error is 0.15″; this is solely a counterfactual stress test. The angle-dependent results in the machine-readable receipt are **not the direction of a jet, gas cloud, motion or emission**.

## Extended gas could also evade a small aperture

Convolve the analytic two-Gaussian point-spread profile with a *Gaussian intrinsic source of increasing sky size* (all test cases retain the same total integrated 5.5 nJy):

| Intrinsic source FWHM | Fraction recovered by this aperture | Recovered source contribution |
|---:|---:|---:|
| Unresolved | 76.8% | 4.22 nJy |
| 0.071″ | 74.2% | 4.08 nJy |
| 0.141″ | 59.5% | 3.27 nJy |
| 0.235″ | 35.3% | 1.94 nJy |
| 0.353″ | 16.4% | 0.90 nJy |
| 0.471″ | 8.1% | 0.45 nJy |

Thus the experiment rejects neither **faint extended line-emitting gas** nor **mislocalized/shifted emission**; the absence in a 0.12″ *point-source* aperture is not a valid spatially integrated upper limit for the proposed AGN-illuminated cloud.

## Scientific interpretation for source / distance / direction / reason

1. **Compact persistent emitter:** disfavored **conditionally** if (i) line → F150W total 5.5 nJy and (ii) image/source WCS agree to within approximately 0.05–0.075″ and (iii) actual PSF broadly follows the two STScI-anchored families.
2. **Transient event:** compatible, e.g. Type IIn or TDE, but **not confirmed** because the NIRSpec photometric calibration, observing epochs and sky-to-slit alignment have not been jointly solved with calibrated Level-2 spectra. The missing Paγ companion remains a major counterargument to the Paβ(z≈0.249) identity.
3. **Extended gas illuminated by another galaxy's AGN:** remains viable under this small-aperture test. A multi-aperture or spatially complete matched-filter search and real NIRSpec IFU tomography are required. Projected neighboring galaxy positions alone cannot establish a causal direction or distance.
4. **Spectral-background / MOS overlap:** also not definitively excluded by a broad-band imaging null; use the independent 16 NIRSpec MSA mask and calibrated GWCS spectrum workflow to test line provenance.

**We cannot yet claim the physical class, redshift, 3D propagation vector or maximum extent of an emitting filament.** The present independently repeated, falsifiable result is narrower: a *persistent correctly registered unresolved source* at the published flux has difficulty hiding in these four NIRCam images, but shifted/extended emitters can.

## Reproduce and audit

```bash
# Python needs numpy, scipy and astropy.
python3 tools/test_nircam_f150w_analytic_injection.py
python3 tools/measure_nircam_matched_covariance.py --data-dir /tmp/janus_f150w_time_domain --output /tmp/f150w_covariance_repeat.json
python3 tools/inject_nircam_f150w_empirical_psf_brackets.py --data-dir /tmp/janus_f150w_time_domain --output /tmp/f150w_psf_injection_repeat.json
```

The last script checks all four original MAST SCI-strip SHA256 values against `NIRCAM_F150W_TIME_DOMAIN_RECEIPT_2026-10-16.json`. Results with all control locations, all offsets, both PSF parameters and each epoch's injected flux are preserved in `NIRCAM_F150W_PSF_INJECTION_OBSERVED_2026-10-19.json`. The existing `measure_nircam_matched_covariance.py` and `recheck_nircam_time_domain.py` supply the trusted original byte ranges and sky-coordinate machinery.

**Primary references**

- STScI, *NIRCam Point Spread Functions*, resampled `i2d` F150W EE50/EE80 table: https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-point-spread-functions
- MOCSS HUDF 16004 Å original preprint, https://arxiv.org/html/2610.12553v1 (its F150W predicted 5.5 nJy remains an **assumption**, not our independent Level-2 calibration).
