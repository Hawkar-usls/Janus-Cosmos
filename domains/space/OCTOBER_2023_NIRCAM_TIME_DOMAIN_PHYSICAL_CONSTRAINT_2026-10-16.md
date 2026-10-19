# JANUS COSMOS: independent October 2023 time-domain test of the HUDF 16004 Å emitter

**Research state, 2026-10-16:** This is new independent *imaging photometry* from actual JWST/NIRCam `i2d` SCI pixels, not an independent detection of the NIRSpec line. It sharply constrains **persistent** emission but does **not** identify the emitting species, its distance, or a gas trail. All source-identification gates remain open.

## What was found, not just proposed

The line was reported by Huberty et al. (2026) with flux `(2.5 ± 0.3) × 10^-18 erg s^-1 cm^-2` on **23–24 August 2023** and **19 November 2023**, at `RA=53.13127°, Dec=-27.78908°` and `16004 ± 1 Å`. The authors' `synphot` model predicts that the **line alone**, if present at the reported flux, produces **5.5 nJy total F150W flux**. This is a prediction conditional on the line interpretation/calibration, not an independent broad-band measurement.

The JADES imaging emphasized in the original paper was taken in **October 2022**, about a year before spectroscopy. We independently found and processed **four additional F150W imaging observations in October 2023**, precisely *between* the two spectroscopy epochs. This closes that timing loophole for a **continuously bright** source.

| Actual NIRCam F150W date (UTC) | Distinct GO-1180 i2d mosaic | Empirical background-subtracted 0.12″ aperture (nJy) |
| --- | --- | ---: |
| 2023-10-01 | `jw01180-o014_t023_nircam_clear-f150w_i2d.fits` | `+0.453 ± 1.078` |
| 2023-10-02 | `jw01180-o016_t013_nircam_clear-f150w_i2d.fits` | `+0.341 ± 0.995` |
| 2023-10-03 | `jw01180-o013_t011_nircam_clear-f150w_i2d.fits` | `+0.504 ± 0.905` |
| 2023-10-09 | `jw01180-o024_t020_nircam_clear-f150w_i2d.fits` | `−0.546 ± 1.298` |
| **Inverse-variance combination** | Four independently dated GO-1180 observations | **`+0.278 ± 0.521`** |

All four are *forced* photometry at the **frozen** sky position, with a circular radius of `0.12″`, local background annulus `0.24–0.38″`, calibrated SCI pixel units `MJy/sr`, exact FITS WCS and 350 local blank-sky apertures for a robust per-mosaic noise estimate. These are not PSF-aperture-corrected or a vetted transient survey light curve. The third observation cross-checks an earlier independent 901-row cutout analysis (`+0.554 ± 0.966 nJy`) to within its expected difference from sky-control selection.

### Positive controls and astrometric stress test

- In the **same four image cutouts**, the three faint JADES neighboring galaxies are recovered at combined 0.12″-aperture fluxes **NW `3.030 ± 0.521`, East `3.256 ± 0.521`, South `3.268 ± 0.521` nJy**; compare the *total* F150W fluxes `4.75`, `4.69`, `3.57 nJy` in Huberty et al. Table 3. Extended-galaxy morphology and different apertures mean these are sensitivity/WCS controls, not exact aperture-correction standards.
- A preregisterable **113-position grid** in a radius `0.15″` around the putative source, spacing `0.025″`, yields a *maximum* weighted-aperture `S/N=0.884` (uncorrected for grid trials). At the frozen target position the grid gives `0.633σ`. No hidden point source is apparent anywhere in this limited positional range.
- Both facts make a gross WCS mistake or local dead patch a poorer explanation of the faint October result. **Larger** positional systematics, source-confusion/deblending, non-Gaussian mosaic covariance, and the original spectroscopic source position still need formal independent review.

### Physically interpretable quantitative bound, conditional

STScI's NIRCam reference PSF has F150W 50% encircled-energy radius `0.051″` and 80% radius `0.134″` on resampled images. A `0.12″` aperture therefore contains a substantial but not complete fraction of a centered point source. For an **illustrative 65–80% recovery fraction**, the paper's 5.5 nJy line-only prediction becomes **`3.58–4.40 nJy` in the measured aperture**, versus **`0.278 ± 0.521 nJy`** actually found. Under fixed-source position, independent-Gaussian-noise, no line-flux change and correct photometric mapping assumptions, that is a **nominal ~6.3–7.9 standard deviation discrepancy**. **This is not a final calibrated exclusion significance:** the recovery fraction, per-pointing weight maps, correlated noise, contamination and source-position uncertainty were not determined by injected-source tests.

At a representative `70%` recovery, a diagnostic **3σ upper bound** on the mean October aperture is `0.278 + 3×0.521 ≈ 1.84 nJy`. Converted through the *assumed* PSF fraction and author's spectral-to-F150W calibration this is a line-flux proxy **`≲1.2×10^-18 erg s^-1 cm^-2`** during October, versus the quoted **`2.5×10^-18`** in August/November. The contrast is ≥factor 2 in this **conditional** diagnostic, not a formally published line-flux limit.

The actual cadence implies (i) line high on Aug 23/24, (ii) F150W broad-band low from Oct 1–9, (iii) line high again Nov 19, *if* the NIRSpec flux and inferred sky position are accurate and the flux variation is intrinsic. A typical single-monotonic SN/TDE light curve or a spatially extended, continuously lit narrow-line gas cloud **cannot simply be assumed** from these dates. An intermittent compact event, rebrightening interaction, varying nuclear source, or a systematic spectral-position/flux issue requires testing.

If the change within approximately `38–41 observer-frame days` is intrinsic and **the entire emitting region varies coherently**, the simple light-crossing scale is `R ≲ c Δt_obs / (1+z) ≈ 0.033 pc/(1+z)`. Examples: `≲0.026 pc` for Paβ at z≈0.248, `≲0.014 pc` for Hα at z≈1.438, `≲0.010 pc` for [O III] at z≈2.196. This **conditional upper scale applies only to the variable emitting component**; it does not prove the gas has that size or constrain a faint, constant surrounding trail. A kpc-scale homogeneous gas structure cannot vary coherently on this timescale under ordinary light-travel-time arguments, although a small changing knot within a larger cloud could.

## Independent spectroscopic host/distance crossmatch

Read the official **JADES NIRSpec DR4 spectroscopic catalogue** (`Combined_DR4_external_v1.2.1.fits`, 5,190 target rows, SHA-256 `4b4e7fd4306208ee801406c7ed73fdab895a2e67184af7f639a91e5f9c8ab7d4`). The **nearest selected NIRSpec source center is 5.282″ away** and has `z_spec=-1`/flag E (no spectroscopic redshift). **Zero spectroscopically targeted JADES DR4 sources are within 5″** of the orphan, including all three galaxies within 1.1″. The next selected source at 5.859″ has only flag-D tentative `z=1.84051`, not a secure distance to the orphan. This fully prevents claiming that a neighbor has a spectroscopic redshift matching the mysterious line based on DR4.

An independent official **MAST CAOM 1.008″ cone** (3,288 observation metadata rows, 275 JWST rows, source query SHA-256 `29f75c440f4e05865a4b45a202bb36410cc9fd75da6c67eb822a73d800a70bef`) lists archived **NIRSpec/MSA prism footprint-covering rows** from GO-1210, 1286 and 3215. They are **not** proof that the particular orphan position fell in any active shutter; the original GO-3290 MOCSS products are often indexed as **images**, not spectra. Later NIRCam/GRISM F277W/F356W/F444W data do not cover the `1.6004 µm` line. No independently established spectroscopic identification or new source-centered IFU map was found in this archive crossmatch.

## What the evidence now favors, and what it does not

1. **Most discriminating positive result:** If the published line is truly from the frozen sky point, it was **not continuously present** at the August/November flux during the interleaved October 2023 images. This raises the priority of intermittent compact emission and of a rigorous spectroscopic origin systematic check.
2. **Paβ at z≈0.248 (dusty SN IIn / interaction):** nearby galaxies have compatible *but highly uncertain photometric* redshifts; Paγ at ~13657 Å is not convincingly detected. Historical 2014–15 MUSE Hα absence does not exclude a 2023 transient.
3. **Hα at z≈1.438 (dusty transient/TDE):** no expected optical pair establishes the redshift, and the timing is unconventional for a single smooth transient; physically still open.
4. **[O III] 5008 at z≈2.196 (compact transient AGN / ionized knot):** mandatory [O III] 4960 companion is **15850.79 Å**, *not 15864 Å* as quoted by the original paper. The predicted companion flux (~`0.83×10^-18`) is below or comparable to the paper's local `3σ` limit (`~0.94×10^-18`), so its apparent absence alone cannot decisively reject [O III]. The proper local-background-subtracted 15850.8Å extraction is still unperformed.
5. **Directed gas trail:** there is **no established emission-line map or trail direction**. The quoted coordinates define the direction *to* the sky point, not an outgoing physical beam. Neither the three galaxies' position angles nor the detector spectral dispersion supplies a gas-flow direction; no physical length or distance can be assigned without an identified transition/redshift and resolved line map.

**Net answer at current evidence level:** The line is a reported astronomical emission feature; its October 2023 absence in **four independently checked F150W image epochs** now strongly disfavors a *continuously luminous compact line source* at that sky position. It remains **scientifically impossible to specify its exact physical source, distance, and gas-trail endpoint from currently extracted lines alone**. A number would be fabricated. An independent 15850.8Å/Paγ JWST spectral extraction or new simultaneous IFU spectroscopy and imaging is the missing identifying measurement.

## Reproducibility and original sources

- Primary discovery: https://arxiv.org/html/2610.12553v1 (published 2026; observation timeline and 5.5 nJy `synphot` calculation).
- Official MAST archive API: https://mast.stsci.edu/api/v0/invoke ; GET byte-range science FITS by `mast:JWST/product/` URIs.
- NIRCam resampled encircled energy reference: https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-point-spread-functions (Table 3).
- JADES DR4 official spectroscopic catalog: https://jades.herts.ac.uk/DR4/Combined_DR4_external_v1.2.1.fits ; readme https://jades.herts.ac.uk/DR4/Readme_DR4_catalogues.md.
- Analysis scripts: `tools/recheck_nircam_time_domain.py`, `tools/summarize_nircam_time_domain.py`, `tools/crossmatch_jades_dr4_mast_cone.py`.
- Immutable summary receipts: `domains/space/NIRCAM_F150W_TIME_DOMAIN_RECEIPT_2026-10-16.json`, `domains/space/NIRCAM_F150W_ASTROMETRY_RECEIPT_2026-10-16.json`, `domains/space/JADES_DR4_MAST_COVERAGE_RECEIPT_2026-10-16.json`.
- Raw external FITS and MAST response JSON are not committed because of size; exact archive product names, checked SCI byte ranges, retrieved input SHA-256, official URLs, and results are recorded for independent replay.

**Scientific status gates:** `NEW_TIME_DOMAIN_CONSTRAINT_VALIDATED_WITH_LOCAL_CONTROLS`; `SPECTRAL_LINE_IDENTIFICATION_UNRESOLVED`; `DISTANCE_UNRESOLVED`; `TRAIL_UNDETECTED_AND_UNBOUNDED`.
