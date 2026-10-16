# JANUS COSMOS / TOPA — 16-exposure MOCSS geometry, line identity, and origin

**Research state: `PARTIAL_INPUT_VERIFICATION`, 2026-10-16.** We have not
independently reconstructed the 16004 Å line in detector-coordinate spectra,
nor mapped its spatial extent. The source type and distance remain unknown.

## Directly verified new archival facts

The official JWST GO-3290 MAST observation list has **16** long MOCSS G140M/F100LP
exposures (each about 8315.667 seconds), with **two detector RATE files per
exposure** (32 detector records). Using exact 0–32767-byte FITS header ranges for
30 frames and two previously downloaded full RATE frames, we verified the
`MSAMETFL` and `MSAMETID` primary-header metadata for every exposure. The
16 source FITS headers and detector mappings are preserved in
`MOCSS_16_EXPOSURE_MANIFEST_2026-10-16.json`.

The NIRSpec MSA archive metadata is unusually important in this all-open-shutter
technique. Two MAST **original, complete 80-MB MSA FITS files** (one per
orientation) each contain `SHUTTER_INFO` rows from **all eight dithers** in that
orientation, so the `msa_metadata_id` entries for all 16 exposures can be
cross-validated. **16/16 IDs matched.** A reproducibility pitfall is that
`SHUTTER_INFO.msa_metadata_id` is encoded as signed 16-bit, whereas the
original RATE header IDs extend past 65535 (e.g. 66129 → 593). The verifier
normalizes these as signed `int16` after a modulo-65536 transform. Comparing
the raw integers would incorrectly reject legitimate shutters.

For representative RATE observations from **both observing orientations**,
MAST listed RATE, RATEINTS, UNCAL and MSA products but did **not** list
`cal`, `s2d` or `x1d` calibrated products. This is not an exhaustive search of
all higher-level science-product archives; it documents what the selected
MAST observation product lists offered when queried.

**Crucial limit:** The 32 FITS *headers* and all 16 MSA associations were checked,
but we did **not** download/reduce all 32 complete detector science arrays.
Only two full RATE frames had previously been inspected at detector level.
The published line is not independently re-extracted. Source sky-to-detector
mapping needs the stage-2 `assign_wcs` transformation, pinned JWST pipeline
v1.19 / CRDS `jwst_1413.pmap`, per-exposure shutter culling, correct IDs,
row-wise sky subtraction, dither alignment, and line-injection recovery.

## Physically useful new checks

1. **Published [O III] wavelength inconsistency.** For an observed primary
   [O III] 5008.239 Å line at 16004 Å, the companion [O III] 4960.295 Å
   should appear at **15850.79 Å**. The discovery paper explicitly states
   **15864 Å**, **13.21 Å** away, while quoting a nondetection limit evaluated
   at the latter position. This is large compared to the reported 13.2 Å
   FWHM of the primary. A negative at the wrong wavelength does not test the
   right wavelength. Moreover, the expected companion flux is only about
   **0.83 × 10⁻¹⁸ erg/s/cm²**, below the paper's quoted **0.94 × 10⁻¹⁸**
   (3σ) local limit. The exact 15850.8 Å neighborhood must be independently
   measured before inferring anything from its nondetection. **We are not
   claiming a companion detection.**

2. **No significant four-stack flux variability yet.** Using the paper's own
   A/B/C/D fluxes and uncertainties (2.3±0.5, 3.5±0.6, 2.9±0.7, 2.9±1.2,
   in units of 10⁻¹⁸ erg/s/cm²), the inverse-variance weighted constant-flux
   model has χ²=2.3835 on 3 degrees of freedom (survival probability ≈0.497).
   A transient lasting 88 days is not excluded; a strong change between the
   August and November spectroscopy is not established. The 2023-10-03
   NIRCam F150W nondetection in our earlier archive replay, *between* the two
   spectral epochs, makes a simple continuous narrow-line source worth
   scrutinizing—but its small uncorrected aperture and positional/PSF
   uncertainties are not a certified flux limit.

3. **Possible AGN cloud illumination directions, conditional on a common
   redshift.** Applying the paper's JADES Table 3 sky coordinates and a
   z≈0.249 proper angular scale of 4.02 kpc/arcsec:

   | Proposed host | Relative sky bearing from line (E of N) | Angular separation | Projected separation at z≈0.249 | Minimum one-way illumination time if truly associated |
   | --- | ---: | ---: | ---: | ---: |
   | NW, z_phot≈0.22 | 328.23° | 0.6352″ | 2.553 kpc | ~8,328 yr |
   | S, z_phot≈0.19 | 174.54° | 1.0379″ | 4.172 kpc | ~13,608 yr |
   | E, z_phot≈6.03 | 95.14° | 1.0041″ | 4.036 kpc | ~13,165 yr, but likely unrelated |

   These separations are *projected*, and the physical numbers presume the
   Paβ redshift for **both** the cloud and the neighboring galaxy. They are
   not distances to Earth, not measured beam lengths, and not evidence for
   an AGN: all three associations await spectroscopy. The NW/S photometric
   uncertainties do not by themselves establish z≈0.249.

4. **One NIRSpec IFU footprint covers the immediate neighborhood.** Official
   STScI documentation gives **3″ × 3″**, sampled at 0.1″. Centered on the
   line-emitter, it can include all three candidate neighbors at ≤1.04″.
   A genuinely blind ±5″ off-source map needs a larger mosaic/older MOCSS
   reduction. Test spatial coincidence of line emission with nuclei and
   photoionized gas, not just alignment of raw detector spectral traces.

## Explicit hypothesis gates

| Candidate | In its favor | Needs independent confirmation |
| --- | --- | --- |
| Type IIn SN / Paβ at z≈0.249 | Narrow line, weak possible nearby host, delayed continuum | Paγ at ~1.365 μm (dust-corrected flux ratio), real host z; light curve including October 2023 |
| Dusty TDE / Hα | Transient continuum and line ratios can be unusual | Local host/black hole, other lines, IFU kinematics, multi-epoch line profile |
| Gas cloud illuminated by neighbor AGN | Nearby galaxies; narrow, off-nuclear line and extreme equivalent-width geometry | Shared spectroscopic z, two-dimensional 16004 Å morphology, AGN diagnostics, line ratios |
| [O III] 5008 at z≈2.196 | Narrow line; incorrect published 4960 companion check motivates correct one | 15850.8 Å spectrum and sky position, Hβ near 15538.9 Å, host redshift |

## Reproduce and next action

```sh
python3 -m pip install astropy scipy numpy
python3 tools/build_mocss_16_exposure_manifest.py --archive /tmp/janus_jwst_3290_observations.json --data-dir /tmp
python3 tools/analyze_mocss_geometry_and_hypotheses.py --data-dir /tmp
```

The second command requires archived MSA files
`jw03290001001_05_msa.fits`, `jw03290002001_01_msa.fits` and representative
MAST products-query JSONs as described in the verification script. The
script makes **no automatic claim of a spectral detection**.

Before a physical classification or a spatial extent can be accepted, a
follow-up experiment must project all 16 MOCSS RATE observations through
the right per-shutter WCS, independently find the 16004±10 Å line, measure
correct 15850.8 Å / Paγ windows, perform off-source tests and an
injection/recovery calibration. If that is inconclusive, request fresh
same-epoch NIRSpec IFU + NIRCam observations.

### Primary sources

- Huberty et al., source paper, https://arxiv.org/html/2610.12553v1
- Huberty et al., MOCSS methods Appendices C–E, https://arxiv.org/html/2609.36003
- JWST GO-3290 public MAST data, https://doi.org/10.17909/tct8-jr92
- STScI NIRSpec IFU FoV, https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-observing-modes/nirspec-ifu-spectroscopy
