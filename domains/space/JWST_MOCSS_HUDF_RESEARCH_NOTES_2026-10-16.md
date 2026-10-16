# JANUS COSMOS × TOPA: HUDF 16004 Å line — archival follow-up

**Status (2026-10-16):** Real published JWST spectral-line detection; physical emitter unknown. This addendum includes a limited independent archive replay and a pre-registered spatial-origin protocol. No independent 1-D spectral re-extraction or full survey sky cube has been claimed.

## Primary materials

- Huberty et al. (2026), *A Single Line-Emitter Identified in an Unbiased Search with JWST/NIRSpec*, DOI `10.3847/1538-4357/ae9a94`, https://arxiv.org/html/2610.12553v1.
- Huberty et al. (2026), *MOCSS: Multi-Object Coronagraphy for JWST/NIRSpec Slitless Spectroscopy*, https://arxiv.org/html/2609.36003v1; Table 1 is the 24-source catalog (23 known emitters plus orphan source ID 24).
- MAST archive JWST GO-3290, https://doi.org/10.17909/tct8-jr92.

The target lies at **RA 53.13127°, Dec −27.78908° (J2000)**. The paper reports λ=16004±1 Å, integrated flux (2.5±0.3)×10⁻¹⁸ erg s⁻¹ cm⁻², and a line visible in four stacks from 2023-08-23/24 and 2023-11-19. Those spectral values are **author-reported**, not independently remeasured here.

## What was independently processed

1. Queried MAST GO-3290 (80 observation-detector records, including 32 long detector-science records) and the 7.2″ sky cone.
2. Downloaded and parsed **two original NIRSpec NRS2 rate FITS** (2023-08-24 and 2023-11-19, each 2048×2048, six integrations, 8315.667 s exposure). Measured detector quality/flag summaries. These images were **not** projected into a sky-and-wavelength 1-D spectrum because the custom open-MSA MOCSS geometry and reductions are still required.
3. Fetched three public NIRCam **F150W** source catalogs and exact rows of calibrated `SCI` FITS mosaics. In each epoch, transformed the frozen source RA/Dec to the proper FITS WCS pixel, summed a 0.12″-radius aperture after local 0.24–0.38″ annulus background, and estimated a robust per-aperture noise from 250 blank-sky controls.
4. Queried catalogued pulsars only as a sanity check: ATNF mirror had zero entries within 1°, positive-control Crab pulsar returned a nonempty catalog stream. This does not address unknown pulsars and does not identify this emitter.

| NIRCam F150W epoch | Diagnostic background-subtracted 0.12″ aperture flux | Empirical SNR |
| --- | ---: | ---: |
| **2023-10-03** | **0.554 ± 0.966 nJy** | +0.64 |
| 2025-10-06 | −0.358 ± 1.479 nJy | −0.12 |
| 2026-10-10 | −1.143 ± 1.679 nJy | −0.74 |

The 2023 October image is **between** the August and November positive spectroscopic measurements, making its nondetection an important challenge to a simply continuous line-emission model. The paper predicted **approximately 5.5 nJy total F150W equivalent** for the 2023 line flux alone, but this number is not the expected flux *inside our smaller uncorrected aperture*. Before imposing a quantitative exclusion, measure the F150W encircled-energy correction, source astrometric uncertainty, mosaicked inverse variance/PSF and injection-recovery completeness.

All three independent catalogues likewise have **zero catalog entries within 2″**. This catalog-level negative result is weaker than calibrated forced photometry and does not exclude faint objects (the paper specifically discusses three close JADES faint sources below or around detection limits).

## Does the line appear somewhere else, and which way does it go?

- **Known sky direction:** the RA/Dec above, corresponding approximately to galactic longitude 223.54°, latitude −54.42°. This identifies the incoming line of sight, **not** a directed beam axis.
- **Neighboring galaxies** reported in the original paper: northwest at approximately 0.64″, east at 1.00″, and south at 1.04″. They may in principle host a nearby ionizing source, but have not been connected to the spectral line by a measured common redshift.
- **Published same-wavelength catalog screen:** using Table 1 of the MOCSS method paper, 23 other sources and 41 representative listed transitions were checked by λobs=λrest(1+z). **Zero** predicted λobs within ±100 Å of 16004 Å. Closest representative match: Hα of MOCSSHUDF1 at approximately **15722 Å**, position **4.8″** from the target. This is **not an exhaustive wavelength-bin search** through all raw MSA shutters or the larger Universe.
- **Spatial extent:** not measured. Slitless detector dispersion is not a projected light trail. The original paper's two dispersion angles constrain the inferred source position; they cannot alone determine how far line-emitting gas extends across space.

### Conditional distances if the transition can be identified

For approximate flat ΛCDM (H₀=67.66 km s⁻¹ Mpc⁻¹, Ωₘ=0.3111), these are scenarios, **not** determined distances:

| Proposed transition | Redshift | Lookback time | Luminosity distance | Proper scale at the source |
| --- | ---: | ---: | ---: | ---: |
| Paβ (possibly SN IIn) | 0.249 | 3.02 Gyr | 1.29 Gpc | 4.02 kpc/″ |
| Hα | 1.438 | 9.36 Gyr | 10.63 Gpc | 8.67 kpc/″ |
| [O III] 5008 | 2.196 | 10.80 Gyr | 17.84 Gpc | 8.47 kpc/″ |
| [O II] 3730 | 3.291 | 11.86 Gyr | 29.07 Gpc | 7.65 kpc/″ |
| Lyα | 12.165 | 13.43 Gyr | 132.32 Gpc | 3.70 kpc/″ |

Luminosity distance is a standard flux-distance definition and **must not be reported as the geometric length of a beam**, especially in an expanding Universe. The corresponding lookback times describe how long the photons traveled. The redshift and physical size of this source remain unknown.

## Reproduce the archive subset

In the Janus-Cosmos checkout (Python 3.11+):

```sh
python3 -m pip install numpy astropy
python3 tools/fetch_jwst_mocss_mast_subset.py --output-dir /tmp
JANUS_DATA_ROOT=/tmp JANUS_RECEIPT_OUTPUT=/tmp/janus_mast_replay.json python3 tools/replay_jwst_mocss_mast_subset.py
python3 tools/verify_jwst_mocss_archive_receipt.py
```

The fetcher downloads exact MAST `SCI` byte ranges rather than the full multi-GB mosaic; filenames, inclusive ranges, byte counts and SHA-256 of observed inputs are recorded in the committed receipt. The offline verifier checks the **committed** archive receipt and science gates, not a network rerun. Compare the newly generated `/tmp/janus_mast_replay.json` with the committed record and account for any subsequent archive reprocessing.

## Next decisive tests

1. Reconstruct all **16 MOCSS spectral exposures** in two sky orientations with pinned CRDS/pipeline, blind off-source controls and line injection/recovery. This must produce an independent λ=16004Å image, a 2-D map and line-flux/line-ratio limits.
2. Search fixed sky apertures across at least **±5″** for **16004±10 Å** separately in both orientations. Correct for the spatial trials factor, wavelength calibration and overlapping dispersed traces.
3. New **NIRSpec IFU** map plus same-epoch F150W imaging at the orphan and neighboring galaxies; test Hα, [O III] 4960 near **15850.8 Å**, Hβ near **15538.9 Å** under the [O III] model, and Paγ around **1.365 µm** under the Paβ/z≈0.249 scenario.
4. Compare light curves and common redshifts for a Type IIn SN, dusty TDE, nearby AGN-illuminated gas cloud, persistent metal-poor galaxy and alternative pipeline artifacts. **No physical class gets promoted without a decisive independent diagnostic.**
