# Real MAST Oct-2023 NIRCam source-position search: local maxima and sky-control look-elsewhere check

**2026-10-19 UTC. New independent GitHub Actions result:** [NIRCam four-epoch original-pixel spatial offset workflow 39256147894](https://github.com/Hawkar-usls/Janus-Cosmos/actions/runs/39256147894), `success`, scientific artifact ID **12186849693**; artifact ZIP SHA256 **98ed0db2e0fb10b9774e13d0f2a79ce6dc620e432c474f04f5d831a42aac502f**. An earlier attempt (run 39255670376) failed before *any* observation was analyzed because its fresh temporary MAST byte-range cache directory did not yet exist. The replay explicitly fixed the directory creation and succeeded, rather than silently treating infrastructure failure as scientific evidence.

## Measured, not simulated

Independently download and verify the four original JW01180 F150W `i2d` SCI byte-range SHA256s acquired **2023 Oct 1, 2, 3, 9**. This scan uses the *same* 0.12″ aperture, 0.24–0.38″ local median-background annulus, same four fixed epoch weights and same 95 jointly clipped matched-sky controls as the previous covariance receipt. There is **no synthetic point source injected in this experiment**.

Scan **113 offset positions**, separated by 0.075″ on a regular two-dimensional grid, within 0.45″ of RA=53.13127°, Dec=−27.78908°. Repeat the *entire maximum-over-grid search* around each of **95 matched-sky control centers**. This explicitly guards against interpreting the largest value among many searched positions as a single pre-specified measurement. The controls include static faint galaxies and spatial correlations; their empirical exceedance fractions are **not independent-trial formal p-values**.

| Search radius | Maximum target aperture brightness (nJy) | Median of 95 control-field maxima (nJy) | Controls with equal/higher maxima | Target maximum offset (east,north), arcsec |
|---|---:|---:|---:|---:|
| 0.15″ | **0.502** | 0.714 | **65/95** | (0,−0.150) |
| 0.225″ | **0.502** | 1.009 | **92/95** | (0,−0.150) |
| 0.30″ | **0.502** | 1.138 | **94/95** | (0,−0.150) |
| 0.45″ | **0.950** | 1.469 | **86/95** | (−0.225,+0.375) |

The central forced-aperture flux is ~0.371 nJy on this 95-control median convention, consistent with the prior same-frame analysis. **No unusual nearby pointlike F150W source has been detected by this search.** The field is not a clean set of independent detector noise trials, so the control counts must not be converted to detection or rejection `p` values.

## Does the weak outer maximum point to the north-west galaxy?

The nearest cataloged weak galaxy is at RA 53.131165°, Dec −27.788930°, offset approximately (−0.3344,+0.5400) arcsec from our reference; projected separation ~0.6352″ and position angle **328.23°** (east of north). The target field's maximized aperture location (−0.225,+0.375) lies at a position angle ~**329.04°**, a **0.80°** angular alignment with the nearby galaxy's direction.

However, inspecting the full stored **113-point map** finds a very similar opposite-side candidate at (+0.225,−0.375): **0.875 nJy** compared to **0.950 nJy** toward the north-west. The NW–SE difference is only **0.075 nJy**, much smaller than the ~0.663 nJy matched-sky covariance SD (which is itself not a calibrated error on this directional contrast). A separate off-axis peak toward the south occurs inside radius 0.15″. This is *not* compelling evidence for a one-directional filament. A background/neighbor brightness gradient and aperture-annulus effects are entirely plausible.

## Relation to the PSF-injection experiment

The separate COSMOS [original-pixel analytic PSF injection receipt](NIRCAM_F150W_PSF_INJECTION_OBSERVED_2026-10-19.json) showed that a correctly registered unresolved 5.5 nJy source would contribute around 4.2 nJy in these apertures. It also showed that deliberately offsetting that source 0.15″ can cause a *single fixed-coordinate* aperture to miss most of the signal. Our new spatial search covers offsets within 0.45″, so an **ordinary pointlike persistent counterpart** of the published flux is even less plausible in October if the original NIRSpec→F150W 5.5 nJy conversion and analytic image profile are reliable. Grid spacing and unrelated sky contaminants add uncertainty; we do not report a formal exclusion confidence.

Crucially, an **extended gas cloud** with FWHM ~0.35″ could yield under 1 nJy in a 0.12″ aperture at the same integrated total flux (from the *model*, not observed gas). The present small-aperture search does **not** exclude this gas-cloud scenario or set a spatially integrated F150W upper limit.

The broad F150W imaging data cannot isolate the **16004 Å line**; there is no spectral-channel spatial direction, line centroid velocity, unique host, redshift or physical filament length measured here. Next is a larger-aperture/segmentation-controlled spatial test and, most decisively, an independent *NIRSpec* two-orientation common-sky-source/line-centroid solution followed by calibrated Level-2+IFU data.

### Audit and replay

The new GitHub script `tools/scan_nircam_f150w_local_offset_candidates.py` and its [workflow](../../.github/workflows/jwst-nircam-2023-local-offset-search.yml) are published in the COSMOS Draft PR. The successful Actions artifact has been frozen as `NIRCAM_F150W_LOCAL_OFFSETS_OBSERVED_2026-10-19.json`, containing every target-grid aperture, the fixed weights, the original 4 hashes, and the radius-specific 95-control *field max* summaries. All four MAST original SCI byte subsets have independent SHA-256 pins in `NIRCAM_F150W_TIME_DOMAIN_RECEIPT_2026-10-16.json`.

**Scientific status:** `OFFSET_POINT_SOURCE_BRIGHT_COUNTERPART = NOT_FOUND_ABOVE_MATCHED_SKY_FIELD_MAX_CONTROLS`; `DIRECTED_16004A_GAS_TRAIL = NOT_DETECTED`; `COMPACT_PERSISTENT_SOURCE_TENSION = CONDITIONAL`; `EXTENDED_GAS = OPEN`; `HOST_CLASS_REDSHIFT_DISTANCE = OPEN`.
