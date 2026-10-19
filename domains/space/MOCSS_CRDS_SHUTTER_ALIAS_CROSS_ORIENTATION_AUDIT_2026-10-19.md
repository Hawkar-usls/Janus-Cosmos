# MOCSS 16004 Å — independent CRDS WCS slitless spectral-alias audit

**UTC 2026-10-19. Status: verified geometry, not confirmed spectral contamination or a new host identification.**

We inspected all four official WCS GitHub Actions artifacts from successful [COSMOS workflow run 39251950005](https://github.com/Hawkar-usls/Janus-Cosmos/actions/runs/39251950005), preserving original ZIP SHA256s. Each artifact contains four exposure-specific detector-footprint receipts, for **16 distinct SHA-identified original NIRSpec RATE exposures** and original MSA OPEN-state masks. The geometry used **JWST 1.19.1, CRDS `jwst_1413.pmap`**. The WCS is calculated on a *synthetic 9-slit background-source probe* with physically real 3×3 MSA shutter coordinates; it is **not** full detector flux calibration.

## A measured geometric ambiguity (not a measured second line)

Compare the virtual WCS detector centers of three adjacent, actually OPEN MSA rows for a hypothetical line. Evaluate NIRSpec CRDS at 15850.7933 Å and 16004 Å; use the local wavelength-dispersion slope to estimate which neighbor-row photon wavelength would map to the detector X of the central 16004 Å line.

| Slit model | Source line wavelength producing central 16004 Å pixel (first-order WCS approximation) | Across 16 exposures |
|---|---:|---:|
| Source at central row | 16004.00 Å | fixed reference |
| Source at row −1 | **15987.38 Å** | 15987.345–15987.442 Å |
| Source at row +1 | **16020.62 Å** | 16020.557–16020.655 Å |

The two virtual adjacent-row positions differ from the central trace by ≈2.6 detector X pixels at the *same* 16004 Å wavelength. Shifting the **hypothetical** neighbor wavelength by ≈16.6 Å realigns their detector X coordinates; the linear-model Y residual is below **0.016 detector pixel** across all 32 neighbor tests. This demonstrates a finite *single-orientation* NIRSpec slitless sky-position/wavelength degeneracy; it does **not** show photons at 15987.38 or 16020.62 Å.

The original CRDS probe artifacts additionally report **48/48 finite on-detector WCS centers for each of the three line windows** — 13656.748 Å (Paγ under the Paβ identification), 15850.793 Å (the corrected [OIII] companion) and 16004 Å (main line), for three sample shutters × 16 exposures. This establishes only that the *geometric sampled centers* are on detector; it does **not** establish throughput, noise sensitivity, a detected companion line, or the absence of contamination at those locations.

## Essential two-orientation correction: the original authors already knew this

The original MOCSS [paper Section III.2](https://arxiv.org/html/2610.12553v1) explicitly describes that a detector pixel can mean different wavelengths for sources at different sky positions. This motivated their two position angles separated by approximately **89°**. They enumerated candidate host shutters in both orientations and required **a common sky host plus wavelength agreement within 6 Å**. The paper's 16004 ±1 Å sky origin was reported after this cross-orientation matching and 2D-spectrum WCS verification.

Our independent first-order same-signed alias wavelengths agree **by wavelength alone** between orientations, within about 0.08 Å:

| Hypothetical row alias in both orientations | First orientation (A/B), median | Second orientation (C/D), median |
|---|---:|---:|
| row −1 | 15987.347 Å | 15987.427 Å |
| row +1 | 16020.653 Å | 16020.571 Å |

**This is insufficient to propose a new source.** These rows refer to **different sky directions after an ~89° rotation**: we have not established that any single actual astronomical object lies in the relevant shutters at *both* orientations. The authors' common-sky-host condition is an additional nontrivial filter. Our 16-WCS geometry audit has **not replicated** this complete two-orientation condition; it does **not falsify** their published localization.

STScI explicitly documents NIRSpec MOS overlapping/dispersed-source contamination and source-centering constraints. However, an instrumental *wavelength alias* is not a measured Doppler shift or a directed gas flow.

## Consequences for physical source, distance, direction and why

1. It remains possible in principle that light from some displaced source contributes to a measured detector spectral pixel; identifying one requires a common **sky** origin across both near-orthogonal orientations, not merely a matching apparent detector line.
2. The three nearby faint HUDF galaxies are **not** assigned to these ±1-row WCS traces by this experiment. Without projecting their real RA/Dec through all 16 true MSA masks one must not identify any of them as the emitter.
3. The NIRCam October 2023 lack of a compact persistent counterpart at the reported source position remains an independently measured, **conditional** tension. A different sky origin could in principle resolve it, but this audit has **not found** such an origin.
4. Paβ at z≈0.249, an AGN-illuminated cloud, a Type IIn/TDE transient and spectral/background contamination all retain untested alternatives. A true redshift, luminosity distance, source identity, propagation angle and physical trail length are **not determined** here.

## Reproducibility

```bash
python3 tools/analyze_mocss_crds_adjacent_shutter_aliases.py \
  --zip-dir /tmp --output /tmp/mocss_alias_check.json
```

The four exact ZIP filenames required in `/tmp` are `janus_mocss_neighbor_a.zip` through `_d.zip`; archive artifact IDs and SHA256 checksums are hard-coded in the standalone script and recorded in `MOCSS_CRDS_16RATE_ADJACENT_SHUTTER_ALIASES_2026-10-19.json`. Rerunning the script fails on any changed artifact. The published machine-readable receipt includes both orientation-specific wavelength medians and the explicit *sky host not verified* gate.

**Next decisive physical test:** independently repeat the *paper's* full two-orientation common-sky-host and ≤6 Å selection using every original MSA OPEN shutter, the actual sky coordinates of all neighboring galaxies, and the original line centroids. Only survivors should go to physical STPSF cross-contamination fitting on real calibrated detector SCI/ERR/DQ images.

Primary sources:

- Huberty et al., [A Single Line-Emitter Identified in an Unbiased Search with JWST/NIRSpec](https://arxiv.org/html/2610.12553v1), in particular Section III.2–III.3; original G140M/F100LP MOCSS intentionally uses almost all shutters open.
- STScI, [NIRSpec Multi-Object Spectroscopy](https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-observing-modes/nirspec-multi-object-spectroscopy) and [MOS Recommended Strategies](https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-observing-strategies/nirspec-mos-recommended-strategies).
