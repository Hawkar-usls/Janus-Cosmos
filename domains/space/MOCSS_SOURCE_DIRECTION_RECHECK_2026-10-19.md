# JANUS COSMOS / TOPA: What was emitting at 16004 A, and in which direction?

**Status 2026-10-19: independent WCS-neighbor null controls measured; source type,
redshift, true emission/illumination axis and trail still unidentified.** We do
not promote a noisy detector residual to an astrophysical source.

## What we can actually test today

The published HuDF MOCSS line is 16004+-1 A at RA 53.13127, Dec -27.78908.
The spectra were taken on August 23-24 and November 19, 2023, in four groups
A-D across two nearly perpendicular dispersion orientations. The original
paper independently identifies a real candidate. A celestial location is
not a measured direction of propagation *away from* that location.

Our upstream COSMOS run 39230522837 reconstructed twelve (4 groups x 3 line
windows) 2D cutouts from all sixteen original JWST RATE frames, SHA-checked
against the original MAST inputs. We downloaded artifact 12176970014 and
independently matched its exact 2D NPZ SHA256
`0b8974e9f761e111abce46d16e0db322f368fd67306b7e8920cfc92fefbe8190`.
Recompute from the frozen archive with
`tools/scan_mocss_four_stack_sky_azimuths.py`.

The new test uses existing independently constructed local sky-coordinate WCS
Jacobians. For each nearby published galaxy, we evaluate the *same 16004-A*
line-template response at that separate sky position in every A/B/C/D stack,
and compare it to 72 matched-radius azimuthal directions at intervals of 5 deg.
These 72 measurements are correlated and are **not** a 72-trial p-value or SNR.

| Galaxy direction from line | Angular separation | A | B | C | D |
| --- | ---: | ---: | ---: | ---: | ---: |
| NW (PA 328.2 deg) | 0.635 arcsec | -0.18 | +1.86 | -0.19 | -1.85 |
| S (PA 174.5 deg) | 1.038 arcsec | +1.24 | -0.00 | -2.34 | -0.02 |
| E (PA 95.1 deg) | 1.004 arcsec | -0.23 | +0.49 | -2.10 | +0.89 |

Entries are *signed deviations from the local azimuth-ring median in robust
MAD units*, not Gaussian sigma. **None of the three galaxy directions remains
positive in all four stacks.** Exploratory samples on all three bearings at
0.35, 0.60, 0.85 and 1.10 arcsec likewise have no stable direction-by-radius
coherence across A-D. The full per-radius numbers are frozen in
`MOCSS_FOUR_STACK_NEIGHBOR_AZIMUTH_NULLS_2026-10-19.json`.

**Crucial falsification limit:** These native RATE proxy stacks do NOT apply
the physical shutter-open fraction to *each detector pixel*. They lack MOCSS
row-sky illumination, properly dispersed PSFs, real line-flux calibration,
sky-background covariance and a trials-corrected significance threshold.
They therefore do not establish that the neighboring positions contain or
do not contain the true line, nor do they give an upper limit on trail length.
There is no evidence of a real extended beam out to 1.1 arcsec *in this test*,
but 1.1 arcsec is a test boundary, not a measured beam boundary.

## Conditional geometry: from whom, and toward where?

If the 16004 A line is Pa-beta at z=0.24856, the source luminosity distance is
approximately 1.293 Gpc and measured *flux* implies an observed Pa-beta line
luminosity near 5.0e38 erg/s. If a neighboring galaxy's AGN illuminated a
separate gas cloud at the orphan position, the assumed ionizing propagation
directions (bearings clockwise from celestial north) are:

| Hypothetical illuminating galaxy | Galaxy -> cloud PA | Projected path | At least this many years for light to reach cloud |
| --- | ---: | ---: | ---: |
| NW | 148.23 deg (SE) | 2.55 kpc | 8,328 |
| S | 354.54 deg (N) | 4.17 kpc | 13,608 |
| E | 275.14 deg (W) | 4.04 kpc | 13,165, but its zphot~6.03 makes this association improbable |

Only the NW and S neighbors have broad photometric redshift distributions
possibly including z=0.249. None has a confirmed spectroscopic association.
These are projected *minimum direct path lengths*, **not** a measured time
delay of the echo as seen by Earth, which depends on the unknown 3D depth.
Emission-line photons from the gas are not intrinsically a directed signal.

## Why would the line be isolated? A quantitative dust test

Huberty et al.'s Case-B expectation Pa-beta/Pa-gamma=1.6 would predict a
Pa-gamma observed flux of 1.56e-18 erg/s/cm^2 absent extra attenuation.
Using the *published broad 3-sigma line-limit range* 0.82-1.06e-18 as a
**hypothetical** Pa-gamma-specific limit and CCM89 Rv=3.1 at the two *rest*
wavelengths, the required differential near-IR extinction corresponds to
**Av >= 5.35-8.88 mag**. The observed Pa-beta luminosity would need dust
correction to at least 1.9-4.6e39 erg/s. The exact Pa-gamma spectral noise
has not been independently calibrated; the shocked SN gas need not obey
Case B or a uniform CCM foreground screen. Therefore **high dust is a
conditional requirement of one model, not measured extinction**.

## Why the 2023 October images matter

The *discovery paper's* principal F150W comparison was October **2022**, but
COSMOS separately found four F150W archival mosaics dated October **1-9,
2023**, **between** the August and November spectroscopic detections.
The diagnostic matched-sky covariance flux at the line position is
0.2785+-0.663 nJy in a 0.12-arcsec radius aperture. This is *not*
PSF-aperture-corrected or a final absolute-flux science measurement; it must
not be compared directly as if the published 5.5 nJy *total* line-only model
were measured in that small aperture.

If a real turn-off and reappearance is eventually demonstrated, each happens
within ~39-41 observer-frame days, or ~31-33 rest-frame days if z=0.249.
For a simple recombining H gas region at T=10^4 K to respond on <=31.2 days,
alpha_B=2.6e-13 cm^3/s implies **ne >=1.43e6 cm^-3**. Typical extended
ionized gas with ne=100 cm^-3 has recombination time ~1219 years; its
brightness cannot track a sharp one-month ionizing-source dip under this
simple model. Dense shock gas, dust variability and geometrically thin
echo surfaces have different physics. *The October dip itself is not proven.*

## Epistemic result, not a claim of solution

**Source identity: unknown. Distance: conditional on line identification.
Incoming photons: from the published HUDF direction. Ionizing source -> cloud
path: only vectors under counterfactual galaxy association. Trail length:
not measured. Why the line is seen: constrained but not identified.**

Frozen reproducibility:

- `domains/space/MOCSS_FOUR_STACK_NEIGHBOR_AZIMUTH_NULLS_2026-10-19.json`
- `domains/space/MOCSS_NEIGHBOR_ILLUMINATION_DUST_RECOMBINATION_2026-10-19.json`
- `tools/scan_mocss_four_stack_sky_azimuths.py`
- `tools/model_mocss_illumination_recombination_and_pabeta_dust.py`
- `tools/verify_mocss_direction_2026_10_19.py`

Original paper: https://arxiv.org/html/2610.12553v1
CCM dust law: https://adsabs.harvard.edu/pdf/1989ApJ...345..245C
AGN cloud geometry: https://arxiv.org/abs/1110.6921
Official MAST data: https://doi.org/10.17909/tct8-jr92
