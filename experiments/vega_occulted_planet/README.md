# Vega Occulted / Unresolved Inner Planet Arm

Status: **HYPOTHESIS ARM · NOT A DISCOVERY CLAIM**

This branch tests a narrow question: could a low/moderate-mass planet in the inner Vega system remain unresolved or projected close to the star while still being compatible with the unusually smooth debris disk?

## Scientific motivation

Published JWST/MIRI imaging finds Vega's debris disk to be remarkably smooth and axisymmetric, with a broad outer belt at roughly 80–170 au, a shallow dip/gap near 60 au, and an inferred warm-dust inner edge near roughly 3–5 au. The same work argues against Saturn-mass planets outside about 10 au under the embedded-perturber scenario, while noting that a modest-mass / Neptune-size planet could plausibly shepherd the inner edge of the warm debris.

JWST/NIRCam coronagraphy reaches model-dependent limits of about <3 Mjup at 1 arcsec (~7.7 au), <2 Mjup at 2 arcsec (~15 au), and ~0.5 Mjup beyond 5 arcsec (~38 au) for the adopted system age. Those limits do not exclude terrestrial/Neptune-mass planets in the inner several au.

Primary references:

- Su et al. 2024, ApJ 977, 277, DOI `10.3847/1538-4357/ad8cde`.
- Beichman et al. 2025, AJ 169, 17, DOI `10.3847/1538-3881/ad890d`.

## Frozen hypothesis

`H1_OCCULTED_UNRESOLVED_INNER_PLANET`

A planet with an orbit predominantly inside the direct-imaging sensitivity regime may exist without producing a large, easily visible perturbation in the outer debris disk.

The experiment does **not** treat alignment with the stellar PSF/saturation zone as evidence for a planet. It only asks whether a region of mass–semimajor-axis parameter space remains observationally admissible.

## Gate order

```text
PUBLISHED CONSTRAINTS
  -> generate mass/orbit grid
  -> reject direct-imaging-excluded regions
  -> reject disk-dynamics-excluded regions
  -> retain inner low-mass admissible cells
  -> emit TOPA hypothesis queue
  -> emit Spider acquisition targets
  -> freeze receipt
```

No OCR, face search, cipher search, semantic pattern matching, or post-hoc threshold tuning is part of this arm.

## TOPA role

TOPA receives explicit competing hypotheses rather than a preferred story:

- H0: no planet is required; dust transport is drag-dominated.
- H1: unresolved inner low/moderate-mass planet.
- H2: planet associated with the ~60 au dip.
- H3: reduction / PSF / coronagraph systematics.
- H4: alternative disk dynamics without a planet.

TOPA is allowed to rank *tests*, not to convert missing data into evidence.

## Spider role

The Spider queue requests data that could actually discriminate H0/H1:

- public JWST/NIRCam F444W coronagraphic products and contrast curves;
- JWST/MIRI F1550C/F2300C/F2550W disk products;
- HST/STIS coronagraphic products and reference-PSF metadata;
- ALMA 1.34 mm outer-belt products;
- published or archival radial-velocity time series where available;
- Gaia / Hipparcos astrometric acceleration constraints where scientifically valid;
- multi-epoch direct-imaging astrometry for reported point-like/extended sources.

## Claim ceiling

A surviving grid cell means only:

> `NOT EXCLUDED BY THE ENCODED CONSTRAINTS`

It does not mean a planet has been detected.
