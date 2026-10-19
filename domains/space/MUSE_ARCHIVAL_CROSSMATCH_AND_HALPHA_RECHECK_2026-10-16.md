# JWST MOCSS 16004 Å: independent MUSE/Chandra crossmatch and Hα spectral slice

**Status (2026-10-16): physical source, line identification, redshift, distance,
and a directed gas trail remain UNCONFIRMED.** This recheck adds falsifiable
observational constraints, not a claimed physical discovery. In particular, a
detector dispersion direction is not a sky-plane propagation direction.

## Question and sky position

- Orphan line reported in Huberty et al. (2026, ApJ 1009, 117):
  RA **53.13127°**, Dec **−27.78908°**; observed vacuum **16004 ± 1 Å**.
- Four NIRSpec MOCSS stacks August/November 2023: combined flux
  **(2.5 ± 0.3) × 10⁻¹⁸ erg s⁻¹ cm⁻²**. These are not new data here.
- Previous JANUS reprocessing used independent WCS across 16 MOCSS exposures
  but found no demonstrable coherent line extension between the orphan and
  any of the three near neighbors. Its high-pass RATE contrasts are *not*
  calibrated line-surface-brightness maps.

## Direct MUSE HUDF DR2 cone query (new)

Catalog: [Bacon et al. 2023, MUSE UDF DR2](https://vizier.cds.unistra.fr/viz-bin/cat/J/A+A/670/A4),
table `J/A+A/670/A4/dr2-main`, 3″ cone around the published orphan. Query:

`https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J%2FA%2BA%2F670%2FA4%2Fdr2-main&-c=53.13127+-27.78908&-c.rs=3&-out=ID%2C_r%2CDataSet%2CDepth%2CFrom%2Czconf%2CMconf%2CIflag%2Czsys%2Ce_zsys%2Cr_z%2Cz%2Ce_z%2CRAJ2000%2CDEJ2000%2CCENTER%2CIN-HST%2CSNRMax%2CLyaEmiFlux%2CLyaEmiSNR`

Exactly one catalogued MUSE emitter within 3″:

| Field | Catalogued value |
| --- | --- |
| MUSE ID | **7974** |
| Catalog RA, Dec | 53.1318°, −27.7893° (catalog only prints 4 decimals) |
| Separation/bearing from orphan | **1.87″, 115.1° east of north** |
| Identified emission | Lyα, proposed **z = 5.980558 ± 0.001299** |
| Significance | **SNR = 3.592** |
| Redshift-confidence flag | **1 out of 3 (LOW)** |
| Matching-confidence flag | 2 out of 3 |
| MUSE dataset | MOSAIC, nominal depth 10 h, pre-JWST epoch |

This is an **independent candidate in the environment**, not a measured redshift
of the JWST orphan, nor a validated match to the nearby eastern faint galaxy
(which has photo-z≈6.03). Projecting the orphan's λ16004 into this candidate
redshift yields rest **2292.7 Å** without demonstrating a consistent line
identification. There is **no MUSE catalogue object centered within 1″** of
the orphan, and no robust MUSE z≈0.248 near-neighbor spectroscopy from this cone.

## Independent Chandra 7 Ms cone query (new)

Catalog: [Luo et al. 2017, Chandra Deep Field-South 7 Ms](https://vizier.cds.unistra.fr/viz-bin/cat/J/ApJS/228/2),
`J/ApJS/228/2/sources`, cone 10″. Exactly two catalogued X-ray AGNs:

| AGN sky position | Offset from orphan | Catalog z | Direction from orphan |
| --- | ---: | ---: | ---: |
| RA 53.130646°, Dec −27.790217° | **4.55″** | secure spectroscopic **0.666** | **205.9°** |
| RA 53.130445°, Dec −27.791120° | **7.80″** | adopted photo-z **3.240** | **199.7°** |

Neither AGN coincides with the line, and neither redshift matches the Paβ
z≈0.248 or [O III] z≈2.196 hypothesis. A directed illumination cone from
either AGN to the orphan is **not detected**. Angular proximity is not a 3-D
physical association; two-dimensional bearings are not gas velocities.

## Actual spatial/wavelength slice from the MUSE DR2 FITS cube (new)

Official public file:
`https://amused.univ-lyon1.fr/data/UDF/HUDF/download/DATACUBE_UDF-MOSAIC.fits`.
MUSE is **air wavelength** (`CTYPE3=AWAV`), with 4750 Å + 1.25 Å × channel.
The FITS primary header and first SCI extension show `NAXIS1=945`,
`NAXIS2=947`, `NAXIS3=3681`, float32, rotated TAN WCS. We directly fetched
HTTP `206 multipart/byteranges` for spectral channels **2740–2767**, only
spatial FITS rows **108–162** of each 945-pixel-wide plane, totaling
**5,824,700 bytes on wire**. This avoids the approximately 26-GB full cube.
The source falls at 0-indexed MUSE pixel (**803.28**, **134.34**).

SHA-256 of downloaded 64-KB FITS header:
`2ad7cce7b475cfb991eed30d8952bdec015245aedef1d2b90fb2ba161ec18d74`.
SHA-256 of the 28-plane HTTP multipart body:
`4c269e53562e4a7d851b91e21d757cf83cde5ebe11ea43d8a23a15e980de2424`.

If the orphan's λ16004 Å is H I Paβ (rest vacuum 12821.6 Å), then
`z = 16004/12821.6 − 1 = 0.24821`. Rest **air Hα 6562.8 Å** should appear
near **8191.73 Å** in the MUSE cube, assuming common gas kinematics. It is
inside our 8175–8208.75 Å slice.

We applied an **uncorrected 0.6″-radius aperture** and local sky annulus
0.95–1.65″, integrating an approximately 6.25 Å line window at 8191.73 Å.
No significant *spatially localized* Hα excess stands out relative to 17
valid nearby blank apertures. Indicative numbers:

| Aperture | Locally background-subtracted window integral (erg/s/cm²) |
| --- | ---: |
| Orphan | **+3.37 × 10⁻¹⁹** |
| NW candidate | +3.61 × 10⁻¹⁹ |
| Southern candidate | −2.95 × 10⁻¹⁹ |
| Eastern candidate | +0.43 × 10⁻¹⁹ |
| Median, 17 blank controls | +1.48 × 10⁻¹⁹ |
| Control MAD-to-sigma estimate | **2.05 × 10⁻¹⁹** |

The orphan is only **+0.92 local control-MAD units** above the sky-control
median; this is **not a formal significance**. The apparent uncorrected
integral contains alternating emission/absorption due to local sky structure;
there is **no coherent narrow spike centered exactly at 8191.73 Å**.

Under standard unreddened hydrogen case-B recombination conditions, a Paβ
flux of 2.5×10⁻¹⁸ in 2023 would imply Hα/Paβ of order **17–18**, or
approximately **4.4 × 10⁻¹⁷ erg/s/cm²** in a persistent *same-state* object.
That greatly exceeds the 2014–2015 local MUSE diagnostic. Thus we find no
evidence for an **unobscured, long-lived, steady Paβ gas source at identical
flux extending back to 2014–2015**. This does **not** eliminate a new 2023
transient, recently illuminated nebula, strong differential extinction, a
spatially offset extended source, or systematic aperture/spectral noise.

This is intentionally not a rigorously calibrated Hα flux limit: we have not
reprocessed the MUSE variance/exposure cube, correlated sky noise, PSF
encircled energy, or completeness. The MUSE data preceded the 2023 JWST
event by ~8 years. Avoid assigning a sigma or a Bayes factor for source
class from this look.

Offline replay: `tools/check_muse_dr2_halpha_slice.py` plus companion
`domains/space/MUSE_DR2_HALPHA_SLICE_RECEIPT_2026-10-16.json`. The latter
contains channel-by-channel sky-subtracted spectra and blank controls.

## Distance is conditional on identifying the rest-frame line

For transparent comparison we use a **flat ΛCDM** model with `H0=70
km/s/Mpc`, `Ωm=0.3`, `ΩΛ=0.7`. These are **model alternatives, not
measurements**:

| Proposed line ID | z calculated from vacuum λ | Luminosity distance | kpc per arcsecond |
| --- | ---: | ---: | ---: |
| H I Paβ | 0.24821 | **1.25 Gpc** | 3.89 |
| Hα | 1.43792 | **10.36 Gpc** | 8.45 |
| [O III] 5008 | 2.19553 | **17.42 Gpc** | 8.27 |
| Lyα | 12.16476 | **129.61 Gpc** | 3.63 |

The luminosity distance can be greater than the current age of the universe
expressed in light years; it is **not** a photon path length or direct
light-travel time. Without a securely identified rest wavelength, no
unique distance, physical extent, or luminosity can be inferred.

## Source/direction conclusion

1. **Physical source:** not identified. The data confirm a sky-projected
   emission-line signal with nearby galaxies and unrelated or not-yet-related
   spectroscopic/AGN sources. Paβ/gas-cloud and [O III]/faint-galaxy
   scenarios remain different models. Claiming SN IIn, TDE, AGN echo, or
   an instrument artifact as the **established** source would be fabrication.
2. **Distance:** no measured source redshift. Values above are explicitly
   conditional alternatives. The nearby MUSE z≈5.98 is low-confidence **for
   another sky position**; the X-ray AGNs have distinct redshifts.
3. **Directed gas trail:** **not detected** in the existing approximate
   MOCSS sky-rate 3″ field. No position angle, velocity field, lower/upper
   physical trail length, or 3-D trajectory has been measured.

The problem cannot be honestly declared physically resolved from one line
without calibrated second-line identification or contemporary spatial IFU
data. The only honest closed result is that the proposed source/distance/trail
association remains **scientifically non-identifiable** under these
observations, with several hypotheses still consistent with the evidence.

## Primary literature and official provenance

- [Huberty et al. 2026, original MOCSS report](https://arxiv.org/html/2610.12553v1).
- [Bacon et al. 2023, MUSE HUDF DR2](https://arxiv.org/abs/2211.08493),
  [DR2 catalogs](https://vizier.cds.unistra.fr/viz-bin/cat/J/A+A/670/A4),
  [AMUSED full cube](https://amused.univ-lyon1.fr/project/UDF/download).
- [Luo et al. 2017, 7 Ms CDF-S](https://vizier.cds.unistra.fr/viz-bin/cat/J/ApJS/228/2).
- [Hummer & Storey 1987, hydrogen recombination case B](https://academic.oup.com/mnras/article/224/3/801/980922).

