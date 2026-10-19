# JANUS Cosmos V5 — real Fermi-LAT calibration tables and GBM TTE pulse checks

**Status:** instrument-provenance and pulse-diagnostic **PASS**; joint observational LIV hypothesis test **BLOCKED**. Not a detection of discrete spacetime, AMPS firewall, or Black Hole Boom.

## Scope and lineage

V5 extends the frozen V4 three-burst selection, without altering V1–V4. GRB 090510 (`bn090510016`, z=0.903), GRB 090902B (`bn090902462`, z=1.822), and GRB 130427A (`bn130427324`, z=0.34) retain 175 real V4 LAT diagnostic SOURCE-region candidates. These are *not* 175 proven source photons. The V4 event FITS and spacecraft/GTI filter remains the provenance for each LAT candidate.

New inputs are nine `P8R3_SOURCE_V3` LAT instrument response FITS from the official Fermitools 2.5.8 distribution (three response kinds × `FB`, `PSF`, `EDISP` partitions) and, for each GRB, two NaI and one BGO detector's real TTE, CSPEC and time-dependent RSP2: **27 GBM FITS files**. Detector choice is a fixed pilot cohort, not an optimized GBM spectral selection.

| GRB | NaI + BGO | GBM exact `TRIGTIME`, mission MET seconds | Prompt window relative to trigger |
| --- | --- | ---: | --- |
| 090510 | n6, n7, b1 | 263607781.971090 | 0–2 s, 32 ms bins |
| 090902B | n0, n1, b0 | 273582310.312714 | 0–40 s, 256 ms bins |
| 130427A | n6, n7, b1 | 388741629.420158 | 0–140 s, 256 ms bins |

Inputs have per-file SHA-256, size, original archive URL and package provenance in `data/cosmos/COSMOS-ALI-V5-FERMI-LAT-GBM-CALDB-MANIFEST-2026-10-19.json`. Source data are **not** bundled into Git; they total approximately 77 MB for GBM, while the Fermitools distribution package is ~252 MB. The nine CALDB FITS can be extracted from that exact pinned package. The old HEASARC CALDB endpoint can return `P8R3_SOURCE_V2`; it is *not* substituted for `V3`.

## Real checks implemented

`workspace/janus_cosmos_fermi_calibration_v5.py` validates all input SHA-256 hashes and sizes before reading, exact selected detector identities, three-way `TRIGTIME` consistency, 128 GBM channel edges, FITS data extents and EVENT PHA/time ordering, relative-time GBM GTIs, preset prompt and off-source intervals, and presence/shape of LAT Aeff/PSF/EDISP grids. It also verifies GBM RSP2 time interval order and counts response versions.

Some *actual* products use different EBOUNDS channel calibrations. For GRB 130427A, CSPEC `v01` and TTE `v00` differ by up to ~2.15% in a selected detector's channel edge; other RSP2–TTE discrepancies are retained for each detector. In GRB 090902B RSP2, one terminal extension has zero duration. The audit records it but does not count it as usable response exposure. These conditions are **not** quietly corrected through channel relabelling.

GBM TTE event `TIME` and `GTI` are in **seconds relative to trigger**, while FITS primary `TSTART`/`TSTOP` and RSP2 response spans are mission MET seconds. Mixing them creates an empty or incorrect prompt histogram. The replay asserts the time-scale distinction.

Reported light curves use NaI channel-edge overlap with 25–100 and 100–300 keV, or BGO overlap with 300–3000 keV. Count bands are channel-proxy diagnostics, not incident-photon spectral fits. A stationary pre-burst [-20,-5] s count rate is a background *illustration*, not a modeled time-dependent background. PSF, energy dispersion and photon source association are **not** evaluated. Nearest-bin FRONT/BACK effective area samples (m²) ignore incidence azimuth, livetime, energy migration, and celestial exposure and therefore are **not** proper LAT event weights or fluxes.

The 9 selected GBM TTE files contain **6,962,621 individual events**. Counts and per-bin time series reside in `data/cosmos/COSMOS-ALI-V5-CALDB-GBM-LAT-INSTRUMENT-PULSE-2026-10-19.json`.

## Joint propagation/intrinsic emission null gate

For a toy leading-order propagation coefficient (K_j(E,z_j)), the burst-specific timing model is

\[
 t_{ji} = \mu_j + [\eta + \beta_j]K_j(E_{ji},z_j) + \epsilon_{ji}.
\]

For unrestricted source-emission slopes \(\beta_j\), the transformation \(\eta\to\eta+\delta\), \(\beta_j\to\beta_j-\delta\) leaves predictions unchanged *exactly*. The real GBM pulse peaks by themselves do not establish an energy-resolved **GeV** emission-time prior that would break this degeneracy. The V5 toy suite compares *explicitly hypothetical* intrinsic-source Gaussian priors with widths 0.1, 0.5 and 2.0 in its toy units, and source-only mimic injections. This sensitivity study is **not** an observed test, significance or bound.

`workspace/janus_cosmos_fermi_joint_v5.py` never submits real LAT events or GBM counts to its synthetic Gaussian fit. Its outputs are frozen separately in `data/cosmos/COSMOS-ALI-V5-SOURCE-LAG-TOY-SENSITIVITY-2026-10-19.json`.

## Reproduction

1. Obtain the exact official Fermitools package at the manifest's `lat_caldb.package_url`; verify the *entire package* SHA-256; extract `share/fermitools/data/caldb/data/glast/lat/bcf/{ea,psf,edisp}/*P8R3_SOURCE_V3*` under a private directory `CALDB_ROOT`. Do not use `V2` as a replacement.
2. Download the manifest's 27 GBM `origin_url` files into `GBM_ROOT` under their exact file names, checking the declared SHA-256 and byte count for each. The repository deliberately doesn't fetch large archival sources during CI.
3. From the repository root, run:

```bash
python3 workspace/janus_cosmos_fermi_calibration_v5.py \
  --caldb-dir CALDB_ROOT \
  --gbm-dir GBM_ROOT \
  --output data/cosmos/COSMOS-ALI-V5-CALDB-GBM-LAT-INSTRUMENT-PULSE-2026-10-19.json

python3 workspace/janus_cosmos_fermi_joint_v5.py \
  --instrument-audit data/cosmos/COSMOS-ALI-V5-CALDB-GBM-LAT-INSTRUMENT-PULSE-2026-10-19.json \
  --output data/cosmos/COSMOS-ALI-V5-SOURCE-LAG-TOY-SENSITIVITY-2026-10-19.json

python3 -m unittest discover -s tests -q
```

The shipped JSON receipts and offline integrity/adversarial tests require **only the Python standard library**. The full instrument replay additionally requires the byte-pinned external FITS. No network access is needed for offline CI.

## Boundaries and next physically testable step (V6)

A genuine observed LIV test additionally requires: validated event-type-resolved LAT PSF and energy-dispersion evaluation (plus effective exposure, incidence/livetime corrections, diffuse/transient background and event association), GBM time-variable background with properly regridded response matrices, time-dependent source spectral evolution, and an independently justified hierarchical GeV–keV/MeV source-intrinsic lag model. It should then predeclare null calibrations, injection/coverage checks, a burst selection/selection-function analysis, and a trials budget **before** reporting any observation-based statistics.

There remains no observed evidence here for a superluminal photon, Lorentz violation, granular spacetime, or a quantum black-hole firewall. The AMPS firewall is a theoretical debate and BHB Boom remains an unconfirmed project hypothesis. Gamma-ray bursts far from a black-hole horizon do not test horizon firewall emission directly.

### Primary instrument references

- [NASA LAT IRF overview](https://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_LAT_IRFs/IRF_overview.html)
- [NASA LAT effective area](https://fermi.gsfc.nasa.gov/ssc/data/analysis/lat_irfs/irf_effective_area.html)
- [NASA GBM data and response documentation](https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm/DOCUMENTATION.html)
- [NASA GBM GRB analysis tutorial](https://fermi.gsfc.nasa.gov/ssc/data/p7v6/analysis/scitools/gbm_grb_analysis.html)
- [Official NASA HEASARC GBM trigger archive](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/)
