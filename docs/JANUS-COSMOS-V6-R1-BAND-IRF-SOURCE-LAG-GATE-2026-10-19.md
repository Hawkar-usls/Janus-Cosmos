# JANUS Cosmos V6-R1: real GBM Band+PL pilot, LAT CALDB audit, source-lag falsification

Date (UTC): 2026-10-19. Status: **PARTIAL SCIENCE GATE** — no full LAT/GBM
likelihood, no observed LIV p-value, no horizon/firewall/BHB Boom detection.

This step follows the frozen [V6-R0 protocol](JANUS-COSMOS-LAT-GBM-JOINT-V6-R0-2026-10-19.md).
The V5 LAT/GBM manifest remains the only authority for the 42 raw NASA
FITS/CALDB sources: 27 GBM TTE/CSPEC/RSP2, 9 P8R3_SOURCE_V3 IRF tables,
and six LAT EVENT/SC files. Source binaries are rehydrated with strict SHA256
checks and are **not** committed to GitHub.

## Real GBM spectral model

The measured counts in detector channel `c` over time interval `t` are modeled
using each detector's actual OGIP RSP2 matrix, time-overlap weighted for
attitude-dependent matrices:

`S_c(t) = duration(t) * A * integral[dE R_c(E,t) * (Band(E;alpha,beta,Epeak) + q*Band(1000 keV)*(E/1000 keV)^(-1.8))]`.

`A` is one shared photon normalization for the two NaI and one BGO detector.
`q` is the ratio of the extra PL to Band at 1 MeV. A **fixed** grid is searched
over α = (-1.4,-1.0,-0.6), β = (-3,-2.4), Epeak = (150,400,1000,2500) keV,
q = (0,0.5,2,8); there are 96 shapes per segment. Band-only (`q=0`) is nested
in the *same* grid, but the model selection remains descriptive; 96-grid search
statistics have not been calibrated by null simulations.

The actual response is evaluated using 4-subinterval Simpson true-energy
quadrature rather than V6-R0 midpoint integration. Source spectral shapes and
normalizations can evolve independently across 17 preregistered coarse windows.
The windows are **not statistically independent**, since they share off-source
counts, backgrounds, and instrument nuisance parameters.

### Time-variable background and independent checks

For each GBM detector, two **pre-trigger** off-source windows train a common
log-linear rate slope; a third, later pre-trigger window is held out. All
windows end at least 5 s before the instrument trigger. Counts in each
detector channel are handled with the exact Poisson ON/OFF profile, using the
time-integrated exponential slope to convert background exposure between OFF
and ON windows. Slope uncertainty, shared-OFF covariance and possible
precursor contamination are **not** integrated into a global likelihood.

For each of the 17 intervals we also independently fit **NaI-only** models
and score the BGO on-source counts as a held-out detector. This is a real
falsifiable cross-detector stress test, not a calibrated goodness-of-fit p.
The detector EBOUNDS mismatch and relative detector normalization have not
been formally reconciled. Extrapolating a two-window pre-trigger drift over
long bright bursts is particularly weak; an extrapolation warning is recorded.

### Frozen results: descriptions, not detections

| GRB and relative interval | Three-detector Band→Band+PL deviance improvement | NaI-trained held-out BGO deviance | Interpretation |
| --- | ---: | --- | --- |
| 090510, all six 0–2 s windows | 0 each | identical | No preference for PL on this coarse grid |
| 090902B, 20–30 s | 109.476 | 628.120→460.797 | Improvement in both fit and BGO holdout; neither is calibrated significance |
| 130427A, 20–40 s | 250.290 | 127.873→130.236 | Apparent fit improvement fails this independent BGO check |
| 130427A, 40–70 s | 69.824 | 201.853→139.119 | Partial cross-detector improvement; long-background extrapolation warning |

The early, bright 130427A BGO window 0–10 s has a heldout deviation of
`13043.242`. The source/background/response model is plainly insufficient
there: do not use this as a valid astrophysical parameter estimator. The
additional nonthermal spectral component in GRB 090902B is already known
in published Fermi work; this reproduction does **not** discover it.

Primary reference: Abdo et al., *Fermi Observations of GRB 090902B: A Distinct
Spectral Component in the Prompt and Delayed Emission*,
https://arxiv.org/abs/0909.2470 . Official response/spectral workflow:
https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm/gbm_data_tools/gdt-docs/notebooks/SpectralFitting.html .

## LAT IRF physical model readiness: honest BLOCK

We rehydrated and parsed all **nine** NASA Pass8 P8R3_SOURCE_V3 FITS table
families for effective area, PSF and energy dispersion in the FRONT/BACK,
PSF0–3 and EDISP0–3 partitions (30 instrument-type tables in total).

- Validate PSF Moffat core/tail parameters, positivity and energy scaling
  from `PSF_SCALING_PARAMS`; report their actual 1 GeV/30° values.
- Validate EDISP mixture fractions and asymmetric scale/skew/bias shape
  coefficients, as well as `EDISP_SCALING_PARAMS` for every event type.
- Validate azimuthal dependence and livetime-efficiency *parameter tables*;
  no approximate correction is passed off as an applied efficiency.

This is a **real CALDB parameter audit**. We intentionally do **not** claim
an event-level normalized PSF/EDISP source likelihood: we still need full
IRF interpolation, angular and energy forward convolution, satellite
azimuth/efficiency live-time weighting, diffuse Galactic/isotropic templates,
unbinned likelihood, and benchmark comparison with official Fermitools.

NASA instrument documentation:

- https://fermi.gsfc.nasa.gov/ssc/data/analysis/lat_irfs/irf_overview.html
- https://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_LAT_IRFs/IRF_PSF.html
- https://fermi.gsfc.nasa.gov/ssc/data/analysis/lat_irfs/irf_energy_dispersion.html
- https://fermi.gsfc.nasa.gov/ssc/data/access/lat/BackgroundModels.html

## Source-intrinsic delay and synthetic injection coverage

For burst `j`, the illustrative delay observable obeys

`observed_delay_j = K_j * eta + beta_j + Gaussian_measurement_error_j`.

Without an **independently measured** population prior for intrinsic
`beta_j`, the propagation coefficient `eta` is exactly unidentifiable:
for any candidate `eta'`, `beta'_j=beta_j+(eta-eta')*K_j` leaves all
predicted delays identical. No quantity inferred from the same observed
delay sample can independently establish that prior.

The synthetic-only demonstration uses deterministic PRNG seed 20261019,
2000 replicates per scenario and arbitrary dimensionless propagation
weights, **not cosmological distances or Fermi observations**:

| Synthetic scenario | Naive 95% CI coverage, assuming every intrinsic lag is zero |
| --- | ---: |
| Zero delay, correct source assumption | 94.7% |
| Propagation injection, correct source assumption | 95.35% |
| No propagation, intrinsic delay designed to mimic it | **0%** |

This is a severe negative control for unjustified observed LIV limits and
false p-values. It is **not** instrument-level Poisson injection coverage.
See Vasileiou et al., *Constraints on Lorentz invariance violation from
Fermi-LAT observations of gamma-ray bursts*, Phys. Rev. D 87, 122001
(https://doi.org/10.1103/PhysRevD.87.122001) for the importance of
conservative treatment of intrinsic spectral evolution.

## Reproduction

Use Python 3.11+, `numpy`, `scipy` for the GBM **real-data optimization**.
Test suite, pure analytic Band/source-negative-controls, and frozen output
integrity checks do not need SciPy or NumPy. Strict download manifests are
checked before decoding. From the JANUS Cosmos repository root:

```bash
python3 scripts/fetch_fermi_v6_inputs.py \
  --gbm-dir /tmp/janus_v6r1_gbm --caldb-dir /tmp/janus_v6r1_caldb \
  --package /tmp/janus_v6r1_fermitools.conda
python3 scripts/fetch_fermi_lat_v4.py --raw-dir /tmp/janus_v6r1_lat
PYTHONPATH=workspace python3 workspace/janus_cosmos_fermi_spectral_v6r1.py \
  --gbm-dir /tmp/janus_v6r1_gbm
PYTHONPATH=workspace python3 workspace/janus_cosmos_fermi_lat_irf_v6r1.py \
  --caldb-dir /tmp/janus_v6r1_caldb
PYTHONPATH=workspace python3 workspace/janus_cosmos_fermi_source_lag_v6r1.py
python3 -m unittest discover -s tests -q
```

Three frozen scientific outputs:

1. `data/cosmos/COSMOS-ALI-V6R1-GBM-BAND-POWERLAW-DRIFT-2026-10-19.json`
2. `data/cosmos/COSMOS-ALI-V6R1-LAT-FULL-IRF-READINESS-2026-10-19.json`
3. `data/cosmos/COSMOS-ALI-V6R1-SOURCE-LAG-TOY-INJECTION-COVERAGE-2026-10-19.json`

## Unfinished V6-R1 gates

1. Fit realistic background evolution with sufficient OFF controls,
   shared-background covariance and quantified drift uncertainty; calibrate
   TTE/CSPEC/RSP2 EBOUNDS and cross-detector normalization systematics.
2. Implement/verify full official LAT PSF and EDISP convolution, phi and
   efficiency corrections, Galactic + isotropic diffuse models, and
   event-type likelihoods against Fermitools output.
3. Build and independently justify an intrinsic GeV emission-lag
   hierarchical prior using separate source data or conservative bounds.
4. Validate instrument-level injection/recovery under realistic sources and
   backgrounds, pre-register coverage and false-positive gates, then unblind
   a real multi-GRB LRT.

Until **all** four pass: `observed_liv_p_value=null`,
`discrete_spacetime_detected=false`, `AMPS_FIREWALL=OPEN_UNCONFIRMED`,
`BHB_BOOM=NOT_ESTABLISHED`. Black-hole horizon/fall-in radiation is not
probed directly by these three GRB emission-delay experiments.
