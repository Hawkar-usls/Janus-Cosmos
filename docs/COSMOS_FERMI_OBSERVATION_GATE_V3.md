# JANUS Cosmos V3 — Fermi-LAT real published photons, source-delay integrity gate

**Experiment:** `JANUS-COSMOS-ALI-V3-FERMI-PHOTON-TOF-2026-10-18`  
**Status:** verified deterministic analysis of **published astronomical table values**, without full FITS photon event data or source/instrument likelihood. The selected discrete-space dispersion **has not** been supported or excluded by this experiment.

## Why this is a new scientific stage

V1 audited the video-inspired claims; V2 created twelve completely synthetic photon time-of-flight benchmarks for one specified nearest-neighbor lattice law. V3 uses **real, published Fermi-LAT photon energy and arrival-time numbers** while explicitly refusing an invalid Lorentz-violation conclusion.

Source: Ajello et al. / Fermi-LAT Collaboration, *A Decade of Gamma-Ray Bursts Observed by Fermi-LAT: The Second GRB Catalog*, [arXiv:1906.11403v1, Table 7](https://arxiv.org/html/1906.11403v1). The primary table lists gamma-ray photons above 10 GeV with burst-association probability above 90%, their energies in GeV, their arrival times in seconds **after the Fermi GBM trigger**, and redshift where known. The 11 known-redshift bursts in that table are recorded in `data/cosmos/COSMOS-ALI-V3-FERMI-PUBLISHED-HE-PHOTONS-2026-10-18.json` with published numeric precision, not fabricated extra digits.

## Evidence levels that MUST remain distinct

| Input or inference | Actual status | Admissible use |
| --- | --- | --- |
| 38 high-energy photon energies, times after GBM trigger; 11 bursts with redshift | `REAL_PUBLISHED_ROUNDED_CATALOG_VALUES` | Conditional prediction overlay and basic timing-order diagnostic |
| Raw Fermi-LAT extended event FITS, spacecraft orientation, IRF/PSF, individual probabilities | `NOT_ACQUIRED` | **Not** claimed as imported or analyzed |
| Actual low-energy matched photons at a common emission epoch | `NOT_AVAILABLE_IN_SELECTED_TABLE` | 1 GeV reference is purely hypothetical |
| Intrinsic energy-dependent source-lag model, detector time/energy response | `NOT_ESTABLISHED` | Cannot attribute catalog arrival differences to travel time |
| Validated LIV statistic or bound on lattice spacing | `BLOCKED` | Requires further work |

For an individual photon, the reported time is `t_arrival - t_GBM_trigger`, **not** `t_high - t_low` from two photons known to have been emitted simultaneously. The truth of a common emission time is not provided by Table 7. Moreover the full event selection and instrument response cannot be reconstructed from the rounded published event values alone. [NASA Fermi LAT data columns](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/lat_data_columns.html) document that proper event products contain energy, time, direction and reconstruction diagnostics; FSSC specifically [recommends extended data for GRBs](https://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/extract_latdata.html).

## Predeclared physics projection, not a fit

Reuse **unchanged** V2 dispersion:

\[
\omega(k)=\frac{2c}{a}\sin(ka/2),\quad v_g/c=\sqrt{1-q^2},\quad q(z)=\frac{aE_{\rm obs}(1+z)}{2\hbar c}.
\]

\[
\Delta t_{\rm lattice}\simeq\int_0^{z_s}\frac{\delta[q(E_h,z)]-\delta[q(E_l,z)]}{H(z)}dz,\quad
\delta(q)=\frac{q^2}{1+\sqrt{1-q^2}}.
\]

The **unobserved** 1 GeV reference (E_l), the purely illustrative `a=4e-27 m` and `a=1.616255e-35 m`, and V2 flat illustrative cosmology (`H0=70`, `Omega_m=0.3`) were fixed before reading the Table 7 delays. For each catalogued photon, the program calculates the *model-only* shift relative to that reference. It never subtracts an actual measured 1 GeV photon, never estimates `a`, and never calls the model delay an observed propagation delay.

Example: for [GRB 090510](https://arxiv.org/html/1906.11403v1), Table 7 reports a 30.0 GeV photon at **0.83 s after trigger**, redshift **0.90**. The stipulated toy `a=4e-27 m` predicts about **0.0299 s** relative to a hypothetical 1 GeV photon, and a Planck-like choice predicts about **4.88×10⁻¹⁹ s**. **Neither is an observational measurement, fitted bound or signal**; the 0.83 s timestamp is not an emission-lag-free travel-time excess.

## Empirical timing-order negative control

The selected 38 photons yield 162 within-burst comparisons of different high photon energies:

- 61 times the higher-energy photon arrived **earlier**;
- 101 times the higher-energy photon arrived **later**.

One decisive *methodological* example is **GRB 090902B**: a 40.0 GeV photon arrived at 81.7 s after the trigger, while a 21.7 GeV photon arrived at 331.9 s. Their observed difference is −250.2 s. A dispersion-only story with simultaneous emission would give the wrong explanation. But these are not simultaneous-emission pairs, so neither sign identifies nor falsifies the true photon propagation law. Instrument exposure/selection and source evolution must also be modeled.

These counts are **not** a statistical test with independent trials: many photon pairs share the same events, and the 11-burst subset is a selected group of redshift-known LAT detections.

## Frozen replay

```bash
python3 workspace/janus_cosmos_fermi_observation_gate.py \
  --batch data/cosmos/COSMOS-ALI-V3-FERMI-PUBLISHED-HE-PHOTONS-2026-10-18.json \
  > /tmp/janus-cosmos-v3-fermi.json
cmp /tmp/janus-cosmos-v3-fermi.json \
  data/cosmos/COSMOS-ALI-V3-FERMI-PUBLISHED-TOF-METHOD-RESULT-2026-10-18.json
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

The source transcription has its exact SHA-256 pinned in the engine, and V2's revealed receipt plus V1's audit remain checksum-locked. The score reports 38 actual published entries, 162 order pairs, a **conditional** delay per photon for both fixed spacings, and explicit science gates preventing a model fit to the trigger-relative timestamps. Tests reject a tampered catalog, forged probabilities, falsely claiming raw FITS/known source lags, missing events, or promoting AMPS/LIV conclusions.

## Next permitted experiment: V4 source and instrument forward model

1. Retrieve and SHA-pin NASA Fermi-LAT **Extended** photon event FITS and matching spacecraft pointing file for preregistered GRBs. Record exact search region, temporal and energy windows, event class/type, zenith cuts, IRF version and requested retrieval job ID. This is *not yet done*.
2. Retrieve compatible GBM prompt event data and source redshifts, verify MET origin, trigger epoch and leap-second conversions, and avoid energy/time precision upgrades from rounded tables.
3. Apply detector response/PSF/exposure and background selection, then construct a source emission model including energy-dependent intrinsic lags or population-level priors.
4. Fit a joint multi-burst likelihood with null (Lorentz-invariant) and *exact selected* toy dispersion; preregister likelihood tests, source nuisance marginalization, systematic variations and holdout bursts.
5. Publish source/selection reproducibility, calibration injections, trials-aware model-specific limits and independent checks against published Fermi and MAGIC work. Only then can model-specific observational exclusion or a candidate anomaly even be discussed.

**Scientific non-promotion:** `V3_PUBLISHED_PHOTONS != PHOTON_PROPAGATION_LIV_DETECTION`; `AMPS_FIREWALL = OPEN_UNCONFIRMED`; `BHB_EXPLOSION = NOT_ESTABLISHED`; `SIMULATED_UNIVERSE = NOT_ESTABLISHED`; `OSIRIS_PROOF_AUTHORITY = NO_CHANGE`.
