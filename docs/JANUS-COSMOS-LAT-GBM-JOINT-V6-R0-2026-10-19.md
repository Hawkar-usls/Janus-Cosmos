# JANUS Cosmos V6-R0 — real GBM spectral forward-fold pilot and LAT exposure/control diagnostics

**Scientific status:** source-verified instrument and real count-analysis progression, **not a complete calibrated joint LAT/GBM likelihood**. Consequently, **no observed Lorentz-invariance-violation p-value, lattice-spacing limit, AMPS firewall, or Black Hole Boom detection** follows from V6-R0.

## Verified legacy and official data

- The 9 original `P8R3_SOURCE_V3` CALDB Aeff/PSF/EDISP FITS, all 27 GBM TTE/CSPEC/RSP2 inputs and all six LAT Extended EVENTS/spacecraft files for GRBs **090510, 090902B, 130427A** were rehydrated using their existing archived SHA-256 pins. Every byte matched the frozen V4/V5 manifests.
- The V5 instrument receipt reproduced **byte-for-byte**: SHA-256 `53c4e8bacb55c64914ae263c0295aad53055d342cdac8fa95bb46e7974d38424`. It accounts for 6,962,621 GBM individual TTE events and 175 selected LAT 5° ROI candidates. The older V1–V5 result files were not modified.
- For the LAT instrument response, the official FSSC recommends event-class- and event-type-matched P8R3 IRFs with separate effective area, PSF and energy-dispersion factors. A simplistic effective-area lookup is **not** exposure, nor does it establish photon source identity.

## Actual advance 1: OGIP GBM RSP2 spectral forward fold

`workspace/janus_cosmos_fermi_response_v6.py` implements a bounded reader for the OGIP `SPECRESP MATRIX` binary heaps. It handles `PI(1)`/`PI(128)` channel-group arrays and `PE(128)` response arrays; importantly, `F_CHAN=1..128` in these NASA RSP2 files is **one-based**, while GBM TTE PHA channels use **zero-based** indexing. Every matrix validates 140 true-energy bins, 128 PHA channels, nonnegative finite response values, heap pointer extent, true-energy ordering and response time ranges. It excludes the zero-duration terminal RSP2 extension in the 090902B archive.

For each preregistered time window it mixes attitude-specific response matrices by **their overlap exposure fraction**. The resulting forward fold is

\[
s_{dct} = \Delta t_t\,A_t\sum_e R_{dce}^{(t)} \int_{E_{e,\min}}^{E_{e,\max}} \left(\frac{E}{100\,\mathrm{keV}}\right)^{\alpha_t}\!\exp\!\left[-\frac{(2+\alpha_t)E}{E_{\mathrm{peak},t}}\right]dE.
\]

The pilot evaluates each energy-bin integral by a declared midpoint rule. `R` has units cm²; `A` has units photons cm⁻² s⁻¹ keV⁻¹ at 100 keV. There are two NaI plus one BGO detector per GRB and **17** time-resolved spectral intervals in total. Two NaI detectors share the 8–900 keV selected PHA-overlap band; BGO uses 300–30,000 keV.

Observed on-source counts and the frozen pre-trigger window `[-20,-5] s` are fitted using a per-channel profiled on/off Poisson nuisance model:

\[
n_{dc}\sim\mathrm{Poisson}(s_{dc}+b_{dc}), \quad m_{dc}\sim\mathrm{Poisson}(\tau b_{dc}),\quad \tau=\Delta t_{\mathrm{off}}/\Delta t_{\mathrm{on}}.
\]

There is no Gaussian approximation for low-count channels. A bounded 1D profile of spectral normalization is evaluated on a predeclared shape grid \(\alpha\in\{-1.6,-1.2,-0.8,-0.4\}\), \(E_{\mathrm{peak}}\in\{150,300,600,1200,2400,4800\}\) keV. The same offsource counts are reused between time bins, so **the individual interval likelihood values must not be multiplied or summed as independent evidence** without a single shared background model.

The spectral-mismatch descriptive stress threshold is Poisson deviance per fitted channel above 5 (not a calibrated goodness-of-fit p-value). It identified particularly poor single-component fits:

| GRB | Stress windows (relative to trigger) | Interpretation |
| --- | --- | --- |
| 090510 | None above descriptive threshold | Does **not** prove a single cutoff power law is complete |
| 090902B | 5–10 s | Single cutoff law misses structure; additional spectral components may be required |
| 130427A | 0–10 s; 10–20 s | Substantial model residuals during bright emission |

These real findings prohibit using these spectral pilots as a validated source photon-time prior, and make a Band-plus-power-law fit, more robust background, and predictive adequacy checks the next mandatory model improvements. The full per-window model parameters, observed count totals, response weights and *diagnostic* deviances are retained in `data/cosmos/COSMOS-ALI-V6-GBM-RSP2-ONOFF-SPECTRAL-2026-10-19.json`.

## Actual advance 2: LAT control region and approximate effective exposure

`workspace/janus_cosmos_fermi_lat_diagnostics_v6.py` re-parses NASA Extended `EVENTS`, `GTI`, and `SC_DATA`. It exactly reconstructs **all 175** frozen LAT >1 GeV SOURCE-class 5° source-region candidates, then independently applies the same cuts out to 10° and creates a **5–10° sky annulus control**. Pass8 conversion/PSF/EDISP event-type bit partitions are validated for every selected LAT event.

The nominal time-interval response proxy integrates spacecraft livetime with the nearest FRONT+BACK effective-area bin evaluated at each SC_DATA boresight-source inclination and 1/10 GeV. It accounts for real good-time cuts, sky pointing and nominal livetimes, but not continuous incidence integration, azimuth dependence, livetime-efficiency corrections, PSF/EDISP, or a measured diffuse background. It is **not** a LAT flux, calibrated exposure, or likelihood.

Using exact GBM triggers, the diagnostic prompt 5° LAT ROI counts were 19 (090510, 0–2 s), 42 (090902B, 0–40 s), and 47 (130427A, 0–140 s). Prompt 5–10° control-annulus counts were 0, 0, and 1 respectively. These are **candidate events**, not background-subtracted source photon memberships; the off-pulse control can include GRB afterglow, and the annulus is not an adequate diffuse model.

## Cross-instrument V6 lock

`workspace/janus_cosmos_fermi_joint_v6.py` combines both real instrument diagnostics and the frozen GBM soft-band peak context, and explicitly refuses to promote a synthetic source-lag toy into a real observation. With the propagation time coefficient \(K_j(E,z_j)\),

\[
t_{ji}=\mu_j+(\eta+\beta_j)K_j(E_{ji},z_j)+\epsilon_{ji}
\]

has a structural degeneracy \(\eta\to\eta+\delta\), \(\beta_j\to\beta_j-\delta\) when burst-specific intrinsic delays are unconstrained. Even accurate GBM keV timing does not independently fix the emission clock of LAT GeV photons. No observed joint LIV coefficient, confidence interval or p-value is computed. This is appropriate: the astrophysical source-lag premise is not yet established, and LAT/GBM full response models remain incomplete.

## Reproduction

After obtaining the frozen NASA/Fermitools input binaries at the paths below, run from the Cosmos repository root:

```bash
# Downloads ~330 MB of byte-pinned upstream FITS/official Fermitools sources.
python3 scripts/fetch_fermi_v6_inputs.py \
  --gbm-dir /tmp/janus_v6_gbm --caldb-dir /tmp/janus_v6_caldb \
  --package /tmp/janus_v6_fermitools.conda
python3 scripts/fetch_fermi_lat_v4.py --raw-dir /tmp/janus_v6_lat_v4_raw

python3 workspace/janus_cosmos_fermi_response_v6.py --gbm-dir /tmp/janus_v6_gbm
python3 workspace/janus_cosmos_fermi_lat_diagnostics_v6.py \
  --caldb-dir /tmp/janus_v6_caldb --raw-dir /tmp/janus_v6_lat_v4_raw
python3 workspace/janus_cosmos_fermi_joint_v6.py
python3 -m unittest discover -s tests -q
```

Offline tests operate on the committed JSON receipts and synthetic FITS/OGIP fixtures, with Python standard library only. The upstream files remain referenced by exact SHA-256 rather than embedded in Git.

## V6-R1 — exact next falsification gate

1. Fit and falsify more complete time-resolved **Band + extra hard power-law** GBM source spectra with time-variable background polynomials and explicit CSPEC/RSP2 channel-energy calibration nuisance terms, comparing held-out times/detectors.
2. Integrate LAT PSF, energy dispersion, full effective exposure/azimuth/efficiency and Galactic+isotropic background with time/energy-dependent source membership, ideally through official ScienceTools or cross-checked against it.
3. Adopt a defensible independently evidenced hierarchical GeV emission-delay prior; include burst-to-burst variation, selection effects and instrumental clock uncertainty.
4. Preregister injection–recovery and coverage simulations, null/permutation controls and detection thresholds. Only after the entire end-to-end physical likelihood is validated may any *observed* LIV p-value or limit be released.

The AMPS quantum firewall and BHB Boom still lack an explicit calibrated microscopic horizon-to-exterior emission prediction; Fermi GRB data do not directly establish that phenomenon.

### Instrument primary sources

- [NASA LAT P8R3 instrument response](https://fermi.gsfc.nasa.gov/ssc/data/analysis/lat_irfs/irf_overview.html)
- [NASA LAT background models](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/BackgroundModels.html)
- [NASA GBM response, spectral folding and time-dependent DRM](https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm/gbm_data_tools/gdt-docs/notebooks/Responses.html)
- [NASA GBM source/background spectral analysis](https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm/gbm_data_tools/gdt-docs/notebooks/PhaExport.html)
- [Vasileiou et al. (2013), GRB LIV constraints and intrinsic spectral evolution](https://arxiv.org/abs/1305.3463)
