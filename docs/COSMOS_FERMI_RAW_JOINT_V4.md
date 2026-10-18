# JANUS Cosmos V4 — individual Fermi-LAT Extended FITS and joint-GRB model

Date: 2026-10-18. Experiment: `JANUS-COSMOS-ALI-V4-FERMI-LAT-RAW-JOINT-2026-10-18`.

## Observed data: exactly what has been acquired

Three targeted requests to NASA FSSC `LATDataQuery.cgi` were completed using
**Extended** Pass 8 event data (not the FSSC `Photon` stream) and spacecraft
information. Sky/time/energy boxes, query identifiers, direct FSSC result
URLs, binary sizes and SHA256 sums for each **six original downloaded FITS
files** are recorded in the [immutable-input manifest](../data/cosmos/COSMOS-ALI-V4-FERMI-RAW-FITS-MANIFEST-2026-10-18.json).
The binaries are ~1.5 MB collectively and are **not checked into Git**;
the included download/verification script retrieves their exact bytes when
NASA retains the generated query products. When a source URL expires, users
must repeat the identical FSSC Extended query and audit the new hashes, not
substitute another stream or loosen a checksum.

| GRB | Adopted redshift (approx.) | FSSC extended ROI rows | Fixed diagnostic cuts, retained | Nominal spacecraft livetime (s) |
|---|---:|---:|---:|---:|
| 090510 | 0.903 | 271 | 32 | 554.347038 |
| 090902B | 1.822 | 1,226 | 67 | 1,409.424064 |
| 130427A | 0.34 | 5,023 | 76 | 4,269.321754 |
| **Total** | — | **6,520** | **175** | — |

The exported 175 rows include Fermi MET, provisional offset from the GBM
trigger, measured MeV energy, angular distance, incident angle, zenith angle,
`EVENT_CLASS`/`EVENT_TYPE` masks, event and run identifiers. They are real
**LAT photons in a sky aperture**, *not* 175 proven GRB-origin photons. Each
source's actual temporal and spatial background, PSF, and exposure matter.

## Fixed diagnostic cut / instrumentation

1. Both files identify FSSC LAT; the EVENT HDU contains 32-bit Pass8 event
   masks, and the query `POS`, `TIME`, `ENERGY` FITS keywords reproduce the
   request. The event and spacecraft FITS query windows must agree.
2. Events must lie inside a valid LAT `GTI` and one spacecraft interval with
   `DATA_QUAL=1`, `LAT_CONFIG=1`, `IN_SAA=false`.
3. Select events with the Pass8 **SOURCE** mask bit 128 and recorded energy
   1–300 GeV; reject `ZENITH_ANGLE>90°`, `THETA>65°` and angular separation
   from the nominated GRB position `>5°`. These are *identical, intentionally
   simple diagnostic cuts* across the three targets and not a scientifically
   optimized LAT GRB analysis.
4. Preserve individual event/run identities and source-file SHA256, plus
   `PROC_VER`, `TIMESYS`, `TIMEREF`, GTI counts, spacecraft time intervals and
   nominal livetime. File content must match exact recorded bytes or fail.

Instrument quantities still missing for a physical analysis: Pass8 `P8R3`
effective area vs energy/incidence/event type, point-spread distribution,
energy dispersion, varying exposure, diffuse and transient source background,
probabilistic GRB association, systematic uncertainties, and a scientifically
validated GBM–LAT cross-instrument clock correction. Reading `THETA`, livetime
and data-quality bits is **not the same as forward-folding the IRF**.

Primary references:

- [NASA LAT data extraction and Extended selection](https://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/extract_latdata.html)
- [NASA LAT event columns](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/lat_data_columns.html)
- [NASA Pass8 event classes and IRF guidance](https://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Pass8_usage.html)
- [NASA FSSC query format](https://fermi.gsfc.nasa.gov/ssc/LATDataQuery_help.html)
- [Published 2FLGC (Ajello et al., 2019)](https://arxiv.org/abs/1906.11403) — separate rounded V3 benchmark; the V4 175-photon extraction is **not** Table 7.

## Explicit multi-GRB source + propagation toy likelihood

For burst \(j\), synthetic photon \(i\), toy Gaussian measurement noise
\(\epsilon_{ji}\), source epoch \(\mu_j\) and free intrinsic energy-dependent
source slope \(\beta_j\):

\[
t_{ji}=\mu_j+(\eta+\beta_j)K_j(E_{ji},z_j)+\epsilon_{ji},
\quad \epsilon_{ji}\sim N(0,\sigma_j^2).
\]

Here \(K_j\) is the V2 illustrative lattice/FLRW leading-order propagation
delay relative to a *hypothetical* 1 GeV photon, evaluated at an illustrative
spacing \(a_{\rm ref}=4\times10^{-27}\,\mathrm m\). The synthetic test
parameter \(\eta=(a/a_{\rm ref})^2\) is constrained nonnegative, with its
toy profile evaluated on \([0,2]\). For a fixed energy-redshift profile,
each burst independently observes \(\eta+\beta_j\), **not** \(\eta\).

An unconstrained \(\beta_j\) permits the transformation
\((\eta,\beta_j)\mapsto(\eta+\delta,\beta_j-\delta)\), leaving every
predicted event time unchanged. **Different GRB redshifts do not remove this
degeneracy if each intrinsic slope is free.** The engine profiles each
\(\mu_j\) and \(\beta_j\) analytically, then compares:

- **Known-zero source slope**, an intentionally idealized conditional
  assumption: the seeded \(\eta=0.8\) is recovered numerically at ~0.799.
- **Gaussian source-slope prior**, another declared assumption: \(\eta\)
  becomes conditionally identifiable, but its result depends on the prior
  width and must not be called an empirical inference.
- **Independent unbounded slopes**: the profile is flat and the output is
  `eta_mle=null`, `TS=0`, even with three redshifts.
- **Adversarial injection**: setting true propagation \(\eta=0\) but
  intrinsic \(\beta_j=0.8\) falsely reproduces ~0.799 if slopes are
  incorrectly fixed at zero. The correctly free-slope model cannot infer
  propagation.

The deterministic 160-replicate toy null exercise calibrates only this simple
Gaussian detector-free likelihood, with a parameter on the \(\eta=0\)
boundary. The output's Monte Carlo probability is explicitly
`not_an_observed_p_value=true`. Its precision, trials, burst distribution and
assumed noise are **inapplicable to real Fermi data**.

## Reproduce

Python 3.11+ standard library only:

```bash
# Download if NASA query URLs are still available, otherwise requery and audit.
python3 scripts/fetch_fermi_lat_v4.py --raw-dir /tmp/janus-v4-lat
python3 scripts/fetch_fermi_lat_v4.py --raw-dir /tmp/janus-v4-lat --verify-only

# Byte-exact observed event extraction after checking all six SHA256 values.
python3 workspace/janus_cosmos_fermi_raw_v4.py \
  --raw-dir /tmp/janus-v4-lat \
  --output /tmp/v4-lat-observed.json
cmp /tmp/v4-lat-observed.json data/cosmos/COSMOS-ALI-V4-FERMI-RAW-EVENT-DIAGNOSTICS-2026-10-18.json

# Offline synthetic multiburst likelihood and its negative controls.
python3 workspace/janus_cosmos_fermi_joint_v4.py --output /tmp/v4-joint-toy.json
cmp /tmp/v4-joint-toy.json data/cosmos/COSMOS-ALI-V4-JOINT-TOY-LIKELIHOOD-2026-10-18.json
python3 -m unittest discover -s tests -p 'test_cosmos_fermi_raw_v4.py' -v
```

CI verifies the manifest/provenance assertions, frozen event-level result,
synthetic repeatability and all previous V1–V3 gates offline. It does not
re-download archive FITS on every PR because generated NASA URLs can expire.

## Next scientific gate

Before attempting any *observed* multi-burst likelihood, require all of:

1. NASA calibration IRFs (effective area, PSF, energy dispersion) and
   exposure/background templates for each event class and time interval.
2. Correlated GBM prompt light curves/spectra with clock synchronization,
   a parameterized source-energy and pulse-evolution emission model.
3. Realistic, preregistered priors or independently constrained source lags,
   and injection/recovery plus mock catalog coverage **through the IRF**.
4. A control sample and source-population/selection model, blind test
   statistic, trials accounting and independent numerical replication.

Until these are present, `ASTROPHYSICAL_MULTIBURST_LRT=BLOCKED`,
`LIV_OBSERVED=NOT_ESTABLISHED`, `LATTICE_SPACING_CONSTRAINT=NOT_DERIVED`,
`AMPS_FIREWALL=OPEN_UNCONFIRMED`, `BHB_BOOM=NOT_ESTABLISHED`.
