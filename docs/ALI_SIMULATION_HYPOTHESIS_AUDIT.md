# JANUS Cosmos V1 — ALI video, simulation hypothesis, black holes, horizons

## Scope and source reliability

User-provided video: [Что-то странное происходит с миром | ALI](https://www.youtube.com/watch?v=QBpuorUxuSM), first published 2025-03-15, approximately 1 h 14 m.

YouTube presented an automated-traffic verification page when we tried to inspect it. We found public timestamped [Video Highlight](https://videohighlight.com/v/QBpuorUxuSM) and a [Russian EchoNote summary](https://www.echonote.ru/s/tRQhUfrOnGOglhuOy9os-PRGbar_JrwL1aSaAV7vTlo). Neither is a verified verbatim transcript. All 18 thematic assertions in `data/cosmos/COSMOS-ALI-SIMULATION-CLAIMS-2026-10-18.json` are **our paraphrases of third-party summaries**, not authenticated quotations by ALI. Timestamps are summary-derived. They can contain errors and should be re-audited against the original if a transcript is supplied.

The video explores simulation metaphors with physically real phenomena (relativity, horizons, black-hole entropy) and open questions (AMPS, cosmological topology, quantum spacetime, consciousness). JANUS keeps these categories separate. A mathematical construction is neither an observation nor a proof that somebody is running a computer simulation.

## Replay and provenance

```bash
python3 workspace/janus_cosmos_simulation_audit.py --batch data/cosmos/COSMOS-ALI-SIMULATION-CLAIMS-2026-10-18.json
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

The first command's exact UTF-8 bytes match the frozen `data/cosmos/COSMOS-ALI-SIMULATION-AUDIT-RESULT-2026-10-18.json`. The audit checks the SHA-256 of the earlier BHB-R3A result **before** it is allowed to reuse exterior Schwarzschild kinematics. All R0–R3A and FW1 receipts remain unchanged.

## Six executable physics / logic controls

| Mode | Well-defined mathematical control | What cannot be concluded |
|---|---|---|
| `lorentz` | `gamma=1/sqrt(1-beta²)`, `L/L0=1/gamma`, `Δtau/Δt=1/gamma`, velocity addition `(u+beta)/(1+u*beta)` | Time dilation is CPU lag; `c` is universal video FPS |
| `schwarzschild` | Exterior lapse `alpha=sqrt(1-r_s/r)`, radial infaller `v/c=sqrt(1-alpha²/e²)`, local photon `c`, outgoing redshift `E_inf/E_loc=alpha` | Massive object reaches `c`; physical AMPS firewall; radiation flash |
| `bh_area_entropy` | `A=4*pi*r_s²`, `S_BH/kB=A/(4*lP²)`, bits `S/(kB*ln 2)` | Pixels on a literal display, universal quantum computer |
| `flat_torus` | Periodically identified Euclidean cuboid `T³`; intrinsic curvature zero, finite volume, no boundary | Our universe is demonstrably finite or has wraparound images |
| `lattice_wave` | Selected 1D nearest-neighbor model `omega=2c*sin(ka/2)/a`; `v_g/c=cos(ka/2)` | Observed Planck lattice, Lorentz-invariance violation or simulator |
| `conditional_observer_count` | `N_sim/(N_sim+N_base)` for stipulated equally weighted observer populations | Bayesian probability that our world is simulated |

The lattice is explicitly a *toy* with a preferred rest frame and no calibrated physical spacing `a`. Its high-wave-number deviations do not follow from every discrete or simulated theory. The entropy law is semiclassical and this benchmark simply evaluates it for hypothetical radii measured in Planck lengths. A flat compact 3-torus is a geometric counterexample to the claim that flatness necessarily entails infinite volume; observational cosmic topology remains unconstrained by this construction.

## Representative stable outputs

- At `v=0.5c`, `gamma=1.154700538...`, length ratio `0.866025403...`, transformed photon speed `c`.
- For radial infall at `r/r_s=1+10^-24`, speed relative to a static exterior tetrad remains strictly below `c`; the earlier R3A outer-photon energy and negative observation gate are preserved, not replaced.
- Doubling hypothetical Schwarzschild radius from 10 to 20 Planck lengths increases Bekenstein-Hawking entropy by exactly 4.
- A flat periodic cuboid with side lengths 2, 3 and 5 has volume 30, finite, without spatial boundary.
- Selected lattice nearest-neighbor group speed tends to `c` for `ka→0`, but tends to zero as `ka→pi`; this is a diagnostic of **that discretization**, not evidence of cosmic frames.
- With stipulated observer counts `(N_base,N_sim)=(1,0), (1,1), (1,999)`, illustrative fractions are 0, 0.5 and 0.999. They carry **no measured prior or likelihood**.

## Hypothesis ledger and status gates

The detailed input records 18 time-coded claims with scientific domain, paraphrase, evidence class, independent source links, negative control, linked numerical case IDs, and an explicit epistemic ceiling. Results include the distribution of assessments and reproducible individual case digests. The meta registry stores a corresponding extensive structured receipt with SHA-256 evidence pointers and no retroactive change to earlier BHB records.

Particularly important distinctions:

1. **AMPS firewall ≠ ordinary Hawking radiation ≠ real detected membrane.** AMPS [2012](https://arxiv.org/abs/1207.3123) shows tension among three assumptions about unitarity, effective field theory and smooth horizon crossing. It does not supply a confirmed material shell, interaction cross section or distant photon flux.
2. **Holography ≠ simulation hypothesis.** An entropy-area law can be physically powerful without identifying an external computer or renderer. A detailed microscopic entropy derivation exists only in selected cases, e.g. [Strominger and Vafa](https://arxiv.org/abs/hep-th/9601029).
3. **Bostrom trilemma ≠ a measured 99.9% or billion-to-one probability.** The [original argument](https://simulation-argument.com/simulation.pdf) is conditional on disputed premises and reference classes.
4. **Discrete time ≠ quantum theory's established result.** Specific lattice choices may induce phenomenology as explored by [Beane et al.](https://arxiv.org/abs/1210.1847); the observed universe is not proven to have this grid.
5. **Speed of light in SI.** The decimal figure [299,792,458 m/s](https://www.nist.gov/si-redefinition/definitions-si-base-units) is fixed exactly by the metre definition; it is not a measurement of hypothetical CPU ticks.
6. **Flatness ≠ infinite volume.** The [cosmic-topology literature](https://arxiv.org/abs/1604.02179) treats observational searches for globally identified spaces separately from local curvature.
7. **Quantum delayed-choice ≠ human consciousness sending information backwards in time.** Retrocausal signaling would need a separate no-signaling-violating dynamical prediction; none is implemented here.

The root `SIMULATION_HYPOTHESIS_OBSERVED` remains `NOT_ESTABLISHED`. `AMPS_FIREWALL_PHYSICAL_EXISTENCE` is `OPEN_UNCONFIRMED`. `BHB_EXPLOSION_OBSERVED` is `NOT_ESTABLISHED`. `OSIRIS_PROOF_AUTHORITY` stays `NO_CHANGE`. Future work would need an original transcript and preregistered distinct empirical predictions before any promotion.
