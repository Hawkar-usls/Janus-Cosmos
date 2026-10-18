# JANUS Cosmos V2 — distinguish a specified lattice photon model from Lorentz-invariant propagation

**Experiment:** `JANUS-COSMOS-ALI-V2-LATTICE-TOF-2026-10-18`  
**Status:** deterministic **synthetic** injection-and-recovery demonstration; **no astronomical observation or Lorentz-invariance violation detected**.  
**Lineage:** video topic [`QBpuorUxuSM`](https://www.youtube.com/watch?v=QBpuorUxuSM) → [ALI V1 audit](ALI_SIMULATION_HYPOTHESIS_AUDIT.md) → V2 fixed forward prediction. Source video transcript was not independently verified; this construction is a separately specified physical hypothesis, not a purported equation taken from the video.

## Hypotheses and falsifiable discriminant

- **H0 (Lorentz-invariant reference):** photon propagation group velocity is the same vacuum speed `c`, independent of photon energy. The *propagation* contribution to the high-minus-low energy arrival-time difference is exactly zero. Intrinsic source emission timing and instrument systematics need separate modeling.
- **H1 (a particular discretization toy):** for a one-dimensional nearest-neighbor field lattice, posit the dispersion `omega(k)=(2*c/a)*sin(k*a/2)`, with `E=hbar*omega`. Interpret a single radial propagation direction with a preferred comoving lattice frame. It follows that `v_group/c=cos(k*a/2)=sqrt(1-q^2)` with `q=a*E/(2*hbar*c)`.

This **one-dimensional dispersion embedded as an isotropized radial phenomenological ansatz is not an established three-dimensional discretized spacetime or a complete quantum-gravity theory**. The model is testable precisely because its assumptions are fixed. Limits on this toy cannot exclude *all* discrete spacetime hypotheses or establish a simulation.

Use the explicit, conditional redshift prescription `E(z)=E_observed*(1+z)`, flat illustrative Lambda-CDM `H(z)=H0*sqrt(Omega_m*(1+z)^3+Omega_lambda)`, and leading-order small group-speed defect:

\[
\delta(E,z)=1-v_g/c=\frac{q(E,z)^2}{1+\sqrt{1-q(E,z)^2}},\qquad q(E,z)=\frac{aE_{\rm obs}(1+z)}{2\hbar c}.
\]

\[
\Delta t_{\rm propagation}^{\rm H1}(E_h,E_l,z_s)\simeq
\int_0^{z_s}\frac{\delta(E_h,z)-\delta(E_l,z)}{H(z)}\,dz,
\qquad \Delta t_{\rm propagation}^{\rm H0}=0.
\]

These are leading-order **stipulated** FLRW propagation assumptions; a microscopic expanding lattice and cosmological photon-energy conservation are not derived. Require `q <= 0.01` along the entire trajectory. The deterministic evaluator uses 70-digit `decimal.Decimal` arithmetic and fixed 384-panel Simpson quadrature; constants `c`, `h`, `eV` are exact SI, while `Mpc`, cosmology parameters, source redshifts, lattice spacings and timing uncertainties are **declared model inputs**, not claimed measurements.

Observed pair lag is modeled as

\[
\Delta t_{\rm measured}=(1+z_s)\Delta t_{\rm intrinsic,source}
+\Delta t_{\rm propagation}+\epsilon_{\rm timing}.
\]

The source-lag term is essential. One two-energy event *cannot* determine both an unknown source emission lag and a propagation lag. `FREE_UNBOUNDED` explicitly returns `UNIDENTIFIABLE_FREE_INTRINSIC_LAG`, not a model preference. `KNOWN_ZERO`, `KNOWN_VALUE`, and `BOUNDED_SYMMETRIC` are assumptions whose independent justification would be necessary for genuine data.

## Synthetic preregistration and blind replay

In `data/cosmos/`, five frozen V2 JSON assets are versioned:

| JSON name stem | Purpose | Read by scoring? |
| --- | --- | --- |
| `...-PROTOCOL-2026-10-18.json` | Constants, cosmology, model, fixed policies, 12 case IDs, heuristic decision rule | Yes |
| `...-SEALED-TRUTH-2026-10-18.json` | Synthetic underlying H0/H1 labels, source lags and noise | **No** |
| `...-OBSERVATIONS-2026-10-18.json` | Generated energy-pair lags and per-case uncertainty, with truth-file SHA commitment but no truth labels | Yes |
| `...-BLIND-SCORE-2026-10-18.json` | Model predictions, profiled chi-squares, heuristic decisions | Produced here |
| `...-REVEAL-2026-10-18.json` | SHA-validated holdout reopening, recovery and false-positive counter | Only in reveal |

This is **process blinding**, not cryptographic secrecy or independent human blinding: the truth file is intentionally committed for public reproducibility. `abs(chi2_H0-chi2_H1) >= 25` is an explicitly preregistered **engineering decision heuristic**. It is neither a calibrated five-sigma significance statement nor a likelihood-ratio inference with justified trials/systematics.

Run from the JANUS Cosmos repository root:

```bash
python3 workspace/janus_cosmos_lattice_tof.py --prepare \
  --protocol data/cosmos/COSMOS-ALI-V2-LATTICE-TOF-PROTOCOL-2026-10-18.json \
  --truth data/cosmos/COSMOS-ALI-V2-LATTICE-TOF-SEALED-TRUTH-2026-10-18.json \
  > /tmp/janus-v2-observations.json
cmp /tmp/janus-v2-observations.json data/cosmos/COSMOS-ALI-V2-LATTICE-TOF-OBSERVATIONS-2026-10-18.json

python3 workspace/janus_cosmos_lattice_tof.py --score \
  --protocol data/cosmos/COSMOS-ALI-V2-LATTICE-TOF-PROTOCOL-2026-10-18.json \
  --observations data/cosmos/COSMOS-ALI-V2-LATTICE-TOF-OBSERVATIONS-2026-10-18.json \
  > /tmp/janus-v2-blind-score.json
cmp /tmp/janus-v2-blind-score.json data/cosmos/COSMOS-ALI-V2-LATTICE-TOF-BLIND-SCORE-2026-10-18.json

python3 workspace/janus_cosmos_lattice_tof.py --reveal \
  --protocol data/cosmos/COSMOS-ALI-V2-LATTICE-TOF-PROTOCOL-2026-10-18.json \
  --observations data/cosmos/COSMOS-ALI-V2-LATTICE-TOF-OBSERVATIONS-2026-10-18.json \
  --blind-score data/cosmos/COSMOS-ALI-V2-LATTICE-TOF-BLIND-SCORE-2026-10-18.json \
  --truth data/cosmos/COSMOS-ALI-V2-LATTICE-TOF-SEALED-TRUTH-2026-10-18.json \
  > /tmp/janus-v2-reveal.json
cmp /tmp/janus-v2-reveal.json data/cosmos/COSMOS-ALI-V2-LATTICE-TOF-REVEAL-2026-10-18.json

python3 -m unittest discover -s tests -p 'test_*.py' -v
```

## What the 12 tests mean

An intentionally **nonphysical demonstration spacing** `a=4e-27 m`, `z=0.5`, `E_low=1 GeV`, `E_high=1 TeV`, yields about **15.52845 s** of propagation delay in this exact ansatz; this is not a fitted astronomical effect. The specified Planck-length-like benchmark `a=1.616255e-35 m` yields about **2.54e-16 s**, far below the illustrative `0.5 s` timing noise. Fixing source lag as independently known zero allows the synthetic cases to prefer their respective injected models, but permitting a free source lag makes both models degenerate.

Most important negative control: an H0/Lorentz-invariant synthetic event **with an intrinsic source lag deliberately equal to the H1 predicted propagation delay** spuriously prefers H1 if one *incorrectly* assumes zero intrinsic lag. Reveal explicitly reports one such false positive. With `FREE_UNBOUNDED` on the same event, the scorer correctly refuses identification.

Frozen reveal counts: **4 correct preferences, 1 deliberately induced false preference, 5 inconclusive and 2 unidentifiable**. These counts measure test harness behavior, not a discovery rate, an observational bound, or evidence of discreteness.

## Experimental promotion gate: what real measurement would require

1. Preselect an actual transient sample and instrument responses, energy/time calibration, observatory clock cross-check and photon association rules before inspecting arrival-time residuals.
2. Quantify source-intrinsic spectral evolution, unresolved overlapping pulses, energy-dependent source lags, redshift uncertainty, EBL absorption/selection bias and pulse-window choice; an unconstrained emission lag breaks this one-event test.
3. Fit H0 and the **predeclared** H1 model across *independent* bursts at different redshifts and energies using a joint likelihood with astrophysical nuisance parameters, population scatter and null injections.
4. Calibrate real-data uncertainties, trials and selection functions with simulated instrument/source data and control samples; derive exclusion bands for `a` without recycling this toy `Delta-chi2=25` as a p-value.
5. Require independent reproduction and consistency with already published Lorentz-invariance constraints before claiming a new limit. Absence of observed energy-delay scaling constrains only selected model assumptions. An apparent positive scaling is **not evidence of a virtual/simulated universe**, quantum horizon firewall or BHB explosion.

## Scientific references and inherited gates

- MAGIC Collaboration, *Bounds on Lorentz Invariance Violation from MAGIC Observation of GRB 190114C*, [arXiv:2001.09728](https://arxiv.org/abs/2001.09728). This is real scientific literature used for comparison of propagation delays and intrinsic time-lag confounding; **no GRB 190114C photon data were used here**.
- Beane et al., *Constraints on the Universe as a Numerical Simulation*, [arXiv:1210.1847](https://arxiv.org/abs/1210.1847). Tests depend on particular lattice/discretization assumptions.
- [NIST SI definitions](https://www.nist.gov/si-redefinition/definitions-si-base-units).
- Pinned V1 Cosmos audit SHA-256: `3b78e75eb93933e25755f501b45db1d78c4385356e0778e251055fd4353cba49`. V1 already pins BHB R3A. **AMPS quantum firewall: unconfirmed. BHB Boom: not established. OSIRIS proof authority: unchanged.**
