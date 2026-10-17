# JANUS Cosmos — BHB (Black Hole Boom)

`BHB-R1` is a computational research lane for **near-horizon local collision
kinematics**, not a discovered physical explosion, not a high-energy observatory,
and not a numerical black-hole spacetime solver.

## Repository integration

- Executable engine: `workspace/janus_bhb.py` (Python standard library, Decimal 90-digit precision).
- Frozen upstream evidence: `data/bhb/BHB-R0-UPSTREAM-2026-10-17.json`, SHA-256
  `3fa8b5027edd79a9a9e5d1d7561d7ca029debc7f6656142a34167842ef76f067`.
- Upstream Meta Registry commit: `5b41a09c0e7885a649610a9e3b3ffe2a1c4778f2`, original blob `10ca85cab017a936ef5d3511d932dcc85cbb7696`.
- Local input: `data/bhb/BHB-COSMOS-R1-CASES-2026-10-17.json`.
- Frozen local output: `data/bhb/BHB-COSMOS-R1-RESULT-2026-10-17.json`.
- Negative controls and replay: `tests/test_bhb_cosmos.py`.
- CI: `.github/workflows/bhb-schwarzschild-gate.yml`.

## How to run

```bash
python3 workspace/janus_bhb.py --replay-upstream
python3 workspace/janus_bhb.py --batch data/bhb/BHB-COSMOS-R1-CASES-2026-10-17.json
python3 -m unittest discover -s tests -p 'test_bhb_cosmos.py' -v
```

A new batch file contains a `cases` array. Each case has `r_over_rs` (decimal
string greater than 1), `physical_context` (`collapsing_black_hole` or
`eternal_schwarzschild`), and two particle objects. Each particle specifies
`direction` (`in` or `out`), `origin` (`infinity`, `white_hole`,
`local_production`, `unverified`), `rest_energy_GeV` (optional, nominal
0.938272), and `specific_energy` (positive conserved energy per rest mass,
default 1). `spacetime` must be `schwarzschild` when present. `r/r_s=1` is
rejected: the static tetrad is undefined on the horizon itself.

The runner emits machine-readable JSON to standard output; it does **not**
change weights, registry authority, observational ledgers or canonical OSIRIS
proof state. The CI job compares the generated R1 JSON byte-for-byte with its
frozen receipt and retains a run receipt as a workflow artifact.

## Physical kernel

Use geometrized units for the radial Schwarzschild background:

\[
f=1-\frac{r_s}{r},\qquad
\gamma_i=\frac{e_i}{\sqrt f},\qquad
\beta_i=\sigma_i\sqrt{1-\frac{f}{e_i^2}}.
\]

For particle rest energies \(m_1,m_2\) (here expressed in GeV),

\[
E_\mathrm{cm}^2=m_1^2+m_2^2+2m_1m_2\gamma_{12},\qquad
\gamma_{12}=\gamma_1\gamma_2(1-\beta_1\beta_2).
\]

Same-sign near-horizon velocities require rationalized evaluation to avoid
catastrophic cancellation. An independent local Minkowski tetrad invariant
checks the computed result. Both values are strictly local to a collision
event and do not represent energy detectable at infinity.

## Frozen result interpretation

- Ingoing + ingoing, both with \(e=1\): local \(E_\mathrm{cm}=1.876544\) GeV for two illustrative nominal protons, including very near the horizon.
- Ingoing + ingoing, \(e_1=1,e_2=2\): finite, approximately 1.99038 GeV at \(r=1.000001r_s\).
- Idealized ingoing + outgoing \(e=1\): 1876.54494 GeV at the same radius, but outgoing radial free fall is only conditionally modeled in the **eternal white-hole extension**.
- Applying that white-hole worldline to a **collapsing** Schwarzschild black hole fails the trajectory-origin gate.
- Locally produced outgoing matter retains **source and energy-budget debt**; kinematic possibility does not establish a production mechanism.

### Scientific gates

| Gate | Status |
| --- | --- |
| `SCHWARZSCHILD_LOCAL_INVARIANT` | `PASS` (analytic model and numeric checks) |
| `ALL_TRAJECTORIES_ADMISSIBLE` | `NOT_ESTABLISHED` for mixed case matrix |
| `ESCAPING_SIGNAL` | `OPEN` (no emission/ray tracing) |
| `KERR_ROTATION` | `NOT_MODELED` |
| `REAL_BHB_DETECTION` | `NOT_ESTABLISHED` |

This module cannot be wired into the existing Newtonian velocity-Verlet
`janus_celestial_dynamics.py` as though the latter described event horizons.
They are explicitly separate model classes. Neither code path modifies
`PROJECT_STATUS.json` or asserts new astronomical evidence.

## Research lineage

- Bañados, Silk & West (2009), `arXiv:0909.0169`: Kerr horizon collision effects, outside the scope of this Schwarzschild R1.
- Toporensky & Zaslavskii, `arXiv:2411.01989`: near-horizon Schwarzschild head-on collisions and white-hole provenance.
- Toporensky et al., `arXiv:2609.14701`: classification, physical realizability and kinematic censorship.

The next scientific gate is a **separate** rotating Kerr geodesic model with
checked accessibility of critical trajectories, local collision energy,
collision-product transport, and energy measured at infinity. None of these
gates is silently passed by the R1 radial result.
