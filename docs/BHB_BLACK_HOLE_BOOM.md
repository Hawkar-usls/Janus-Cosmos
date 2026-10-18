# JANUS Cosmos — BHB (Black Hole Boom)

`BHB-R1/R2` is a computational research lane for **near-horizon local collision
kinematics**, not a discovered physical explosion, not a high-energy observatory,
and not a numerical black-hole spacetime solver.

## Repository integration

- Executable engine: `workspace/janus_bhb.py` (Python standard library, Decimal 90-digit precision).
- Kerr engine: `workspace/janus_bhb_kerr.py` (equatorial Kerr, Decimal 90-digit precision).
- Frozen upstream evidence: `data/bhb/BHB-R0-UPSTREAM-2026-10-17.json`, SHA-256
  `3fa8b5027edd79a9a9e5d1d7561d7ca029debc7f6656142a34167842ef76f067`.
- Upstream Meta Registry commit: `5b41a09c0e7885a649610a9e3b3ffe2a1c4778f2`, original blob `10ca85cab017a936ef5d3511d932dcc85cbb7696`.
- Local input: `data/bhb/BHB-COSMOS-R1-CASES-2026-10-17.json`.
- Frozen local output: `data/bhb/BHB-COSMOS-R1-RESULT-2026-10-17.json`.
- Kerr case matrix and frozen output: `data/bhb/BHB-COSMOS-R2-CASES-2026-10-18.json`, `data/bhb/BHB-COSMOS-R2-RESULT-2026-10-18.json`.
- Negative controls and replay: `tests/test_bhb_cosmos.py`.
- Kerr null-ray, radial barrier and invariance controls: `tests/test_bhb_kerr.py`.
- CI: `.github/workflows/bhb-schwarzschild-gate.yml`.

## How to run

```bash
python3 workspace/janus_bhb.py --replay-upstream
python3 workspace/janus_bhb.py --batch data/bhb/BHB-COSMOS-R1-CASES-2026-10-17.json
python3 workspace/janus_bhb_kerr.py --batch data/bhb/BHB-COSMOS-R2-CASES-2026-10-18.json
python3 -m unittest discover -s tests -p 'test_bhb_cosmos.py' -v
python3 -m unittest discover -s tests -p 'test_bhb_kerr.py' -v
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
| `ESCAPING_SIGNAL` | `OPEN` (no emission/ray tracing in R1) |
| `KERR_ROTATION` | `NOT_MODELED` (historical R1 gate; see R2 below) |
| `REAL_BHB_DETECTION` | `NOT_ESTABLISHED` |

This module cannot be wired into the existing Newtonian velocity-Verlet
`janus_celestial_dynamics.py` as though the latter described event horizons.
They are explicitly separate model classes. Neither code path modifies
`PROJECT_STATUS.json` or asserts new astronomical evidence.

## Research lineage

- Bañados, Silk & West (2009), `arXiv:0909.0169`: Kerr horizon collision effects, outside the scope of this Schwarzschild R1.
- Toporensky & Zaslavskii, `arXiv:2411.01989`: near-horizon Schwarzschild head-on collisions and white-hole provenance.
- Toporensky et al., `arXiv:2609.14701`: classification, physical realizability and kinematic censorship.

## R2 — rotating Kerr, equatorial test geodesics

All radii are `r/M`, where geometrized `G=M=c=1`; `spin_a_over_M` is signed
and satisfies `|a/M|<=1`. The code rejects `r<=r_+=1+sqrt(1-a²)` and does not
interpret the horizon itself as a place where a massive particle attains c.

For an equatorial test particle with specific conserved `E` and `Lz/M`,
`Delta=r²-2r+a²`, `P=E(r²+a²)-aL`, the radial polynomial is

\[
R=P^2-\Delta\,[r^2+(L-aE)^2].
\]

The program computes the local invariant two-particle `E_cm` from the ZAMO
orthonormal frame, and independently validates it using the coordinate Kerr
metric contraction `-g(u_1,u_2)`. For particles claimed to be **infalling
from infinity**, it checks `E>=1`, future-directedness, and the *global*
minimum of the cubic `R/r` on `[r_event,infinity)`. A forbidden interval or
double-root asymptotic separatrix fails that provenance gate even if `R` is
locally nonnegative. The incoming critical orbit near an **extremal** Kerr
hole is admitted only at strictly exterior finite `r>r_+`; no finite collision
at the exact horizon is asserted.

An independent hypothetical **outward equatorial photon** is specified by its
ZAMO-frame azimuthal sine `s`, with radial cosine `sqrt(1-s²)>0` (no photon
production mechanism assumed). Its conserved energy transfer is

\[
\frac{E_\infty}{E_{\rm ZAMO}}=\alpha+\omega\sqrt{g_{\phi\phi}}s,\qquad
b=\frac{\sqrt{g_{\phi\phi}}s}{E_\infty/E_{\rm ZAMO}}.
\]

For positive `E_infinity`, the null radial polynomial is
`R_null/r = r³+(a²-b²)r+2(b-a)²`. Its minimum on
`[r_emission,infinity)` analytically determines whether that outward
**equatorial** null geodesic is clear of exterior radial barriers. The
frame-local radial null norm is independently rechecked. A positive escape
gate is for **that given ray**, not for photons emitted in the collision.

### Reproducible R2 result snapshot (nominal illustrative protons)

| Scenario | Local `E_cm` | Trajectory gate | Independent outgoing photon |
| --- | ---: | --- | --- |
| Schwarzschild limit `a=0,r=5` | 1.876544 GeV | admissible | radial ray escapes; transfer `0.77460` |
| Generic Kerr `a=0.9,r=1.6,L₂=2` | 2.429815 GeV | admissible | radial ray escapes; transfer `0.19695` |
| Kerr `a=0.9,r=1.5,L₂=2.7` | locally computable, not admitted from infinity | radial barrier blocks | independent radial ray escapes |
| Extremal Kerr limit `a=1,r=1.0001,L₂=2` | 143.633128 GeV | exterior idealized geodesics admitted | radial ray escapes; transfer `0.00005000` |
| Schwarzschild `r=2.5,s=0.99` | 1.876544 GeV | admissible | ray has exterior barrier: **blocked** |

Results assume eternal stationary Kerr geometry, noninteracting test particles
until the stated local event and exact equatorial symmetry. They ignore
backreaction, finite-size particles, electromagnetic interactions, scattering
matrix elements, photon spectra, detector response and event rates. The photon
is prescribed independently of particle momentum conservation. In particular,
even when the `E_cm` grows, there is **no established BHB explosion, escaped
collision luminosity or telescope detection**.

### R2 scientific gates

| Gate | R2 verdict |
| --- | --- |
| `KERR_EQUATORIAL_INVARIANT` | `PASS` for valid evaluated local cases, coordinate/ZAMO crosscheck |
| `RADIAL_TRAJECTORY_ACCESS` | checked case-by-case; some rejected |
| `HYPOTHETICAL_EQUATORIAL_PHOTON_TRANSPORT` | checked per prescribed ray, escaped and blocked controls |
| `ESCAPED_COLLISION_SIGNAL` | `OPEN_NO_PRODUCT_MODEL` |
| `REAL_BHB_DETECTION` | `NOT_ESTABLISHED` |

Primary references: Bañados, Silk & West `arXiv:0909.0169`; Jacobson &
Sotiriou `arXiv:0911.3363` (physical limitations); Yan, Guo & Chen
`arXiv:2104.07889` (Kerr photon escape and energy shifts). R2 does not claim
to reproduce their full emission probability calculations.

## FW1 — the quantum black-hole firewall (AMPS), not the classical horizon

`FW1` is a separate **hypothesis overlay**, not a modification to the R2
classical Kerr equations. The original Almheiri–Marolf–Polchinski–Sully
argument (`arXiv:1207.3123`) makes three statements mutually incompatible
under its assumptions for an old, evaporating black hole:

1. Hawking evaporation is unitary (radiation eventually has a pure state).
2. Low-energy semiclassical field theory describes the exterior sufficiently
   close to the horizon.
3. An infalling observer encounters a locally ordinary, smooth horizon
   (*no drama*).

The hypothetical firewall is one **proposed resolution**, not a universally
accepted or directly observed physical membrane. Entanglement structure,
post-Page-time state, quantum backreaction, possible nonlocality and
microstate-based alternatives require theory well beyond classical R2
geodesics. See also AMPS et al. `arXiv:1304.6483` and Mathur
`arXiv:0909.1038`. Neither ``the event horizon'' nor ``the AMPS firewall''
means that a massive infaller achieves the local speed of light.

The JSON contract is in
`data/bhb/BHB-FIREWALL-FW1-HYPOTHESES-2026-10-18.json`; the machine-generated
frozen output is in
`data/bhb/BHB-FIREWALL-FW1-RESULT-2026-10-18.json`.

```bash
python3 workspace/janus_bhb_firewall.py
python3 -m unittest discover -s tests -p 'test_bhb_firewall.py' -v
```

FW1 binds the **byte-exact** R2 frozen JSON SHA-256 and its internal receipt
hash; it rejects modified baselines, altered quantum premises and any attempt
to set `physical_firewall_presence=true` or claim a detection. For each of
the seven R2 cases it records the coordinate gap `r-r+` solely as geometry,
not an assumed firewall thickness or proper distance. It preserves local
classical `E_cm` only as a **conditional** baseline, with quantum-firewall
collision energy, deposited energy, emitted spectrum and escaped luminosity
all `null` until a defensible dynamic model exists.

| FW1 scientific gate | State |
| --- | --- |
| `CLASSICAL_R2_RECEIPT_INTEGRITY` | `PASS` |
| `AMPS_LOGICAL_PREMISES` | `REGISTERED_CONDITIONAL_CONFLICT` |
| `FIREWALL_PHYSICAL_EXISTENCE` | `OPEN_UNCONFIRMED` |
| `BLACK_HOLE_AGE_POST_PAGE_TIME` | `UNKNOWN_NO_QUANTUM_STATE_DATA` |
| `QUANTUM_FIREWALL_DYNAMICS` | `NOT_MODELED` |
| `ESCAPED_COLLISION_SIGNAL` | `OPEN_NO_PRODUCT_MODEL` |
| `REAL_BHB_OR_FIREWALL_DETECTION` | `NOT_ESTABLISHED` |

The physical firewall must not be confused with the **JANUS science
firewall**, which is only a software evidence boundary preventing promotion
of unverified conjectures into scientific claims.
