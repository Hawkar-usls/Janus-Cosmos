# JANUS COSMOS — SPACE · EARTH · OCEAN · OSIRIS / S𓂸ḥ

**Janus-Cosmos is the planetary-scale observation and falsification domain of JANUS. It is not limited to outer space.**

Its scope includes:

- outer space and astronomical observations;
- near-Earth and atmospheric phenomena;
- Earth-surface and geophysical data;
- oceans and seas;
- seafloor bathymetry, multibeam, sidescan and hydroacoustics;
- navigation, sound-speed and instrument-calibration problems;
- cross-domain cases where a signal can only be understood by joining sky, atmosphere, land, ocean and instrument state.

The core rule is simple:

```text
RAW DATA > STORY
MEASUREMENT > INTERPRETATION
REPLICATION > PRETTY PATTERN
NEGATIVE RESULT = FIRST-CLASS RESULT
```

Janus-Cosmos therefore keeps positive, negative, ambiguous and blocked results; freezes coordinates and hypotheses before unblinding where possible; separates instrument artefacts from physical features; and requires independent replication before strong claims.

In short:

```text
COSMOS = SPACE + EARTH + OCEAN + ATMOSPHERE + THEIR COUPLINGS
```

## Ocean / Cousteau domain

The ocean branch lives on `janus-echo-cousteau` and extends JANUS into marine observation and calibration: multibeam, sidescan, bathymetry, hydroacoustics, navigation, sound-speed profiles, seafloor morphology, instrument-failure controls, synthetic replay and blind target discrimination.

The current Hannah/BODC lane uses a deliberately bidirectional analysis:

```text
START → changes → target interval ← changes ← EOF
```

followed by:

```text
TIME convergence
→ exact survey line / timestamp
→ SPACE convergence
→ ping / beam / footprint
→ fixed seafloor coordinate
→ independent repeat-pass test
→ artefact-vs-relief test
→ EA600 H0/H1/H2
→ blind synthetic replay
→ Echo-Pyramid comparison
```

A critical retained hypothesis is `H2_REAL_TERRAIN_TRIGGERED_FAILURE`: a single-beam bottom tracker can genuinely fail because real steep or rough terrain causes loss of bottom. Therefore an instrument failure is neither automatically discarded as noise nor promoted to morphology; it must be compared with independently georeferenced EM122/TOBI morphology and repeated passes.

Bathymetric geometry and acoustic resonance remain separate gates. EM122/TOBI can support morphology. A resonance claim requires suitable raw/complex acoustic data with phase, Q/linewidth, decay/ringdown, multi-aspect response, instrument response and environmental correction.

## OSIRIS / S𓂸ḥ computational domain

Dedicated proof-carrying SAT/complexity superproject for **OSIRIS** and the modern JANUS gate alias **S𓂸ḥ**.

```text
CNF
  → CNF0 structural witness
  → OSIRIS v2 technical SAT mechanics
  → OSIRIS v2.1 / PR192 verified pair-product quotient
  → S𓂸ḥ/0 verified K=1 articulation decomposition
  → S𓂸ḥ/1 verified K=2 decomposition
  → S𓂸ḥ/2 automatic minimum-separator search k=1..4 + preregistered holdout
  → OSIRIS v3 ORIGIN_PRIME spiral compute
  → generic OSIRIS fallback / UNKNOWN_BUDGET where applicable
```

Current state: **S𓂸ḥ/0 PASS · S𓂸ḥ/1 PASS · S𓂸ḥ/2 PASS · OSIRIS v3 spiral compute PASS**.

**P_VS_NP = OPEN.** Neither `P=NP` nor `P!=NP` is established.

## The state model is now a spiral, not a circle

Legacy state intuition:

```text
A → B → C → A
```

Current runtime law:

```text
ORIGIN_n
  → COMPUTE / EXPERIENCE_n
  → VERIFIED_RETURN_n
  → ORIGIN_PRIME_(n+1)
```

The important distinction is:

```text
POSITION may repeat
STATE must advance
RETURN != RESET
```

`ORIGIN_PRIME` is computationally active. Exact formula-bound verified experience from one generation may be reused by the next generation, but memory itself never receives verdict authority.

### Safe reuse implemented in v3

**SAT:** a prior model may skip repeated solving only after the stored assignment is rechecked against the exact current canonical CNF.

**UNSAT:** a prior `UNSAT` label is never reused as a verdict shortcut. A stored separator route may avoid repeated exhaustive root separator discovery only after the separator and component partition are revalidated on the exact current residual CNF; current boundary/component closure still executes.

Frozen K=4 self-test observation:

```text
SAT generation 1 root discovery candidates:   793
SAT generation 2 exact-input candidates:         0  (reverified witness)

UNSAT generation 1 root discovery candidates: 793
UNSAT generation 2 exact-input candidates:       0  (revalidated route; current closure rerun)

negative controls: 8/8 REJECT
```

This is a repeated exact-input reuse result, **not** a general complexity result.

## Clone once, get the whole project

```bash
git clone --recurse-submodules https://github.com/Hawkar-usls/Janus-Cosmos.git
cd Janus-Cosmos
python run_osiris.py
```

If already cloned:

```bash
git submodule update --init --recursive
python run_osiris.py
```

`python run_osiris.py` executes the current spiral self-test.

The inherited frozen S𓂸ḥ/2 baseline remains available:

```bash
python run_osiris.py --legacy
```

To solve a JSON CNF with persistent ORIGIN_PRIME state:

```bash
python run_osiris.py --cnf formula.json --budget 50000 --state ~/.janus/osiris_spiral_state.json
```

Accepted CNF input forms:

```json
[[1, -2], [2, 3]]
```

or:

```json
{
  "formula": [[1, -2], [2, 3]],
  "budget": 50000
}
```

## BHB — Black Hole Boom: relativistic collision research

`Janus-Cosmos` has a **separate general-relativistic kinematics lane** for
near-horizon radial Schwarzschild (R1) and equatorial rotating Kerr (R2)
test-particle collisions. It replays the seven-point BHB R0 benchmark from
JANUS Meta Registry and freezes provenance-aware collision case matrices.

```bash
python3 workspace/janus_bhb.py --replay-upstream
python3 workspace/janus_bhb.py --batch data/bhb/BHB-COSMOS-R1-CASES-2026-10-17.json
python3 workspace/janus_bhb_kerr.py --batch data/bhb/BHB-COSMOS-R2-CASES-2026-10-18.json
python3 workspace/janus_bhb_firewall.py
python3 workspace/janus_bhb_collision_radiation.py --batch data/bhb/BHB-R3A-QED-COLLISION-CASES-2026-10-18.json
python3 -m unittest discover -s tests -p 'test_bhb*.py' -v
```

The scientific gates distinguish invariant local collision energy from a
physically available outgoing trajectory, and both from escaping radiation.
R2 adds Kerr angular momentum, analytic radial-barrier tests for infinity
infall, and hypothetical equatorial photon escape with a ZAMO redshift
factor. **It does not compute photons actually made in a collision.** This is
not measured black-hole emission or an explosion. See the
[BHB operating guide](docs/BHB_BLACK_HOLE_BOOM.md).

The additional **FW1 AMPS quantum-firewall hypothesis** is registered and
checked separately via `workspace/janus_bhb_firewall.py`; it audits the AMPS
information-paradox premises, replays the frozen R2 receipt, and deliberately
produces **no physical firewall spectrum, collision-energy prediction or
observation claim**. It does not modify the classical Kerr model.

**R3A** adds an explicitly specified *conditional* collision-product
benchmark: infalling electron plus externally supported static positron
annihilating into two radially opposed photons outside a Schwarzschild
horizon. Four-momentum and Killing energy are conserved; the outward product
is checked with the R2 null-geodesic kernel. The static target's required
support acceleration diverges in the horizon limit. Neither a natural
positron shell, an AMPS firewall, an event rate nor an observed gamma-ray
flash is inferred from this synthetic reaction.

## Video-inspired Cosmos V1: ALI simulation hypothesis audit

The [ALI video audit guide](docs/ALI_SIMULATION_HYPOTHESIS_AUDIT.md)
maps 18 time-coded topics from a **third-party summary** of
`Что-то странное происходит с миром | ALI` into six auditable models:
Lorentz kinematics, Schwarzschild horizon kinematics, area-law entropy,
compact flat spatial topology, conditional lattice-wave dispersion and
conditional observer counting. Direct video transcript access was blocked;
the summaries are **not verified quotes**. None proves a simulated universe.

```bash
python3 workspace/janus_cosmos_simulation_audit.py --batch data/cosmos/COSMOS-ALI-SIMULATION-CLAIMS-2026-10-18.json
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

This additive lane locks the existing BHB R3A SHA-256 receipt, checks 16
cases, retains 18 claim provenance records and refuses to convert an AMPS
firewall, a discretization toy or philosophical counting into astrophysical
evidence or BHB radiation. The exact output is frozen under `data/cosmos/`.

## Cosmos V2: falsifiable photon-lattice time-of-flight benchmark

The [V2 photon TOF experiment](docs/COSMOS_LATTICE_TOF_V2.md) defines a
specific nearest-neighbor lattice dispersion relation and compares its
energy-dependent photon arrival delays against Lorentz-invariant propagation
in an illustrative expanding universe. Twelve **synthetic** holdout cases
cover source-intrinsic lag degeneracy, false-positive controls, a near-Planck
spacing limit, timing uncertainties and zero-baseline comparisons.

```bash
python3 workspace/janus_cosmos_lattice_tof.py --score \
  --protocol data/cosmos/COSMOS-ALI-V2-LATTICE-TOF-PROTOCOL-2026-10-18.json \
  --observations data/cosmos/COSMOS-ALI-V2-LATTICE-TOF-OBSERVATIONS-2026-10-18.json
python3 -m unittest discover -s tests -p 'test_cosmos_lattice_tof.py' -v
```

The blind score alone cannot reveal truth injection labels; separate `--prepare`
and `--reveal` phases reproduce and audit them. **No real GRB photons, observed
Lorentz violation, spacetime lattice, simulation evidence, black-hole quantum
firewall, or BHB explosion is claimed.**

## Repository layout

- `vendor/Janus-Fundamentum/` — exact inherited project/source snapshot at the S𓂸ḥ/2 head.
- `vendor/janus-meta-registry/` — exact inherited receipt/source registry snapshot.
- `experiments/direct` — convenience symlink into the inherited canonical experiment directory.
- `workspace/` — native OSIRIS v3+ development and spiral runtime.
- `workspace/janus_bhb.py` — independent relativistic Schwarzschild collision engine.
- `workspace/janus_bhb_kerr.py` — Kerr equatorial collision and null-geodesic transport engine.
- `workspace/janus_bhb_firewall.py` — conditional quantum-firewall hypothesis audit over pinned R2.
- `workspace/janus_bhb_collision_radiation.py` — R3A QED pair-annihilation products, conserved energy and radial escape.
- `workspace/janus_cosmos_simulation_audit.py` — video-inspired six-model physics / interpretation gate.
- `workspace/janus_cosmos_lattice_tof.py` — V2 fixed photon-dispersion model, process-blind synthetic scoring and label reveal.
- `data/bhb/` — pinned external BHB reference and frozen Cosmos experiment receipts.
- `data/cosmos/` — versioned ALI topic matrix and deterministic physics-audit result.
- `tests/test_bhb_cosmos.py` — near-horizon replay and physical-provenance negative controls.
- `tests/test_bhb_kerr.py` — Kerr radial-barrier and photon escape negative controls.
- `tests/test_bhb_firewall.py` — AMPS premises, unchanged-source and non-promotion controls.
- `tests/test_bhb_collision_radiation.py` — R3A local/infinity conservation, escape and support negative controls.
- `tests/test_cosmos_simulation_audit.py` — Lorentz, Schwarzschild, entropy, topology, lattice and provenance gates.
- `tests/test_cosmos_lattice_tof.py` — 12-case TOF replay, source-lag degeneracy and false-positive controls.
- `receipts` — convenience symlink into the pinned meta-registry data directory.
- `data/cousteau/` — marine/ocean research artifacts on the `janus-echo-cousteau` branch.
- `docs/` — architecture, provenance, lineage, scope and scientific boundaries.
- `archive/` — pointer to the preserved pre-OSIRIS astronomy repository state.
- `run_osiris.py` — single canonical OSIRIS entrypoint.

## Historical/source firewall

`S𓂸ḥ` is a modern project alias/overlay. It is not claimed to be an ancient Egyptian spelling or SAT terminology. Historical Pyramid Text material was heuristic semantic/operator inspiration only. Correctness and authority come from executable CNF mechanics, exact replay, current verification and certificates.

The frozen spiral semantic specification explicitly keeps:

```text
ANCIENT_TEXT != MODERN_ALGORITHM
STRUCTURAL_PARALLEL != HISTORICAL_INTENT
PT_NUMBER_ORDER != PHYSICAL_WALL_ORDER
```

The same firewall applies across Janus-Cosmos: an unusual sonar, astronomical or geophysical pattern is not promoted to an exotic explanation until ordinary instrument, propagation, biological, geological and processing controls have been tested.

## Current canonical source pins

- `Janus-Fundamentum@00b09d778a6d57fc7f905df9b6235fb30e29c5a3`
- `janus-meta-registry@be1c7538932a823907b805516efaca035029ea8b`
- inherited S𓂸ḥ/2 result integrity: `9dc0748a974b890a934192144d6c806d2aa455ac4b96d744ff59e3bae598e94d`
- OSIRIS v3 spiral result integrity: `8324d217751d3d087973447c7af0990d0da654b9c571864832c705ba355f9a87`

The previous HST/MAST/Palomar/LUCI Janus-Cosmos is preserved intact on branch `legacy-cosmos-pre-osiris-2026-08-19`.
