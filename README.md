# JANUS COSMOS — OSIRIS / S𓂸ḥ

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

## Repository layout

- `vendor/Janus-Fundamentum/` — exact inherited project/source snapshot at the S𓂸ḥ/2 head.
- `vendor/janus-meta-registry/` — exact inherited receipt/source registry snapshot.
- `experiments/direct` — convenience symlink into the inherited canonical experiment directory.
- `workspace/` — native OSIRIS v3+ development and spiral runtime.
- `receipts` — convenience symlink into the pinned meta-registry data directory.
- `docs/` — architecture, provenance, lineage and scientific boundaries.
- `archive/` — pointer to the preserved pre-OSIRIS astronomy repository state.
- `run_osiris.py` — single canonical entrypoint.

## Historical/source firewall

`S𓂸ḥ` is a modern project alias/overlay. It is not claimed to be an ancient Egyptian spelling or SAT terminology. Historical Pyramid Text material was heuristic semantic/operator inspiration only. Correctness and authority come from executable CNF mechanics, exact replay, current verification and certificates.

The frozen spiral semantic specification explicitly keeps:

```text
ANCIENT_TEXT != MODERN_ALGORITHM
STRUCTURAL_PARALLEL != HISTORICAL_INTENT
PT_NUMBER_ORDER != PHYSICAL_WALL_ORDER
```

## Current canonical source pins

- `Janus-Fundamentum@00b09d778a6d57fc7f905df9b6235fb30e29c5a3`
- `janus-meta-registry@be1c7538932a823907b805516efaca035029ea8b`
- inherited S𓂸ḥ/2 result integrity: `9dc0748a974b890a934192144d6c806d2aa455ac4b96d744ff59e3bae598e94d`
- OSIRIS v3 spiral result integrity: `8324d217751d3d087973447c7af0990d0da654b9c571864832c705ba355f9a87`

The previous HST/MAST/Palomar/LUCI Janus-Cosmos is preserved intact on branch `legacy-cosmos-pre-osiris-2026-08-19`.
