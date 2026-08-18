# JANUS COSMOS — OSIRIS / S𓂸ḥ

Dedicated proof-carrying SAT/complexity superproject for **OSIRIS** and the modern JANUS gate alias **S𓂸ḥ**.

```text
CNF
  → CNF0 structural witness
  → OSIRIS v2 technical SAT mechanics
  → OSIRIS v2.1 / PR192 verified pair-product quotient
  → S𓂸ḥ/0 verified K=1 articulation decomposition
  → S𓂸ḥ/1 verified K=2 decomposition + ORIGIN_PRIME ribbon state
  → S𓂸ḥ/2 automatic minimum-separator search k=1..4 + preregistered holdout
  → generic OSIRIS fallback / UNKNOWN_BUDGET
```

Current frozen finite-gate state: **S𓂸ḥ/0 PASS · S𓂸ḥ/1 PASS · S𓂸ḥ/2 PASS**.

**P_VS_NP = OPEN.** Neither `P=NP` nor `P!=NP` is established.

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

## Repository layout

- `vendor/Janus-Fundamentum/` — exact project/source snapshot at the S𓂸ḥ/2 head.
- `vendor/janus-meta-registry/` — exact receipt/source registry snapshot.
- `experiments/direct` — convenience symlink into the canonical runnable experiment directory.
- `receipts` — convenience symlink into the complete meta-registry data directory.
- `docs/` — architecture, provenance, lineage and scientific boundaries.
- `archive/` — pointer to the preserved pre-OSIRIS astronomy repository state.
- `run_osiris.py` — single canonical entrypoint.

## Historical/source firewall

`S𓂸ḥ` is a modern project alias/overlay. It is not claimed to be an ancient Egyptian spelling or SAT terminology. Historical Pyramid Text material was heuristic operator inspiration only. Correctness and authority come from executable CNF mechanics, exact replay and certificates.

## Current canonical source pins

- `Janus-Fundamentum@00b09d778a6d57fc7f905df9b6235fb30e29c5a3`
- `janus-meta-registry@be1c7538932a823907b805516efaca035029ea8b`
- S𓂸ḥ/2 result integrity: `9dc0748a974b890a934192144d6c806d2aa455ac4b96d744ff59e3bae598e94d`

The previous HST/MAST/Palomar/LUCI Janus-Cosmos is preserved intact on branch `legacy-cosmos-pre-osiris-2026-08-19`.
