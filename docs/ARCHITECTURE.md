# Architecture

Canonical inherited solver priority:

1. `CNF0_STRUCTURAL_WITNESS`
2. `PR192_VERIFIED_PAIR_PRODUCT_PREBIRTH_QUOTIENT`
3. `S𓂸ḥ/2` automatic verified minimum separator search for `k=1..4`
4. generic OSIRIS v2 exact residual engine or fail-closed `UNKNOWN_BUDGET`

S𓂸ḥ/2 exhausts all candidate separators below and at the first successful `k`, proves minimality in its frozen finite scope, evaluates all `2^k` boundary rows and all components, and reports discovery and boundary work separately.

## OSIRIS v3: active spiral state

The state model is now computationally active:

```text
ORIGIN_n
  -> COMPUTE_EXPERIENCE_n
  -> VERIFIED_RETURN_n
  -> ORIGIN_PRIME_(n+1)
```

A repeated technical POSITION may be identical while the state lineage advances:

```text
POSITION_(n+1) == POSITION_n   may hold
STATE_(n+1) != STATE_n         must hold
RETURN != RESET
```

The previous ribbon implementation proved that history could be bound into a new state. OSIRIS v3 additionally allows exact formula-bound verified experience to influence the next search.

### Experience classes

`SAT_WITNESS`:
- exact formula/budget/provider match required;
- experience commitment must verify;
- stored assignment must satisfy the current canonical CNF again;
- only then may the next traversal bypass repeated solving.

`SEPARATOR_ROUTE`:
- exact formula/residual/provider match required;
- separator and claimed component partition are revalidated on the current residual CNF;
- a valid route may bypass repeated exhaustive root minimum-separator discovery;
- fresh minimality is not claimed on the reuse traversal;
- current boundary/component solve still executes.

`UNSAT_MEMORY`:
- never a verdict shortcut in v3;
- current closure is mandatory unless a future separately frozen contract introduces an independently checkable reusable UNSAT proof object.

### Authority firewall

```text
MEMORY != TRUTH
HISTORY != CERTIFICATE
EXPERIENCE_REUSE != VERDICT_AUTHORITY
CACHED_UNSAT != CURRENT_UNSAT_PROOF
```

Experience may change search/repeated-discovery cost. It may not change the meaning of SAT/UNSAT or promote `UNKNOWN_BUDGET`.

## Frozen v3 observation

On the preregistered K=4 repeated-input controls:

```text
SAT:   root discovery 793 -> 0 on generation 2
UNSAT: root discovery 793 -> 0 on generation 2
```

SAT generation 2 reverified the stored witness against the exact current CNF. UNSAT generation 2 revalidated the stored separator route and reran current closure. Eight tamper/mismatch controls rejected.

These counts are not combined with heterogeneous boundary/verification work and do not establish a general complexity improvement.

## Source and development boundary

The exact inherited solver code is pinned under `vendor/Janus-Fundamentum`. The active v3 runtime, semantic specification and frozen contract live natively under `workspace/`.

Historical Pyramid Text semantics are an inspiration layer only:

```text
ANCIENT_TEXT != MODERN_ALGORITHM
STRUCTURAL_PARALLEL != HISTORICAL_INTENT
PT_NUMBER_ORDER != PHYSICAL_WALL_ORDER
```

`P_VS_NP` remains `OPEN`.
