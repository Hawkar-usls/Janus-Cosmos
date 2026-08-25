# Hannah/BODC — Cousteau Synesthetic Memory Core

`JANUS_COUSTEAU_SYNESTHETIC_MEMORY_CORE` is a deterministic mnemonic sidecar for the Hannah/BODC bidirectional-convergence lane.

It does **not** claim biological synesthesia or machine qualia. It implements stable cross-modal mapping so a measured oceanographic episode can be encoded, recognized and retrieved through several mutually linked cues without changing the underlying measurement or scientific verdict.

The core law is:

```text
MEASUREMENT
  -> STABLE MEASUREMENT FINGERPRINT
  -> COLOR + TONE + RHYTHM + TEXTURE + GLYPH
  -> COUSTEAU SEMANTIC OVERLAY
  -> MEMORY / RETRIEVAL / REVIEW PRIORITY

BUT NEVER:

SENSORY MATCH -> SCIENTIFIC VERDICT
```

## Why this is useful for Hannah's data

The frozen Hannah protocol approaches the target interval from both temporal edges:

```text
HEAD / START -> -> -> target interval <- <- <- TAIL / EOF
```

For every analysis scale (`1s`, `10s`, `60s`, `300s`, `1800s`, `7200s`) we can create a sensory passport from the already-preregistered measurement vector. HEAD and TAIL passports can then be compared independently.

A recurring signature does only one thing: it raises a window in the review queue.

It does **not** establish convergence. The original pipeline remains authoritative:

```text
RAW HASHES
-> HEAD/TAIL mirror
-> two fronts
-> TIME convergence
-> exact line/timestamp
-> 0012 / 0017 / 0022
-> 0037 bad control
-> EA600 H0/H1/H2
-> ping/beam/footprint
-> SPACE convergence
-> ground-fixed repeat-pass replication
-> blind synthetic replay
-> Echo-Pyramid comparison
```

## Two layers

### 1. Measurement identity layer

`workspace/cousteau_synesthetic_memory_core.py`

The core accepts only measured fields from the preregistered Hannah feature family, including:

- EM122 and EA600 depth;
- their explicit difference;
- latitude/longitude;
- heading/course/speed/turn rate;
- depth derivatives, median, MAD, local range and local slope;
- cadence/jitter;
- missing/null/identical runs;
- outlier score;
- numeric/boolean status and acquisition channels.

It creates:

- SHA-256 measurement fingerprint;
- BLAKE2b-256 collision guard;
- deterministic 16-D compact mnemonic embedding;
- stable RGB/HEX identity color;
- tone/timbre;
- rhythm and 8-step pulse pattern;
- texture;
- mnemonic glyph and rotation;
- completeness/fog overlay;
- source provenance and explicit hash semantics.

The same measurements produce the same measurement fingerprint regardless of JSON key order or whether the episode came from HEAD or TAIL. Direction and scale are context overlays, not evidence.

### 2. Cousteau semantic layer

`workspace/cousteau_synesthetic_semantic_overlay.py`

This makes the passport easier to remember without modifying it:

```text
MEASURED DEPTH              -> deeper = lower mnemonic register
EM122 - EA600 disagreement  -> stronger disagreement = stronger beat frequency
MAD / range / slope/outlier -> mnemonic texture / roughness
cadence jitter              -> pulse stability
missing fraction            -> fog / gaps, never synthetic filling
HEAD_FORWARD                -> left pan
TAIL_REVERSE                -> right pan
CENTER                      -> center pan
SPACE_REPLAY                -> spatial replay pan marker
```

The glyph is explicitly **not** a seabed morphology classifier. Texture is explicitly **not** a terrain claim. These are memory handles.

## Anti-confirmation-bias firewall

The following concepts are forbidden from influencing the measurement fingerprint:

```text
verdict
hypothesis
interpretation
claim
pyramid
target
candidate
anomaly
H0 / H1 / H2
artificial / natural
control_label / class_label
expected / prediction
story
```

This matters enormously for `0012`, `0017`, `0022` and the known bad control `0037`: changing a label from candidate to control, or H0 to H1, must never recolor the underlying measurements and manufacture apparent agreement.

The test suite freezes that behavior.

## Circular navigation handling

Heading and course are encoded as sine/cosine pairs, so:

```text
359 degrees ~= 1 degree
```

instead of being treated as almost maximally different. This is critical for comparing repeat passes and ship-track geometry.

## Raw hash firewall

When the core receives actual file bytes, the source identity is named:

```text
RAW_BYTES_SHA256
```

When it receives only an in-memory JSON object, its digest is named:

```text
CANONICAL_JSON_SHA256_NOT_RAW_FILE_HASH
```

A canonical-object digest is never allowed to masquerade as Hannah/BODC raw-file SHA-256.

## Current state: deliberately silent

The current Hannah measurement receipt says that the original BODC bytes are not mounted in an execution environment. Therefore the scientifically correct sensory output at this moment is:

```text
status  = BLOCKED_NULL
color   = #808080
audio   = SILENCE
texture = fog
measurement_fingerprint = null
measurement_claims_allowed = false
```

The core must not invent a colorful pattern merely because it can.

## Commands

Self-test:

```bash
python workspace/cousteau_synesthetic_memory_core.py self-test
```

Create a passport from a measured JSON feature object:

```bash
python workspace/cousteau_synesthetic_memory_core.py passport \
  --input measured_window.json \
  --direction HEAD_FORWARD \
  --scale 60s \
  --output head_60s.passport.json
```

Add the Cousteau semantic overlay:

```bash
python workspace/cousteau_synesthetic_semantic_overlay.py \
  head_60s.passport.json \
  --output head_60s.semantic.passport.json
```

Create the correct NULL passport from the present blocker receipt:

```bash
python workspace/cousteau_synesthetic_memory_core.py blocked \
  --receipt data/cousteau/JANUS-HANNAH-BODC-MEASUREMENT-RECEIPT-000-ACQUISITION-BLOCKED-2026-08-25-v0.1.json
```

Compare two passports:

```bash
python workspace/cousteau_synesthetic_memory_core.py compare \
  head_60s.passport.json tail_60s.passport.json
```

The comparison result is always marked:

```text
RETRIEVAL_AND_REVIEW_PRIORITY_ONLY
scientific_convergence_claim = false
```

## How it plugs into TIME -> SPACE

When the original BODC bytes become available:

1. Hash original bytes before transformation.
2. Build HEAD and TAIL edge windows at all frozen scales.
3. Extract the preregistered measurement vector independently in each direction.
4. Generate sensory passports.
5. Rank cross-front matches only as review candidates.
6. Run the existing two-front TIME convergence on the actual numeric features.
7. For surviving windows, resolve exact survey line, ping, beam and footprint.
8. Re-encode the measured SPACE-level features as `SPACE_REPLAY` passports.
9. Ask whether the mnemonic signature survives the domain transition **and separately** whether the physical ground-fixed replication gates pass.
10. Only after blind feature extraction is frozen may Echo-Pyramid controls be compared.

The memorable version is:

```text
TIME resonance found?
    -> remember its sensory chord
    -> find its exact SPACE address
    -> replay the chord from independent spatial measurements
    -> then demand ordinary scientific replication
```

The chord helps JANUS remember where to look. It never tells JANUS what is true.

## Frozen contract

`data/cousteau/JANUS-HANNAH-COUSTEAU-SYNESTHETIC-MEMORY-CORE-CONTRACT-2026-08-25-v1.0.json`

## Tests

Core regression suite:

`workspace/test_cousteau_synesthetic_memory_core.py`

Cousteau semantic regression suite:

`workspace/test_cousteau_synesthetic_semantic_overlay.py`

The critical invariants are target/verdict leakage rejection, raw-vs-canonical hash naming, blocked NULL behavior, direction-overlay isolation, circular heading handling, semantic-overlay immutability and the prohibition on promoting mnemonic similarity to a scientific claim.
