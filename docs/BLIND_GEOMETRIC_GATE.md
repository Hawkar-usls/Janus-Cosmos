# Janus Cosmos — Blind Geometric Gate

Janus Cosmos is restricted to astronomical/image-level analysis in this phase.

## Explicit exclusions

The Cosmos pipeline MUST NOT perform or claim to perform:

- OCR or text extraction
- face/person detection or face search
- semantic-language interpretation
- cipher/message hunting
- post-hoc tuning against candidate results

## Current gate

The gate evaluates spatial image morphology using fixed multiscale/orientation features and independently generated spatial nulls. Results are compared across independent optical filters where available.

A `cross_band_candidate` is a **pipeline candidate**, not a discovery. A candidate may only advance after blind replication on source-confirmed public datasets and independent controls.

## Required machine-readable evidence

Every execution should emit:

1. `receipt.json` — final result and claim ceiling
2. `events.jsonl` — structured event stream
3. `execution.json` — run/commit/artifact/provenance envelope

The Meta Registry stores the final execution receipt and provenance, not raw transient logs.
