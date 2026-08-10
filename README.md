<div align="center">

# JANUS COSMOS × FRACTALGPT
### Blind multiscale anomaly search in public astronomical data

![Status](https://img.shields.io/badge/status-research%20pilot-1f6feb)
![Method](https://img.shields.io/badge/method-blind%20%2B%20null--tested-8957e5)
![Semantics](https://img.shields.io/badge/semantic%20search-OFF-6e7681)

</div>

## Mission

Janus Cosmos is an experimental research repository for applying the project's FractalGPT blind-planning methodology to public astronomical observations.

The first target is **Hubble Space Telescope multi-filter imagery**, with MAST as the authoritative archive and JWST as a later confirmation corpus.

The goal is not to search images for faces, words, messages, or visually suggestive shapes. The primary gate asks a stricter question:

> **Can a multiscale spatial structure be detected reproducibly, survive matched null models, persist across relevant observations/bands, and replicate independently?**

## Pipeline

```text
Public HST data / MAST
        ↓
source-provenance freeze
        ↓
image normalization
        ↓
FractalGPT blind multiscale proposals
        ↓
spatial null models
        ↓
candidate ranking
        ↓
cross-band persistence
        ↓
independent replication
        ↓
scientific review
```

## FractalGPT boundary

FractalGPT is a **blind planner**, not an oracle or discovery claim generator.

During the primary anomaly gate:

- OCR is disabled.
- Face detection/search is disabled.
- Cipher/message search is disabled.
- Semantic interpretation is disabled.
- Post-hoc parameter retuning is forbidden.
- A single visually striking image is never sufficient.

A candidate is interesting only when it survives the pre-registered statistical and replication gates.

## Expected discoveries

The prior expectation is ordinary astrophysics and imaging structure, including:

- galaxy mergers and tidal tails;
- shells, arcs and rings;
- gravitational-lensing morphologies;
- jets and cavities;
- faint filaments and asymmetric structures;
- instrumental or processing artifacts;
- rare but astrophysically ordinary morphology.

A genuinely interesting result would be a reproducible, previously under-characterized structure that survives the controls.

## First corpus

The project uses a small, source-confirmed HST pilot before scaling to larger archives. Source identifiers, filters, observation metadata and provenance are preserved in machine-readable manifests and receipts.

## Repository layout

```text
experiments/
  fractalgpt/
    run_search.py
  data/
    hst_source_manifest.json
schemas/
  candidate_receipt.schema.json
.github/workflows/
  janus-cosmos.yml
```

## Claim ceiling

This project can establish an image/data-level anomaly candidate. It cannot, by itself, establish extraterrestrial intelligence, hidden communication, unknown physics, or any semantic interpretation.

## Status

**RESEARCH PILOT — NO ASTRONOMICAL DISCOVERY CLAIMED.**

The repository is intentionally conservative: failed hypotheses, null results, and methodological corrections are retained as part of the experiment's provenance.
