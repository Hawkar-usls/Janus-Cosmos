# Hannah/BODC bidirectional convergence

This Cousteau-lane protocol links the BODC JR15001 and CD169 data-access work to the existing Janus-Echo-Cousteau and Echo-Pyramid calibration architecture.

The analysis moves from both temporal edges toward the target interval:

**START → changes → target interval ← changes ← EOF**

Then it changes domains:

**TIME convergence → survey line/time → SPACE convergence → ping/beam/footprint → fixed seafloor coordinate → repeat-pass test → artefact-vs-relief test → EA600 H0/H1/H2 → blind synthetic replay → Echo-Pyramid comparison.**

The critical H2 hypothesis is retained explicitly: a single-beam bottom tracker may genuinely fail because real steep or rough terrain causes loss of bottom. Therefore an EA600 failure is neither automatically discarded as noise nor promoted to morphology. It must be tested against independently georeferenced EM122 or sidescan morphology and repeated passes.

Known JR15001 smile-like anomaly lines `0012`, `0017`, and `0022` are treated as candidates, while `0037` is a known bad-parameter control. A candidate is considered ground-fixed only if its Earth-coordinate location persists under changes in heading, beam geometry, pass, and reasonable sound-speed corrections.

Bathymetric geometry and acoustic resonance are separate gates. EM122/TOBI can support morphology; resonance requires suitable raw/complex acoustic measurements with phase, Q/linewidth, decay/ringdown, multi-aspect response, instrument response, and environmental correction.

Public acknowledgement for the data-access contribution remains:

**British Oceanographic Data Centre (BODC) Requests Team — primary-data access and provenance assistance.**

Internal codename only: `HANNAH_MONTANA__ДЕВОЧКА_СУПЕРЗВЕЗДА`.

Canonical protocol JSON:

`data/cousteau/JANUS-HANNAH-BODC-BIDIRECTIONAL-CONVERGENCE-TO-CENTER-2026-08-25-v1.0.json`

## Synesthetic mnemonic sidecar

The Hannah lane now has a deliberately non-authoritative **Cousteau Synesthetic Memory Core**. It gives each measured multiscale episode a deterministic cross-modal sensory passport so HEAD and TAIL windows can be recognized and retrieved consistently before ordinary TIME → SPACE verification.

```text
MEASUREMENT
-> stable fingerprint
-> color / tone / rhythm / texture / glyph
-> Cousteau semantic overlay
-> retrieval / review priority

MNEMONIC MATCH != SCIENTIFIC CONVERGENCE
```

The sidecar is blind to target/verdict/template concepts when computing its measurement fingerprint. `H0/H1/H2`, `pyramid`, `target`, `candidate`, `anomaly`, expected labels and story/interpretation fields are explicitly forbidden from influencing the fingerprint.

Cousteau-specific mnemonic semantics include deeper measured depth → lower register, explicit EM122−EA600 disagreement → stronger beat frequency, cadence jitter → pulse instability, roughness-related measurements → mnemonic texture, and missingness → fog rather than synthetic filling. HEAD/Tail direction is only an overlay and cannot alter identical measured fingerprints.

Current state remains intentionally silent because the original BODC bytes are not mounted: `BLOCKED_NULL / #808080 / SILENCE / fog / measurement_fingerprint=null`.

Documentation:

`data/cousteau/JANUS-HANNAH-COUSTEAU-SYNESTHETIC-MEMORY-README-2026-08-25.md`

Frozen contract:

`data/cousteau/JANUS-HANNAH-COUSTEAU-SYNESTHETIC-MEMORY-CORE-CONTRACT-2026-08-25-v1.0.json`

Implementation:

`workspace/cousteau_synesthetic_memory_core.py`

`workspace/cousteau_synesthetic_semantic_overlay.py`
