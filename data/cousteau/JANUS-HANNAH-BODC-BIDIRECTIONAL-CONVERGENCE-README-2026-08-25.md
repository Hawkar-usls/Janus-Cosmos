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
