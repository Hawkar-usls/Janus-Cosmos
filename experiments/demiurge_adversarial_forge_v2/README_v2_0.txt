JANUS COSMOS v2.0 — DEMIURGE-ADVERSARIAL DETECTOR FORGE
========================================================

PATCH RELEASE v2.0.2
--------------------
v2.0.1 still failed on Windows with Python 3.14, NumPy 2.5.2 and SciPy 1.18.0:
rounding raw synthetic metrics to 12 decimal places was not a sound definition
of detector identity.

v2.0.2 separates two proof layers:

  PORTABLE DETECTOR IDENTITY
    genome + normalized source hashes + manifests + validation decision

  PLATFORM NUMERICAL EVIDENCE
    exact raw metrics + exact metrics_sha256 + a source-bound 2e-6 tolerance gate

The scientific validation gate continues to use raw metrics. The portable
freeze never hides the numerical evidence; it gives that evidence its own exact
hash and checks it against a narrow registered conformance envelope. Material
metric drift, genome drift, source drift, manifest drift or a failed validation
decision still blocks the sky run.

MISSION
-------
The task is to CHECK THE COSMOS, not to manufacture a spectacular p-value.

v2.0 therefore inserts a hard scientific firewall between detector design and
astronomical evaluation:

  SYNTHETIC ADVERSARIAL FORGE
        -> validate against held-out synthetic controls
        -> FREEZE genome + SHA-256
        -> BLIND WALL
        -> unrelated real-sky specificity controls
        -> Orion DSS2 + 2MASS validation
        -> NGC1425 HST validation

The detector may not adapt after the freeze.

WHY THIS VERSION EXISTS
-----------------------
v1.6 reproduced the Orion response in DSS2 and independently in 2MASS, and the
DSS upper-right polygon was not sufficient to explain the response. That was an
important result, but it still left a dangerous alternative: perhaps the JANUS
morphology score marks many ordinary sky fields as extreme.

v2.0 attacks exactly that possibility.

DEMIURGE BRIDGE
---------------
Janus-Demiurge is used as an architectural ancestor, not as scientific evidence.
The bridge reuses only bounded optimization ideas:

- architecture/genome mutation;
- population/evolutionary selection;
- paired counterfactual evaluation;
- immutable experiment ledger;
- freeze-before-target evaluation.

Legacy filter_37, digital-root, "tachyonic", resonance, mood and similar
heuristics are explicitly forbidden from the Cosmos gate.

FORGE FITNESS
-------------
The forge never sees Orion, NGC1425 or the blind sky controls.

It receives synthetic sky-like base fields plus paired interventions:

  SIGNAL intervention:
    coherent arcs, symmetric nodes, filaments or nested ellipses

  ARTIFACT intervention:
    polygon/no-data wedge, plate step, saturation cross, edge truncation,
    mosaic gain blocks or a broad halo

Only the six feature-group weights are allowed to evolve. The preprocessing,
feature definitions, phase/IAAFT null, block-shuffle null and their core
parameters stay fixed.

The selected genome must improve the held-out signal/artifact sensitivity ratio
relative to the canonical equal-weight detector, reduce artifact sensitivity,
and retain planted-signal sensitivity. Validation is pass/fail; target data are
not used to rank or repair the genome.

PACKAGED FROZEN DETECTOR
------------------------
The deterministic forge packaged with this release produced:

  genome_sha256:
    ebec8aaf0b623e8f805615130362570738ddacf4c5f4a90cc347652209fd20d6

  freeze_sha256:
    55bdf9f19820623a392b6b7a859408d3beffb68f28dd4a32001bb1903025403a

The runtime re-forges the detector and refuses the portable identity check if
the detector changes. It separately verifies the exact platform metric receipt
and the registered numerical conformance envelope.

Held-out synthetic validation improved versus the equal-weight canonical detector:

  signal/artifact sensitivity ratio: 0.35108 -> 0.42829  (+21.99%)
  paired artifact delta:              1.79144 -> 1.68276  (-6.07%)
  paired planted-signal delta:        0.62894 -> 0.72071  (+14.59%)

These are forge validation metrics only. They are NOT astronomical evidence.

REAL-SKY SPECIFICITY GATE
-------------------------
Four unrelated sky fields are fixed before the run from one SHA-256 seed. Their
coordinates are stored in SKY_MANIFEST_v2_0.json.

Each field is downloaded in:

  DSS2 Red
  DSS2 Blue
  2MASS J
  2MASS K

Each field is tested in two geometries:

  WHOLE
  PSEUDO_BELT_CORRIDOR

The control family is corrected separately for whole-field and corridor tests.
A control field is a false positive only when a robust result spans BOTH survey
families (at least one DSS2 band and at least one 2MASS band, with both primary
null models passing).

The gate allows at most:

  1 / 4 whole-field false-positive fields
  1 / 4 corridor false-positive fields

If this gate fails, Orion and NGC1425 target claims are BLOCKED even if their
individual p-values look impressive.

ORION VALIDATION
----------------
Orion is tested in four bands / two independent survey families:

  DSS2 Red
  DSS2 Blue
  2MASS J
  2MASS K

Pre-registered variants:

  WHOLE
  BELT_CORRIDOR

Family size:

  4 bands x 2 variants x 2 null models = 16 tests
  alpha = 0.05 / 16 = 0.003125

The strong status SKY_FIXED_MORPHOLOGY_CANDIDATE requires:

1. the unrelated real-sky specificity gate to PASS;
2. the Belt corridor to be robust in both DSS2 bands;
3. the Belt corridor to be robust in both 2MASS bands.

This is still a validation of a known target, NOT a blind astronomical discovery.

NGC1425 VALIDATION
------------------
The package downloads the official HST/WFPC2 F555W and F814W mosaics.

The historical parent family correction is preserved:

  alpha = 0.05 / 28 = 0.0017857142857142859

The HST candidate is admitted only if both filters pass both primary null models
AND the real-sky specificity gate passes.

MONTE CARLO
-----------
Full run:

  test nulls/model = 768
  calibration nulls/model = 96
  seed chunks = 3

Minimum empirical p = 1/769 ~= 0.00130039, which can resolve every pre-registered
alpha used in this package.

PRIMARY NULLS
-------------
- morphology-preserving phase/IAAFT
- local block shuffle

No OCR, face search, semantic analysis, cipher search or post-hoc target tuning
enters the statistical gate.

QUICK START
-----------
Double-click:

  run_janus_cosmos_v2_0.bat

The script will:

  [1] check/install dependencies
  [2] run offline scientific/self-consistency tests
  [3] deterministically re-forge and verify the frozen detector
  [4] download 22 astronomy FITS products
  [5] run the frozen real-sky experiment

Checkpoints are saved per field / variant / null model. If the long analysis is
interrupted, run the BAT again and valid completed checkpoints are reused.

EXPECTED RUNTIME / DOWNLOAD
---------------------------
Runtime depends heavily on CPU. The powered analysis performs roughly tens of
thousands of 128x128 surrogate evaluations, so a desktop run can take many
minutes. HST NGC1425 mosaics are the largest downloads; HiPS control/Orion files
are smaller.

FILES TO SEND BACK
------------------
After the full run, send:

  results_v2_0\janus-cosmos-v2.0-report.json
  results_v2_0\janus-cosmos-v2.0-events.jsonl
  results_v2_0\terminal.log
  results_v2_0\SUMMARY_v2.0.txt

CLAIM CEILING
-------------
The strongest allowed result is:

  IMAGE-LEVEL SKY-FIXED MORPHOLOGY CANDIDATE

This package does NOT establish an artificial structure, hidden message,
censorship/concealment, alien technology, or new physics.
