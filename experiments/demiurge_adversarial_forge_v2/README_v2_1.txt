JANUS COSMOS v2.1.1 — PARALLEL DETECTOR SPECIFICITY REPAIR
==========================================================

PURPOSE
-------
v2.0.2 completed correctly but failed detector specificity: ordinary blind sky
fields repeatedly saturated the synthetic-null p-value floor. That run is now
preserved under evidence/v2_0_2_specificity_negative as a hash-bound
VALID_NEGATIVE_DETECTOR_SPECIFICITY certificate.

v2.1 does not retrain or tune the detector after seeing Orion or NGC1425. The
v2.0.2 genome and portable freeze identity remain unchanged. v2.1 freezes a new
admission protocol that asks whether a target is more extreme than comparable
real observations.

v2.1.1 is an execution-only repair. It does not change the v2.1.0 scientific
protocol, detector, null models, null counts, seeds, control cohorts, scores or
admission rules. It schedules up to ten independent fields in Windows-safe
spawned processes and then restores the frozen protocol order before writing
the report.

WHAT CHANGED
------------
1. Orion is ranked against 20 fresh deterministic DSS2/2MASS sky fields. The
   four already observed v2.0.2 controls are excluded from admission.
2. Every corridor is calibrated against 255 deterministic random positions and
   orientations inside the same image.
3. Whole-image and corridor morphology must agree across DSS2 and 2MASS after
   rank photometric standardization, Gaussian PSF matching and 64x64 resolution
   matching. Agreement itself is ranked against the 20 real sky controls.
4. NGC1425 no longer uses the WFPC2 mosaic silhouette. F555/F814 are read from
   the precommitted WF3 chip, cropped to common positive weight-map support and
   compared with 20 other real MAST SGAL WFPC2 fields.
5. Real-field ties block admission. All 20 controls are mandatory. Synthetic
   p-values remain diagnostics and cannot admit a candidate by themselves.
6. Downloads use a bounded thread pool; CPU-heavy independent field evaluations
   use a bounded ten-process pool. Internal BLAS pools are capped to one thread
   to avoid CPU oversubscription.
7. The console shows an animated progress bar, completed model count, active
   workers, elapsed time and ETA. A machine-readable progress snapshot is
   written to results_v2_1/progress_v2.1.json.
8. Model and whole-field checkpoints use atomic replacement. Partial reports
   are also written after every completed field, and prior reports/logs are
   rotated into results_v2_1/history instead of being overwritten.

RUN ON WINDOWS
--------------
Double-click:

    run_janus_cosmos_v2_1.bat

The BAT verifies dependencies, the v2.0.2 negative certificate, the frozen
v2.1 protocol and the unchanged parent detector. It then downloads 168 frozen
FITS products and runs the evaluation. Completed downloads and v2.1.0 model
checkpoints remain compatible, so the interrupted original v2.1 run can be
resumed by launching the same BAT again.

The expanded control cohort remains intentionally larger than v2.0.2. The
parallel scheduler reduces wall-clock time without reducing that workload. Do
not reduce the cohort or manually edit the protocol/expected files; either
change fails closed or disables admission.

OUTPUTS
-------
    results_v2_1/janus-cosmos-v2.1-report.json
    results_v2_1/janus-cosmos-v2.1-events.jsonl
    results_v2_1/SUMMARY_v2.1.txt
    results_v2_1/terminal_v2.1.log
    results_v2_1/progress_v2.1.json

INTERPRETATION
--------------
PASS means that the frozen protocol completed without an integrity/runtime
error. A target is admitted only when its own status is explicitly
SKY_FIXED_MORPHOLOGY_CANDIDATE or HST_CROSS_FILTER_MORPHOLOGY_CANDIDATE.

Claim ceiling:
Image-level sky-fixed morphology candidate only after all frozen real-field
specificity gates. No artificial-structure, hidden-message, censorship,
astronomical-discovery, or new-physics claim.
