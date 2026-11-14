# JANUS COSMOS — FERMI-V7 event evidence lab

FERMI-V7 is an interactive, reproducible, read-only examination of **existing archival Fermi-LAT Pass 8 event diagnostics**. It is not a new telescope observation, scientific breakthrough, raw FITS re-download, LIV model fitting result, or externally independent validation.

Website after successful deploy:
https://hawkar-usls.github.io/Janus-Cosmos/fermi/

## Sources, not identical populations

- FERMI-V3: 38 greater-than-10-GeV Table 7 published/rounded photons from 11 GRBs for overlay/methodology validation.
- FERMI-V4 manifest: 3 NASA FSSC Extended Pass8 queries; three EV00 plus three SC00 FITS files with byte SHA-256. Original binaries are not in Git.
- FERMI-V4 event diagnostics: 6,520 events in the three query ROIs; 175 diagnostic SOURCE-class candidates surviving preselected quality, spacecraft, GTI, energy and angular cuts.
- FERMI-V6: instrument and physics gate history. No complete physical multi-burst LAT+GBM LIV likelihood.
- FERMI-V7: static inspector of the archived V4 individual events with source-bound provenance and safeguards. V7 does not create independent astronomical observations.

The 175 selected candidates are **NOT** 175 confirmed photons from the GRBs; the source-class and position cuts do not provide per-event association probabilities. Unlike Table 7, these individual V4 event energies may be as low as 1 GeV and are not identical to published high-energy photon association rows.

## Build and independent implementation cross-check

    python tools/build_cosmos_fermi_v7.py --out docs/fermi/data/v7-events.json
    python tools/verify_fermi_v7_independent.py --view docs/fermi/data/v7-events.json
    python -m unittest discover -s tests -p 'test_cosmos_fermi_v7.py' -v

The primary generator checks the original V4 JSON-manifest digest, frozen cohort, all six FSSC provenance hashes and sizes, event ordering and unique identifiers, SOURCE class mask, ROI/zenith/incidence bounds, estimated GBM-trigger arithmetic, and event/rejection accounting. It emits SHA-256 of input JSON bytes and canonical event rows.

A **second independently implemented JSON comparator** audits the complete event multiset against V4, source file bytes, source counts and scientific no-claims; it does not import the builder's event validator. Deliberate corruption tests target events, time, class bits, provenance and physics-status promotion. These are internal software controls on the same archived source, not an outside astrophysical replication.

## Optional FITS-byte replay

Normal Pages builds do not re-download or re-verify original FITS binaries. The page consequently displays:

    RAW_FITS_NOT_RECHECKED_THIS_BUILD

For exact binary replay using original pinned FSSC query outputs:

    python scripts/fetch_fermi_lat_v4.py --raw-dir /path/to/fermi
    python tools/build_cosmos_fermi_v7.py --raw-dir /path/to/fermi --out /tmp/fermi-v7-replay.json
    python tools/verify_fermi_v7_independent.py --view /tmp/fermi-v7-replay.json

All six FITS files must match archived hashes and sizes, and the old V4 FITS ingestion must reproduce every V4 event row. That is **byte re-verification plus internal same-parser replay**, not independent parser validation. Missing or changed FSSC query outputs fail closed. Query URLs may expire; new files must be stored as new distinct evidence with provenance, never silently substituted.

## Scientific gates

Time after the estimated Fermi GBM trigger is **not the time the photon was emitted at the source**. A SOURCE-class Pass8 event and a source-centered ROI are not an event-by-event probability of origin. Without a calibrated source-emission-lag model, GBM response, LAT effective area/PSF/EDISP and diffuse backgrounds, we cannot infer physical propagation delay.

- No LIV p-value / astrophysical LIV bound.
- No inferred lattice spacing or faster-than-light/retrocausality claim.
- No complete joint LAT+GBM source likelihood.
- No external independent astrophysical replication.

NASA FSSC documentation:
- https://fermi.gsfc.nasa.gov/ssc/data/access/
- https://fermi.gsfc.nasa.gov/ssc/LATDataQuery_help.html
- https://fermi.gsfc.nasa.gov/ssc/data/access/lat/lat_data_columns.html
- https://fermi.gsfc.nasa.gov/ssc/data/analysis/likelihood/likelihood_overview.html

NEXT: exact-byte six-file replay (or source-controlled replacement query if ephemeral links expire), a separate implementation or external Fermitools cross-check, event-level source/background likelihood and physical source-lag constraints, followed by an independent pre-registered holdout.

STATUS: INTERNAL SOFTWARE AUDIT; NOT INDEPENDENTLY REPLICATED; NOT EXTERNALLY VALIDATED.
