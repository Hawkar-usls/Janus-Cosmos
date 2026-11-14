# JANUS COSMOS Volume Observatory v1

This standalone read-only WebGL2 (Canvas fallback) research companion to the Earth/Space Observer is inspired by the interaction concept at https://now.intellectronica.net/hydrogen-milky-way/ . That reference host was unavailable in the implementation environment. No original code or datasets were copied, and no exact pixel or scientific reproduction is claimed.

Deployment URL after merge and passing GitHub Pages checks:
https://hawkar-usls.github.io/Janus-Cosmos/galaxy/

## Three lenses

1. MODEL (default): seeded procedural spiral arms and luminous gas — *illustration*, not a measured 3D map of H I in the Milky Way. Every particle has a model provenance flag. No original dataset or source distances are inferred.
2. FERMI: data/cosmos/COSMOS-ALI-V3-FERMI-PUBLISHED-HE-PHOTONS-2026-10-18.json is copied unchanged into the site at data/fermi.json by GitHub Actions. The viewer rejects wrong schemas and maps 38 published high-energy photons from 11 GRBs into feature space: X = redshift z, Y = log10(seconds after GBM trigger), Z = energy in GeV. This is NOT RA/Dec, physical distance, source-emission time, true inter-photon time of flight, LIV constraint, or proof of retrocausality. Clicking a point shows source values and a commit-pinned link.
3. HI ANGULAR: manual import of a decimated, source-provenance-bearing HI table, displayed on a fixed-radius *angular sky shell*. Import must explicitly identify units, provenance and valid Galactic longitude, latitude and positive column densities; it never constructs unknown distance.

## HI angular import schema

The importer requires an object such as:

    {
      "schema": "janus.cosmos.hi.angular_samples.v1",
      "source_url": "https://example.org/original-observation",
      "nhi_unit": "cm-2",
      "samples": [
        {"l_deg": 17, "b_deg": -4, "nhi_cm2": 1.2e21},
        {"l_deg": 31, "b_deg": 6, "nhi_cm2": 7.1e20}
      ]
    }

These numerical rows are **synthetic format examples only**, NOT HI4PI measurements. Genuine rows must be transcribed or computed from a verified source under reproducible processing, with a pinned original file, checksums, pixelization, calibration and license.

Useful published HI data products:
- https://doi.org/10.1051/0004-6361/201629178
- https://lambda.gsfc.nasa.gov/product/foreground/fg_hi4pi_info.html
- https://cade.irap.omp.eu/doku.php?id=hi4pi

The H I integrated column-density field does NOT uniquely determine radial distances of the gas. An angular sky projection or arbitrary-radius shell is NOT a reconstructed 3D gas map.

## Boundaries and operation

- Original source JSON is never modified.
- The model is deterministic and never published as an observational claim.
- No hidden telemetry, credentials, external JavaScript libraries, trained model changes or auto-adjudication.
- The existing Earth/Space Observer is retained; a navigation link is the only change to its UI.
- Scientific ceilings: NO_LIV_BOUND, NO_RETROCAUSALITY_DETECTION, NO_SOURCE_LAG_CALIBRATION.
- Site deployment only after merge and successful Pages workflow, not merely after creation of a branch.

Validation:
    node --check docs/galaxy/model.mjs
    node --check docs/galaxy/view.mjs
    node --test tests/test_cosmos_volume.mjs

The PR checks these contracts. A live browser GPU performance benchmark and cross-browser visual fidelity test are separate, not yet established.

Next scientific engineering step: a checksum-verified and downsampled HI4PI angular sample package generated deterministically from FITS, with preserved provenance. A genuinely 3D Galactic H I reconstruction would need independent distance measurements/inference, physically justified uncertainty and external validation.
