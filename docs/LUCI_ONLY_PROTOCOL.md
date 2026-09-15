# JANUS COSMOS — LUCI/LUCIFER-only protocol v1

## Scope

This branch isolates observations from the Large Binocular Telescope near-infrared instruments LUCI1/LUCI2, historically named LUCIFER. It is intentionally separate from the HST/MAST canonical corpus.

LUCI is a ground-based LBT instrument, not a space telescope. The branch preserves the historical `LUCIFER` alias only for provenance matching.

## Instrument boundary

A frame is admitted only when the FITS header proves `INSTRUME` is LUCI/LUCIFER. If a telescope header is present, it must identify the Large Binocular Telescope. Non-LUCI files fail closed.

The geometry pipeline currently admits **imaging only**. Spectroscopic detector frames are rejected because their 2-D geometry is not the same statistical object as a sky image.

## Spectral scope

Near infrared, approximately 0.9–2.5 micrometres (exact limits depend on LUCI1/LUCI2 configuration).

## Data source

Primary source: official LBT Archive / IA2. The helper

```bash
python experiments/luci/build_luci_archive_manifest.py --metadata-only
```

checks the official TAP service and discovers its ObsCore schema. A live public two-band manifest can be built with:

```bash
python experiments/luci/build_luci_archive_manifest.py --limit 2000 --max-targets 5
```

## Run

```bash
python -m janus_cosmos.luci_pipeline \
  --manifest data/runtime/luci/luci_archive_manifest.json \
  --nulls 1024 \
  --seeds 20260815,20260816,20260817
```

For infrastructure smoke only:

```bash
python -m janus_cosmos.luci_pipeline \
  --manifest <manifest.json> \
  --nulls 32 \
  --seeds 1,2 \
  --allow-underpowered
```

## Required scientific gates

1. LUCI/LUCIFER FITS provenance passes.
2. Imaging mode only.
3. Source bytes are hashed and cached.
4. Existing morphology-preserving JANUS nulls are reused unchanged.
5. Bonferroni power guard remains active.
6. A target requires at least two independently tested LUCI bands to become a cross-filter image candidate.
7. No HST/JWST/WISE/Spitzer or other instrument may contribute to the LUCI-only result.

## Interpretation ceiling

A passing result is only a LUCI near-infrared image-level geometric candidate requiring localization, reduction-quality checks, injection-recovery sensitivity tests, matched controls, and independent replication. It is not evidence of unknown physics, artificial structures, hidden communication, or extraterrestrial intelligence.
