# Replay instructions: GO-3290 MOCSS sky coordinates and raw detector controls

**Purpose:** independent verification of 16 exposure sky-to-detector wavelength geometry and exploratory source/background controls. This is **not** a flux-calibrated spectroscopy pipeline or evidence of a physical beam.

Requirements: Python 3.11–3.12; `jwst==1.19.1` (with AstroPy/NumPy/SciPy transitively installed), `matplotlib` only for PNG quick-looks; network to public MAST and CRDS for the initial downloads. The CI validator needs only Python's standard library and the committed receipt.

From repository root, create a temporary analysis folder and set the **official CRDS context**:

```bash
mkdir -p janus_jwst
export CRDS_PATH=/tmp/janus_jwst_crds
export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
export CRDS_CONTEXT=jwst_1413.pmap
python -m pip install 'jwst==1.19.1' matplotlib
```

Download the 32 RATE FITS and two MSA FITS (~2.85 GB total); the downloader resumes files already present in `/tmp`. No user credentials are required for public MAST GO-3290 products.

```bash
python tools/fetch_mocss_all_16_rate.py --detector nrs1
python tools/fetch_mocss_all_16_rate.py --detector nrs2
```

Calibrate a 3×3 coarse shutter grid from two anchor NRS2 RATE frames; these virtual shutter objects have **NO detection/line-flux status**. Fit two local affine shutter proposals from those WCS-recovered centers:

```bash
python tools/prepare_mocss_wcs_probe.py --detector nrs2
python tools/calibrate_mocss_wcs_grid.py
python tools/locate_mocss_target_shutters.py
```

Generate and process **all 16 per-detector WCS-probe copies**, using exposure-specific RA_REF/DEC_REF and MSA IDs. The script explicitly handles unsigned-16-bit overflow in TEST copies only.

```bash
python tools/prepare_all_mocss_target_wcs.py --detector nrs1 --half-width 4
python tools/prepare_all_mocss_target_wcs.py --detector nrs2 --half-width 4 --output-dir /tmp/janus_mocss_wcs_all16_nrs2
python tools/trace_mocss_all_16.py
python tools/trace_mocss_all_16.py --inputs /tmp/janus_mocss_wcs_all16_nrs2/all16_input_receipt.json --output janus_jwst/mocss_16_exposure_nrs2_wcs_trace_2026-10-16.json
```

Run the uncalibrated RATE-pixel checks, four spectral coadds with background-method sensitivity, and near-source sky projection and descriptive controls. **Do not cite the output z-like numbers as formal significance.**

```bash
python tools/probe_mocss_rate_line_contrasts.py
python tools/stack_mocss_wcs_registered_rate.py
python tools/stack_mocss_wcs_registered_rate.py --background-method local_median_41 --output janus_jwst/mocss_wcs_registered_rate_highpass_2026-10-16.json --image-output /tmp/mocss_rate_registered_highpass_four_groups.npz
python tools/stack_mocss_wcs_registered_rate.py --background-method local_median_41 --median-width 31 --output janus_jwst/mocss_wcs_registered_rate_highpass_31px_2026-10-16.json --image-output /tmp/mocss_rate_highpass_31px.npz
python tools/stack_mocss_wcs_registered_rate.py --background-method local_median_41 --median-width 61 --output janus_jwst/mocss_wcs_registered_rate_highpass_61px_2026-10-16.json --image-output /tmp/mocss_rate_highpass_61px.npz
python tools/map_mocss_local_sky.py
python tools/assess_mocss_local_sky_coherence.py
python tools/compile_mocss_sky_registered_receipt.py
```

Use `python tools/verify_mocss_sky_registered_rate_receipt.py` to verify the committed evidence receipt offline. The machine receipt contains SHA-256 of every full original detector RATE/representative MSA file, exposing reprocessing changes. Never replace original archival FITS with edited TEST products. To test flux significance/extension or Paβ–Paγ, complete the independent detector-background + flat-field + STScI Stage 2 spectrophotometric pipeline and injected-line recovery documented in `domains/space/MOCSS_16_SKY_REGISTERED_RATE_RECHECK_2026-10-16.md`.
