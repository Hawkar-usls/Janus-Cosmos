# LUCI-PALOMAR-JPFM-2F-E — Staged Exact-WCS Protocol

## Why 2F-E exists

2F-D completed the exhaustive archive-first crossmatch of all 122,820 frozen POSS-I S0 coordinates, but its preregistered 250-FITS safety cap blocked all pixel inspection after the frozen metadata-WCS set expanded to 918 pairs across 798 LUCI files.

2F-E does **not** weaken that cap or rewrite 2F-D. It uses the already frozen parent artifact as an immutable input and separates two operations that were previously coupled:

1. deterministic staged FITS **header/WCS** replay with no image-array access;
2. a much smaller, pre-pixel-frozen representative set for R1 injection/recovery and counterpart morphology.

## Parent binding

The workflow downloads artifact `9247843315` from workflow run `31888129868` and requires:

- artifact ZIP SHA-256 `1240a600dfb0189243dfb3188ab53dcc8ad6f7b270236c56756b49c2e4fc6184`;
- frozen metadata-WCS CSV SHA-256 `aa57e64deca11bb2afd09364fa8b72837e92fc3c2187357c5041b273023754c6`;
- 918 pairs, 64 Palomar sources, 798 LUCI files.

If any binding changes, 2F-E fails before replay.

## Exposure-admissibility firewall

Before any new FITS access, rows whose literal archive `filters` string contains `blind` or `PV lens` (case-insensitive) are excluded from **counterpart inference only**. The parent exhaustive coverage record remains untouched.

The bound parent set deterministically yields 443 exposure-admissible pairs across 42 sources and 403 files.

## Staged exact WCS

The 403 admissible filenames are sorted and partitioned into deterministic contiguous header-only shards of at most 225 files. Each shard stays below the inherited 250-file safety ceiling.

For every file:

- the complete FITS file is downloaded and SHA-256 recorded;
- only FITS headers/WCS are inspected;
- image arrays are not read;
- exact celestial-WCS containment is evaluated with the existing `exact_wcs_pixel` implementation.

Only after both shards complete is the exact FITS-WCS pair set written and hashed.

## Representative selection before pixels

For each Palomar source with one or more exact FITS-WCS overlaps:

1. prefer exact rows whose archive filter configuration matches `clear z/J/H/K/Ks`;
2. if none exist, use all exposure-admissible exact rows;
3. sort unique files by `date_obs`, then filename;
4. freeze the earliest, median, and latest files, deduplicated.

This yields no more than three representatives per exact source. The representative CSV and SHA are materialized before the first image array is read.

## Pixel gate

Only frozen representatives may enter image analysis.

Every selected LUCI file must independently pass the R1 PSF-relative gate:

- all-star recovery >= 0.80;
- recovery at SNR >= 8 >= 0.90;
- hot-pixel acceptance <= 0.05.

At the exact Palomar coordinate:

- a detected R1 PSF source may enter same-frame matched morphology;
- a non-detection is admitted only if local coordinate injections at SNR 8 and 12 are both recovered using the frame-derived R1 PSF width;
- otherwise the no-source result is blocked by local sensitivity.

Matched controls remain same-frame, within 300 px, peak-SNR ratio 0.5–2.0, minimum 8 and maximum 20.

## Claim ceiling

2F-E is a representative exact-overlap near-IR counterpart test. It is **not** a complete photometric reprocessing or stacked-depth survey of every LUCI exposure. No temporal/UAP labels are used in admission or scoring, and no anomaly, artificial-origin, UAP-origin, ETI, or causal claim is permitted.
