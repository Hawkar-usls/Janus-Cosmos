# LUCI–PALOMAR JPFM-2F-C protocol

`LUCI-PALOMAR-JPFM-2F-C` is the first JANUS COSMOS gate intended to permit a direct, independent, decades-later near-infrared counterpart test for frozen Palomar/POSS candidate coordinates.

It is intentionally stricter than the earlier LUCI transfer pilot. The earlier raw morphology PASS was voided after a corrective audit showed that unresolved/single-pixel structures could enter the source catalogue. 2F-C therefore requires a PSF-aware multi-pixel detector and injection–recovery on real LUCI backgrounds before any morphology inference is admitted.

## Frozen Palomar corpus

The experiment expands the 64-source pilot to 640 coordinates without consulting LUCI outcomes. For each of the 16 already-frozen structural clusters it selects 20 median-near `anomaly_score` rows as `typical`, then 20 remaining highest-`anomaly_score` rows as `unusual`, with `src_id` tie-breaking. The resulting CSV is written and SHA-256 hashed before the first LUCI archive-overlap query.

This is a selection/preregistration device, not a statement that the selected POSS-I detections are astronomical anomalies.

## PSF-aware source admission

A LUCI source candidate must be a resolved multi-pixel object. The detector smooths only for peak discovery; shape admission is measured on the unsmoothed background-subtracted image. It rejects components smaller than five pixels, minor-axis FWHM below 0.8 px, major-axis FWHM above 12 px, and elongation above 4.

The detector must then pass injection–recovery on real LUCI frames. Gaussian PSFs are injected across a preregistered FWHM/SNR grid and compared against single-pixel impulses of matched peak amplitude. The gate requires at least 80% recovery overall, 90% recovery for injections with peak SNR >= 8, and <=5% hot-pixel acceptance. All preregistered LUCI validation frames must pass before the scientific chain is admitted.

## Direct overlap chain

The execution order is fixed:

`Palomar coordinate -> coarse LUCI archive opportunity -> exact downloaded-FITS WCS containment -> overlap-frame injection/recovery -> IR source/no-source -> PSF morphology -> matched local controls`

A coarse archive-center match is never considered an overlap. The exact downloaded FITS file must contain the Palomar coordinate under its celestial WCS. If an exact overlap exists, that exact frame is separately subjected to injection–recovery before source presence is inspected. A detected counterpart is compared only against local PSF sources from the same LUCI frame; at least five local controls are required for morphology comparison.

## Firewalls and claim ceiling

No UAP, nuclear-window, witness, temporal-event, or later outcome labels enter the expanded Palomar selection. CI success only means the protocol executed; a scientific `BLOCKED` receipt remains blocked.

The maximum allowed interpretation is an independent near-IR counterpart test at frozen historical sky coordinates. It does not establish an anomaly, artificial origin, UAP origin, extraterrestrial intelligence, or causality.
