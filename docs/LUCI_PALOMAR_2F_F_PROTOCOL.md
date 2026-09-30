# LUCI-PALOMAR-JPFM-2F-F

Targeted unresolved-sensitivity recovery only. The 92 qualified 2F-E negative pair results are immutable inputs and are not re-read at pixel level.

The gate freezes the eight 2F-E problem sources, then tests the 100 exact-WCS exposures that were never selected for 2F-E pixel replay. The inherited detector and R1 frame gate remain unchanged. A narrowly scoped edge repair is allowed only when the old local sensitivity block was caused by the fixed 12-pixel margin: the modeled Gaussian PSF must retain at least 99.5% of its flux on detector, and both SNR 8 and 12 injections must still be recovered by the unchanged `measure_psf_at` detector.

A source-level PASS means every one of the 42 frozen Palomar sources has at least one sensitivity-qualified LUCI no-counterpart epoch and no LUCI counterpart candidate was produced. This is not photometric completeness and carries no origin/causality claim.
