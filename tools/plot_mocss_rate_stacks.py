#!/usr/bin/env python3
"""Quick-look raw DN/s coadd images; NOT flux-calibrated detection maps."""

import argparse
import numpy as np
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='/tmp/mocss_rate_registered_four_groups.npz')
    ap.add_argument('--output', default='janus_jwst/mocss_raw_rate_4stacks_quicklook_2026-10-16.png')
    args = ap.parse_args()
    d = np.load(args.input)
    fig, axs = plt.subplots(2, 4, figsize=(17, 8), constrained_layout=True)
    for j, group in enumerate('ABCD'):
        for i, label in enumerate(['16004A', '15850p8A']):
            full = d[group + '_' + label]
            data = full[26 - 14:26 + 15, 76 - 25:76 + 26]
            ax = axs[i, j]
            vmin, vmax = np.nanpercentile(data, [8, 92])
            ax.imshow(data, aspect='auto', origin='lower', cmap='coolwarm',
                      extent=[-25.5, 25.5, -14.5, 14.5],
                      vmin=vmin, vmax=vmax)
            ax.plot(0, 0, color='lime', marker='+', markersize=12, markeredgewidth=1.5)
            ax.set_title(f'Group {group}: {label} / raw DN/s')
            ax.set_xlabel('Δ dispersion pixel')
            ax.set_ylabel('Δ cross-dispersion pixel')
    fig.suptitle('JWST/MAST GO-3290: WCS-aligned raw-rate cutouts — exploratory, NOT science-calibrated',
                 fontsize=13)
    fig.savefig(args.output, dpi=140)
    print('WROTE', args.output)


if __name__ == '__main__':
    main()
