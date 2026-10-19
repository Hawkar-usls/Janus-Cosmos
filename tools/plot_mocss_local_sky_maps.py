#!/usr/bin/env python3
"""Plot exploratory per-orientation 16004/15850-A sky-projected RATE controls."""

import argparse
import json
import numpy as np
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', default='/tmp/mocss_16004_local_sky_maps.npz')
    p.add_argument('--receipt', default='janus_jwst/mocss_16004A_local_sky_projection_2026-10-16.json')
    p.add_argument('--output', default='janus_jwst/mocss_local_sky_quicklook_2026-10-16.png')
    args = p.parse_args()
    maps = np.load(args.source)
    receipt = json.load(open(args.receipt))
    fig, axs = plt.subplots(2, 4, figsize=(17, 8), constrained_layout=True)
    for j, group in enumerate('ABCD'):
        for i, label in enumerate(['16004A', '15850p8A']):
            ax = axs[i, j]
            data = maps[group + '_' + label]
            scale = np.nanpercentile(np.abs(data), 90)
            ax.imshow(data, origin='lower', extent=[-1.55, 1.55, -1.55, 1.55],
                      cmap='coolwarm', vmin=-scale, vmax=scale)
            ax.plot(0, 0, '+', color='lime', markersize=12, markeredgewidth=2)
            for name, (e,n) in ((name, (r['east_arcsec'], r['north_arcsec']))
                                for name,r in receipt['neighbors'].items()):
                ax.plot(e, n, marker='o', markerfacecolor='none',
                        markeredgecolor='yellow', markersize=9, linewidth=1)
                ax.annotate(name.split('_')[0], (e,n), color='yellow',
                            xytext=(4,5), textcoords='offset points', fontsize=8)
            ax.set_title(f'{group} / {label} / noncalibrated rate')
            ax.set_xlabel('East of source (arcsec)')
            ax.set_ylabel('North of source (arcsec)')
            ax.set_xlim(-1.5, 1.5)
            ax.set_ylim(-1.5, 1.5)
    fig.suptitle('MOCSS local sky projection ±1.5″ | WCS+raw-rate high-pass exploratory controls',
                 fontsize=13)
    fig.savefig(args.output, dpi=140)
    print('WROTE', args.output)


if __name__ == '__main__':
    main()
