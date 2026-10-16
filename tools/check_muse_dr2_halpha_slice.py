#!/usr/bin/env python3
"""Independent MUSE DR2 H-alpha quick-look from a byte-range slice.

The UDF-MOSAIC MUSE cube is ~26 GB.  A contiguous HTTP byte range for 28
spectral planes (indices 2740..2767 inclusive) from its first SCI extension
requires ~100 MB.  This script operates offline on that downloaded range.

This is a *diagnostic*, sky-subtracted, uncorrected-aperture test.  It cannot
exclude a source which was transient in 2023, eight years after MUSE exposures.
"""

import argparse
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


RA0 = 53.1631122222
DEC0 = -27.8244911111
CD = np.array([
    [-4.1285823637633e-05, 3.7173922575492e-05],
    [3.7173922575492e-05, 4.1285823637633e-05],
])
NX = 945
NY = 947
K0 = 2740
K1 = 2767
FIRST_SCI_BYTE = 28800
WAVE0 = 4750.0
DELTA_LAMBDA = 1.25
AIR_HALPHA_ANGSTROM = 6562.8
VACUUM_PABETA_ANGSTROM = 12821.6
ORPHAN_VACUUM_WAVELENGTH = 16004.0

TARGETS = {
    "orphan": (53.13127, -27.78908),
    "NW_candidate": (53.131165, -27.788930),
    "S_candidate": (53.131301, -27.789367),
    "E_candidate": (53.131584, -27.789105),
    "MUSE_7974_Lyalpha_low_confidence": (53.1318, -27.7893),
}


def world_to_pixel(ra_deg, dec_deg):
    ra, dec, ra0, dec0 = map(math.radians, (ra_deg, dec_deg, RA0, DEC0))
    dra = ra - ra0
    cosc = math.sin(dec0) * math.sin(dec) + math.cos(dec0) * math.cos(dec) * math.cos(dra)
    xi = math.cos(dec) * math.sin(dra) / cosc
    eta = (math.cos(dec0) * math.sin(dec) - math.sin(dec0) * math.cos(dec) * math.cos(dra)) / cosc
    # FITS CRPIX1=CRPIX2=1 => zero-index pixels are CD^{-1} * TAN offsets.
    x, y = np.linalg.solve(CD, np.array([xi, eta]) * 180.0 / math.pi)
    return float(x), float(y)


def circular_mask(x, y, r_arcsec, yy, xx):
    # These are 0.2 arcsec MUSE spaxels to within ~0.1%; the CD matrix includes rotation.
    pix_radius = r_arcsec / 0.2
    return (xx - x) ** 2 + (yy - y) ** 2 <= pix_radius ** 2


def robust_spectrum(cube, x, y, radius_arcsec=0.6, annulus=(0.95, 1.65)):
    xlo, xhi = max(0, int(x - 12)), min(NX, int(x + 13))
    ylo, yhi = max(0, int(y - 12)), min(NY, int(y + 13))
    yy, xx = np.mgrid[ylo:yhi, xlo:xhi]
    src = circular_mask(x, y, radius_arcsec, yy, xx)
    ann = circular_mask(x, y, annulus[1], yy, xx) & ~circular_mask(x, y, annulus[0], yy, xx)
    data = cube[:, ylo:yhi, xlo:xhi]
    spec = []
    frac = []
    for plane in data:
        good = np.isfinite(plane)
        src_good = src & good
        ann_good = ann & good
        frac.append(round(int(src_good.sum()) / int(src.sum()), 3))
        if src_good.sum() < src.sum() * 0.7 or ann_good.sum() < 20:
            spec.append(float("nan"))
            continue
        # The SCI cube contains flux per spaxel per Angstrom in units of 1e-20.
        sky_flux_per_spaxel = np.median(plane[ann_good])
        spec.append(float(np.sum(plane[src_good] - sky_flux_per_spaxel)))
    return np.array(spec), frac, int(src.sum()), int(ann.sum())


def continuum_subtract(spec):
    wave = WAVE0 + DELTA_LAMBDA * np.arange(K0, K1 + 1)
    expected = AIR_HALPHA_ANGSTROM * ORPHAN_VACUUM_WAVELENGTH / VACUUM_PABETA_ANGSTROM
    baseline = (np.abs(wave - expected) >= 7.0) & (np.abs(wave - expected) <= 17.5)
    cont = float(np.nanmedian(spec[baseline]))
    # Independent sideband scatter, not a fully modeled uncertainty.
    median_abs_dev = float(np.nanmedian(np.abs(spec[baseline] - cont)))
    sigma = 1.4826 * median_abs_dev
    window = np.abs(wave - expected) <= 3.125
    line_flux = float(np.nansum((spec[window] - cont) * DELTA_LAMBDA) * 1e-20)
    noise_approx = sigma * DELTA_LAMBDA * math.sqrt(int(window.sum())) * 1e-20
    return {
        "candidate_Halpha_flux_unaperturecorrected_cgs": line_flux,
        "continuum_median_in_cube_flux_units": cont,
        "sideband_MAD_sigma_in_cube_flux_units": sigma,
        "heuristic_MAD_noise_cgs": noise_approx,
        "diagnostic_SNR_not_formal": line_flux / noise_approx if noise_approx > 0 else None,
        "window_channel_indices": (np.where(window)[0] + K0).tolist(),
        "baseline_channel_indices": (np.where(baseline)[0] + K0).tolist(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("slice_file", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--multipart", action="store_true", help="Parse HTTP 206 multipart byte-range response")
    parser.add_argument("--row-start", type=int, default=108)
    parser.add_argument("--row-end", type=int, default=162)
    args = parser.parse_args()
    observed_bytes = args.slice_file.stat().st_size
    if args.multipart:
        raw = args.slice_file.read_bytes().lstrip(b"\r\n")
        first_line = raw.split(b"\r\n", 1)[0]
        if not first_line.startswith(b"--"):
            raise ValueError("Invalid multipart/byteranges MIME boundary")
        boundary = first_line
        cube = np.full((K1 - K0 + 1, NY, NX), np.nan, dtype=np.float32)
        seen = set()
        cursor = 0
        while True:
            begin = raw.find(boundary, cursor)
            if begin < 0 or raw[begin + len(boundary):begin + len(boundary) + 2] == b"--":
                break
            meta_start = begin + len(boundary) + 2
            header_end = raw.find(b"\r\n\r\n", meta_start)
            if header_end < 0:
                raise ValueError("Missing multipart header delimiter")
            header = raw[meta_start:header_end].decode("ascii")
            match = re.search(r"Content-Range:\s*bytes\s+(\d+)-(\d+)/(\d+)", header, re.I)
            if not match:
                raise ValueError(f"Missing Content-Range: {header!r}")
            start, end, total = map(int, match.groups())
            if total != 26353494720:
                raise ValueError("Unexpected upstream MUSE FITS size")
            expected_width = (args.row_end - args.row_start + 1) * NX * 4
            if end - start + 1 != expected_width:
                raise ValueError("Wrong HTTP byte-range payload width")
            offset = start - FIRST_SCI_BYTE - args.row_start * NX * 4
            if offset % (NX * NY * 4):
                raise ValueError("Byte-range does not align to expected spectral plane")
            plane_id = offset // (NX * NY * 4)
            if not K0 <= plane_id <= K1 or plane_id in seen:
                raise ValueError(f"Unexpected or duplicate spectral plane: {plane_id}")
            payload_start = header_end + 4
            payload_end = payload_start + expected_width
            data = np.frombuffer(raw[payload_start:payload_end], dtype=">f4").reshape(args.row_end - args.row_start + 1, NX)
            cube[plane_id - K0, args.row_start:args.row_end + 1, :] = data
            seen.add(plane_id)
            cursor = payload_end
        if seen != set(range(K0, K1 + 1)):
            raise ValueError(f"Missing spectral planes: {set(range(K0, K1+1)) - seen}")
    else:
        expected_bytes = (K1 - K0 + 1) * NY * NX * 4
        if observed_bytes != expected_bytes:
            raise ValueError(f"Expected {expected_bytes} bytes; got {observed_bytes}")
        cube = np.memmap(args.slice_file, dtype=">f4", mode="r", shape=(K1 - K0 + 1, NY, NX))
    results = {}
    all_specs = {}
    for name, coordinates in TARGETS.items():
        x, y = world_to_pixel(*coordinates)
        spec, fractions, pixels_src, pixels_ann = robust_spectrum(cube, x, y)
        all_specs[name] = spec
        results[name] = {
            "ra_deg": coordinates[0], "dec_deg": coordinates[1],
            "pixel_x0": x, "pixel_y0": y,
            "aperture_pixels": pixels_src, "annulus_pixels": pixels_ann,
            "valid_source_fraction_min": min(fractions),
            **continuum_subtract(spec),
            "sky_subtracted_spectrum_1e-20_per_A": [float(v) if np.isfinite(v) else None for v in spec],
        }

    center = TARGETS["orphan"]
    # Reproducible blank apertures near the target: offsets > 4 arcsec.
    controls = []
    for east in (-8, -6, -4, 4, 6, 8):
        for north in (-8, -6, -4, 4, 6, 8):
            ra = center[0] + east / (3600 * math.cos(math.radians(center[1])))
            dec = center[1] + north / 3600
            x, y = world_to_pixel(ra, dec)
            if x < 15 or y < 15 or x > NX - 15 or y > NY - 15:
                continue
            spec, frac, _, _ = robust_spectrum(cube, x, y)
            if min(frac) < 0.7 or not np.all(np.isfinite(spec)):
                continue
            cont = continuum_subtract(spec)
            if min(frac) >= 0.7 and math.isfinite(cont["candidate_Halpha_flux_unaperturecorrected_cgs"]):
                controls.append(cont["candidate_Halpha_flux_unaperturecorrected_cgs"])
    ctrl = np.array(controls)
    ctrl_median = float(np.median(ctrl)) if ctrl.size else None
    ctrl_sigma = float(1.4826 * np.median(np.abs(ctrl - np.median(ctrl)))) if ctrl.size else None
    if ctrl_sigma and ctrl_sigma > 0:
        results["orphan"]["blank_aperture_SNR_diagnostic"] = (
            results["orphan"]["candidate_Halpha_flux_unaperturecorrected_cgs"] - ctrl_median
        ) / ctrl_sigma
    receipt = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_url": "https://amused.univ-lyon1.fr/data/UDF/HUDF/download/DATACUBE_UDF-MOSAIC.fits",
        "muse_observation_epochs": "2014-2015 (non-contemporaneous with JWST Aug/Nov 2023)",
        "dataset": "MUSE UDF-MOSAIC DR2 10h", "science_hdu": 1,
        "source_cube_byte_start": FIRST_SCI_BYTE,
        "subset_k_start": K0, "subset_k_end_inclusive": K1,
        "byte_start_in_original": FIRST_SCI_BYTE + K0 * NX * NY * 4,
        "byte_end_inclusive_in_original": FIRST_SCI_BYTE + (K1 + 1) * NX * NY * 4 - 1,
        "extracted_rows_in_2d_plane": [args.row_start, args.row_end] if args.multipart else [0, NY - 1],
        "actual_http_payload_file_size_bytes": observed_bytes,
        "wavelength_type": "air wavelengths (CTYPE3=AWAV)",
        "wavelength_start_A": WAVE0 + K0 * DELTA_LAMBDA,
        "wavelength_step_A": DELTA_LAMBDA,
        "expected_Halpha_air_wavelength_A_if_Pabeta": AIR_HALPHA_ANGSTROM * ORPHAN_VACUUM_WAVELENGTH / VACUUM_PABETA_ANGSTROM,
        "line_model": "Pa-beta vacuum 12821.6 A observed vacuum 16004 A; rest Halpha air 6562.8 A",
        "flux_units_per_spaxel_per_A": "1e-20 erg s-1 cm-2 A-1",
        "aperture_arcsec_radius": 0.6,
        "background_annulus_arcsec_radius": [0.95, 1.65],
        "analyses": results,
        "blank_control_count": len(controls),
        "blank_control_fluxes_cgs": controls,
        "blank_control_median_cgs": ctrl_median,
        "blank_control_sigma_MAD_cgs": ctrl_sigma,
        "limitations": [
            "Only MUSE SCI extension, not the propagated variance or exposure cube",
            "Aperture flux is uncorrected for PSF encircled energy",
            "Simple local sky annulus may oversubtract extended emission",
            "Only 35 A of local continuum available in this thin spectral slice",
            "MUSE epochs 2014-2015 cannot constrain a transient specific to 2023",
            "OH sky residuals and correlated noise can invalidate simplistic SNR",
        ],
    }
    if args.output:
        args.output.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    for name, d in results.items():
        print(f"{name:40} Halpha={d['candidate_Halpha_flux_unaperturecorrected_cgs']:.3e} "
              f"sidebandMADsnr={d['diagnostic_SNR_not_formal']:.2f} ",
              f"pixel=({d['pixel_x0']:.2f},{d['pixel_y0']:.2f})")
    print(f"blank_count={len(ctrl)} blank_median={ctrl_median} blank_sigma_MAD={ctrl_sigma}")


if __name__ == "__main__":
    main()
