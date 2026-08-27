#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
DATA = OUT / "miri_final"
PROV = OUT / "vega_miri_final_provenance.json"
REPORT = OUT / "vega_miri_resolvability_gate.json"
PROFILES = OUT / "vega_miri_radial_profiles.json"

DISTANCE_PC = 7.68
JWST_DIAMETER_M = 6.5
TARGET_EDGE_AU = (3.0, 5.0)


def choose_image(path: Path) -> tuple[np.ndarray, fits.Header, int]:
    with fits.open(path, memmap=False) as hdul:
        candidates = []
        for idx, hdu in enumerate(hdul):
            if hdu.data is None:
                continue
            a = np.asarray(hdu.data)
            b = np.squeeze(a)
            if b.ndim == 2 and min(b.shape) >= 16:
                candidates.append((b.size, idx, b.astype(float, copy=True), hdu.header.copy()))
        if not candidates:
            raise RuntimeError(f"no usable 2D FITS image found in {path.name}")
        _, idx, image, header = max(candidates, key=lambda q: q[0])
        return image, header, idx


def pixel_scale_arcsec(header: fits.Header) -> tuple[float | None, str]:
    try:
        w = WCS(header)
        if w.has_celestial:
            scales_deg = proj_plane_pixel_scales(w.celestial)
            scale = float(np.sqrt(abs(scales_deg[0] * scales_deg[1])) * 3600.0)
            if math.isfinite(scale) and 0.001 < scale < 10.0:
                return scale, "WCS_PROJ_PLANE_GEOMETRIC_MEAN"
    except Exception:
        pass
    for key in ("PIXELSCL", "PIXSCALE", "PIXSCAL1"):
        try:
            value = float(header[key])
            if math.isfinite(value) and 0.001 < value < 10.0:
                return value, f"HEADER_{key}"
        except Exception:
            pass
    return None, "UNAVAILABLE"


def center_pixel(header: fits.Header, shape: tuple[int, int]) -> tuple[float, float, str]:
    # FITS CRPIX is 1-indexed. The final reductions are expected to be registered
    # on Vega, but we record the convention explicitly rather than silently
    # claiming an independent astrometric recentering.
    try:
        x = float(header["CRPIX1"]) - 1.0
        y = float(header["CRPIX2"]) - 1.0
        if -0.5 <= x < shape[1] - 0.5 and -0.5 <= y < shape[0] - 0.5:
            return x, y, "FITS_CRPIX"
    except Exception:
        pass
    return (shape[1] - 1.0) / 2.0, (shape[0] - 1.0) / 2.0, "ARRAY_GEOMETRIC_CENTER_FALLBACK"


def diffraction_scales_arcsec(wavelength_micron: float) -> dict:
    lam = wavelength_micron * 1e-6
    rad_to_arcsec = 206264.806247
    return {
        "fwhm_proxy_1p03_lambda_over_D_arcsec": 1.03 * lam / JWST_DIAMETER_M * rad_to_arcsec,
        "rayleigh_1p22_lambda_over_D_arcsec": 1.22 * lam / JWST_DIAMETER_M * rad_to_arcsec,
    }


def robust_profile(
    image: np.ndarray,
    cx: float,
    cy: float,
    pix_arcsec: float,
) -> list[dict]:
    yy, xx = np.indices(image.shape, dtype=float)
    r_pix = np.hypot(xx - cx, yy - cy)
    r_au = r_pix * pix_arcsec * DISTANCE_PC
    finite = np.isfinite(image)
    if not np.any(finite):
        return []
    rmax = float(np.nanpercentile(r_au[finite], 98.0))
    # 1 AU bins are deliberately descriptive, not independent resolution elements.
    edges = np.arange(0.0, max(2.0, math.floor(rmax) + 1.0), 1.0)
    rows: list[dict] = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = finite & (r_au >= lo) & (r_au < hi)
        n = int(np.count_nonzero(m))
        if n < 8:
            continue
        vals = image[m]
        med = float(np.nanmedian(vals))
        mad = float(1.4826 * np.nanmedian(np.abs(vals - med)))
        rows.append(
            {
                "r_lo_au": float(lo),
                "r_hi_au": float(hi),
                "r_mid_au": float(0.5 * (lo + hi)),
                "n_finite_pixels": n,
                "median_flux": med,
                "robust_sigma_flux": mad,
            }
        )
    return rows


def finite_support_fraction(
    image: np.ndarray,
    cx: float,
    cy: float,
    pix_arcsec: float,
    rlo_au: float,
    rhi_au: float,
) -> float | None:
    yy, xx = np.indices(image.shape, dtype=float)
    r_au = np.hypot(xx - cx, yy - cy) * pix_arcsec * DISTANCE_PC
    ann = (r_au >= rlo_au) & (r_au < rhi_au)
    n = int(np.count_nonzero(ann))
    if n == 0:
        return None
    return float(np.count_nonzero(np.isfinite(image) & ann) / n)


def main() -> int:
    if not PROV.exists():
        raise SystemExit("missing vega_miri_final_provenance.json; run fetch_vega_miri_final.py first")
    provenance = json.loads(PROV.read_text(encoding="utf-8"))
    filters = []
    profiles = {}

    edge_arcsec = [TARGET_EDGE_AU[0] / DISTANCE_PC, TARGET_EDGE_AU[1] / DISTANCE_PC]
    for product in provenance["products"]:
        path = DATA / product["filename"]
        image, header, hdu_index = choose_image(path)
        scale, scale_source = pixel_scale_arcsec(header)
        cx, cy, center_source = center_pixel(header, image.shape)
        diff = diffraction_scales_arcsec(float(product["wavelength_micron"]))
        fwhm = diff["fwhm_proxy_1p03_lambda_over_D_arcsec"]
        rayleigh = diff["rayleigh_1p22_lambda_over_D_arcsec"]

        if edge_arcsec[1] <= fwhm:
            resolvability = "ENTIRE_3_TO_5_AU_EDGE_BELOW_FWHM_PROXY"
        elif edge_arcsec[0] < fwhm < edge_arcsec[1]:
            resolvability = "3_TO_5_AU_EDGE_STRADDLES_FWHM_PROXY"
        else:
            resolvability = "NOMINALLY_ABOVE_FWHM_PROXY_DIFFRACTION_ONLY"

        info = {
            "filter": product["filter"],
            "filename": product["filename"],
            "sha256": product["sha256"],
            "wavelength_micron": product["wavelength_micron"],
            "hdu_index": hdu_index,
            "shape_yx": list(image.shape),
            "bunit": str(header.get("BUNIT", "")),
            "pixel_scale_arcsec": scale,
            "pixel_scale_source": scale_source,
            "center_xy_zero_indexed": [cx, cy],
            "center_source": center_source,
            "finite_fraction_full_image": float(np.mean(np.isfinite(image))),
            "target_edge_au": list(TARGET_EDGE_AU),
            "target_edge_arcsec": edge_arcsec,
            "diffraction": diff,
            "resolvability": resolvability,
            "coronagraph_iwa_included": False,
        }
        if scale is not None:
            info["target_edge_pixels"] = [edge_arcsec[0] / scale, edge_arcsec[1] / scale]
            info["finite_support_fraction_3_to_5_au"] = finite_support_fraction(
                image, cx, cy, scale, *TARGET_EDGE_AU
            )
            info["finite_support_fraction_5_to_10_au"] = finite_support_fraction(
                image, cx, cy, scale, 5.0, 10.0
            )
            profiles[product["filter"]] = {
                "units": str(header.get("BUNIT", "")),
                "bin_width_au": 1.0,
                "warning": "1 AU radial bins are descriptive oversampling; adjacent bins are not independent at MIRI diffraction resolution.",
                "rows": robust_profile(image, cx, cy, scale),
            }
        else:
            info["target_edge_pixels"] = None
            info["finite_support_fraction_3_to_5_au"] = None
            profiles[product["filter"]] = {
                "status": "NO_PIXEL_SCALE_NO_PHYSICAL_RADIAL_PROFILE"
            }
        filters.append(info)

    all_fully = all(f["resolvability"] == "ENTIRE_3_TO_5_AU_EDGE_BELOW_FWHM_PROXY" for f in filters)
    any_straddle = any(f["resolvability"] == "3_TO_5_AU_EDGE_STRADDLES_FWHM_PROXY" for f in filters)
    if all_fully:
        gate = "DIRECT_3_TO_5_AU_EDGE_NOT_DIFFRACTION_RESOLVED"
    elif any_straddle:
        gate = "DIRECT_3_TO_5_AU_EDGE_ONLY_PARTIALLY_RESOLVED_AT_SHORTEST_MIRI_BAND"
    else:
        gate = "DIFFRACTION_GATE_ALONE_DOES_NOT_BLOCK_EDGE_RESOLUTION"

    report = {
        "schema": "janus.cosmos.vega.miri_resolvability_gate.v1.4",
        "distance_pc": DISTANCE_PC,
        "jwst_primary_diameter_m": JWST_DIAMETER_M,
        "target_warm_disk_inner_edge_au": list(TARGET_EDGE_AU),
        "target_warm_disk_inner_edge_arcsec": edge_arcsec,
        "filters": filters,
        "gate": gate,
        "interpretation": {
            "status": "REAL_MIRI_PRODUCTS_INGESTED_DIRECT_EDGE_RECOVERY_LIMITED_BY_RESOLUTION",
            "next_valid_test": "Use the final MIRI products for radial-profile/SED forward modeling and broader disk morphology; do not treat a 3-5 AU localized image residual as a resolved planet/shepherd detection.",
            "why": "At Vega's distance the 3-5 AU edge subtends only about 0.39-0.65 arcsec, comparable to or below JWST/MIRI diffraction scales. Coronagraphic inner-working-angle effects are not even included in this optimistic gate.",
        },
        "claim_ceiling": "RESOLVABILITY_AND_PROFILE_DIAGNOSTIC_ONLY_NO_PLANET_DETECTION",
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    PROFILES.write_text(
        json.dumps(
            {
                "schema": "janus.cosmos.vega.miri_radial_profiles.v1.4",
                "distance_pc": DISTANCE_PC,
                "profiles": profiles,
                "claim_firewall": "Profiles are descriptive reductions of the authors' final images; they are not an inverse solution for a planet mass or orbit.",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("VEGA MIRI RESOLVABILITY GATE PASS")
    print("gate =", gate)
    print("target_edge_arcsec =", edge_arcsec)
    for f in filters:
        print(
            f["filter"],
            "shape=", f["shape_yx"],
            "pix_arcsec=", f["pixel_scale_arcsec"],
            "fwhm_arcsec=", f["diffraction"]["fwhm_proxy_1p03_lambda_over_D_arcsec"],
            "rayleigh_arcsec=", f["diffraction"]["rayleigh_1p22_lambda_over_D_arcsec"],
            "edge_pixels=", f["target_edge_pixels"],
            "support_3_5=", f["finite_support_fraction_3_to_5_au"],
            "resolvability=", f["resolvability"],
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
