from __future__ import annotations

import hashlib
import io
import json
import math
import tarfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
from astropy.io import fits
from astropy.wcs import WCS

RA_DEG = 223.415064157
DEC_DEG = 33.979315670
CUTOUT_SIZE = 64
VERSION = "neo2"
BANDS = (1, 2)
ANNULUS_RMIN = 6.0
ANNULUS_RMAX = 14.0
STRONG_Z = 5.0
WEAK_Z = 3.0
URL = (
    "https://unwise.me/cutout_fits"
    f"?version={VERSION}&ra={RA_DEG:.12f}&dec={DEC_DEG:.12f}"
    f"&size={CUTOUT_SIZE}&bands=12&file_img_m=on&file_invvar_m=on"
)
OUT = Path("data/stargate/STARGATE-ABYDOS-RAW-UNWISE-FROZEN-CENTER-v1-LATEST-RECEIPT.json")


def finite(v):
    if v is None:
        return None
    try:
        x = float(v)
    except Exception:
        return None
    return x if math.isfinite(x) else None


def robust_center_sigma(values: np.ndarray) -> tuple[float | None, float | None, int]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 8:
        return None, None, int(x.size)
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    sigma = 1.4826 * mad
    if not np.isfinite(sigma) or sigma <= 0:
        q25, q75 = np.percentile(x, [25, 75])
        sigma = float((q75 - q25) / 1.349) if q75 > q25 else None
    return med, sigma, int(x.size)


def aperture3(data: np.ndarray, ix: int, iy: int) -> float | None:
    if iy - 1 < 0 or iy + 1 >= data.shape[0] or ix - 1 < 0 or ix + 1 >= data.shape[1]:
        return None
    a = data[iy - 1:iy + 2, ix - 1:ix + 2]
    if not np.all(np.isfinite(a)):
        return None
    return float(np.sum(a))


def aperture3_variance(invvar: np.ndarray | None, ix: int, iy: int) -> float | None:
    if invvar is None:
        return None
    if iy - 1 < 0 or iy + 1 >= invvar.shape[0] or ix - 1 < 0 or ix + 1 >= invvar.shape[1]:
        return None
    a = np.asarray(invvar[iy - 1:iy + 2, ix - 1:ix + 2], dtype=float)
    good = np.isfinite(a) & (a > 0)
    if not np.all(good):
        return None
    return float(np.sum(1.0 / a))


def measure_pair(img_path: Path, ivar_path: Path | None, band: int) -> dict:
    with fits.open(img_path, memmap=False) as hdul:
        data = np.asarray(hdul[0].data, dtype=float)
        header = hdul[0].header.copy()
    invvar = None
    if ivar_path is not None and ivar_path.exists():
        with fits.open(ivar_path, memmap=False) as hdul:
            invvar = np.asarray(hdul[0].data, dtype=float)

    wcs = WCS(header)
    x, y = wcs.world_to_pixel_values(RA_DEG, DEC_DEG)
    ix, iy = int(round(float(x))), int(round(float(y)))
    yy, xx = np.indices(data.shape)
    rr = np.hypot(xx - float(x), yy - float(y))
    ann = data[(rr >= ANNULUS_RMIN) & (rr <= ANNULUS_RMAX) & np.isfinite(data)]
    bg, sigma_pix, n_bg = robust_center_sigma(ann)

    center_pixel = None
    if 0 <= iy < data.shape[0] and 0 <= ix < data.shape[1] and np.isfinite(data[iy, ix]):
        center_pixel = float(data[iy, ix])

    pixel_excess = None if center_pixel is None or bg is None else center_pixel - bg
    pixel_emp_z = None if pixel_excess is None or not sigma_pix else pixel_excess / sigma_pix
    pixel_ivar_z = None
    center_invvar = None
    if invvar is not None and 0 <= iy < invvar.shape[0] and 0 <= ix < invvar.shape[1]:
        center_invvar = finite(invvar[iy, ix])
        if center_invvar is not None and center_invvar > 0 and pixel_excess is not None:
            pixel_ivar_z = pixel_excess * math.sqrt(center_invvar)

    target_ap = aperture3(data, ix, iy)
    target_ap_excess = None if target_ap is None or bg is None else target_ap - 9.0 * bg
    ap_var = aperture3_variance(invvar, ix, iy)
    ap_ivar_z = None
    if target_ap_excess is not None and ap_var is not None and ap_var > 0:
        ap_ivar_z = target_ap_excess / math.sqrt(ap_var)

    comp_aps = []
    # Frozen comparison centers: all integer centers in the 6..14 pixel annulus,
    # excluding positions whose 3x3 aperture would leave the image.
    for cy in range(1, data.shape[0] - 1):
        for cx in range(1, data.shape[1] - 1):
            r = math.hypot(cx - float(x), cy - float(y))
            if not (ANNULUS_RMIN <= r <= ANNULUS_RMAX):
                continue
            a = aperture3(data, cx, cy)
            if a is not None:
                comp_aps.append(a)
    ap_bg, ap_sigma, n_comp = robust_center_sigma(np.asarray(comp_aps, dtype=float))
    ap_emp_z = None
    if target_ap is not None and ap_bg is not None and ap_sigma:
        ap_emp_z = (target_ap - ap_bg) / ap_sigma

    # Empirical aperture z is primary because coadd inverse-variance pixels are not
    # assumed independent after resampling. Ivar z is preserved as a diagnostic only.
    primary_z = ap_emp_z

    return {
        "band": f"W{band}",
        "image_file": img_path.name,
        "image_sha256": hashlib.sha256(img_path.read_bytes()).hexdigest(),
        "invvar_file": None if ivar_path is None else ivar_path.name,
        "invvar_sha256": None if ivar_path is None else hashlib.sha256(ivar_path.read_bytes()).hexdigest(),
        "shape": [int(data.shape[0]), int(data.shape[1])],
        "target_pixel_float_zero_indexed": {"x": float(x), "y": float(y)},
        "target_pixel_nearest_zero_indexed": {"x": ix, "y": iy},
        "background_annulus_pixels": [ANNULUS_RMIN, ANNULUS_RMAX],
        "background": {"median_pixel": bg, "robust_sigma_pixel": sigma_pix, "sample_count": n_bg},
        "central_pixel": {
            "value": center_pixel,
            "background_subtracted": pixel_excess,
            "empirical_z": pixel_emp_z,
            "inverse_variance": center_invvar,
            "ivar_diagnostic_z": pixel_ivar_z,
        },
        "central_3x3": {
            "sum": target_ap,
            "background_subtracted_using_pixel_median": target_ap_excess,
            "ivar_sum_variance": ap_var,
            "ivar_diagnostic_z": ap_ivar_z,
            "comparison_aperture_median": ap_bg,
            "comparison_aperture_robust_sigma": ap_sigma,
            "comparison_aperture_count": n_comp,
            "empirical_z_primary": ap_emp_z,
        },
        "primary_empirical_z": primary_z,
    }


def main() -> None:
    payload = {
        "schema": "janus.cosmos.stargate_abydos.raw_unwise_center.receipt.v1",
        "experiment_id": "STARGATE-ABYDOS-RAW-UNWISE-FROZEN-CENTER-v1",
        "status": "INCOMPLETE",
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "RAW_COADD_BRIGHTNESS_AT_EXACT_FROZEN_COORDINATE",
        "anomaly_scoring_used": False,
        "frozen_target": {"frame": "ICRS", "ra_deg": RA_DEG, "dec_deg": DEC_DEG},
        "data_source": {"service": "unWISE cutout_fits", "version": VERSION, "url": URL, "bands": ["W1", "W2"]},
        "known_nearest_catalogued_ir_group_arcsec": 8.537689265953338,
        "measurements": [],
        "service_errors": [],
    }

    work = Path(".janus_unwise_raw")
    work.mkdir(exist_ok=True)
    try:
        response = requests.get(URL, timeout=180)
        response.raise_for_status()
        payload["download"] = {
            "http_status": int(response.status_code),
            "bytes": len(response.content),
            "sha256": hashlib.sha256(response.content).hexdigest(),
            "content_type": response.headers.get("content-type"),
        }
        with tarfile.open(fileobj=io.BytesIO(response.content), mode="r:gz") as tf:
            members = [m for m in tf.getmembers() if m.isfile()]
            safe = []
            for m in members:
                p = Path(m.name)
                if p.is_absolute() or ".." in p.parts:
                    continue
                safe.append(m)
            tf.extractall(work, members=safe, filter="data")

        for band in BANDS:
            imgs = sorted(work.rglob(f"*w{band}-img-m.fits"))
            if not imgs:
                payload["service_errors"].append({"band": f"W{band}", "error": "NO_IMAGE_FILE_IN_CUTOUT_ARCHIVE"})
                continue
            for img in imgs:
                ivar_candidate = img.with_name(img.name.replace("-img-m.fits", "-invvar-m.fits"))
                ivar = ivar_candidate if ivar_candidate.exists() else None
                try:
                    payload["measurements"].append(measure_pair(img, ivar, band))
                except Exception as exc:
                    payload["service_errors"].append({"band": f"W{band}", "image": img.name, "error": f"{type(exc).__name__}: {exc}"})
    except Exception as exc:
        payload["service_errors"].append({"service": "unWISE cutout_fits", "error": f"{type(exc).__name__}: {exc}"})

    by_band = {}
    for band in ("W1", "W2"):
        zs = [m["primary_empirical_z"] for m in payload["measurements"] if m["band"] == band and m["primary_empirical_z"] is not None and math.isfinite(m["primary_empirical_z"])]
        by_band[band] = {
            "measurement_count": len(zs),
            "median_primary_empirical_z": float(np.median(zs)) if zs else None,
            "min_primary_empirical_z": float(np.min(zs)) if zs else None,
            "max_primary_empirical_z": float(np.max(zs)) if zs else None,
        }

    z1 = by_band["W1"]["median_primary_empirical_z"]
    z2 = by_band["W2"]["median_primary_empirical_z"]
    complete = z1 is not None and z2 is not None
    if not complete:
        verdict = "INCONCLUSIVE_RAW_COADD_GATE"
        plain = "Raw unWISE gate did not recover usable W1 and W2 coadd measurements for the exact frozen coordinate."
    elif z1 >= STRONG_Z and z2 >= STRONG_Z:
        verdict = "STRONG_RAW_COADD_EXCESS_BOTH_W1_W2__BLEND_IDENTITY_UNRESOLVED"
        plain = "The exact frozen coordinate has a >=5 sigma empirical 3x3 coadd brightness excess in both W1 and W2. Because a catalogued IR group lies only 8.54 arcsec away, this cannot be promoted to an independent centered source without PSF deblending."
    elif z1 >= WEAK_Z or z2 >= WEAK_Z:
        verdict = "WEAK_OR_SINGLE_BAND_RAW_COADD_EXCESS__BLEND_IDENTITY_UNRESOLVED"
        plain = "The exact frozen coordinate shows at least a 3 sigma empirical coadd brightness excess in W1 or W2, but not a frozen strong two-band detection. Nearby-source blending remains unresolved."
    else:
        verdict = "NO_SIGNIFICANT_RAW_COADD_EXCESS_AT_FROZEN_CENTER"
        plain = "The exact frozen coordinate does not show a >=3 sigma empirical 3x3 coadd brightness excess in either W1 or W2 under the frozen raw-coadd test."

    payload["status"] = "COMPLETE" if complete else "INCONCLUSIVE"
    payload["summary"] = {
        "per_band": by_band,
        "verdict": verdict,
        "plain_en": plain,
        "thresholds": {"weak_sigma": WEAK_Z, "strong_sigma": STRONG_Z, "strong_requires_both_bands": True},
    }
    payload["firewall"] = {
        "raw_coadd_excess_is_independent_centered_source": False,
        "raw_coadd_excess_is_real_Abydos": False,
        "raw_coadd_excess_is_planet": False,
        "raw_coadd_excess_is_anomaly": False,
        "nearby_8p54arcsec_source_blending_unresolved": True,
        "coadd_inverse_variance_z_is_diagnostic_not_primary": True,
        "claim_ceiling": "RAW_COADD_BRIGHTNESS_EXCESS_AT_EXACT_FROZEN_COORDINATE_WITHOUT_SOURCE_IDENTITY",
    }
    payload["next_gate"] = "PSF_DEBLEND_FROZEN_CENTER_VS_KNOWN_8P54ARCSEC_IR_NEIGHBOR_IF_RAW_EXCESS_PRESENT"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload["summary"], ensure_ascii=False, indent=2, sort_keys=True))
    if not complete:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
