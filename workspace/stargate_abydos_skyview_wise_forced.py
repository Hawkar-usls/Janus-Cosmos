from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astroquery.skyview import SkyView

OUT = Path("data/stargate/STARGATE-ABYDOS-SKYVIEW-WISE-FORCED-CENTER-v1-LATEST-RECEIPT.json")
TARGET = SkyCoord(223.415064157 * u.deg, 33.979315670 * u.deg, frame="icrs")
SURVEYS = ["WISE 3.4", "WISE 4.6"]
APERTURE_ARCSEC = 4.0
CONTROL_RING_ARCSEC = 45.0
CONTROL_COUNT = 32
THRESHOLD_SIGMA = 5.0


def aperture_flux(data: np.ndarray, x: float, y: float, r_pix: float) -> tuple[float | None, dict[str, Any]]:
    yy, xx = np.indices(data.shape)
    rr = np.hypot(xx - x, yy - y)
    ap = (rr <= r_pix) & np.isfinite(data)
    ann = (rr >= 2.0 * r_pix) & (rr <= 3.0 * r_pix) & np.isfinite(data)
    if ap.sum() < 3 or ann.sum() < 16:
        return None, {"aperture_pixels": int(ap.sum()), "annulus_pixels": int(ann.sum())}
    bg = float(np.nanmedian(data[ann]))
    flux = float(np.nansum(data[ap] - bg))
    return flux, {"aperture_pixels": int(ap.sum()), "annulus_pixels": int(ann.sum()), "background_median": bg}


def measure(hdu, survey: str) -> dict[str, Any]:
    data = np.asarray(hdu.data, dtype=float)
    while data.ndim > 2:
        data = data[0]
    w = WCS(hdu.header).celestial
    x, y = w.world_to_pixel(TARGET)
    scales = proj_plane_pixel_scales(w) * 3600.0
    pixscale = float(np.sqrt(abs(scales[0] * scales[1])))
    r_pix = APERTURE_ARCSEC / pixscale
    source_flux, meta = aperture_flux(data, x, y, r_pix)

    controls = []
    for k in range(CONTROL_COUNT):
        pa = (360.0 * k / CONTROL_COUNT) * u.deg
        c = TARGET.directional_offset_by(pa, CONTROL_RING_ARCSEC * u.arcsec)
        cx, cy = w.world_to_pixel(c)
        if cx < 3*r_pix or cy < 3*r_pix or cx > data.shape[1]-3*r_pix or cy > data.shape[0]-3*r_pix:
            continue
        f, _ = aperture_flux(data, cx, cy, r_pix)
        if f is not None and math.isfinite(f):
            controls.append(float(f))

    median_control = None
    robust_sigma = None
    z = None
    if source_flux is not None and len(controls) >= 8:
        arr = np.asarray(controls, dtype=float)
        median_control = float(np.median(arr))
        mad = float(np.median(np.abs(arr - median_control)))
        robust_sigma = 1.4826 * mad
        if robust_sigma > 0 and math.isfinite(robust_sigma):
            z = float((source_flux - median_control) / robust_sigma)

    digest = hashlib.sha256(np.nan_to_num(data, nan=0.0).astype(np.float32).tobytes()).hexdigest()
    return {
        "survey": survey,
        "pixel_scale_arcsec": pixscale,
        "aperture_radius_arcsec": APERTURE_ARCSEC,
        "source_pixel": [float(x), float(y)],
        "source_flux_native_units_background_subtracted": source_flux,
        "control_ring_arcsec": CONTROL_RING_ARCSEC,
        "control_count_requested": CONTROL_COUNT,
        "control_count_used": len(controls),
        "control_median_flux": median_control,
        "control_robust_sigma": robust_sigma,
        "forced_control_z": z,
        "passes_frozen_5sigma_gate": bool(z is not None and z >= THRESHOLD_SIGMA),
        "image_sha256_float32_nan0": digest,
        "bunit": hdu.header.get("BUNIT"),
        **meta,
    }


def main() -> None:
    results: dict[str, Any] = {}
    errors = []
    for survey in SURVEYS:
        item: dict[str, Any] = {"survey": survey, "status": "UNKNOWN"}
        try:
            images = SkyView.get_images(
                position=TARGET,
                survey=[survey],
                radius=1.5 * u.arcmin,
                pixels=256,
                cache=False,
                show_progress=False,
            )
            if not images:
                item["status"] = "NO_IMAGE"
            else:
                item["status"] = "OK"
                item["measurement"] = measure(images[0][0], survey)
        except Exception as exc:
            item["status"] = "QUERY_ERROR"
            item["error"] = f"{type(exc).__name__}: {exc}"
            errors.append({"survey": survey, "error": item["error"]})
        results[survey] = item

    usable = [item for item in results.values() if item.get("status") == "OK" and item.get("measurement", {}).get("forced_control_z") is not None]
    passes = [item["survey"] for item in usable if item["measurement"]["passes_frozen_5sigma_gate"]]
    if len(usable) < 2:
        verdict = "RAW_IMAGE_GATE_INCONCLUSIVE"
    elif len(passes) == 2:
        verdict = "RAW_W1_W2_BOTH_PASS_5SIGMA__BLEND_UNRESOLVED"
    elif len(passes) == 1:
        verdict = "RAW_SINGLE_BAND_PASS_5SIGMA__BLEND_UNRESOLVED"
    else:
        verdict = "RAW_NEITHER_BAND_PASS_5SIGMA"

    payload = {
        "schema": "janus.cosmos.stargate_abydos.skyview_wise_forced_center.receipt.v1",
        "experiment_id": "STARGATE-ABYDOS-SKYVIEW-WISE-FORCED-CENTER-v1",
        "status": "COMPLETE" if len(usable) == 2 else "INCONCLUSIVE",
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "IMAGE_LEVEL_FORCED_PHOTOMETRY_WITHOUT_RECENTERING",
        "anomaly_scoring_used": False,
        "frozen_target": {"frame": "ICRS", "ra_deg": float(TARGET.ra.deg), "dec_deg": float(TARGET.dec.deg)},
        "frozen_rules": {
            "surveys": SURVEYS,
            "aperture_radius_arcsec": APERTURE_ARCSEC,
            "control_ring_arcsec": CONTROL_RING_ARCSEC,
            "control_count": CONTROL_COUNT,
            "detection_threshold_sigma": THRESHOLD_SIGMA,
            "recenter_for_detection": False,
            "no_post_hoc_radius_or_threshold_changes": True,
        },
        "known_nearest_catalogued_ir_group_arcsec": 8.537689265953338,
        "results": results,
        "service_errors": errors,
        "summary": {
            "usable_band_count": len(usable),
            "surveys_passing_forced_5sigma": passes,
            "verdict": verdict,
            "per_band_z": {s: results[s].get("measurement", {}).get("forced_control_z") for s in SURVEYS},
        },
        "firewall": {
            "raw_aperture_excess_is_independent_source": False,
            "raw_aperture_excess_is_real_Abydos": False,
            "raw_aperture_excess_is_planet": False,
            "raw_aperture_excess_is_anomaly": False,
            "nearby_8p54arcsec_source_blending_unresolved": True,
            "claim_ceiling": "PIXEL_LEVEL_WISE_PERSISTENCE_AT_FROZEN_COORDINATE_WITH_BLEND_WARNING",
        },
        "next_gate": "FIXED_POSITION_PSF_DEBLEND_CENTER_VS_NEARBY_IR_GROUPS_IF_RAW_EXCESS_PRESENT",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload["summary"], ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
