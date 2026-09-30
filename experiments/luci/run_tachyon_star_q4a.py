#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
import urllib.request
from pathlib import Path

from astropy.io import fits
from astropy.wcs import WCS

FRAMES = (
    ("A_BEFORE", "luci1.20220319.0114.fits.gz", "2022-03-19T11:44:15.5176"),
    ("B_CANDIDATE", "luci1.20220319.0116.fits.gz", "2022-03-19T11:48:08.1003"),
    ("C_AFTER", "luci1.20220319.0120.fits.gz", "2022-03-19T11:55:20.3598"),
)
BASE_URL = "https://archive.lbto.org/files/lbt/"
CANDIDATE_X = 1437.216755721343
CANDIDATE_Y = 2038.7227543099102
PIXEL_ORIGIN = 0
ROUNDTRIP_TOL_PX = 0.001
PATCH_HALF = 8
CONTROL_X = (128.0, 256.0, 384.0, 512.0, 640.0, 768.0, 896.0, 1024.0, 1152.0, 1280.0, 1664.0, 1792.0, 1920.0)
CONTROL_Y = CANDIDATE_Y
CLAIM = "HEADER_WCS_GEOMETRY_FREEZE_ONLY__NO_NEW_PIXEL_MEASUREMENT__NO_ORIGIN_CONCLUSION__NO_TACHYON_OR_FTL_CLAIM"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download_exact(filename: str, cache: Path) -> Path:
    cache.mkdir(parents=True, exist_ok=True)
    out = cache / filename
    if out.exists() and out.stat().st_size > 0:
        return out
    req = urllib.request.Request(BASE_URL + filename, headers={"User-Agent": "Janus-Cosmos-Q4A/1.0"})
    with urllib.request.urlopen(req, timeout=300) as r:
        out.write_bytes(r.read())
    return out


def header_wcs(path: Path) -> tuple[fits.Header, WCS]:
    # Header-only access. Q4A intentionally never dereferences an image HDU .data array.
    hdr = fits.getheader(path, ext=0)
    w = WCS(hdr)
    if not w.has_celestial:
        raise RuntimeError(f"HDU0 WCS is not celestial: {path.name}")
    return hdr, w


def pix_to_sky(w: WCS, x: float, y: float) -> tuple[float, float]:
    world = w.all_pix2world([[float(x), float(y)]], PIXEL_ORIGIN)[0]
    ra, dec = float(world[0]), float(world[1])
    if not (math.isfinite(ra) and math.isfinite(dec)):
        raise RuntimeError("nonfinite pix->sky coordinate")
    return ra, dec


def sky_to_pix(w: WCS, ra: float, dec: float) -> tuple[float, float]:
    pix = w.all_world2pix([[float(ra), float(dec)]], PIXEL_ORIGIN)[0]
    x, y = float(pix[0]), float(pix[1])
    if not (math.isfinite(x) and math.isfinite(y)):
        raise RuntimeError("nonfinite sky->pix coordinate")
    return x, y


def inside_patch(shape: tuple[int, int], x: float, y: float, half: int = PATCH_HALF) -> bool:
    h, w = shape
    return bool(half <= x <= (w - 1 - half) and half <= y <= (h - 1 - half))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", default="results/tachyon_star_q4a")
    ap.add_argument("--cache-dir", default=".cache/tachyon_star_q4a")
    args = ap.parse_args()
    out = Path(args.output_dir)
    cache = Path(args.cache_dir)
    out.mkdir(parents=True, exist_ok=True)

    frame_rows = []
    frame_map = {}
    blocked = []
    for role, filename, frozen_date in FRAMES:
        try:
            path = download_exact(filename, cache)
            hdr, w = header_wcs(path)
            n1, n2 = int(hdr.get("NAXIS1", 0)), int(hdr.get("NAXIS2", 0))
            if n1 <= 0 or n2 <= 0:
                raise RuntimeError("invalid NAXIS1/NAXIS2")
            row = {
                "role": role,
                "filename": filename,
                "url": BASE_URL + filename,
                "fits_gz_sha256": sha256_file(path),
                "header_date_obs": str(hdr.get("DATE-OBS", "")),
                "frozen_date_obs": frozen_date,
                "instrument": str(hdr.get("INSTRUME", "")),
                "filter1": str(hdr.get("FILTER1", hdr.get("FILTER", ""))),
                "filter2": str(hdr.get("FILTER2", "")),
                "exptime": hdr.get("EXPTIME"),
                "naxis1": n1,
                "naxis2": n2,
                "wcs_ctype": [str(q) for q in w.wcs.ctype],
                "image_pixel_array_dereferenced": False,
            }
            frame_rows.append(row)
            frame_map[role] = {"path": path, "hdr": hdr, "wcs": w, "shape": (n2, n1), "row": row}
        except Exception as exc:
            blocked.append({"role": role, "filename": filename, "error": f"{type(exc).__name__}: {exc}"})

    if blocked or len(frame_map) != 3:
        rec = {
            "schema": "janus.cosmos.tachyon_star.q4a.receipt.v1",
            "experiment_id": "JANUS-TACHYON-STAR-Q4A-LUCI-HEADER-WCS-GEOMETRY-FREEZE",
            "status": "BLOCKED_GEOMETRY_INCOMPLETE",
            "blocked": blocked,
            "frames": frame_rows,
            "new_forced_photometry_performed": False,
            "image_pixel_arrays_dereferenced": False,
            "claim_ceiling": CLAIM,
        }
        (out / "receipt.json").write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(rec, indent=2, sort_keys=True))
        return 3

    wb = frame_map["B_CANDIDATE"]["wcs"]
    ra, dec = pix_to_sky(wb, CANDIDATE_X, CANDIDATE_Y)
    bx, by = sky_to_pix(wb, ra, dec)
    roundtrip = math.hypot(bx - CANDIDATE_X, by - CANDIDATE_Y)
    if roundtrip > ROUNDTRIP_TOL_PX:
        blocked.append({"reason": "B_CANDIDATE_WCS_ROUNDTRIP_FAIL", "roundtrip_error_px": roundtrip})

    sky_track = []
    detector_track = []
    for role, filename, _ in FRAMES:
        fm = frame_map[role]
        sx, sy = sky_to_pix(fm["wcs"], ra, dec)
        sky_inside = inside_patch(fm["shape"], sx, sy)
        det_inside = inside_patch(fm["shape"], CANDIDATE_X, CANDIDATE_Y)
        if not sky_inside:
            blocked.append({"role": role, "reason": "SKY_TRACK_PATCH_OUT_OF_BOUNDS", "x": sx, "y": sy})
        if not det_inside:
            blocked.append({"role": role, "reason": "DETECTOR_TRACK_PATCH_OUT_OF_BOUNDS", "x": CANDIDATE_X, "y": CANDIDATE_Y})
        sky_track.append({
            "role": role,
            "filename": filename,
            "x_zero_based": sx,
            "y_zero_based": sy,
            "patch_17x17_inside": sky_inside,
            "separation_from_detector_fixed_px": math.hypot(sx - CANDIDATE_X, sy - CANDIDATE_Y),
        })
        dra, ddec = pix_to_sky(fm["wcs"], CANDIDATE_X, CANDIDATE_Y)
        detector_track.append({
            "role": role,
            "filename": filename,
            "x_zero_based": CANDIDATE_X,
            "y_zero_based": CANDIDATE_Y,
            "ra_deg_at_detector_pixel": dra,
            "dec_deg_at_detector_pixel": ddec,
            "patch_17x17_inside": det_inside,
        })

    controls = []
    shape_b = frame_map["B_CANDIDATE"]["shape"]
    for x in CONTROL_X:
        ok = inside_patch(shape_b, x, CONTROL_Y)
        if not ok:
            blocked.append({"reason": "CONTROL_PATCH_OUT_OF_BOUNDS", "x": x, "y": CONTROL_Y})
        controls.append({
            "x_zero_based": x,
            "y_zero_based": CONTROL_Y,
            "patch_17x17_inside": ok,
            "candidate_distance_px": abs(x - CANDIDATE_X),
        })

    status = "PASS_GEOMETRY_FROZEN" if not blocked else "BLOCKED_GEOMETRY_INCOMPLETE"
    manifest = {
        "schema": "janus.cosmos.tachyon_star.q4a.geometry_manifest.v1",
        "experiment_id": "JANUS-TACHYON-STAR-Q4A-LUCI-HEADER-WCS-GEOMETRY-FREEZE",
        "candidate_frame": "luci1.20220319.0116.fits.gz",
        "candidate_centroid_zero_based": [CANDIDATE_X, CANDIDATE_Y],
        "candidate_sky_coordinate_deg": {"ra": ra, "dec": dec},
        "candidate_roundtrip_error_px": roundtrip,
        "pixel_origin": PIXEL_ORIGIN,
        "patch_size_px": 17,
        "sky_fixed_track": sky_track,
        "detector_fixed_track": detector_track,
        "detector_controls": controls,
        "frame_headers": frame_rows,
        "image_pixel_arrays_dereferenced": False,
        "new_forced_photometry_performed": False,
    }
    manifest_path = out / "q4a_geometry_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    rec = {
        "schema": "janus.cosmos.tachyon_star.q4a.receipt.v1",
        "experiment_id": "JANUS-TACHYON-STAR-Q4A-LUCI-HEADER-WCS-GEOMETRY-FREEZE",
        "status": status,
        "blocked": blocked,
        "frames_retrieved": len(frame_rows),
        "candidate_sky_coordinate_deg": {"ra": ra, "dec": dec},
        "candidate_roundtrip_error_px": roundtrip,
        "control_count": len(controls),
        "sky_detector_separation_px": {q["role"]: q["separation_from_detector_fixed_px"] for q in sky_track},
        "geometry_manifest_sha256": sha256_file(manifest_path),
        "new_forced_photometry_performed": False,
        "image_pixel_arrays_dereferenced": False,
        "claim_ceiling": CLAIM,
    }
    (out / "receipt.json").write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(rec, indent=2, sort_keys=True))
    return 0 if status == "PASS_GEOMETRY_FROZEN" else 3


if __name__ == "__main__":
    raise SystemExit(main())
