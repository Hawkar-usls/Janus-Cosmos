from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, ICRS, SkyCoord, get_body
from astropy.time import Time


def _time_grid(start: str, stop: str, step_seconds: int) -> Time:
    t0 = Time(start, scale="utc")
    t1 = Time(stop, scale="utc")
    span_s = float((t1 - t0).to_value(u.s))
    n = int(math.floor(span_s / step_seconds)) + 1
    return t0 + np.arange(n, dtype=float) * step_seconds * u.s


def _unit_from_altaz_arrays(alt_deg: np.ndarray, az_deg: np.ndarray) -> np.ndarray:
    alt = np.deg2rad(np.asarray(alt_deg, dtype=float))
    az = np.deg2rad(np.asarray(az_deg, dtype=float))
    return np.column_stack([
        np.cos(alt) * np.sin(az),  # East
        np.cos(alt) * np.cos(az),  # North
        np.sin(alt),               # Up
    ])


def _unit_from_altaz(alt_deg: float, az_deg: float) -> np.ndarray:
    return _unit_from_altaz_arrays(np.array([alt_deg]), np.array([az_deg]))[0]


def _altaz_from_unit(vec: np.ndarray) -> tuple[float, float]:
    v = np.asarray(vec, dtype=float)
    v = v / np.linalg.norm(v)
    alt = math.degrees(math.asin(float(np.clip(v[2], -1.0, 1.0))))
    az = math.degrees(math.atan2(float(v[0]), float(v[1]))) % 360.0
    return alt, az


def _angle_deg_rows(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    aa = np.asarray(a, dtype=float)
    bb = np.asarray(b, dtype=float)
    dot = np.sum(aa * bb, axis=1) / (np.linalg.norm(aa, axis=1) * np.linalg.norm(bb, axis=1))
    return np.rad2deg(np.arccos(np.clip(dot, -1.0, 1.0)))


def _radec_strings(coord: SkyCoord) -> dict:
    icrs = coord.icrs
    return {
        "ra_deg": float(icrs.ra.deg),
        "dec_deg": float(icrs.dec.deg),
        "ra_hms": icrs.ra.to_string(unit=u.hourangle, sep=":", precision=3, pad=True),
        "dec_dms": icrs.dec.to_string(unit=u.deg, sep=":", precision=3, alwayssign=True, pad=True),
    }


def _evaluate(times: Time, cfg: dict, location: EarthLocation, love: SkyCoord) -> list[dict]:
    frame = AltAz(obstime=times, location=location, pressure=0.0 * u.hPa)
    love_altaz = love.transform_to(frame)
    moon = get_body("moon", times, location=location)
    moon_altaz = moon.transform_to(frame)

    love_alt = np.asarray(love_altaz.alt.deg, dtype=float)
    love_az = np.asarray(love_altaz.az.deg, dtype=float)
    moon_alt = np.asarray(moon_altaz.alt.deg, dtype=float)
    moon_az = np.asarray(moon_altaz.az.deg, dtype=float)
    love_vec = _unit_from_altaz_arrays(love_alt, love_az)
    moon_vec = _unit_from_altaz_arrays(moon_alt, moon_az)

    normal_alt = float(cfg["khufu_faces"]["face_normal_altitude_deg"])
    out = []
    for face, face_az in cfg["khufu_faces"]["outward_normal_azimuth_deg"].items():
        n = _unit_from_altaz(normal_alt, float(face_az))
        love_dot = love_vec @ n
        # Reciprocal source direction. If an incoming ray from s_edem hits the
        # facet, specular reflection leaves the facet toward s_love.
        edem_vec = 2.0 * love_dot[:, None] * n[None, :] - love_vec
        edem_vec /= np.linalg.norm(edem_vec, axis=1)[:, None]
        edem_alt = np.rad2deg(np.arcsin(np.clip(edem_vec[:, 2], -1.0, 1.0)))
        moon_dot = moon_vec @ n
        edem_dot = edem_vec @ n
        sep = _angle_deg_rows(edem_vec, moon_vec)
        valid = (
            (love_alt > 0.0)
            & (moon_alt > 0.0)
            & (edem_alt > 0.0)
            & (love_dot > 0.0)
            & (moon_dot > 0.0)
            & (edem_dot > 0.0)
        )
        if not np.any(valid):
            out.append({"face": face, "status": "NO_VALID_SAMPLES"})
            continue
        valid_idx = np.flatnonzero(valid)
        j = int(valid_idx[np.argmin(sep[valid_idx])])
        edem_alt_j, edem_az_j = _altaz_from_unit(edem_vec[j])
        local_frame = AltAz(obstime=times[j], location=location, pressure=0.0 * u.hPa)
        edem_icrs = SkyCoord(az=edem_az_j * u.deg, alt=edem_alt_j * u.deg, frame=local_frame).transform_to(ICRS())
        moon_icrs = moon[j].transform_to(ICRS())
        out.append({
            "face": face,
            "status": "OK",
            "index": j,
            "time_utc": times[j].utc.isot,
            "angular_error_to_moon_deg": float(sep[j]),
            "love_alt_deg": float(love_alt[j]),
            "love_az_deg": float(love_az[j]),
            "moon_alt_deg": float(moon_alt[j]),
            "moon_az_deg": float(moon_az[j]),
            "edem_alt_deg": float(edem_alt_j),
            "edem_az_deg": float(edem_az_j),
            "edem_icrs": _radec_strings(edem_icrs),
            "moon_icrs": _radec_strings(moon_icrs),
        })
    return out


def _best(rows: list[dict]) -> dict:
    ok = [r for r in rows if r.get("status") == "OK"]
    if not ok:
        raise RuntimeError("No physically valid Giza reverse-spear samples in search window")
    return min(ok, key=lambda r: r["angular_error_to_moon_deg"])


def run(prereg_path: Path, output_path: Path) -> dict:
    cfg = json.loads(prereg_path.read_text(encoding="utf-8"))
    observer = cfg["observer"]
    location = EarthLocation.from_geodetic(
        lon=float(observer["longitude_deg_east"]) * u.deg,
        lat=float(observer["latitude_deg"]) * u.deg,
        height=float(observer["height_m"]) * u.m,
    )
    target = cfg["love_target"]
    love = SkyCoord(
        ra=float(target["ra_deg_icrs"]) * u.deg,
        dec=float(target["dec_deg_icrs"]) * u.deg,
        distance=float(target["distance_pc"]) * u.pc,
        frame="icrs",
    )

    tw = cfg["time_window"]
    coarse_times = _time_grid(tw["start_utc"], tw["stop_utc"], int(tw["coarse_step_seconds"]))
    coarse_rows = _evaluate(coarse_times, cfg, location, love)
    coarse_best = _best(coarse_rows)

    center = Time(coarse_best["time_utc"], scale="utc")
    half = int(tw["refine_half_window_seconds"])
    step = int(tw["refine_step_seconds"])
    offsets = np.arange(-half, half + step, step, dtype=float) * u.s
    refine_times = center + offsets
    refine_rows = _evaluate(refine_times, cfg, location, love)
    refined_best = _best(refine_rows)

    err = float(refined_best["angular_error_to_moon_deg"])
    gates = cfg["reverse_spear_geometry"]
    if err <= float(gates["exact_center_gate_deg"]):
        gate = "EXACT_CENTER_GATE_PASS"
    elif err <= float(gates["strong_gate_deg"]):
        gate = "STRONG_GATE_PASS"
    elif err <= float(gates["primary_gate_deg"]):
        gate = "PRIMARY_GATE_PASS"
    else:
        gate = "GEOMETRY_GATE_FAIL"

    result = {
        "schema": "janus.cosmos.edem.giza_reverse_spear.result.v1",
        "experiment_id": cfg["experiment_id"],
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_label": cfg["analysis_label"],
        "observer": observer,
        "frame_policy": {
            "all_reflection_geometry_observer": "GIZA_GREAT_PYRAMID_APPROX_CENTER",
            "love_projected_topocentrically_at_finite_distance": True,
            "moon_projected_topocentrically_from_giza": True,
            "reported_celestial_coordinate_frame": "ICRS",
            "telescope_locations_used_for_geometry": False,
        },
        "coarse_best": coarse_best,
        "refined_best": refined_best,
        "edem_geometry_direction_candidate": refined_best["edem_icrs"],
        "gate_status": gate,
        "coordinate_interpretation": (
            "This is the Giza-anchored reciprocal sky direction implied by the frozen Love vector, "
            "Khufu face geometry, and the preregistered Moon selector. It is a directional coordinate, "
            "not an identification of a physical object named Eden/Edem."
        ),
        "edem_identity_confirmed": False,
        "claim_ceiling": cfg["claim_ceiling"],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "experiment_id": result["experiment_id"],
        "gate_status": result["gate_status"],
        "face": refined_best["face"],
        "time_utc": refined_best["time_utc"],
        "error_deg": refined_best["angular_error_to_moon_deg"],
        "edem_ra_deg": refined_best["edem_icrs"]["ra_deg"],
        "edem_dec_deg": refined_best["edem_icrs"]["dec_deg"],
        "edem_ra_hms": refined_best["edem_icrs"]["ra_hms"],
        "edem_dec_dms": refined_best["edem_icrs"]["dec_dms"],
    }, indent=2))
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description="Resolve a Giza-anchored reciprocal Edem sky direction from the frozen Love vector")
    ap.add_argument("--prereg", default="data/love/EDEM_GIZA_REVERSE_SPEAR_PREREG.json")
    ap.add_argument("--output", default="results/edem_giza_reverse_spear/edem-giza-result.json")
    args = ap.parse_args()
    run(Path(args.prereg), Path(args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
