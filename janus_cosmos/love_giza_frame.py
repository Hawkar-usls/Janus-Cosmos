from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, GCRS, SkyCoord
from astropy.time import Time


def _unit_from_altaz(alt_deg: float, az_deg: float) -> np.ndarray:
    alt = math.radians(float(alt_deg))
    az = math.radians(float(az_deg))
    # Local ENU: azimuth 0=N, 90=E.
    return np.array([
        math.cos(alt) * math.sin(az),
        math.cos(alt) * math.cos(az),
        math.sin(alt),
    ], dtype=float)


def _angle_deg(a: np.ndarray, b: np.ndarray) -> float:
    x = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    return math.degrees(math.acos(max(-1.0, min(1.0, x))))


def _time_grid(start: str, stop: str, step_seconds: int) -> Time:
    t0 = Time(start, scale="utc")
    t1 = Time(stop, scale="utc")
    span_s = float((t1 - t0).to_value(u.s))
    n = int(math.floor(span_s / step_seconds)) + 1
    return t0 + np.arange(n, dtype=float) * step_seconds * u.s


def _surface_parallax_mas(target: SkyCoord, location: EarthLocation, times: Time) -> tuple[np.ndarray, float]:
    """Return Giza-vs-geocenter angular shift in a common GCRS vector space.

    The tiny angle is evaluated with atan2(|a x b|, a dot b), which remains
    stable at microarcsecond scale where arccos of a normalized dot product
    can round to zero in double precision.
    """
    site_pos, _ = location.get_gcrs_posvel(times)
    geocentric = target.transform_to(GCRS(obstime=times))
    star_xyz = np.moveaxis(geocentric.cartesian.xyz.to_value(u.m), 0, -1)
    site_xyz = np.moveaxis(site_pos.xyz.to_value(u.m), 0, -1)
    topo_xyz = star_xyz - site_xyz

    cross_norm = np.linalg.norm(np.cross(star_xyz, topo_xyz), axis=1)
    dot_raw = np.sum(star_xyz * topo_xyz, axis=1)
    angle_rad = np.arctan2(cross_norm, dot_raw)
    mas = np.degrees(angle_rad) * 3600.0 * 1000.0

    baseline_m = float(np.max(np.linalg.norm(site_xyz, axis=1)))
    distance_m = float(target.distance.to_value(u.m))
    upper_bound_mas = math.degrees(math.asin(min(1.0, baseline_m / distance_m))) * 3600.0 * 1000.0
    return mas, upper_bound_mas


def run(prereg_path: Path, output_path: Path) -> dict:
    prereg = json.loads(prereg_path.read_text(encoding="utf-8"))
    target_cfg = prereg["target"]
    obs_cfg = prereg["observer"]
    window = prereg["time_window"]
    atm = prereg["atmosphere"]
    pyramid = prereg["pyramid_idealized_geometry"]

    location = EarthLocation.from_geodetic(
        lon=float(obs_cfg["longitude_deg_east"]) * u.deg,
        lat=float(obs_cfg["latitude_deg"]) * u.deg,
        height=float(obs_cfg["height_m"]) * u.m,
    )
    target = SkyCoord(
        ra=float(target_cfg["ra_deg_icrs"]) * u.deg,
        dec=float(target_cfg["dec_deg_icrs"]) * u.deg,
        distance=float(target_cfg["distance_pc"]) * u.pc,
        frame="icrs",
    )
    times = _time_grid(window["start_utc"], window["stop_utc"], int(window["step_seconds"]))
    altaz = target.transform_to(AltAz(
        obstime=times,
        location=location,
        pressure=float(atm["pressure_hpa"]) * u.hPa,
    ))
    alt = np.asarray(altaz.alt.deg, dtype=float)
    az = np.asarray(altaz.az.deg, dtype=float)

    giza_vs_geocenter_mas, theoretical_upper_bound_mas = _surface_parallax_mas(target, location, times)

    visible = alt > 0.0
    imax = int(np.nanargmax(alt))
    visibility = {
        "visible_sample_count": int(np.count_nonzero(visible)),
        "first_above_horizon_utc": times[np.flatnonzero(visible)[0]].utc.isot if np.any(visible) else None,
        "last_above_horizon_utc": times[np.flatnonzero(visible)[-1]].utc.isot if np.any(visible) else None,
        "maximum_altitude_deg": float(alt[imax]),
        "maximum_altitude_azimuth_deg": float(az[imax]),
        "maximum_altitude_time_utc": times[imax].utc.isot,
    }

    normal_alt = float(pyramid["face_outward_normal_altitude_deg"])
    face_results = {}
    for face, face_az in pyramid["face_outward_azimuths_deg"].items():
        normal = _unit_from_altaz(normal_alt, float(face_az))
        angles = np.full(len(times), np.nan, dtype=float)
        front = np.zeros(len(times), dtype=bool)
        for i, (a, z) in enumerate(zip(alt, az)):
            if a <= 0:
                continue
            source_vec = _unit_from_altaz(float(a), float(z))
            dot = float(np.dot(source_vec, normal))
            front[i] = dot > 0.0
            if front[i]:
                angles[i] = _angle_deg(source_vec, normal)
        if np.any(np.isfinite(angles)):
            j = int(np.nanargmin(angles))
            face_results[face] = {
                "minimum_source_to_outward_normal_angle_deg": float(angles[j]),
                "time_utc": times[j].utc.isot,
                "source_altitude_deg": float(alt[j]),
                "source_azimuth_deg": float(az[j]),
                "face_outward_normal_altitude_deg": normal_alt,
                "face_outward_normal_azimuth_deg": float(face_az),
                "front_hemisphere": True,
            }
        else:
            face_results[face] = {
                "minimum_source_to_outward_normal_angle_deg": None,
                "time_utc": None,
                "front_hemisphere": False,
            }

    ranked_faces = sorted(
        (
            (face, data["minimum_source_to_outward_normal_angle_deg"])
            for face, data in face_results.items()
            if data["minimum_source_to_outward_normal_angle_deg"] is not None
        ),
        key=lambda x: x[1],
    )

    result = {
        "schema": "janus.cosmos.love.giza_frame.result.v1.2",
        "experiment_id": prereg["experiment_id"],
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_label": target_cfg["analysis_label"],
        "frame_policy": prereg["frame_policy"],
        "observer": obs_cfg,
        "target_icrs": {
            "ra_deg": float(target_cfg["ra_deg_icrs"]),
            "dec_deg": float(target_cfg["dec_deg_icrs"]),
            "distance_pc": float(target_cfg["distance_pc"]),
        },
        "time_window": window,
        "giza_topocentric_visibility": visibility,
        "observer_shift": {
            "meaning": "Angular difference between Earth-geocentric and Giza-surface lines of sight, computed by subtracting the Giza site vector in common GCRS axes.",
            "minimum_mas": float(np.nanmin(giza_vs_geocenter_mas)),
            "maximum_mas": float(np.nanmax(giza_vs_geocenter_mas)),
            "median_mas": float(np.nanmedian(giza_vs_geocenter_mas)),
            "theoretical_surface_baseline_upper_bound_mas": float(theoretical_upper_bound_mas),
            "cannot_explain_arcsecond_or_degree_scale_target_changes": bool(theoretical_upper_bound_mas < 1.0),
        },
        "idealized_pyramid_face_alignment": {
            "face_inclination_deg_above_horizontal": float(pyramid["face_inclination_deg_above_horizontal"]),
            "face_outward_normal_altitude_deg": normal_alt,
            "faces": face_results,
            "best_face": ranked_faces[0][0] if ranked_faces else None,
            "best_face_minimum_angle_deg": float(ranked_faces[0][1]) if ranked_faces else None,
            "interpretation": "Pure topocentric geometry only. This is not evidence of an ancient optical mechanism or a planetary discovery.",
        },
        "love_gate": {
            "reserved_codename": target_cfg["reserved_post_gate_codename"],
            "candidate_activated": False,
            "reason": "A coordinate-frame or pyramid-face alignment cannot activate the planetary codename gate.",
            "next_gate": prereg["admission"]["next_gate"],
        },
        "claim_ceiling": prereg["admission"]["geometry_claim_ceiling"],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "experiment_id": result["experiment_id"],
        "observer": result["observer"]["frame_name"],
        "max_alt_deg": result["giza_topocentric_visibility"]["maximum_altitude_deg"],
        "max_surface_shift_mas": result["observer_shift"]["maximum_mas"],
        "surface_shift_upper_bound_mas": result["observer_shift"]["theoretical_surface_baseline_upper_bound_mas"],
        "best_face": result["idealized_pyramid_face_alignment"]["best_face"],
        "best_face_angle_deg": result["idealized_pyramid_face_alignment"]["best_face_minimum_angle_deg"],
        "love_candidate_activated": result["love_gate"]["candidate_activated"],
    }, indent=2))
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description="Project frozen TARGET_A ICRS direction into a Giza pyramid-centric topocentric frame.")
    ap.add_argument("--prereg", default="data/love/LAVE_TO_LOVE_GIZA_FRAME_PREREG.json")
    ap.add_argument("--output", default="results/love_giza_frame/love-giza-frame-result.json")
    args = ap.parse_args()
    run(Path(args.prereg), Path(args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
