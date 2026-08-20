#!/usr/bin/env python3
"""Frozen-orientation Kamil virtual-pyramid gate.

No free azimuth rotation is allowed:
  SHAFT_44_5  = az 180 deg, alt 44.5 deg
  FACE_51_84  = az 180 deg, alt 51.84 deg
  APEX_BEAM   = zenith, alt 90 deg

The script evaluates the two inherited January 2026 anchor timestamps,
then an any-time-in-January limit and matched isotropic random-sky controls.
It is a geometry/control experiment only; it does not assert an ancient device.
"""

import json
import math
from datetime import datetime, timezone
import numpy as np

LAT_DEG = 22.0183333333
LON_DEG = 26.0877777778
MC_N = 500_000
MC_SEED = 20260820

TARGETS_J2000 = {
    "LOVE": (204.30267916666668, -36.78240527777778),
    "EDEM": (139.22409686590188, 30.26038779947318),
    "ORION_BELT_CENTROID": (84.08137638816612, -1.1480040013137136),
}
OPERATORS = {
    "SHAFT_44_5_SOUTH": {"az_deg": 180.0, "alt_deg": 44.5},
    "FACE_51_84_SOUTH": {"az_deg": 180.0, "alt_deg": 51.84},
    "APEX_BEAM_ZENITH": {"az_deg": 0.0, "alt_deg": 90.0},
}
ANCHORS = [
    datetime(2026, 1, 5, 2, 3, 34, tzinfo=timezone.utc),
    datetime(2026, 1, 10, 3, 47, 22, tzinfo=timezone.utc),
]


def jd(dt):
    return dt.timestamp() / 86400.0 + 2440587.5


def gmst_deg(dt):
    j = jd(dt)
    t = (j - 2451545.0) / 36525.0
    return (280.46061837 + 360.98564736629 * (j - 2451545.0)
            + 0.000387933 * t * t - t * t * t / 38710000.0) % 360.0


def precess_j2000_to_mean_date(ra_deg, dec_deg, dt):
    t = (jd(dt) - 2451545.0) / 36525.0
    zeta = (2306.2181*t + 0.30188*t*t + 0.017998*t*t*t) / 3600.0
    z = (2306.2181*t + 1.09468*t*t + 0.018203*t*t*t) / 3600.0
    theta = (2004.3109*t - 0.42665*t*t - 0.041833*t*t*t) / 3600.0
    a, d = map(math.radians, (ra_deg, dec_deg))
    zr, zz, tr = map(math.radians, (zeta, z, theta))
    A = math.cos(d) * math.sin(a + zr)
    B = math.cos(tr)*math.cos(d)*math.cos(a + zr) - math.sin(tr)*math.sin(d)
    C = math.sin(tr)*math.cos(d)*math.cos(a + zr) + math.cos(tr)*math.sin(d)
    return ((math.degrees(math.atan2(A, B)) + z) % 360.0,
            math.degrees(math.asin(C)))


def fixed_altaz_to_mean_date_radec(az_deg, alt_deg, dt):
    phi, A, h = map(math.radians, (LAT_DEG, az_deg, alt_deg))
    sd = math.sin(phi)*math.sin(h) + math.cos(phi)*math.cos(h)*math.cos(A)
    d = math.asin(max(-1.0, min(1.0, sd)))
    cd = max(1e-15, abs(math.cos(d)))
    sh = -math.sin(A)*math.cos(h) / cd
    ch = ((math.sin(h) - math.sin(phi)*math.sin(d)) /
          (math.cos(phi)*cd))
    H = math.atan2(sh, ch)
    lst = (gmst_deg(dt) + LON_DEG) % 360.0
    return ((lst - math.degrees(H)) % 360.0, math.degrees(d))


def sep_deg(ra1, dec1, ra2, dec2):
    a1, d1, a2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    dot = math.sin(d1)*math.sin(d2) + math.cos(d1)*math.cos(d2)*math.cos(a1-a2)
    return math.degrees(math.acos(max(-1.0, min(1.0, dot))))


def sep_vec_deg(ra, dec, ra2, dec2):
    r1, d1 = np.radians(ra), np.radians(dec)
    r2, d2 = math.radians(ra2), math.radians(dec2)
    dot = np.sin(d1)*math.sin(d2) + np.cos(d1)*math.cos(d2)*np.cos(r1-r2)
    return np.degrees(np.arccos(np.clip(dot, -1.0, 1.0)))


rng = np.random.default_rng(MC_SEED)
rand_ra = rng.uniform(0.0, 360.0, MC_N)
rand_dec = np.degrees(np.arcsin(rng.uniform(-1.0, 1.0, MC_N)))

out = {
    "schema": "janus.cosmos.kamil.frozen_orientation_gate.v1",
    "site": {"lat_deg": LAT_DEG, "lon_deg_east": LON_DEG},
    "mc": {"n": MC_N, "seed": MC_SEED, "distribution": "isotropic_sphere"},
    "operators": {},
}

for op, geom in OPERATORS.items():
    entry = {"geometry": geom, "anchors": [], "targets": {}}
    beams = []
    for dt in ANCHORS:
        bra, bdec = fixed_altaz_to_mean_date_radec(geom["az_deg"], geom["alt_deg"], dt)
        beams.append((bra, bdec))
        anchor = {"time_utc": dt.isoformat(), "beam_ra_mean_date_deg": bra,
                  "beam_dec_mean_date_deg": bdec, "separations_deg": {}}
        for name, (ra0, dec0) in TARGETS_J2000.items():
            tra, tdec = precess_j2000_to_mean_date(ra0, dec0, dt)
            anchor["separations_deg"][name] = sep_deg(bra, bdec, tra, tdec)
        entry["anchors"].append(anchor)

    random_score = np.minimum.reduce([
        sep_vec_deg(rand_ra, rand_dec, bra, bdec) for bra, bdec in beams
    ])
    mid = datetime(2026, 1, 15, tzinfo=timezone.utc)
    _, beam_dec_mid = fixed_altaz_to_mean_date_radec(geom["az_deg"], geom["alt_deg"], mid)

    for name, (ra0, dec0) in TARGETS_J2000.items():
        obs = min(a["separations_deg"][name] for a in entry["anchors"])
        p_mc = (np.count_nonzero(random_score <= obs) + 1) / (MC_N + 1)
        _, tdec_mid = precess_j2000_to_mean_date(ra0, dec0, mid)
        d = abs(tdec_mid - beam_dec_mid)
        lo, hi = max(-90.0, beam_dec_mid-d), min(90.0, beam_dec_mid+d)
        p_full = 0.5*(math.sin(math.radians(hi))-math.sin(math.radians(lo)))
        entry["targets"][name] = {
            "best_two_anchor_separation_deg": obs,
            "two_anchor_isotropic_mc_fraction_leq": p_mc,
            "best_possible_full_january_sidereal_sweep_deg": d,
            "full_january_isotropic_fraction_leq": p_full,
        }
    out["operators"][op] = entry

out["firewall"] = {
    "mean_of_date_approximation": "IAU-1976 precession + mean sidereal time; no nutation/refraction",
    "posthoc_anchor_warning": "January timestamps came from prior LOVE/EDEM project searches and are not independent significance tests.",
    "physical_device_claimed": False,
    "causal_melt_link_claimed": False,
}
print(json.dumps(out, indent=2, sort_keys=True))
