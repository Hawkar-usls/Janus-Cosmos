#!/usr/bin/env python3
"""
JANUS COSMOS CELESTIAL DYNAMICS

Scientific two-body trajectory propagation for planets, moons, asteroids,
comets and spacecraft ephemeris experiments. This module uses Cartesian
state vectors and a velocity-Verlet integrator, then can expand uncertainty
in the initial state into a reproducible ensemble of trajectories.

It is an astronomy/planetary-science module. Model output is a prediction,
not an observation or ephemeris authority.
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
from dataclasses import asdict, dataclass
from typing import Dict, List, Sequence

from janus_world_matrix import DimensionSpec, JanusWorldMatrix, _quantile, _sha256_obj

SCHEMA = "janus.cosmos.celestial-dynamics.v1"
MU_SUN_KM3_S2 = 132_712_440_018.0
AU_KM = 149_597_870.7


@dataclass(frozen=True)
class CartesianState:
    x_km: float
    y_km: float
    z_km: float
    vx_km_s: float
    vy_km_s: float
    vz_km_s: float

    def position(self) -> tuple[float, float, float]:
        return (self.x_km, self.y_km, self.z_km)

    def velocity(self) -> tuple[float, float, float]:
        return (self.vx_km_s, self.vy_km_s, self.vz_km_s)


def _norm3(v: Sequence[float]) -> float:
    return math.sqrt(sum(float(x) * float(x) for x in v))


def acceleration_km_s2(position_km: Sequence[float], mu_km3_s2: float) -> tuple[float, float, float]:
    x, y, z = (float(v) for v in position_km)
    mu = float(mu_km3_s2)
    if not math.isfinite(mu) or mu <= 0:
        raise ValueError("mu_km3_s2 must be finite and > 0")
    r2 = x*x + y*y + z*z
    if not math.isfinite(r2) or r2 <= 0:
        raise ValueError("position radius must be finite and > 0")
    inv_r3 = 1.0 / (r2 * math.sqrt(r2))
    s = -mu * inv_r3
    return (s*x, s*y, s*z)


def step_velocity_verlet(state: CartesianState, mu_km3_s2: float, dt_s: float) -> CartesianState:
    dt = float(dt_s)
    if not math.isfinite(dt) or dt <= 0:
        raise ValueError("dt_s must be finite and > 0")
    r0 = state.position()
    v0 = state.velocity()
    a0 = acceleration_km_s2(r0, mu_km3_s2)
    r1 = tuple(r0[i] + v0[i]*dt + 0.5*a0[i]*dt*dt for i in range(3))
    a1 = acceleration_km_s2(r1, mu_km3_s2)
    v1 = tuple(v0[i] + 0.5*(a0[i] + a1[i])*dt for i in range(3))
    return CartesianState(*r1, *v1)


def specific_orbital_energy(state: CartesianState, mu_km3_s2: float) -> float:
    r = _norm3(state.position())
    v = _norm3(state.velocity())
    return 0.5*v*v - float(mu_km3_s2)/r


def propagate(initial: CartesianState, mu_km3_s2: float, duration_s: float, dt_s: float,
              preserve_trajectory: bool = True) -> Dict[str, object]:
    duration = float(duration_s)
    dt = float(dt_s)
    if duration <= 0 or dt <= 0:
        raise ValueError("duration_s and dt_s must be > 0")
    steps = int(round(duration / dt))
    if steps < 1 or abs(steps*dt - duration) > 1e-9 * max(1.0, duration):
        raise ValueError("duration_s must be an integer multiple of dt_s")
    s = initial
    e0 = specific_orbital_energy(s, mu_km3_s2)
    traj: List[Dict[str, float]] = []
    if preserve_trajectory:
        traj.append({"t_s": 0.0, **asdict(s)})
    for i in range(steps):
        s = step_velocity_verlet(s, mu_km3_s2, dt)
        if preserve_trajectory:
            traj.append({"t_s": (i+1)*dt, **asdict(s)})
    e1 = specific_orbital_energy(s, mu_km3_s2)
    result: Dict[str, object] = {
        "schema": SCHEMA,
        "status": "MODEL_CONDITIONED_CELESTIAL_TRAJECTORY",
        "integrator": "VELOCITY_VERLET_TWO_BODY",
        "mu_km3_s2": float(mu_km3_s2),
        "duration_s": duration,
        "dt_s": dt,
        "steps": steps,
        "initial_state": asdict(initial),
        "final_state": asdict(s),
        "initial_radius_km": _norm3(initial.position()),
        "final_radius_km": _norm3(s.position()),
        "specific_energy_initial_km2_s2": e0,
        "specific_energy_final_km2_s2": e1,
        "specific_energy_relative_drift": abs(e1-e0)/max(abs(e0), 1e-30),
        "prediction_is_observation": False,
        "claim_ceiling": "SCIENTIFIC_TWO_BODY_PROPAGATION__VALIDATE_AGAINST_TRUSTED_EPHEMERIS_FOR_REAL_OBJECTS",
    }
    if preserve_trajectory:
        result["trajectory"] = traj
    result["result_sha256"] = _sha256_obj(result)
    return result


def uncertainty_cloud(base: CartesianState, mu_km3_s2: float, sigmas: Dict[str, float],
                      duration_s: float, dt_s: float, samples: int = 256,
                      seed: int | None = None,
                      namespace: str = "JANUS_COSMOS_CELESTIAL") -> Dict[str, object]:
    dims = (
        DimensionSpec("x_offset_km", 0.0, float(sigmas.get("x_km", 0.0))),
        DimensionSpec("y_offset_km", 0.0, float(sigmas.get("y_km", 0.0))),
        DimensionSpec("z_offset_km", 0.0, float(sigmas.get("z_km", 0.0))),
        DimensionSpec("vx_offset_km_s", 0.0, float(sigmas.get("vx_km_s", 0.0))),
        DimensionSpec("vy_offset_km_s", 0.0, float(sigmas.get("vy_km_s", 0.0))),
        DimensionSpec("vz_offset_km_s", 0.0, float(sigmas.get("vz_km_s", 0.0))),
        DimensionSpec("mu_offset_km3_s2", 0.0, float(sigmas.get("mu_km3_s2", 0.0))),
    )
    cloud = JanusWorldMatrix(dims, samples=samples, seed=seed, namespace=namespace,
                             sampling="latin_hypercube").generate()
    radii: List[float] = []
    rows: List[Dict[str, object]] = []
    for idx, draw in enumerate(cloud["world_matrix"]):
        dx,dy,dz,dvx,dvy,dvz,dmu = (float(x) for x in draw)
        s0 = CartesianState(base.x_km+dx, base.y_km+dy, base.z_km+dz,
                            base.vx_km_s+dvx, base.vy_km_s+dvy, base.vz_km_s+dvz)
        mu = float(mu_km3_s2) + dmu
        if mu <= 0:
            raise ValueError("uncertainty draw produced non-positive mu")
        run = propagate(s0, mu, duration_s, dt_s, preserve_trajectory=False)
        sf = CartesianState(**run["final_state"])
        radius = _norm3(sf.position())
        radii.append(radius)
        rows.append({
            "scenario_id": idx,
            "parameter_draw": {
                "x_offset_km": dx, "y_offset_km":dy, "z_offset_km":dz,
                "vx_offset_km_s":dvx, "vy_offset_km_s":dvy, "vz_offset_km_s":dvz,
                "mu_offset_km3_s2":dmu,
            },
            "final_state": asdict(sf),
            "final_radius_km": radius,
        })
    def stats(vals: Sequence[float]) -> Dict[str,float]:
        return {
            "mean": statistics.fmean(vals), "std": statistics.pstdev(vals),
            "q05": _quantile(vals,.05), "q50":_quantile(vals,.50), "q95":_quantile(vals,.95),
            "min":min(vals), "max":max(vals),
        }
    output: Dict[str, object] = {
        "schema": "janus.cosmos.celestial-dynamics.uncertainty-cloud.v1",
        "status": "MODEL_CONDITIONED_TRAJECTORY_ENSEMBLE",
        "formula": "RESPICIENS_ET_PROSPICIENS",
        "sampling": "DETERMINISTIC_LATIN_HYPERCUBE",
        "seed": cloud["seed"],
        "samples": samples,
        "base_state": asdict(base),
        "mu_km3_s2": float(mu_km3_s2),
        "sigmas": {k:float(v) for k,v in sigmas.items()},
        "duration_s": float(duration_s),
        "dt_s": float(dt_s),
        "final_radius_km": stats(radii),
        "scenarios": rows,
        "prediction_is_observation": False,
        "simulation_to_reality_gap": "NOT_EVALUATED__REQUIRES_COMPARISON_WITH_TRUSTED_EPHEMERIS",
        "epistemic_status": "MODEL_ONLY",
    }
    output["result_sha256"] = _sha256_obj(output)
    return output


def self_test() -> None:
    r = AU_KM
    v = math.sqrt(MU_SUN_KM3_S2/r)
    s0 = CartesianState(r,0.0,0.0,0.0,v,0.0)
    one_day = propagate(s0, MU_SUN_KM3_S2, 86_400.0, 600.0, preserve_trajectory=False)
    assert one_day["specific_energy_relative_drift"] < 1e-9
    assert abs(one_day["final_radius_km"]-r)/r < 1e-6
    cloud_a = uncertainty_cloud(s0, MU_SUN_KM3_S2,
        {"x_km":10.0,"y_km":10.0,"vx_km_s":1e-5,"vy_km_s":1e-5},
        3600.0, 300.0, samples=32, seed=424242)
    cloud_b = uncertainty_cloud(s0, MU_SUN_KM3_S2,
        {"x_km":10.0,"y_km":10.0,"vx_km_s":1e-5,"vy_km_s":1e-5},
        3600.0, 300.0, samples=32, seed=424242)
    assert cloud_a["result_sha256"] == cloud_b["result_sha256"]
    assert cloud_a["prediction_is_observation"] is False
    print("JANUS_CELESTIAL_DYNAMICS_SELF_TEST=PASS")


def main() -> int:
    parser = argparse.ArgumentParser(description="JANUS Cosmos scientific celestial trajectory propagator")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    r = AU_KM
    v = math.sqrt(MU_SUN_KM3_S2/r)
    obj = propagate(CartesianState(r,0,0,0,v,0), MU_SUN_KM3_S2, 86_400.0, 600.0, False)
    print(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
