#!/usr/bin/env python3
"""
JANUS COSMOS WORLD MATRIX
Respiciens et Prospiciens: a neutral Monte Carlo environment/observation engine.

Purpose:
- expand uncertainty in observation parameters into many possible worlds;
- summarize the ensemble before downstream interpretation;
- keep simulation distinct from evidence.

This module is intentionally domain-neutral and contains no weapon, aiming,
projectile, firing, or harm-optimization logic.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import statistics
from dataclasses import dataclass, asdict
from typing import Dict, List, Sequence


FORBIDDEN_PARAMETER_FRAGMENTS = (
    "ballistic",
    "caliber",
    "calibre",
    "projectile",
    "muzzle",
    "bullet",
    "weapon",
    "rifle",
    "firearm",
    "gun_",
    "firing",
    "shot_placement",
    "impact_point",
    "aiming",
)


def _guard_name(name: str) -> None:
    lowered = name.lower()
    if any(fragment in lowered for fragment in FORBIDDEN_PARAMETER_FRAGMENTS):
        raise ValueError(f"unsupported parameter domain: {name!r}")


def _quantile(values: Sequence[float], q: float) -> float:
    if not values:
        raise ValueError("quantile requires at least one value")
    if not 0.0 <= q <= 1.0:
        raise ValueError("q must be in [0, 1]")
    xs = sorted(float(v) for v in values)
    if len(xs) == 1:
        return xs[0]
    pos = (len(xs) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return xs[lo]
    frac = pos - lo
    return xs[lo] * (1.0 - frac) + xs[hi] * frac


@dataclass(frozen=True)
class DimensionSpec:
    name: str
    mean: float
    sigma: float
    low: float | None = None
    high: float | None = None
    distribution: str = "normal"

    def validate(self) -> None:
        _guard_name(self.name)
        if self.sigma < 0:
            raise ValueError(f"{self.name}: sigma must be >= 0")
        if self.low is not None and self.high is not None and self.low > self.high:
            raise ValueError(f"{self.name}: low must be <= high")
        if self.distribution not in {"normal", "uniform"}:
            raise ValueError(f"{self.name}: unsupported distribution {self.distribution!r}")

    def draw(self, rng: random.Random) -> float:
        self.validate()
        if self.distribution == "uniform":
            half = self.sigma * math.sqrt(3.0)
            value = rng.uniform(self.mean - half, self.mean + half)
        else:
            value = rng.gauss(self.mean, self.sigma)
        if self.low is not None:
            value = max(self.low, value)
        if self.high is not None:
            value = min(self.high, value)
        return float(value)


def deterministic_seed(specs: Sequence[DimensionSpec], samples: int, namespace: str) -> int:
    payload = {
        "namespace": namespace,
        "samples": int(samples),
        "dimensions": [asdict(s) for s in specs],
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return int.from_bytes(hashlib.sha256(canonical).digest()[:8], "big")


def _covariance_matrix(rows: Sequence[Sequence[float]]) -> List[List[float]]:
    if not rows:
        return []
    n_dims = len(rows[0])
    cols = [[float(row[j]) for row in rows] for j in range(n_dims)]
    means = [statistics.fmean(c) for c in cols]
    denom = max(1, len(rows) - 1)
    out: List[List[float]] = []
    for i in range(n_dims):
        r = []
        for j in range(n_dims):
            cov = sum((a - means[i]) * (b - means[j]) for a, b in zip(cols[i], cols[j])) / denom
            r.append(cov)
        out.append(r)
    return out


class JanusWorldMatrix:
    """Generate and summarize an ensemble of plausible observational worlds."""

    schema = "janus.cosmos.world-matrix.v1"
    formula = "RESPICIENS_ET_PROSPICIENS"

    def __init__(
        self,
        dimensions: Sequence[DimensionSpec],
        samples: int = 1000,
        seed: int | None = None,
        namespace: str = "JANUS_COSMOS",
    ) -> None:
        if not dimensions:
            raise ValueError("at least one dimension is required")
        if samples < 2:
            raise ValueError("samples must be >= 2")
        self.dimensions = tuple(dimensions)
        for spec in self.dimensions:
            spec.validate()
        self.samples = int(samples)
        self.namespace = str(namespace)
        self.seed = int(seed) if seed is not None else deterministic_seed(self.dimensions, self.samples, self.namespace)

    def generate(self) -> Dict[str, object]:
        rng = random.Random(self.seed)
        matrix = [[spec.draw(rng) for spec in self.dimensions] for _ in range(self.samples)]
        names = [spec.name for spec in self.dimensions]
        summary: Dict[str, Dict[str, float]] = {}
        for idx, name in enumerate(names):
            col = [row[idx] for row in matrix]
            summary[name] = {
                "mean": statistics.fmean(col),
                "std": statistics.pstdev(col),
                "q05": _quantile(col, 0.05),
                "q25": _quantile(col, 0.25),
                "q50": _quantile(col, 0.50),
                "q75": _quantile(col, 0.75),
                "q95": _quantile(col, 0.95),
                "min": min(col),
                "max": max(col),
            }

        return {
            "schema": self.schema,
            "formula": self.formula,
            "namespace": self.namespace,
            "seed": self.seed,
            "samples": self.samples,
            "dimensions": names,
            "world_matrix": matrix,
            "summary": summary,
            "covariance": _covariance_matrix(matrix),
            "epistemic_status": {
                "simulation_is_evidence": False,
                "simulation_is_discovery": False,
                "requires_observed_data_comparison": True,
                "requires_null_or_control": True,
                "priors_should_be_frozen_before_result_inspection": True
            },
            "domain_boundary": {
                "astronomy_and_planetary_search": True,
                "sensor_uncertainty": True,
                "environment_uncertainty": True,
                "weapon_domain": False,
                "aiming_or_harm_optimization": False
            }
        }


def sky_preset() -> Sequence[DimensionSpec]:
    return (
        DimensionSpec("ra_offset_arcsec", 0.0, 0.7),
        DimensionSpec("dec_offset_arcsec", 0.0, 0.7),
        DimensionSpec("time_offset_s", 0.0, 1.5),
        DimensionSpec("seeing_arcsec", 1.4, 0.35, low=0.2),
        DimensionSpec("sensor_noise_sigma", 1.0, 0.15, low=0.01),
        DimensionSpec("photometric_zero_point_offset_mag", 0.0, 0.03),
    )


def planet_preset() -> Sequence[DimensionSpec]:
    return (
        DimensionSpec("latitude_offset_km", 0.0, 0.25),
        DimensionSpec("longitude_offset_km", 0.0, 0.25),
        DimensionSpec("local_time_offset_s", 0.0, 2.0),
        DimensionSpec("terrain_height_error_m", 0.0, 3.0),
        DimensionSpec("albedo_noise_fraction", 0.0, 0.025),
        DimensionSpec("sensor_noise_sigma", 1.0, 0.15, low=0.01),
    )


def self_test() -> None:
    engine_a = JanusWorldMatrix(sky_preset(), samples=64, seed=424242)
    engine_b = JanusWorldMatrix(sky_preset(), samples=64, seed=424242)
    a = engine_a.generate()
    b = engine_b.generate()
    assert a["world_matrix"] == b["world_matrix"]
    assert len(a["world_matrix"]) == 64
    assert len(a["world_matrix"][0]) == len(sky_preset())
    for stats in a["summary"].values():
        assert stats["q05"] <= stats["q50"] <= stats["q95"]
    try:
        DimensionSpec("ballistic_coefficient", 0.0, 1.0).validate()
    except ValueError:
        pass
    else:
        raise AssertionError("domain guard did not reject a forbidden parameter")
    assert a["epistemic_status"]["simulation_is_evidence"] is False
    print("JANUS_WORLD_MATRIX_SELF_TEST=PASS")


def main() -> int:
    parser = argparse.ArgumentParser(description="JANUS Cosmos neutral uncertainty world-matrix engine")
    parser.add_argument("--preset", choices=("sky", "planet"), default="sky")
    parser.add_argument("--samples", type=int, default=256)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--namespace", default="JANUS_COSMOS")
    parser.add_argument("--output", default=None)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        self_test()
        return 0

    specs = sky_preset() if args.preset == "sky" else planet_preset()
    payload = JanusWorldMatrix(specs, samples=args.samples, seed=args.seed, namespace=args.namespace).generate()
    text = json.dumps(payload, ensure_ascii=False, indent=2)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(text + "\n")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
