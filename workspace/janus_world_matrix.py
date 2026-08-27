#!/usr/bin/env python3
"""
JANUS COSMOS WORLD MATRIX
Respiciens et Prospiciens: reproducible uncertainty expansion for astronomy,
planetary science, remote sensing and scientific sensor models.

The useful Janus idea retained here is a reproducible fan of plausible
physical and observational worlds with an explicit epistemic firewall.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import statistics
from dataclasses import asdict, dataclass
from typing import Dict, List, Sequence

NORMAL = statistics.NormalDist()


def _canonical_json(obj: object) -> bytes:
    return json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _sha256_obj(obj: object) -> str:
    return hashlib.sha256(_canonical_json(obj)).hexdigest()


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
        if not math.isfinite(self.mean):
            raise ValueError(f"{self.name}: mean must be finite")
        if not math.isfinite(self.sigma) or self.sigma < 0:
            raise ValueError(f"{self.name}: sigma must be finite and >= 0")
        if self.low is not None and self.high is not None and self.low > self.high:
            raise ValueError(f"{self.name}: low must be <= high")
        if self.distribution not in {"normal", "uniform"}:
            raise ValueError(f"{self.name}: unsupported distribution {self.distribution!r}")

    def transform_unit(self, u: float) -> float:
        self.validate()
        u = min(1.0 - 1e-12, max(1e-12, float(u)))
        if self.sigma == 0:
            value = self.mean
        elif self.distribution == "uniform":
            value = self.mean + (2.0 * u - 1.0) * self.sigma * math.sqrt(3.0)
        else:
            value = self.mean + self.sigma * NORMAL.inv_cdf(u)
        if self.low is not None:
            value = max(self.low, value)
        if self.high is not None:
            value = min(self.high, value)
        return float(value)

    def draw(self, rng: random.Random) -> float:
        if self.distribution == "uniform":
            return self.transform_unit(rng.random())
        if self.sigma == 0:
            return self.transform_unit(0.5)
        value = rng.gauss(self.mean, self.sigma)
        if self.low is not None:
            value = max(self.low, value)
        if self.high is not None:
            value = min(self.high, value)
        return float(value)


def deterministic_seed(specs: Sequence[DimensionSpec], samples: int, namespace: str, sampling: str) -> int:
    payload = {
        "namespace": str(namespace),
        "samples": int(samples),
        "sampling": str(sampling),
        "dimensions": [asdict(s) for s in specs],
    }
    return int.from_bytes(hashlib.sha256(_canonical_json(payload)).digest()[:8], "big")


def _latin_hypercube_units(samples: int, dims: int, rng: random.Random) -> List[List[float]]:
    columns: List[List[float]] = []
    for _ in range(dims):
        bins = list(range(samples))
        rng.shuffle(bins)
        columns.append([(bins[i] + rng.random()) / samples for i in range(samples)])
    return [[columns[j][i] for j in range(dims)] for i in range(samples)]


def _covariance_matrix(rows: Sequence[Sequence[float]]) -> List[List[float]]:
    if not rows:
        return []
    n_dims = len(rows[0])
    cols = [[float(row[j]) for row in rows] for j in range(n_dims)]
    means = [statistics.fmean(c) for c in cols]
    denom = max(1, len(rows) - 1)
    return [
        [
            sum((a - means[i]) * (b - means[j]) for a, b in zip(cols[i], cols[j])) / denom
            for j in range(n_dims)
        ]
        for i in range(n_dims)
    ]


class JanusWorldMatrix:
    """Generate and summarize a reproducible ensemble of plausible observational worlds."""

    schema = "janus.cosmos.world-matrix.v1.1"
    formula = "RESPICIENS_ET_PROSPICIENS"

    def __init__(
        self,
        dimensions: Sequence[DimensionSpec],
        samples: int = 1000,
        seed: int | None = None,
        namespace: str = "JANUS_COSMOS",
        sampling: str = "latin_hypercube",
    ) -> None:
        if not dimensions:
            raise ValueError("at least one dimension is required")
        if samples < 2:
            raise ValueError("samples must be >= 2")
        if sampling not in {"latin_hypercube", "monte_carlo"}:
            raise ValueError("sampling must be latin_hypercube or monte_carlo")
        self.dimensions = tuple(dimensions)
        for spec in self.dimensions:
            spec.validate()
        self.samples = int(samples)
        self.namespace = str(namespace)
        self.sampling = sampling
        self.seed = int(seed) if seed is not None else deterministic_seed(self.dimensions, self.samples, self.namespace, self.sampling)

    def _matrix(self, rng: random.Random) -> List[List[float]]:
        if self.sampling == "latin_hypercube":
            units = _latin_hypercube_units(self.samples, len(self.dimensions), rng)
            return [[spec.transform_unit(units[i][j]) for j, spec in enumerate(self.dimensions)] for i in range(self.samples)]
        return [[spec.draw(rng) for spec in self.dimensions] for _ in range(self.samples)]

    def generate(self) -> Dict[str, object]:
        rng = random.Random(self.seed)
        matrix = self._matrix(rng)
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

        frozen_priors = [asdict(s) for s in self.dimensions]
        prior_contract = {
            "namespace": self.namespace,
            "sampling": self.sampling,
            "samples": self.samples,
            "seed": self.seed,
            "dimensions": frozen_priors,
        }
        payload: Dict[str, object] = {
            "schema": self.schema,
            "formula": self.formula,
            "namespace": self.namespace,
            "sampling": self.sampling,
            "seed": self.seed,
            "samples": self.samples,
            "dimensions": names,
            "frozen_priors": frozen_priors,
            "input_sha256": _sha256_obj(prior_contract),
            "world_matrix": matrix,
            "summary": summary,
            "covariance": _covariance_matrix(matrix),
            "learned_residual": {
                "applied": False,
                "may_replace_physics_baseline": False,
                "sidecar_only_if_validated": True,
                "training_data_hash_required_if_used": True,
            },
            "epistemic_status": {
                "simulation_is_evidence": False,
                "simulation_is_discovery": False,
                "requires_observed_data_comparison": True,
                "requires_null_or_control": True,
                "priors_should_be_frozen_before_result_inspection": True,
                "raw_inputs_and_parameters_are_immutable": True,
                "failed_predictions_remain_searchable": True,
                "simulation_to_reality_gap_must_be_recorded": True,
                "idk_is_valid_output": True,
                "uncalibrated_confidence_forbidden": True,
            },
            "known_debts": {
                "input_cross_dimension_correlation_model": "NOT_IMPLEMENTED_V1_1",
                "posterior_update": "MUST_USE_SEPARATE_VALIDATION_OR_HOLDOUT_GATE",
            },
            "domain_boundary": {
                "astronomy_and_planetary_search": True,
                "sensor_uncertainty": True,
                "environment_uncertainty": True,
                "orbital_dynamics": "SEPARATE_SAFE_MODULE",
                "scope": "ASTRONOMY_PLANETARY_REMOTE_SENSING_AND_ENVIRONMENTAL_UNCERTAINTY",
            },
            "claim_ceiling": "MODEL_CONDITIONED_UNCERTAINTY_ENSEMBLE_NOT_OBSERVATION",
        }
        payload["result_sha256"] = _sha256_obj(payload)
        return payload


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
    a = JanusWorldMatrix(sky_preset(), samples=64, seed=424242).generate()
    b = JanusWorldMatrix(sky_preset(), samples=64, seed=424242).generate()
    assert a["world_matrix"] == b["world_matrix"]
    assert a["result_sha256"] == b["result_sha256"]
    assert a["sampling"] == "latin_hypercube"
    assert len(a["world_matrix"]) == 64
    assert len(a["world_matrix"][0]) == len(sky_preset())
    for stats in a["summary"].values():
        assert stats["q05"] <= stats["q50"] <= stats["q95"]
    assert a["epistemic_status"]["simulation_is_evidence"] is False
    assert a["epistemic_status"]["idk_is_valid_output"] is True
    print("JANUS_WORLD_MATRIX_SELF_TEST=PASS")


def main() -> int:
    parser = argparse.ArgumentParser(description="JANUS Cosmos reproducible uncertainty world-matrix engine")
    parser.add_argument("--preset", choices=("sky", "planet"), default="sky")
    parser.add_argument("--samples", type=int, default=256)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--namespace", default="JANUS_COSMOS")
    parser.add_argument("--sampling", choices=("latin_hypercube", "monte_carlo"), default="latin_hypercube")
    parser.add_argument("--output", default=None)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        self_test()
        return 0

    specs = sky_preset() if args.preset == "sky" else planet_preset()
    payload = JanusWorldMatrix(
        specs, samples=args.samples, seed=args.seed, namespace=args.namespace, sampling=args.sampling
    ).generate()
    text = json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(text + "\n")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
