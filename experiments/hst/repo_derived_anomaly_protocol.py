from __future__ import annotations

"""Repo-derived research mechanics for Janus Cosmos.

Sources: JANUS Fundamentum (attack/survivor discipline), JANUS Lapis
(gates/rejections/machine summaries), HRain (graph representation),
JANUS distributed swarm (observer-first invariants), and Janus Demiurge
(experiment-memory ideas, without adaptive tuning of the blind gate).
"""

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy import ndimage


@dataclass(frozen=True)
class InvariantSet:
    ocr: bool = False
    face_search: bool = False
    semantic_analysis: bool = False
    cipher_search: bool = False
    post_hoc_tuning: bool = False
    human_label_inference: bool = False


@dataclass(frozen=True)
class ObservationNode:
    target: str
    filter_name: str
    seed: int
    observed_score: float
    feature_fingerprint: str


@dataclass(frozen=True)
class DecisionRecord:
    target: str
    filter_name: str
    seed: int
    hypothesis: str
    gate_name: str
    passed: bool
    reason: str
    evidence_sha256: str


def normalize_unit(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float32)
    lo, hi = np.percentile(x, [1.0, 99.0])
    if hi <= lo:
        return np.zeros_like(x, dtype=np.float32)
    return np.clip((x - lo) / (hi - lo), 0.0, 1.0)


def structural_features(image: np.ndarray) -> dict:
    """Blind geometry descriptors; no semantic/object labels are used."""
    x = normalize_unit(image)
    component_summary = []
    for q in (0.70, 0.80, 0.90, 0.95):
        labels, count = ndimage.label(x >= q)
        sizes = np.bincount(labels.ravel())[1:] if count else np.array([], dtype=np.int64)
        component_summary.append({
            "threshold": q,
            "count": int(count),
            "largest_fraction": float(sizes.max() / x.size) if sizes.size else 0.0,
            "area_entropy": _entropy(sizes.astype(np.float64)) if sizes.size else 0.0,
        })

    rotational = {}
    for angle in (45.0, 90.0, 135.0, 180.0):
        r = ndimage.rotate(x, angle, reshape=False, order=1, mode="reflect")
        rotational[str(int(angle))] = _corrcoef(x, r)

    F = np.fft.fftshift(np.fft.fft2(x))
    power = np.abs(F) ** 2
    yy, xx = np.indices(power.shape)
    cy, cx = (np.array(power.shape) - 1) / 2.0
    radius = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)
    bins = np.linspace(0.0, float(radius.max()), 24)
    radial = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (radius >= lo) & (radius < hi)
        radial.append(float(power[m].mean()) if np.any(m) else 0.0)
    radial = np.asarray(radial, dtype=np.float64)
    radial_prob = radial / max(float(radial.sum()), 1e-12)

    smooth = ndimage.gaussian_filter(x, sigma=8.0, mode="reflect")
    residual = x - smooth
    out = {
        "component_summary": component_summary,
        "rotational_agreement": rotational,
        "fourier_radial_power": radial.tolist(),
        "fourier_radial_entropy": _entropy(radial_prob),
        "smooth_mean": float(smooth.mean()),
        "high_frequency_energy_fraction": float(np.mean(residual ** 2) / max(np.mean(x ** 2), 1e-12)),
    }
    out["feature_fingerprint"] = feature_fingerprint(out)
    return out


def _corrcoef(a: np.ndarray, b: np.ndarray) -> float:
    aa = a.ravel().astype(np.float64); bb = b.ravel().astype(np.float64)
    aa -= aa.mean(); bb -= bb.mean()
    den = float(np.sqrt(np.dot(aa, aa) * np.dot(bb, bb)))
    return float(np.dot(aa, bb) / den) if den else 0.0


def _entropy(values: np.ndarray) -> float:
    p = np.asarray(values, dtype=np.float64)
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    p = p / p.sum()
    return float(-np.sum(p * np.log2(p)))


def feature_fingerprint(features: dict) -> str:
    payload = json.dumps(features, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.sha256(payload).hexdigest()


def make_observation_node(target: str, filter_name: str, seed: int, observed_score: float, features: dict) -> dict:
    return asdict(ObservationNode(target, filter_name, int(seed), float(observed_score), features["feature_fingerprint"]))


def decide_gate(target: str, filter_name: str, seed: int, hypothesis: str, gate_name: str, passed: bool, reason: str, evidence: dict) -> dict:
    blob = json.dumps(evidence, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    return asdict(DecisionRecord(target, filter_name, int(seed), hypothesis, gate_name, bool(passed), reason, hashlib.sha256(blob).hexdigest()))


def build_candidate_graph(observations: Iterable[dict], decisions: Iterable[dict]) -> dict:
    nodes = list(observations)
    edges = []
    by_key = {}
    for n in nodes:
        by_key.setdefault((n["target"], n["filter_name"]), []).append(n)
    for group in by_key.values():
        for a, b in zip(group[:-1], group[1:]):
            edges.append({"type": "replicated_seed", "source": a["feature_fingerprint"], "target": b["feature_fingerprint"]})
    for d in decisions:
        edges.append({"type": "gate_decision", "target": f"{d['target']}:{d['filter_name']}", "gate": d["gate_name"], "passed": d["passed"]})
    return {"nodes": nodes, "edges": edges}


def write_protocol_summary(path: Path, invariant_set: InvariantSet, observations: list[dict], decisions: list[dict], graph: dict, run_metadata: dict) -> str:
    payload = {
        "schema": "janus.cosmos.repo_derived_protocol.v0.1",
        "invariants": asdict(invariant_set),
        "run_metadata": run_metadata,
        "observation_count": len(observations),
        "decision_count": len(decisions),
        "graph": graph,
        "negative_results": [d for d in decisions if not d["passed"]],
    }
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    path.write_text(text, encoding="utf-8")
    return hashlib.sha256(text.encode()).hexdigest()
