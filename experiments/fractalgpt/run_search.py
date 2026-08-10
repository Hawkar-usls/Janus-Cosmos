"""Janus Cosmos v0.1: blind multiscale morphology preflight.

This pilot intentionally operates on extracted feature rows. It does not
interpret images semantically and it never claims an astronomical discovery.
The spatial null shuffles x/y coordinates while preserving object/band labels
and signal values.
"""
from __future__ import annotations
import csv, hashlib, json, math, random, sys
from pathlib import Path

SCALES = (8, 16, 32)
ORIENTATIONS = tuple(range(0, 180, 30))
NULLS = 256
SEED = 20260810


def trajectory(x, y, scale, orientation):
    a = math.radians(orientation)
    u = x * math.cos(a) + y * math.sin(a)
    v = -x * math.sin(a) + y * math.cos(a)
    return abs(math.sin(u * scale * 0.17) * math.cos(v * scale * 0.11))


def score(rows):
    if not rows:
        return 0.0
    total = 0.0
    for r in rows:
        total += sum(trajectory(r['x'], r['y'], s, o) for s in SCALES for o in ORIENTATIONS)
    return total / (len(rows) * len(SCALES) * len(ORIENTATIONS))


def main(csv_path):
    rows = []
    with open(csv_path, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            rows.append({
                'object_id': r['object_id'], 'band': r['band'],
                'x': float(r['x']), 'y': float(r['y']),
                'signal': float(r.get('signal', 1.0)),
            })

    observed = score(rows)
    rng = random.Random(SEED)
    xs = [r['x'] for r in rows]
    ys = [r['y'] for r in rows]
    null_scores = []
    for _ in range(NULLS):
        nx, ny = xs[:], ys[:]
        rng.shuffle(nx); rng.shuffle(ny)
        shuffled = [dict(r, x=x, y=y) for r, x, y in zip(rows, nx, ny)]
        null_scores.append(score(shuffled))

    ge = sum(v >= observed for v in null_scores)
    p = (ge + 1) / (NULLS + 1)
    receipt = {
        'schema': 'janus.cosmos.fractalgpt.receipt.v0.1',
        'status': 'CANDIDATE_ONLY' if p < 0.05 else 'NO_CANDIDATE',
        'astronomical_discovery': False,
        'observed_score': observed,
        'null_min': min(null_scores) if null_scores else None,
        'null_median': sorted(null_scores)[len(null_scores)//2] if null_scores else None,
        'null_max': max(null_scores) if null_scores else None,
        'p_empirical': p,
        'scales': SCALES,
        'orientations': ORIENTATIONS,
        'nulls': NULLS,
        'seed': SEED,
        'semantic_analysis': False,
        'source_sha256': hashlib.sha256(Path(csv_path).read_bytes()).hexdigest(),
        'claim_ceiling': 'A planner score is not a discovery; source-level and independent replication are required.'
    }
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('usage: python experiments/fractalgpt/run_search.py <features.csv>')
    main(sys.argv[1])
