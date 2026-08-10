from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
from pathlib import Path

import numpy as np

import run_blind_corpus_morphology as morph
from repo_derived_anomaly_protocol import structural_features
from run_blind_corpus_variant import download_variant, base

DEFAULT_TARGETS = ("NGC1425", "NGC1637")
DEFAULT_SEEDS = (20260810, 20260811, 20260812, 20260813, 20260814)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="data/hst_expanded_live_manifest.json")
    ap.add_argument("--nulls", type=int, default=4096)
    ap.add_argument("--seeds", default=','.join(map(str, DEFAULT_SEEDS)))
    ap.add_argument("--output-dir", default="results/focus_recheck")
    args = ap.parse_args()

    seeds = tuple(int(x) for x in args.seeds.split(',') if x.strip())
    manifest = json.loads(Path(args.manifest).read_text(encoding='utf-8'))
    targets = [t for t in manifest['targets'] if t['target'] in DEFAULT_TARGETS]
    if len(targets) != len(DEFAULT_TARGETS):
        raise RuntimeError(f'Missing focus targets: expected {DEFAULT_TARGETS}')

    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    base.EVENT_LOG = out_dir / 'janus-cosmos-focus-events.jsonl'
    base.RECEIPT = out_dir / 'janus-cosmos-focus-receipt.json'
    base.EVENT_LOG.unlink(missing_ok=True)
    morph.NULLS = args.nulls

    rows = []
    with tempfile.TemporaryDirectory() as td:
        for seed in seeds:
            rng = np.random.default_rng(seed)
            for target in targets:
                for item in target['filters']:
                    path = Path(td) / f"{target['target']}_{item['filter']}_{seed}.fits"
                    source = download_variant(item['dataURI'], path, 'mast_api')
                    image = base.read_image(path)
                    result = morph.analyze_filter(image, rng, target['target'], item['filter'])
                    features = structural_features(image)
                    rows.append({
                        'target': target['target'],
                        'filter': item['filter'],
                        'seed': seed,
                        'legacy_p': result['legacy_pixel_permutation']['p_empirical'],
                        'phase_p': result['morphology_preserving_phase']['p_empirical'],
                        'block_p': result['local_block_shuffle']['p_empirical'],
                        'robust_candidate': result['robust_candidate'],
                        'observed_score': result['observed_score'],
                        'feature_fingerprint': features['feature_fingerprint'],
                        'structural_features': features,
                        'source_products': source,
                    })

    by_target = {}
    for row in rows:
        by_target.setdefault(row['target'], []).append(row)

    consensus = []
    for target, items in sorted(by_target.items()):
        filters = sorted({x['filter'] for x in items})
        passing_filters = []
        per_filter = {}
        for filt in filters:
            ff = [x for x in items if x['filter'] == filt]
            all_pass = all(x['robust_candidate'] for x in ff)
            passing_filters.append(filt) if all_pass else None
            per_filter[filt] = {
                'seed_count': len(ff),
                'all_seeds_robust': all_pass,
                'phase_p_min': min(x['phase_p'] for x in ff),
                'phase_p_max': max(x['phase_p'] for x in ff),
                'block_p_min': min(x['block_p'] for x in ff),
                'block_p_max': max(x['block_p'] for x in ff),
            }
        consensus.append({
            'target': target,
            'seed_count': len(items) // max(1, len(filters)),
            'passing_filters_all_seeds': passing_filters,
            'focus_cross_band_consensus': len(passing_filters) >= 2,
            'per_filter': per_filter,
        })

    receipt = {
        'schema': 'janus.cosmos.hst.focus_recheck_receipt.v0.1',
        'status': 'FOCUS_RECHECK_ONLY',
        'purpose': 'Secondary, preregistered recheck of prior survivors; does not alter the main corpus gate or candidate threshold.',
        'source': manifest['source'],
        'manifest_sha256': manifest['manifest_sha256'],
        'targets': list(DEFAULT_TARGETS),
        'nulls_per_null_model': args.nulls,
        'seeds': list(seeds),
        'results': rows,
        'consensus': consensus,
        'focus_consensus_candidate_count': sum(1 for x in consensus if x['focus_cross_band_consensus']),
        'blind_gate': {
            'semantic_analysis': False,
            'ocr': False,
            'face_search': False,
            'cipher_search': False,
            'post_hoc_tuning': False,
            'human_label_inference': False,
            'main_gate_unchanged': True,
        },
        'claim_ceiling': 'Secondary geometric recheck only; no astronomical discovery claim without independent replication and scientific review.',
    }
    base.RECEIPT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    receipt['event_log'] = {'sha256': hashlib.sha256(base.EVENT_LOG.read_bytes()).hexdigest()}
    receipt['receipt_sha256'] = hashlib.sha256(base.RECEIPT.read_bytes()).hexdigest()
    base.RECEIPT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps({
        'targets': list(DEFAULT_TARGETS),
        'nulls': args.nulls,
        'seeds': list(seeds),
        'focus_consensus_candidate_count': receipt['focus_consensus_candidate_count'],
        'receipt_sha256': receipt['receipt_sha256'],
    }, indent=2))


if __name__ == '__main__':
    main()
