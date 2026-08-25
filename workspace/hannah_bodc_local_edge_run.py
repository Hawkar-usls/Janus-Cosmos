#!/usr/bin/env python3
"""Run frozen Hannah edge mnemonics on already-acquired local BODC bytes."""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from cousteau_synesthetic_memory_core import build_passport, compare_passports
from cousteau_synesthetic_semantic_overlay import enrich_passport
from hannah_bodc_realdata_probe import (
    EDGE_SECONDS,
    SCALE_BY_SECONDS,
    blake2_file,
    infer_depth_index,
    parse_em122_edges,
    parse_timestamp,
    parse_tpl,
    records_for_window,
    sha256_file,
    split_fields,
    window_payload,
)


def find_one(root: Path, name: str) -> Path | None:
    matches = [p for p in root.rglob('*') if p.is_file() and p.name.lower() == name.lower()]
    if not matches:
        return None
    if len(matches) > 1:
        matches.sort(key=lambda p: (len(p.parts), str(p)))
    return matches[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', required=True, type=Path)
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args()

    em = find_one(args.root, 'em122.ACO')
    tplp = find_one(args.root, 'em122.TPL')
    result = {
        'schema': 'janus.cosmos.cousteau.hannah_bodc.real_edge_mnemonic_run.v1',
        'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'root': str(args.root),
        'scientific_claim': False,
        'raw_bytes_redistributed': False,
        'status': 'STARTED',
        'hard_rules': [
            'RAW_BYTES_OUTRANK_MNEMONIC',
            'SENSORY_MATCH_IS_RETRIEVAL_ONLY',
            'TIME_MATCH_REQUIRES_LATER_SPACE_PROOF',
        ],
    }
    if em is None:
        result['status'] = 'EM122_ACO_NOT_FOUND'
    else:
        result['em122_original'] = {
            'size_bytes': em.stat().st_size,
            'sha256': sha256_file(em),
            'blake2b_256': blake2_file(em),
            'relative_path': str(em.relative_to(args.root)),
        }
        sample_fields = []
        with em.open('rb') as f:
            for raw in f:
                fields = split_fields(raw)
                if parse_timestamp(fields) is not None:
                    sample_fields = fields
                    break
        tpl = parse_tpl(tplp)
        mapping = infer_depth_index(sample_fields, tpl)
        result['em122_schema'] = {
            'sample_column_count': len(sample_fields),
            'tpl': tpl,
            'depth_mapping': mapping,
            'raw_sample_emitted': False,
        }
        if not mapping.get('resolved'):
            result['status'] = 'RAW_HASHED_DEPTH_COLUMN_UNRESOLVED'
        else:
            edge = parse_em122_edges(em, int(mapping['index']))
            result['em122_parse_summary'] = edge['summary']
            result['edge_windows'] = {}
            result['head_tail_comparisons'] = {}
            for sec in EDGE_SECONDS:
                scale = SCALE_BY_SECONDS[sec]
                hrec = records_for_window(edge['head_records'], sec, 'HEAD')
                trec = records_for_window(edge['tail_records'], sec, 'TAIL')
                hp, hraw, hmeta = window_payload(hrec)
                tp, traw, tmeta = window_payload(trec)
                hpass = enrich_passport(build_passport(
                    hp,
                    direction='HEAD_FORWARD',
                    scale=scale,
                    provenance={'source': 'BODCREQ-9408/em122.ACO', 'window': 'HEAD', 'seconds': sec},
                    raw_bytes=hraw,
                ))
                tpass = enrich_passport(build_passport(
                    tp,
                    direction='TAIL_REVERSE',
                    scale=scale,
                    provenance={'source': 'BODCREQ-9408/em122.ACO', 'window': 'TAIL', 'seconds': sec},
                    raw_bytes=traw,
                ))
                result['edge_windows'][scale] = {
                    'HEAD': {'meta': hmeta, 'measurement_payload': hp, 'passport': hpass},
                    'TAIL': {'meta': tmeta, 'measurement_payload': tp, 'passport': tpass},
                }
                result['head_tail_comparisons'][scale] = compare_passports(hpass, tpass)
            result['status'] = 'REAL_EM122_EDGE_MNEMONICS_READY'

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps({
        'status': result.get('status'),
        'em122_sha256': (result.get('em122_original') or {}).get('sha256'),
        'parse_summary': result.get('em122_parse_summary'),
        'scientific_claim': False,
    }, indent=2))
    return 0 if result['status'] == 'REAL_EM122_EDGE_MNEMONICS_READY' else 2


if __name__ == '__main__':
    raise SystemExit(main())
