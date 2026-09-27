#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--earth-ledger', type=Path, required=True)
    ap.add_argument('--palomar-ledger', type=Path, required=True)
    ap.add_argument('--prereg', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()

    earth = json.loads(args.earth_ledger.read_text(encoding='utf-8'))
    pal = json.loads(args.palomar_ledger.read_text(encoding='utf-8'))
    prereg = json.loads(args.prereg.read_text(encoding='utf-8'))
    day = {r['date']: r for r in pal.get('rows', [])}

    rows = []
    primary = []
    provisional = []
    for e in earth.get('rows', []):
        ds = e.get('date')
        observed = day.get(ds)
        if observed is None:
            outcome = 'PALOMAR_NO_OBSERVATION_ROW'
            palomar = None
        elif observed.get('observed_flag') is False:
            outcome = 'PALOMAR_NOT_OBSERVED'
            palomar = {'observed_flag': False}
        else:
            outcome = 'PALOMAR_OBSERVED_ON_EARTH_EVENT_DATE'
            palomar = {
                'observed_flag': True,
                'candidate_count': observed.get('candidate_count'),
                'plate_count': observed.get('plate_count'),
                'tile_count': observed.get('tile_count'),
                'exposure': observed.get('exposure'),
                'artifact_state': observed.get('artifact_state'),
                'plate_ids': observed.get('plate_ids'),
            }
        row = {
            'row_id': e.get('row_id'),
            'date': ds,
            'site': e.get('site'),
            'source_grade': e.get('source_grade'),
            'earth_bridge_eligible': e.get('bridge_eligible') is True,
            'outcome': outcome,
            'palomar': palomar,
        }
        rows.append(row)
        (primary if e.get('bridge_eligible') is True else provisional).append(row)

    def counts(xs):
        out = {}
        for r in xs:
            out[r['outcome']] = out.get(r['outcome'], 0) + 1
        return out

    primary_counts = counts(primary)
    exact_matches = [r for r in primary if r['outcome'] == 'PALOMAR_OBSERVED_ON_EARTH_EVENT_DATE']
    if exact_matches:
        verdict = 'PRIMARY_EXACT_DAY_OBSERVATION_OPPORTUNITY_OVERLAP_PRESENT'
    elif primary and all(r['outcome'] == 'PALOMAR_NO_OBSERVATION_ROW' for r in primary):
        verdict = 'NO_EXACT_DAY_PALOMAR_OBSERVATION_OPPORTUNITY_ON_FROZEN_PRIMARY_OAK_RIDGE_DATES'
    else:
        verdict = 'MIXED_OR_BLOCKED_PRIMARY_EXACT_DAY_RESULT'

    out = {
        'schema': 'janus.cosmos.bridge.oak_ridge_exact_day_result.v1',
        'artifact_id': 'BRIDGE-R2-OAK-RIDGE-EXACT-DAY-RESULT-2026-09-27-v1.0',
        'generated_at_utc': datetime.now(timezone.utc).isoformat(),
        'status': verdict,
        'execution_contract': {
            'prereg_artifact_id': prereg.get('artifact_id'),
            'prereg_status': prereg.get('status'),
            'primary_join': prereg.get('independence_contract', {}).get('primary_join'),
            'window_widening_used': False,
            'palomar_day_rows': len(day),
            'palomar_day_ledger_status': pal.get('status'),
        },
        'source_hashes': {
            'earth_ledger_sha256': sha256(args.earth_ledger),
            'palomar_day_ledger_sha256': sha256(args.palomar_ledger),
            'prereg_sha256': sha256(args.prereg),
        },
        'primary_cohort': {
            'n': len(primary),
            'outcome_counts': primary_counts,
            'rows': primary,
        },
        'provisional_excluded_from_primary': {
            'n': len(provisional),
            'outcome_counts': counts(provisional),
            'rows': provisional,
        },
        'all_frozen_dates': rows,
        'interpretation': {
            'if_no_observation': 'PALOMAR_NO_OBSERVATION_ROW means POSS-I provides no exact-date observation opportunity for that frozen Earth date. It is not an observed zero and is not evidence that no sky phenomenon existed.',
            'if_observed': 'An exact-date Palomar observation would establish only observation opportunity and permit candidate-level comparison; it would not establish common cause or physical coupling.',
            'secondary_windows': 'PLUS/MINUS windows remain NOT_EXECUTED in this result. They may be tested only after this exact-day result is frozen under a new preregistration/null contract.',
        },
        'claim_ceiling': 'THIS RESULT TESTS TEMPORAL OBSERVATION OPPORTUNITY ONLY. IT DOES NOT ESTABLISH OR REFUTE UAP IDENTITY, NONHUMAN ORIGIN, NUCLEAR ATTRACTION, OR EARTH-SPACE CAUSAL COUPLING.',
        'canonical_seal': 'THE EARTH DATES WERE FROZEN FIRST. IF PALOMAR WAS NOT LOOKING, JANUS RECORDS NO OBSERVATION OPPORTUNITY — NEVER A ZERO SKY.'
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'status': verdict, 'primary_n': len(primary), 'primary_counts': primary_counts, 'provisional_n': len(provisional)}, sort_keys=True))


if __name__ == '__main__':
    main()
