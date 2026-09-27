#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from datetime import date, datetime, timedelta, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def parse_date(s: str) -> date:
    return date.fromisoformat(s)


def fraction(num: int, den: int) -> float:
    return num / den if den else float('nan')


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--earth-ledger', type=Path, required=True)
    ap.add_argument('--palomar-ledger', type=Path, required=True)
    ap.add_argument('--exact-day-result', type=Path, required=True)
    ap.add_argument('--prereg', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()

    earth = json.loads(args.earth_ledger.read_text(encoding='utf-8'))
    pal = json.loads(args.palomar_ledger.read_text(encoding='utf-8'))
    r2 = json.loads(args.exact_day_result.read_text(encoding='utf-8'))
    prereg = json.loads(args.prereg.read_text(encoding='utf-8'))

    radii = [int(x) for x in prereg['window_contract']['radii_days']]
    assert radii == [1, 3], radii
    widest = max(radii)

    primary_rows = [r for r in earth.get('rows', []) if r.get('bridge_eligible') is True]
    provisional_rows = [r for r in earth.get('rows', []) if r.get('bridge_eligible') is not True]
    primary_rows = sorted(primary_rows, key=lambda r: r['date'])
    primary_dates = [parse_date(r['date']) for r in primary_rows]
    frozen_dates = [parse_date(x) for x in prereg['primary_cohort_contract']['frozen_primary_dates']]
    assert primary_dates == frozen_dates, (primary_dates, frozen_dates)
    assert len(primary_rows) == prereg['primary_cohort_contract']['expected_primary_n'] == 3

    study_start, study_end = [parse_date(x) for x in pal['study_window']]
    pal_rows_by_date = {
        parse_date(r['date']): r
        for r in pal.get('rows', [])
        if r.get('observed_flag') is True
    }
    pal_dates = set(pal_rows_by_date)

    def evaluate(shift_days: int, radius: int, include_details: bool = False):
        union_dates = set()
        event_details = []
        hit_count = 0
        for src in primary_rows:
            original = parse_date(src['date'])
            center = original + timedelta(days=shift_days)
            ws = center - timedelta(days=radius)
            we = center + timedelta(days=radius)
            hits = sorted(d for d in pal_dates if ws <= d <= we)
            union_dates.update(hits)
            hit = bool(hits)
            hit_count += int(hit)
            if include_details:
                event_details.append({
                    'row_id': src.get('row_id'),
                    'original_earth_date': original.isoformat(),
                    'translated_center_date': center.isoformat(),
                    'window_start': ws.isoformat(),
                    'window_end': we.isoformat(),
                    'observation_opportunity_dates': [d.isoformat() for d in hits],
                    'observation_opportunity_count': len(hits),
                    'window_hit': hit,
                })
        out = {
            'events_with_observation_opportunity': hit_count,
            'unique_observation_opportunity_days': len(union_dates),
        }
        if include_details:
            out['event_rows'] = event_details
            out['union_observation_opportunity_dates'] = [d.isoformat() for d in sorted(union_dates)]
        return out

    # One common null support for both radii: every translated +/-3 window must remain
    # fully inside the frozen Palomar study window. This makes the +/-1 and +/-3
    # comparisons use the same rigidly translated Earth-date patterns.
    min_shift = (study_start + timedelta(days=widest) - min(primary_dates)).days
    max_shift = (study_end - timedelta(days=widest) - max(primary_dates)).days
    shifts = list(range(min_shift, max_shift + 1))
    assert 0 in shifts

    actual = {}
    null = {}
    all_zero = True
    any_overlap = False

    for radius in radii:
        key = f'plus_minus_{radius}_day'
        act = evaluate(0, radius, include_details=True)
        actual[key] = act
        if act['events_with_observation_opportunity'] != 0:
            all_zero = False
            any_overlap = True

        hit_hist = Counter()
        union_hist = Counter()
        max_hits = -1
        max_union = -1
        shifts_at_max_hits = []
        shifts_at_max_union = []

        for s in shifts:
            stat = evaluate(s, radius, include_details=False)
            eh = stat['events_with_observation_opportunity']
            ud = stat['unique_observation_opportunity_days']
            hit_hist[eh] += 1
            union_hist[ud] += 1
            if eh > max_hits:
                max_hits = eh
                shifts_at_max_hits = [s]
            elif eh == max_hits:
                shifts_at_max_hits.append(s)
            if ud > max_union:
                max_union = ud
                shifts_at_max_union = [s]
            elif ud == max_union:
                shifts_at_max_union.append(s)

        n = len(shifts)
        actual_hits = act['events_with_observation_opportunity']
        actual_union = act['unique_observation_opportunity_days']
        ge_hits = sum(v for k, v in hit_hist.items() if k >= actual_hits)
        le_hits = sum(v for k, v in hit_hist.items() if k <= actual_hits)
        ge_union = sum(v for k, v in union_hist.items() if k >= actual_union)
        le_union = sum(v for k, v in union_hist.items() if k <= actual_union)

        null[key] = {
            'feasible_rigid_translations': n,
            'event_hit_histogram': {str(k): hit_hist[k] for k in sorted(hit_hist)},
            'unique_observation_day_histogram': {str(k): union_hist[k] for k in sorted(union_hist)},
            'actual_event_hits': actual_hits,
            'actual_unique_observation_days': actual_union,
            'fraction_translations_event_hits_ge_actual': fraction(ge_hits, n),
            'fraction_translations_event_hits_le_actual': fraction(le_hits, n),
            'fraction_translations_with_zero_event_hits': fraction(hit_hist.get(0, 0), n),
            'fraction_translations_unique_days_ge_actual': fraction(ge_union, n),
            'fraction_translations_unique_days_le_actual': fraction(le_union, n),
            'max_event_hits_in_null': max_hits,
            'translations_at_max_event_hits_n': len(shifts_at_max_hits),
            'example_shifts_at_max_event_hits': shifts_at_max_hits[:20],
            'max_unique_observation_days_in_null': max_union,
            'translations_at_max_unique_days_n': len(shifts_at_max_union),
            'example_shifts_at_max_unique_days': shifts_at_max_union[:20],
        }

    provisional_actual = {}
    for radius in radii:
        key = f'plus_minus_{radius}_day'
        rows = []
        for src in provisional_rows:
            center = parse_date(src['date'])
            ws = center - timedelta(days=radius)
            we = center + timedelta(days=radius)
            hits = sorted(d for d in pal_dates if ws <= d <= we)
            rows.append({
                'row_id': src.get('row_id'),
                'date': src.get('date'),
                'source_grade': src.get('source_grade'),
                'window_start': ws.isoformat(),
                'window_end': we.isoformat(),
                'observation_opportunity_dates': [d.isoformat() for d in hits],
                'window_hit': bool(hits),
            })
        provisional_actual[key] = rows

    if all_zero:
        status = 'NO_PALOMAR_OBSERVATION_OPPORTUNITY_WITHIN_PREREGISTERED_PLUS_MINUS_1_OR_3_DAY_PRIMARY_WINDOWS'
    elif any_overlap:
        status = 'PALOMAR_OBSERVATION_OPPORTUNITY_OVERLAP_PRESENT_IN_AT_LEAST_ONE_PREREGISTERED_WINDOW'
    else:
        status = 'MIXED_OR_BLOCKED_WINDOW_NULL_RESULT'

    pairwise_gaps = []
    for i in range(len(primary_dates)):
        for j in range(i + 1, len(primary_dates)):
            pairwise_gaps.append((primary_dates[j] - primary_dates[i]).days)

    out = {
        'schema': 'janus.cosmos.bridge.oak_ridge_window_null_result.v1',
        'artifact_id': 'BRIDGE-R3-OAK-RIDGE-WINDOW-NULL-RESULT-2026-09-27-v1.0',
        'generated_at_utc': datetime.now(timezone.utc).isoformat(),
        'status': status,
        'execution_contract': {
            'prereg_artifact_id': prereg.get('artifact_id'),
            'prereg_status': prereg.get('status'),
            'exact_day_result_artifact_id': r2.get('artifact_id'),
            'exact_day_status': r2.get('status'),
            'primary_n': len(primary_rows),
            'provisional_n': len(provisional_rows),
            'window_radii_days': radii,
            'palomar_study_window': [study_start.isoformat(), study_end.isoformat()],
            'palomar_observation_day_rows': len(pal_rows_by_date),
            'null_name': prereg['calendar_null_contract']['name'],
            'common_null_support_uses_widest_radius_days': widest,
            'rigid_shift_min_days': min_shift,
            'rigid_shift_max_days': max_shift,
            'feasible_rigid_translations': len(shifts),
            'actual_shift_days': 0,
            'pairwise_primary_date_gaps_days': pairwise_gaps,
        },
        'source_hashes': {
            'earth_ledger_sha256': sha256(args.earth_ledger),
            'palomar_day_ledger_sha256': sha256(args.palomar_ledger),
            'exact_day_result_sha256': sha256(args.exact_day_result),
            'prereg_sha256': sha256(args.prereg),
        },
        'primary_actual_windows': actual,
        'calendar_null': null,
        'provisional_descriptive_only': {
            'excluded_from_primary_null': True,
            'rows_by_window': provisional_actual,
        },
        'interpretation': {
            'observation_opportunity_only': 'A Palomar day inside a window means the POSS-I cohort supplied an observation opportunity near that Earth date. It does not establish that any candidate is related to the Earth event.',
            'zero_window': 'No Palomar observation-day row inside a preregistered window is window-negative for observation opportunity only. It is not an observed zero sky.',
            'null': 'The rigid-translation null moves the entire three-date Earth pattern together and exhaustively scores every feasible translation; pairwise Earth-date separations are therefore unchanged.',
            'candidate_counts': 'Candidate counts are intentionally not the primary R3 statistic because R3 asks whether Palomar was observing near the Earth dates, not whether detector candidates prove a common phenomenon.'
        },
        'claim_ceiling': 'THIS RESULT TESTS TEMPORAL OBSERVATION-OPPORTUNITY COINCIDENCE ONLY. IT DOES NOT ESTABLISH OR REFUTE UAP IDENTITY, NONHUMAN ORIGIN, NUCLEAR ATTRACTION, OR EARTH-SPACE CAUSAL COUPLING.',
        'canonical_seal': 'THE WHOLE EARTH DATE PATTERN MOVES; THE WINDOWS DO NOT. NO OBSERVATION OPPORTUNITY IS NEVER AN OBSERVED ZERO.'
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({
        'status': status,
        'primary_n': len(primary_rows),
        'feasible_rigid_translations': len(shifts),
        'actual': {
            k: {
                'events_with_observation_opportunity': v['events_with_observation_opportunity'],
                'unique_observation_opportunity_days': v['unique_observation_opportunity_days'],
            }
            for k, v in actual.items()
        },
    }, sort_keys=True))


if __name__ == '__main__':
    main()
