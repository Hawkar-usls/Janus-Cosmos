#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Callable


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


def add_years(d: date, years: int) -> date:
    return d.replace(year=d.year + years)


def convolve_histograms(histograms: list[Counter]) -> Counter:
    out = Counter({0: 1})
    for hist in histograms:
        nxt = Counter()
        for a, ca in out.items():
            for b, cb in hist.items():
                nxt[a + b] += ca * cb
        out = nxt
    return out


def sorted_hist(c: Counter) -> dict[str, int]:
    return {str(k): c[k] for k in sorted(c)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--earth-ledger', type=Path, required=True)
    ap.add_argument('--palomar-ledger', type=Path, required=True)
    ap.add_argument('--prereg', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()

    earth = json.loads(args.earth_ledger.read_text(encoding='utf-8'))
    pal = json.loads(args.palomar_ledger.read_text(encoding='utf-8'))
    prereg = json.loads(args.prereg.read_text(encoding='utf-8'))

    assert earth['status'] == 'FROZEN_BEFORE_R4_1_PRIMARY_PALOMAR_MEMBERSHIP_LOOKUP'
    assert prereg['status'] == 'FROZEN_BEFORE_R4_1_PRIMARY_PALOMAR_MEMBERSHIP_EVALUATION'
    assert prereg['earth_cohort_artifact'] == earth['artifact_id']

    frozen = prereg['frozen_primary_holdout']['sites']
    site_ids = list(frozen.keys())
    assert site_ids == ['HOLLOMAN_WHITE_SANDS', 'WALKER_AFB_509TH'], site_ids
    assert prereg['frozen_primary_holdout']['site_n'] == 2
    assert prereg['frozen_primary_holdout']['event_date_n'] == 4

    site_rows: dict[str, list[dict]] = {}
    for site in earth['primary_holdout_sites']:
        sid = site['site_id']
        assert sid in frozen, sid
        rows = [
            r for r in site['rows']
            if r.get('cohort_role') == 'HOLDOUT_PRIMARY_NEW_R4_1'
            and r.get('bridge_eligible') is True
        ]
        rows = sorted(rows, key=lambda r: r['date'])
        assert [r['date'] for r in rows] == frozen[sid], (sid, rows, frozen[sid])
        site_rows[sid] = rows

    assert set(site_rows) == set(site_ids)
    assert sum(len(v) for v in site_rows.values()) == 4

    holdback = earth['sealed_future_holdback']
    assert holdback['site_id'] == 'HANFORD_AEC_PLANT'
    assert holdback['date'] == '1950-07-30'
    assert holdback['bridge_eligible'] is False
    assert holdback['palomar_membership_lookup_forbidden'] is True
    assert 'HANFORD_AEC_PLANT' not in site_rows
    assert '1950-07-30' not in [r['date'] for rows in site_rows.values() for r in rows]

    radii = [int(x) for x in prereg['frozen_windows']['radii_days']]
    assert radii == [0, 1, 3], radii
    widest = max(radii)

    study_start, study_end = [parse_date(x) for x in pal['study_window']]
    assert [study_start.isoformat(), study_end.isoformat()] == prereg['palomar_study_window']
    pal_rows_by_date = {
        parse_date(r['date']): r
        for r in pal.get('rows', [])
        if r.get('observed_flag') is True
    }
    pal_dates = set(pal_rows_by_date)

    def event_window_hits(center: date, radius: int) -> list[date]:
        ws = center - timedelta(days=radius)
        we = center + timedelta(days=radius)
        return sorted(d for d in pal_dates if ws <= d <= we)

    def evaluate_site(
        sid: str,
        radius: int,
        transform: Callable[[date], date],
        include_details: bool = False,
    ) -> dict:
        event_hits = 0
        union_dates: set[date] = set()
        details = []
        for src in site_rows[sid]:
            original = parse_date(src['date'])
            center = transform(original)
            hits = event_window_hits(center, radius)
            hit = bool(hits)
            event_hits += int(hit)
            union_dates.update(hits)
            if include_details:
                details.append({
                    'row_id': src['row_id'],
                    'site_id': sid,
                    'original_earth_date': original.isoformat(),
                    'translated_center_date': center.isoformat(),
                    'window_start': (center - timedelta(days=radius)).isoformat(),
                    'window_end': (center + timedelta(days=radius)).isoformat(),
                    'observation_opportunity_dates': [d.isoformat() for d in hits],
                    'observation_opportunity_count': len(hits),
                    'window_hit': hit,
                })
        out = {
            'event_hits': event_hits,
            'site_has_hit': event_hits > 0,
            'unique_observation_days': len(union_dates),
        }
        if include_details:
            out['event_rows'] = details
            out['union_observation_opportunity_dates'] = [d.isoformat() for d in sorted(union_dates)]
        return out

    actual: dict[str, dict] = {}
    for radius in radii:
        key = f'radius_{radius}_days'
        sites = {}
        T = 0
        S = 0
        union: set[str] = set()
        for sid in site_ids:
            ev = evaluate_site(sid, radius, lambda d: d, include_details=True)
            sites[sid] = ev
            T += ev['event_hits']
            S += int(ev['site_has_hit'])
            union.update(ev['union_observation_opportunity_dates'])
        actual[key] = {
            'T_event_date_hits': T,
            'S_sites_with_at_least_one_hit': S,
            'U_unique_observation_days': len(union),
            'union_observation_opportunity_dates': sorted(union),
            'sites': sites,
        }

    # Primary null: independently rigid-shift each site's complete date cluster.
    primary_support: dict[str, list[int]] = {}
    for sid in site_ids:
        dates = [parse_date(r['date']) for r in site_rows[sid]]
        min_shift = (study_start + timedelta(days=widest) - min(dates)).days
        max_shift = (study_end - timedelta(days=widest) - max(dates)).days
        shifts = list(range(min_shift, max_shift + 1))
        assert shifts and 0 in shifts
        primary_support[sid] = shifts

    primary_null: dict[str, dict] = {}
    for radius in radii:
        key = f'radius_{radius}_days'
        event_hists: list[Counter] = []
        site_hists: list[Counter] = []
        by_site = {}
        for sid in site_ids:
            eh = Counter()
            sh = Counter()
            for shift in primary_support[sid]:
                ev = evaluate_site(sid, radius, lambda d, s=shift: d + timedelta(days=s))
                eh[ev['event_hits']] += 1
                sh[int(ev['site_has_hit'])] += 1
            event_hists.append(eh)
            site_hists.append(sh)
            by_site[sid] = {
                'shift_min_days': min(primary_support[sid]),
                'shift_max_days': max(primary_support[sid]),
                'admissible_shift_count': len(primary_support[sid]),
                'event_hit_histogram': sorted_hist(eh),
                'site_has_hit_histogram': sorted_hist(sh),
            }

        T_hist = convolve_histograms(event_hists)
        S_hist = convolve_histograms(site_hists)
        n_vectors = 1
        for sid in site_ids:
            n_vectors *= len(primary_support[sid])
        assert sum(T_hist.values()) == n_vectors
        assert sum(S_hist.values()) == n_vectors

        actual_T = actual[key]['T_event_date_hits']
        actual_S = actual[key]['S_sites_with_at_least_one_hit']
        ge_T = sum(v for k, v in T_hist.items() if k >= actual_T)
        le_T = sum(v for k, v in T_hist.items() if k <= actual_T)
        ge_S = sum(v for k, v in S_hist.items() if k >= actual_S)
        le_S = sum(v for k, v in S_hist.items() if k <= actual_S)
        primary_null[key] = {
            'null_name': prereg['primary_null']['name'],
            'by_site': by_site,
            'aggregate_shift_vector_count': n_vectors,
            'T_event_date_hit_histogram': sorted_hist(T_hist),
            'S_site_hit_histogram': sorted_hist(S_hist),
            'actual_T': actual_T,
            'actual_S': actual_S,
            'fraction_shift_vectors_T_ge_actual': fraction(ge_T, n_vectors),
            'fraction_shift_vectors_T_le_actual': fraction(le_T, n_vectors),
            'fraction_shift_vectors_S_ge_actual': fraction(ge_S, n_vectors),
            'fraction_shift_vectors_S_le_actual': fraction(le_S, n_vectors),
            'fraction_shift_vectors_T_zero': fraction(T_hist.get(0, 0), n_vectors),
            'max_T_in_null': max(T_hist),
            'max_S_in_null': max(S_hist),
        }

    # Seasonality sensitivity: independently shift each site's cluster by whole years.
    year_support: dict[str, list[int]] = {}
    for sid in site_ids:
        dates = [parse_date(r['date']) for r in site_rows[sid]]
        offsets = []
        for y in range(-100, 101):
            try:
                moved = [add_years(d, y) for d in dates]
            except ValueError:
                continue
            if min(moved) - timedelta(days=widest) < study_start:
                continue
            if max(moved) + timedelta(days=widest) > study_end:
                continue
            offsets.append(y)
        assert offsets and 0 in offsets
        year_support[sid] = offsets

    seasonality_null: dict[str, dict] = {}
    for radius in radii:
        key = f'radius_{radius}_days'
        event_hists: list[Counter] = []
        site_hists: list[Counter] = []
        by_site = {}
        for sid in site_ids:
            eh = Counter()
            sh = Counter()
            for years in year_support[sid]:
                ev = evaluate_site(sid, radius, lambda d, y=years: add_years(d, y))
                eh[ev['event_hits']] += 1
                sh[int(ev['site_has_hit'])] += 1
            event_hists.append(eh)
            site_hists.append(sh)
            by_site[sid] = {
                'admissible_year_offsets': year_support[sid],
                'admissible_year_offset_count': len(year_support[sid]),
                'event_hit_histogram': sorted_hist(eh),
                'site_has_hit_histogram': sorted_hist(sh),
            }

        T_hist = convolve_histograms(event_hists)
        S_hist = convolve_histograms(site_hists)
        n_vectors = 1
        for sid in site_ids:
            n_vectors *= len(year_support[sid])
        assert sum(T_hist.values()) == n_vectors
        assert sum(S_hist.values()) == n_vectors

        actual_T = actual[key]['T_event_date_hits']
        actual_S = actual[key]['S_sites_with_at_least_one_hit']
        ge_T = sum(v for k, v in T_hist.items() if k >= actual_T)
        le_T = sum(v for k, v in T_hist.items() if k <= actual_T)
        ge_S = sum(v for k, v in S_hist.items() if k >= actual_S)
        le_S = sum(v for k, v in S_hist.items() if k <= actual_S)
        seasonality_null[key] = {
            'null_name': prereg['seasonality_sensitivity_null']['name'],
            'by_site': by_site,
            'aggregate_year_shift_vector_count': n_vectors,
            'T_event_date_hit_histogram': sorted_hist(T_hist),
            'S_site_hit_histogram': sorted_hist(S_hist),
            'actual_T': actual_T,
            'actual_S': actual_S,
            'fraction_year_shift_vectors_T_ge_actual': fraction(ge_T, n_vectors),
            'fraction_year_shift_vectors_T_le_actual': fraction(le_T, n_vectors),
            'fraction_year_shift_vectors_S_ge_actual': fraction(ge_S, n_vectors),
            'fraction_year_shift_vectors_S_le_actual': fraction(le_S, n_vectors),
            'fraction_year_shift_vectors_T_zero': fraction(T_hist.get(0, 0), n_vectors),
            'max_T_in_null': max(T_hist),
            'max_S_in_null': max(S_hist),
        }

    actual_Ts = [actual[f'radius_{r}_days']['T_event_date_hits'] for r in radii]
    if all(x == 0 for x in actual_Ts):
        status = 'NO_PALOMAR_OBSERVATION_OPPORTUNITY_IN_ANY_PREREGISTERED_R4_1_PRIMARY_WINDOW'
    else:
        status = 'PALOMAR_OBSERVATION_OPPORTUNITY_PRESENT__INTERPRET_ONLY_AGAINST_R4_1_PREREGISTERED_NULLS'

    out = {
        'schema': 'janus.cosmos.bridge.r4_1_second_holdout_null_result.v1',
        'artifact_id': 'BRIDGE-R4.1-SECOND-HOLDOUT-NULL-RESULT-2026-09-28-v1.0',
        'generated_at_utc': datetime.now(timezone.utc).isoformat(),
        'status': status,
        'execution_contract': {
            'earth_cohort_artifact_id': earth['artifact_id'],
            'earth_cohort_status': earth['status'],
            'prereg_artifact_id': prereg['artifact_id'],
            'prereg_status': prereg['status'],
            'primary_holdout_site_n': len(site_ids),
            'primary_holdout_event_date_n': sum(len(v) for v in site_rows.values()),
            'site_ids': site_ids,
            'radii_days': radii,
            'palomar_study_window': [study_start.isoformat(), study_end.isoformat()],
            'palomar_observation_day_rows': len(pal_rows_by_date),
            'primary_null_name': prereg['primary_null']['name'],
            'seasonality_sensitivity_null_name': prereg['seasonality_sensitivity_null']['name'],
            'hanford_holdback_evaluated': False,
            'hanford_holdback_date_disclosed_to_builder_but_not_queried_as_an_event': '1950-07-30'
        },
        'source_hashes': {
            'earth_ledger_sha256': sha256(args.earth_ledger),
            'palomar_day_ledger_sha256': sha256(args.palomar_ledger),
            'prereg_sha256': sha256(args.prereg),
        },
        'actual_windows': actual,
        'primary_null': primary_null,
        'seasonality_sensitivity_null': seasonality_null,
        'sealed_holdback': {
            'HANFORD_AEC_PLANT_1950_07_30_evaluated': False,
            'reason': 'Primary federal byte/page binding is not closed; R4.1 intentionally preserves its future holdout blindness.'
        },
        'interpretation': {
            'T_r': 'Counts R4.1 frozen holdout site-dates with at least one Palomar observation opportunity. It is not a count of anomalous objects.',
            'S_r': 'Counts R4.1 independent sites with at least one temporal observation-opportunity hit.',
            'primary_null': 'Each new site cluster moves independently as a rigid day pattern; exact aggregate T and S distributions are computed by convolution, not Monte Carlo.',
            'seasonality_null': 'Each new site cluster moves only by whole calendar years, preserving month/day as a sensitivity check for seasonal Palomar scheduling.',
            'sequential_status': 'R4.1 stands as a second holdout. R4 is already known and is not retroactively merged into a fictitious jointly blinded p-value.'
        },
        'claim_ceiling': 'THIS RESULT TESTS TEMPORAL OBSERVATION-OPPORTUNITY COINCIDENCE IN A SECOND INDEPENDENT HOLDOUT ONLY. IT DOES NOT ESTABLISH UAP IDENTITY, NONHUMAN ORIGIN, NUCLEAR ATTRACTION, OR EARTH-SPACE CAUSAL COUPLING.',
        'canonical_seal': 'THE SECOND HOLDOUT WAS FROZEN BEFORE THE SKY CALENDAR OPENED; HANFORD REMAINED CLOSED.'
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({
        'status': status,
        'site_n': len(site_ids),
        'event_date_n': sum(len(v) for v in site_rows.values()),
        'actual': {
            key: {
                'T': val['T_event_date_hits'],
                'S': val['S_sites_with_at_least_one_hit'],
                'U': val['U_unique_observation_days']
            }
            for key, val in actual.items()
        },
        'hanford_evaluated': False
    }, sort_keys=True))


if __name__ == '__main__':
    main()
