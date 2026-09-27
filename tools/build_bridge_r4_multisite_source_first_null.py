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
    # The frozen R4 cohort contains no Feb-29 dates. Refuse rather than silently
    # change month/day if a future cohort introduces one.
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

    assert earth['status'] == 'FROZEN_BEFORE_NEW_HOLDOUT_PALOMAR_MEMBERSHIP_LOOKUP'
    assert prereg['status'] == 'FROZEN_BEFORE_NEW_HOLDOUT_PALOMAR_MEMBERSHIP_EVALUATION'
    assert prereg['earth_cohort_artifact'] == earth['artifact_id']

    radii = [int(x) for x in prereg['frozen_windows']['radii_days']]
    assert radii == [0, 1, 3], radii
    widest = max(radii)

    frozen = prereg['frozen_primary_holdout']['sites']
    site_rows: dict[str, list[dict]] = {}
    for site in earth['primary_holdout_sites']:
        sid = site['site_id']
        rows = [r for r in site['rows'] if r.get('cohort_role') == 'HOLDOUT_PRIMARY_NEW' and r.get('bridge_eligible') is True]
        rows = sorted(rows, key=lambda r: r['date'])
        assert [r['date'] for r in rows] == frozen[sid], (sid, rows, frozen[sid])
        site_rows[sid] = rows

    site_ids = list(frozen.keys())
    assert site_ids == ['LOS_ALAMOS', 'SANDIA_BASE', 'SAVANNAH_RIVER_PLANT'], site_ids
    assert len(site_rows) == prereg['frozen_primary_holdout']['site_n'] == 3
    assert sum(len(v) for v in site_rows.values()) == prereg['frozen_primary_holdout']['event_date_n'] == 6

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
        hit_count = 0
        union_dates: set[date] = set()
        details = []
        for src in site_rows[sid]:
            original = parse_date(src['date'])
            center = transform(original)
            hits = event_window_hits(center, radius)
            hit = bool(hits)
            hit_count += int(hit)
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
            'event_hits': hit_count,
            'site_has_hit': hit_count > 0,
            'unique_observation_days': len(union_dates),
        }
        if include_details:
            out['event_rows'] = details
            out['union_observation_opportunity_dates'] = [d.isoformat() for d in sorted(union_dates)]
        return out

    # Freeze actual results only after both Earth cohort and this prereg exist.
    actual: dict[str, dict] = {}
    for radius in radii:
        key = f'radius_{radius}_days'
        sites = {}
        all_union: set[str] = set()
        T = 0
        S = 0
        for sid in site_ids:
            ev = evaluate_site(sid, radius, lambda d: d, include_details=True)
            sites[sid] = ev
            T += ev['event_hits']
            S += int(ev['site_has_hit'])
            all_union.update(ev['union_observation_opportunity_dates'])
        actual[key] = {
            'T_event_date_hits': T,
            'S_sites_with_at_least_one_hit': S,
            'U_unique_observation_days': len(all_union),
            'union_observation_opportunity_dates': sorted(all_union),
            'sites': sites,
        }

    # Primary null: each site cluster gets an independent integer-day rigid shift.
    primary_support: dict[str, list[int]] = {}
    for sid in site_ids:
        dates = [parse_date(r['date']) for r in site_rows[sid]]
        min_shift = (study_start + timedelta(days=widest) - min(dates)).days
        max_shift = (study_end - timedelta(days=widest) - max(dates)).days
        shifts = list(range(min_shift, max_shift + 1))
        assert 0 in shifts
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
            'shift_vectors_at_max_T': T_hist[max(T_hist)],
            'max_S_in_null': max(S_hist),
            'shift_vectors_at_max_S': S_hist[max(S_hist)],
        }

    # Sensitivity null: preserve calendar month/day by independent whole-year shifts.
    year_support: dict[str, list[int]] = {}
    for sid in site_ids:
        src_dates = [parse_date(r['date']) for r in site_rows[sid]]
        offsets = []
        for y in range(-100, 101):
            try:
                moved = [add_years(d, y) for d in src_dates]
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
            for y in year_support[sid]:
                ev = evaluate_site(sid, radius, lambda d, years=y: add_years(d, years))
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
            'year_shift_vectors_at_max_T': T_hist[max(T_hist)],
            'max_S_in_null': max(S_hist),
            'year_shift_vectors_at_max_S': S_hist[max(S_hist)],
        }

    actual_Ts = [actual[f'radius_{r}_days']['T_event_date_hits'] for r in radii]
    if all(x == 0 for x in actual_Ts):
        status = 'NO_PALOMAR_OBSERVATION_OPPORTUNITY_IN_ANY_PREREGISTERED_R4_HOLDOUT_WINDOW'
    else:
        status = 'PALOMAR_OBSERVATION_OPPORTUNITY_PRESENT__INTERPRET_ONLY_AGAINST_PREREGISTERED_NULLS'

    out = {
        'schema': 'janus.cosmos.bridge.multisite_source_first_null_result.v1',
        'artifact_id': 'BRIDGE-R4-MULTISITE-SOURCE-FIRST-NULL-RESULT-2026-09-28-v1.0',
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
            'known_outcome_calibration_excluded_from_primary': True,
            'primary_null_name': prereg['primary_null']['name'],
            'seasonality_sensitivity_null_name': prereg['seasonality_sensitivity_null']['name'],
        },
        'source_hashes': {
            'earth_ledger_sha256': sha256(args.earth_ledger),
            'palomar_day_ledger_sha256': sha256(args.palomar_ledger),
            'prereg_sha256': sha256(args.prereg),
        },
        'actual_windows': actual,
        'primary_null': primary_null,
        'seasonality_sensitivity_null': seasonality_null,
        'calibration_reuse': {
            'Oak_Ridge_in_primary_statistic': False,
            'reason': 'Oak Ridge R2/R3 outcomes were known before R4 and therefore cannot contribute holdout evidence.'
        },
        'interpretation': {
            'T_r': 'Counts frozen holdout site-dates with at least one Palomar observation opportunity. It is not a count of anomalous objects.',
            'S_r': 'Counts independently sourced sites with at least one temporal observation-opportunity hit and is a site-balance diagnostic.',
            'primary_null': 'Each site cluster moves independently as a rigid day pattern; exact aggregate T and S null distributions are computed by count convolution, not Monte Carlo.',
            'seasonality_null': 'Each site cluster moves only by whole calendar years, preserving month/day as a sensitivity check for the Palomar observing calendar.',
            'no_observation_opportunity': 'Means the materialized POSS-I observation-day ledger contains no observing day in the frozen window; it is not an observed-zero sky.',
            'temporal_overlap': 'Means only that Palomar supplied an observation opportunity near the independently frozen Earth date.'
        },
        'claim_ceiling': 'R4 TESTS TEMPORAL OBSERVATION-OPPORTUNITY COINCIDENCE ONLY. IT DOES NOT ESTABLISH THAT AN EARTH EVENT WAS IMAGED, THAT ANY PALOMAR CANDIDATE IS RELATED, OR THAT NUCLEAR ATTRACTION, NONHUMAN ORIGIN, OR CAUSATION EXISTS.',
        'canonical_seal': 'THE EARTH COHORT WAS FROZEN BEFORE THE SKY CALENDAR WAS OPENED. REPORT EVERY RADIUS AND BOTH NULLS; TEMPORAL OPPORTUNITY IS NOT ORIGIN.'
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({
        'status': status,
        'actual_T': {str(r): actual[f'radius_{r}_days']['T_event_date_hits'] for r in radii},
        'actual_S': {str(r): actual[f'radius_{r}_days']['S_sites_with_at_least_one_hit'] for r in radii},
        'primary_upper_tail': {str(r): primary_null[f'radius_{r}_days']['fraction_shift_vectors_T_ge_actual'] for r in radii},
        'seasonal_upper_tail': {str(r): seasonality_null[f'radius_{r}_days']['fraction_year_shift_vectors_T_ge_actual'] for r in radii},
    }, sort_keys=True))


if __name__ == '__main__':
    main()
