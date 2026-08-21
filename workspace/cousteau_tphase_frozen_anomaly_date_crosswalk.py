#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from collections import Counter
from datetime import timezone
from pathlib import Path

from cousteau_ea_tphase_blind_cluster_v7 import acquire_exact_file, parse_exact
from cousteau_tphase_celestial_spacetime_crossmatch import parse_time_code
from cousteau_ea_tphase_blind_cluster import sha256_bytes

FROZEN_DATES = {
    'PALOMAR_XE325': {'month':4,'day':12,'full_date':'1950-04-12','role':'COSMIC_ANOMALY_DATE_PREEXISTING_IN_PROJECT'},
    'SHAG_HARBOUR_1967': {'month':10,'day':4,'full_date':'1967-10-04','role':'HISTORICAL_UAP_DATE'},
    'NIMITZ_2004': {'month':11,'day':14,'full_date':'2004-11-14','role':'HISTORICAL_UAP_DATE'},
    'OMAHA_2019': {'month':7,'day':15,'full_date':'2019-07-15','role':'HISTORICAL_UAP_DATE'},
}
CROSS = Path('data/cousteau/JANUS-ECHO-COUSTEAU-TPHASE-CELESTIAL-SPACETIME-CROSSMATCH-RUN-001-2026-08-21-v1.0.json')


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True); a=ap.parse_args()
    archive,gz,raw,member,trace=acquire_exact_file(); full,_=parse_exact(raw)
    if len(full)!=5943: raise RuntimeError(f'row count drift {len(full)}')
    dts=[parse_time_code(str(x)) for x in full.source_time_code]
    md=[(x.month,x.day) for x in dts]
    counts=Counter(md)
    years=Counter(x.year for x in dts)
    cross=json.loads(CROSS.read_text(encoding='utf-8'))
    top25=cross['top_matches_by_bisector_subpoint_distance']
    top_md=[(parse_time_code(x['source_time_code']).month,parse_time_code(x['source_time_code']).day) for x in top25]
    rows={}
    for k,v in FROZEN_DATES.items():
        key=(v['month'],v['day']); idx=[i for i,x in enumerate(md) if x==key]
        rows[k]={**v,'catalog_events_same_month_day':len(idx),'catalog_fraction_same_month_day':len(idx)/len(full),
                 'top25_celestial_matches_same_month_day':sum(1 for x in top_md if x==key),
                 'years_present_same_month_day':sorted({dts[i].year for i in idx}),
                 'interpretation':'CALENDAR_DATE_COINCIDENCE_ONLY__NOT_SIDEREAL_OR_CAUSAL_ALIGNMENT'}
    top1=top25[0]
    p=rows['PALOMAR_XE325']
    out={
      'artifact_id':'JANUS-ECHO-COUSTEAU-FROZEN-ANOMALY-DATE-CROSSWALK-RUN-001-2026-08-21-v1.0',
      'source':{'authoritative_rows':len(full),'catalog_ascii_sha256':sha256_bytes(raw),'member':member,'year_distribution':{str(k):v for k,v in sorted(years.items())}},
      'frozen_dates':rows,
      'emergent_top_celestial_match':{'event_time_utc':top1['event_time_utc'],'event_lat':top1['event_lat'],'event_lon':top1['event_lon'],
                                      'same_month_day_as_palomar': top1['event_time_utc'][5:10]=='04-12',
                                      'years_after_palomar':2014-1950,
                                      'palomar_month_day_frequency_in_catalog':p['catalog_fraction_same_month_day']},
      'hard_rules':['FROZEN_HISTORICAL_DATES_ONLY','CALENDAR_MATCH_IS_NOT_SIDEREAL_MATCH','NO_MULTIPLE_DATE_FISHING_BEYOND_LISTED_DATES','TARGET_EVIDENCE_NOT_PROMOTED'],
      'status':'RUN_COMPLETE'
    }
    Path(a.output).write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'status':'SUCCESS','top1':out['emergent_top_celestial_match'],'date_counts':{k:v['catalog_events_same_month_day'] for k,v in rows.items()}},indent=2))
    return 0
if __name__=='__main__': raise SystemExit(main())
