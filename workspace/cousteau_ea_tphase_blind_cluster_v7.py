#!/usr/bin/env python3
from __future__ import annotations

import argparse, gzip, io, json, re, tarfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
from bs4 import BeautifulSoup

from cousteau_ea_tphase_blind_cluster import (
    DBSCAN_GRID, NULL_DOMAIN, NULL_SAMPLES, NULL_SEED,
    blind_cluster, reveal_and_score, sha256_bytes,
)

DATASET_UID='30497'; FILE_UID='2504732'; DOI='10.26022/IEDA/330497'
PAPER_REPORTED_EVENT_COUNT=6843
LANDING=f'https://www.marine-geo.org/tools/files/{DATASET_UID}'
UID_URL=f'https://www.marine-geo.org/tools/search/file_uids.php?data_set_uid={DATASET_UID}'
MODAL_URL='https://www.marine-geo.org/services/download/download_modal.php'
ANCHOR_LAT=-3.865418; ANCHOR_LON=3.854924
ROW_RE=re.compile(r'^\s*(\d{14})\s+(\d+)\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s*$')


def now(): return datetime.now(timezone.utc).isoformat()


def acquire_exact_file():
    s=requests.Session(); s.headers.update({'User-Agent':'Janus-Echo-Cousteau/1.6 exact-format scientific replay','Referer':LANDING})
    trace=[]
    u=s.get(UID_URL,timeout=45); u.raise_for_status(); uids=[str(x) for x in u.json()]
    trace.append({'stage':'authoritative_prebind','dataset_uid':DATASET_UID,'file_uids':uids})
    if FILE_UID not in uids: raise RuntimeError(f'file UID drift: {uids}')
    m=s.post(MODAL_URL,data={'FileDownload':FILE_UID,'data_set_uid':DATASET_UID},timeout=60); m.raise_for_status()
    f=BeautifulSoup(m.text,'html.parser').find('form',id='data_link') or BeautifulSoup(m.text,'html.parser').find('form')
    if not f: raise RuntimeError('download form missing')
    action=f.get('action')
    r=s.post(action,data={'purpose':'Research','client':'DataLink','force_download':'1','data_uids':FILE_UID},timeout=180); r.raise_for_status()
    trace.append({'stage':'download','status':r.status_code,'bytes':len(r.content),'content_type':r.headers.get('content-type'),
                  'content_disposition':r.headers.get('content-disposition'),'archive_sha256':sha256_bytes(r.content)})
    target=None; member_listing=[]
    with tarfile.open(fileobj=io.BytesIO(r.content),mode='r:*') as tf:
        for mem in tf.getmembers():
            if not mem.isfile(): continue
            member_listing.append({'name':mem.name,'bytes':mem.size})
            if mem.name.endswith('EA_CTBTO_catalog_all.dat.gz'):
                fo=tf.extractfile(mem)
                if fo: target=(mem.name,fo.read())
    trace.append({'stage':'tar_members','members':member_listing})
    if target is None: raise RuntimeError('EA_CTBTO_catalog_all.dat.gz not found in authoritative archive')
    member_name,gz=target
    raw=gzip.decompress(gz)
    trace.append({'stage':'nested_gzip','member':member_name,'gzip_bytes':len(gz),'gzip_sha256':sha256_bytes(gz),
                  'uncompressed_bytes':len(raw),'uncompressed_sha256':sha256_bytes(raw)})
    return r.content,gz,raw,member_name,trace


def parse_exact(raw: bytes):
    text=raw.decode('utf-8',errors='strict')
    rows=[]; rejected=[]; header=[]
    for lineno,line in enumerate(text.splitlines(),start=1):
        m=ROW_RE.match(line)
        if not m:
            if line.strip():
                if lineno <= 20: header.append({'line':lineno,'text':line})
                elif len(rejected)<20: rejected.append({'line':lineno,'text':line[:300]})
            continue
        t,n,lat,lon,laterr,lonerr,terr,mag=m.groups()
        row={'source_time_code':t,'n_hydrophones':int(n),'lat':float(lat),'lon':float(lon),
             'lat_error_deg':float(laterr),'lon_error_deg':float(lonerr),'source_time_error_s':float(terr),'source_magnitude_db':float(mag),'source_line':lineno}
        if not (-90 <= row['lat'] <= 90 and -180 <= row['lon'] <= 180):
            raise RuntimeError(f'coordinate range violation line {lineno}')
        rows.append(row)
    if len(rows)<1000: raise RuntimeError(f'exact parser produced only {len(rows)} rows')
    df=pd.DataFrame(rows)
    meta={
      'parser':'EXACT_EA_CTBTO_EIGHT_FIELD_WHITESPACE_V1','row_regex':ROW_RE.pattern,'data_rows':len(rows),
      'paper_reported_event_count':PAPER_REPORTED_EVENT_COUNT,'paper_minus_file_rows':PAPER_REPORTED_EVENT_COUNT-len(rows),
      'paper_vs_file_count_match':len(rows)==PAPER_REPORTED_EVENT_COUNT,'header_nonblank_lines':header,'unexpected_nonmatching_data_lines_sample':rejected,
      'n_hydrophones_distribution':{str(k):int(v) for k,v in sorted(Counter(df.n_hydrophones).items())},
      'lat_range':[float(df.lat.min()),float(df.lat.max())],'lon_range':[float(df.lon.min()),float(df.lon.max())],
      'time_code_first':str(df.iloc[0].source_time_code),'time_code_last':str(df.iloc[-1].source_time_code),
      'columns':['source_time_code','n_hydrophones','lat','lon','lat_error_deg','lon_error_deg','source_time_error_s','source_magnitude_db'],
    }
    return df,meta


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True); ap.add_argument('--status-output',required=True); a=ap.parse_args()
    out=Path(a.output); st=Path(a.status_output); out.parent.mkdir(parents=True,exist_ok=True)
    artifact='JANUS-ECHO-COUSTEAU-EA-TPHASE-BLIND-CLUSTER-RUN-007-2026-08-21-v1.6'; started=now()
    try:
        archive,gz,raw,member,trace=acquire_exact_file(); full,parse_meta=parse_exact(raw)
        coords=full[['lat','lon']].copy()

        # PHASE A: blind clustering. Anchor is not an input to this function.
        blind=blind_cluster(coords); blind_hash=blind['freeze_sha256']

        # PHASE B: reveal frozen anchor only after blind hash exists.
        reveal=reveal_and_score(coords,blind,ANCHOR_LAT,ANCHOR_LON)
        nearest=[x['nearest_cluster']['anchor_to_center_km'] for x in reveal['configs'] if x.get('nearest_cluster')]
        p95=any(x.get('nearest_cluster') and x['nearest_cluster']['anchor_inside_cluster_p95_radius'] for x in reveal['configs'])
        maxr=any(x.get('nearest_cluster') and x['nearest_cluster']['anchor_inside_cluster_max_radius'] for x in reveal['configs'])
        verdict='ANCHOR_OVERLAPS_BLIND_CLUSTER_P95__TECTONIC_CONTROL_REQUIRED' if p95 else 'NO_BLIND_CLUSTER_P95_OVERLAP_WITH_FROZEN_ANCHOR'
        count_mismatch=len(full)!=PAPER_REPORTED_EVENT_COUNT

        result={
          'artifact_id':artifact,'research_branch':'Janus-Echo-Кусто','started_at_utc':started,'completed_at_utc':now(),
          'source':{'doi':DOI,'dataset':'EA_Hydroacoustics','data_set_uid':int(DATASET_UID),'file_uid':FILE_UID,
                    'archive_sha256':sha256_bytes(archive),'archive_bytes':len(archive),'member_name':member,'member_gzip_sha256':sha256_bytes(gz),
                    'catalog_ascii_sha256':sha256_bytes(raw),'catalog_ascii_bytes':len(raw),'raw_committed':False,'parse':parse_meta,'acquisition_trace':trace},
          'count_reconciliation_gate':{
             'paper_reported_events':PAPER_REPORTED_EVENT_COUNT,'authoritative_download_rows':int(len(full)),
             'difference':PAPER_REPORTED_EVENT_COUNT-int(len(full)),'status':'OPEN_MISMATCH_REQUIRES_SOURCE_RECONCILIATION' if count_mismatch else 'MATCH',
             'rule':'DO_NOT_SYNTHESIZE_MISSING_ROWS__CLUSTER_ONLY_AUTHORITATIVE_DOWNLOADED_ROWS',
          },
          'preregistration':{'anchor_hidden_during_clustering':True,'clustering_parameters_frozen_before_anchor_reveal':True,
                             'dbscan_grid':DBSCAN_GRID,'null_domain':NULL_DOMAIN,'null_samples':NULL_SAMPLES,'null_seed':NULL_SEED,
                             'parameter_change_from_original_blind_contract':False},
          'blind_phase':blind,'post_reveal':reveal,
          'summary':{'blind_freeze_sha256':blind_hash,'authoritative_rows':int(len(full)),'paper_reported_rows':PAPER_REPORTED_EVENT_COUNT,
                     'count_mismatch':count_mismatch,'n_hydrophones_distribution':parse_meta['n_hydrophones_distribution'],
                     'nearest_catalog_event_to_anchor_km':reveal['nearest_event']['distance_km'],
                     'nearest_blind_cluster_center_across_grid_km':round(min(nearest),3) if nearest else None,
                     'anchor_inside_any_blind_cluster_p95_radius':p95,'anchor_inside_any_blind_cluster_max_radius':maxr,
                     'verdict':verdict,'semantic_status':'UNCONFIRMED','tectonic_control_required':True},
          'hard_rules':['EXACT_FILE_FORMAT_PARSER_ONLY','NO_SYNTHETIC_ROWS_TO_REACH_PAPER_COUNT','BLIND_CLUSTER_BEFORE_ANCHOR_REVEAL',
                        'NO_PARAMETER_RETUNING_AFTER_REVEAL','MID_ATLANTIC_RIDGE_SEISMICITY_IS_MANDATORY_TECTONIC_CONTROL',
                        'DISTANCE_IS_NOT_CAUSATION','NO_RECENTERING','NO_UNDERWATER_PYRAMID_DETECTED_YET'],
          'status':'BLIND_CLUSTER_RUN_COMPLETE_WITH_COUNT_RECONCILIATION_GATE' if count_mismatch else 'BLIND_CLUSTER_RUN_COMPLETE'}
        out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
        st.write_text(json.dumps({'artifact_id':artifact,'status':'SUCCESS','completed_at_utc':now(),'authoritative_rows':int(len(full)),
                                  'paper_reported_rows':PAPER_REPORTED_EVENT_COUNT,'count_mismatch':count_mismatch,'blind_freeze_sha256':blind_hash,
                                  'verdict':verdict,'result_path':str(out)},indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
        print(json.dumps(result['summary'],indent=2,ensure_ascii=False)); return 0
    except Exception as e:
        payload={'artifact_id':artifact,'status':'BLOCKED_DATA_ACQUISITION_OR_PARSE','started_at_utc':started,'completed_at_utc':now(),
                 'error_type':type(e).__name__,'error':str(e),'scientific_interpretation':'BLOCKER_ONLY__NOT_NEGATIVE_CLUSTER_RESULT'}
        st.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n',encoding='utf-8'); print(json.dumps(payload,indent=2)); return 2

if __name__=='__main__': raise SystemExit(main())
