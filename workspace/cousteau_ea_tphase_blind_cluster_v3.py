#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from datetime import datetime, timezone
from pathlib import Path
import requests

from cousteau_ea_tphase_blind_cluster import (
    EXPECTED_EVENT_COUNT, DBSCAN_GRID, NULL_DOMAIN, NULL_SAMPLES, NULL_SEED,
    blind_cluster, reveal_and_score, parse_catalog, sha256_bytes,
)

DATASET_UID=30497
FILE_UID="2504732"
FILE_UID_ENDPOINT=f"https://www.marine-geo.org/tools/search/file_uids.php?data_set_uid={DATASET_UID}"
FDS="https://www.marine-geo.org/services/FileDownloadServer"
DOI="10.26022/IEDA/330497"


def acquire():
    s=requests.Session(); s.headers['User-Agent']='Janus-Echo-Cousteau/1.2 scientific reproducibility audit'
    checks=[]
    r=s.get(FILE_UID_ENDPOINT,timeout=45); r.raise_for_status()
    uids=r.json(); checks.append({'stage':'file_uid_lookup','url':r.url,'status':r.status_code,'payload':uids})
    if FILE_UID not in [str(x) for x in uids]:
        raise RuntimeError(f'frozen file UID {FILE_UID} no longer listed for dataset {DATASET_UID}: {uids}')
    variants=[
      {'data_uid':FILE_UID},
      {'file_uid':FILE_UID},
      {'uid':FILE_UID},
    ]
    last=[]
    for params in variants:
        rr=s.get(FDS,params=params,timeout=120,allow_redirects=True)
        meta={'stage':'file_download','params':params,'requested_url':rr.request.url,'final_url':rr.url,'status':rr.status_code,'bytes':len(rr.content),'content_type':rr.headers.get('content-type'),'content_disposition':rr.headers.get('content-disposition'),'prefix':rr.text[:250] if len(rr.content)<3000 or 'text' in (rr.headers.get('content-type') or '') else None}
        checks.append(meta)
        if rr.status_code!=200: continue
        try:
            coords,parse_meta=parse_catalog(rr.content)
        except Exception as e:
            last.append(f'{params}: {type(e).__name__}: {e}')
            continue
        if len(coords)>=1000:
            return rr.content,coords,parse_meta,meta,checks
        last.append(f'{params}: only {len(coords)} coordinates')
    raise RuntimeError('direct FileDownloadServer variants failed to yield catalog: '+' | '.join(last))


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True); ap.add_argument('--status-output',required=True); a=ap.parse_args()
    out=Path(a.output); st=Path(a.status_output); out.parent.mkdir(parents=True,exist_ok=True)
    artifact='JANUS-ECHO-COUSTEAU-EA-TPHASE-BLIND-CLUSTER-RUN-003-2026-08-21-v1.2'
    started=datetime.now(timezone.utc).isoformat()
    try:
        raw,coords,parse_meta,dlmeta,checks=acquire()
        # BLIND PHASE: no LOVE-EDEM coordinate is passed into clustering.
        blind=blind_cluster(coords)
        blind_hash=blind['freeze_sha256']
        # REVEAL PHASE starts only after blind result is frozen+hashed.
        reveal=reveal_and_score(coords,blind,-3.865418,3.854924)
        nearest_clusters=[x['nearest_cluster']['anchor_to_center_km'] for x in reveal['configs'] if x.get('nearest_cluster')]
        any_p95=any(x.get('nearest_cluster') and x['nearest_cluster']['anchor_inside_cluster_p95_radius'] for x in reveal['configs'])
        any_max=any(x.get('nearest_cluster') and x['nearest_cluster']['anchor_inside_cluster_max_radius'] for x in reveal['configs'])
        verdict='ANCHOR_OVERLAPS_BLIND_CLUSTER_P95__TECTONIC_CONTROL_REQUIRED' if any_p95 else 'NO_BLIND_CLUSTER_P95_OVERLAP_WITH_FROZEN_ANCHOR'
        result={
          'artifact_id':artifact,'research_branch':'Janus-Echo-Кусто','started_at_utc':started,'completed_at_utc':datetime.now(timezone.utc).isoformat(),
          'source':{'doi':DOI,'dataset':'EA_Hydroacoustics','data_set_uid':DATASET_UID,'file_uid':FILE_UID,'file_uid_endpoint':FILE_UID_ENDPOINT,'raw_sha256':sha256_bytes(raw),'raw_bytes':len(raw),'raw_committed':False,'license':'CC BY-NC-SA 3.0','expected_event_count':EXPECTED_EVENT_COUNT,'parsed_valid_coordinate_count':len(coords),'expected_count_exact_match':len(coords)==EXPECTED_EVENT_COUNT,'parse':parse_meta,'download':dlmeta,'acquisition_checks':checks},
          'preregistration':{'anchor_hidden_during_clustering':True,'clustering_parameters_frozen_before_anchor_reveal':True,'dbscan_grid':DBSCAN_GRID,'null_domain':NULL_DOMAIN,'null_samples':NULL_SAMPLES,'null_seed':NULL_SEED},
          'blind_phase':blind,
          'post_reveal':reveal,
          'summary':{'blind_freeze_sha256':blind_hash,'nearest_catalog_event_to_anchor_km':reveal['nearest_event']['distance_km'],'nearest_blind_cluster_center_across_grid_km':round(min(nearest_clusters),3) if nearest_clusters else None,'anchor_inside_any_blind_cluster_p95_radius':any_p95,'anchor_inside_any_blind_cluster_max_radius':any_max,'verdict':verdict,'semantic_status':'UNCONFIRMED'},
          'hard_rules':['FILE_UID_FROZEN_BEFORE_DATA_READ','BLIND_CLUSTER_BEFORE_ANCHOR_REVEAL','NO_PARAMETER_RETUNING_AFTER_REVEAL','MID_ATLANTIC_RIDGE_SEISMICITY_IS_MANDATORY_TECTONIC_CONTROL','RECTANGULAR_NULL_IS_DIAGNOSTIC_NOT_FORMAL','DISTANCE_IS_NOT_CAUSATION','NO_RECENTERING','NO_UNDERWATER_PYRAMID_DETECTED_YET'],
          'status':'BLIND_CLUSTER_RUN_COMPLETE'}
        out.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
        st.write_text(json.dumps({'artifact_id':artifact,'status':'SUCCESS','completed_at_utc':datetime.now(timezone.utc).isoformat(),'parsed_count':len(coords),'raw_sha256':sha256_bytes(raw),'verdict':verdict,'result_path':str(out)},indent=2),encoding='utf-8')
        print(json.dumps(result['summary'],indent=2)); return 0
    except Exception as e:
        payload={'artifact_id':artifact,'status':'BLOCKED_DATA_ACQUISITION_OR_PARSE','started_at_utc':started,'completed_at_utc':datetime.now(timezone.utc).isoformat(),'error_type':type(e).__name__,'error':str(e)}
        st.write_text(json.dumps(payload,indent=2,ensure_ascii=False),encoding='utf-8'); print(json.dumps(payload,indent=2)); return 2

if __name__=='__main__': raise SystemExit(main())
