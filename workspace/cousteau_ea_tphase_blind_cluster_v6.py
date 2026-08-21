#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import io
import json
import tarfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import requests
from bs4 import BeautifulSoup

from cousteau_ea_tphase_blind_cluster import (
    EXPECTED_EVENT_COUNT,
    DBSCAN_GRID,
    NULL_DOMAIN,
    NULL_SAMPLES,
    NULL_SEED,
    blind_cluster,
    reveal_and_score,
    parse_catalog,
    sha256_bytes,
)

DATASET_UID = "30497"
FILE_UID = "2504732"
DOI = "10.26022/IEDA/330497"
LANDING = f"https://www.marine-geo.org/tools/files/{DATASET_UID}"
UID_URL = f"https://www.marine-geo.org/tools/search/file_uids.php?data_set_uid={DATASET_UID}"
MODAL_URL = "https://www.marine-geo.org/services/download/download_modal.php"
ACCEPT_URL = "https://api.marine-geo.org/services/download/download_accept.php"
ANCHOR_LAT = -3.865418
ANCHOR_LON = 3.854924
MAX_DEPTH = 4
MIN_COORDS = 1000


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def looks_text(raw: bytes) -> bool:
    if not raw:
        return False
    s = raw[:8192]
    if s.count(b"\x00") > max(2, len(s)//100):
        return False
    printable = sum((b in b"\t\n\r") or 32 <= b <= 126 or b >= 128 for b in s)
    return printable / max(1, len(s)) >= 0.72


def walk_payloads(raw: bytes, label: str, trace: list[dict], depth: int = 0, seen=None):
    """Recursively unwrap TAR/ZIP/GZIP members. Transport only; no scientific thresholds change."""
    if seen is None:
        seen = set()
    h = sha256_bytes(raw)
    if h in seen or depth > MAX_DEPTH:
        return
    seen.add(h)
    trace.append({"stage":"payload","depth":depth,"label":label,"bytes":len(raw),"sha256":h,"magic_hex":raw[:24].hex()})
    yield label, raw

    # TAR (including compressed tar streams understood by tarfile)
    try:
        with tarfile.open(fileobj=io.BytesIO(raw), mode="r:*") as tf:
            members=[m for m in tf.getmembers() if m.isfile()]
            trace.append({"stage":"archive","depth":depth,"format":"tar","label":label,
                          "members":[{"name":m.name,"bytes":m.size} for m in members[:200]]})
            for m in members:
                if m.size <= 0 or m.size > 25_000_000:
                    continue
                f=tf.extractfile(m)
                if f is None:
                    continue
                b=f.read()
                yield from walk_payloads(b, f"{label}::tar::{m.name}", trace, depth+1, seen)
            return
    except (tarfile.ReadError, EOFError):
        pass

    # ZIP
    if raw.startswith(b"PK\x03\x04"):
        try:
            with zipfile.ZipFile(io.BytesIO(raw)) as zf:
                infos=[i for i in zf.infolist() if not i.is_dir()]
                trace.append({"stage":"archive","depth":depth,"format":"zip","label":label,
                              "members":[{"name":i.filename,"bytes":i.file_size} for i in infos[:200]]})
                for i in infos:
                    if 0 < i.file_size <= 25_000_000:
                        yield from walk_payloads(zf.read(i), f"{label}::zip::{i.filename}", trace, depth+1, seen)
        except Exception as e:
            trace.append({"stage":"archive_error","format":"zip","label":label,"error":f"{type(e).__name__}: {e}"})
        return

    # Standalone nested GZIP (the observed MGDS TAR member is EA_CTBTO_catalog_all.dat.gz).
    if raw.startswith(b"\x1f\x8b"):
        try:
            b=gzip.decompress(raw)
            trace.append({"stage":"archive","depth":depth,"format":"gzip","label":label,
                          "uncompressed_bytes":len(b),"uncompressed_sha256":sha256_bytes(b)})
            yield from walk_payloads(b, f"{label}::gzip", trace, depth+1, seen)
        except Exception as e:
            trace.append({"stage":"archive_error","format":"gzip","label":label,"error":f"{type(e).__name__}: {e}"})


def acquire():
    s=requests.Session()
    s.headers.update({'User-Agent':'Janus-Echo-Cousteau/1.5 nested-archive scientific reproducibility audit','Referer':LANDING})
    trace=[]

    u=s.get(UID_URL,timeout=45); u.raise_for_status(); uids=[str(x) for x in u.json()]
    trace.append({'stage':'authoritative_prebind','dataset_uid':DATASET_UID,'file_uids':uids,'status':u.status_code,'url':u.url})
    if FILE_UID not in uids:
        raise RuntimeError(f'file UID drift: expected {FILE_UID}, got {uids}')

    m=s.post(MODAL_URL,data={'FileDownload':FILE_UID,'data_set_uid':DATASET_UID},timeout=60,allow_redirects=True); m.raise_for_status()
    soup=BeautifulSoup(m.text,'html.parser'); form=soup.find('form',id='data_link') or soup.find('form')
    if form is None:
        raise RuntimeError('MGDS download form missing')
    action=form.get('action') or ACCEPT_URL
    trace.append({'stage':'download_contract','modal_status':m.status_code,'action':action,
                  'payload':{'purpose':'Research','client':'DataLink','force_download':'1','data_uids':FILE_UID}})

    r=s.post(action,data={'purpose':'Research','client':'DataLink','force_download':'1','data_uids':FILE_UID},timeout=180,allow_redirects=True)
    r.raise_for_status()
    trace.append({'stage':'download','status':r.status_code,'url':r.url,'bytes':len(r.content),
                  'content_type':r.headers.get('content-type'),'content_disposition':r.headers.get('content-disposition'),
                  'archive_sha256':sha256_bytes(r.content),'magic_hex':r.content[:32].hex()})

    candidates=[]
    for label,raw in walk_payloads(r.content,'MGDS_Download',trace):
        if len(raw) < 1000 or not looks_text(raw):
            continue
        low=raw[:500].lower()
        if b'<html' in low or b'<!doctype' in low:
            continue
        try:
            coords,meta=parse_catalog(raw)
            row={'label':label,'bytes':len(raw),'sha256':sha256_bytes(raw),'valid_coordinates':int(len(coords)),
                 'delta_from_expected':abs(int(len(coords))-EXPECTED_EVENT_COUNT),'parse':meta}
            trace.append({'stage':'parse_candidate',**row})
            if len(coords) >= MIN_COORDS:
                candidates.append((abs(len(coords)-EXPECTED_EVENT_COUNT),label,raw,coords,meta))
        except Exception as e:
            trace.append({'stage':'parse_candidate','label':label,'bytes':len(raw),'sha256':sha256_bytes(raw),
                          'error':f'{type(e).__name__}: {e}'})

    if not candidates:
        raise RuntimeError('nested TAR/GZIP traversal found no parseable >=1000-coordinate catalog')
    candidates.sort(key=lambda x:(x[0],x[1]))
    _,label,raw,coords,meta=candidates[0]
    trace.append({'stage':'catalog_selected_blindly','selection_rule':'MIN_ABS_COUNT_DELTA_FROM_PUBLISHED_6843_THEN_LABEL',
                  'label':label,'valid_coordinates':int(len(coords)),'sha256':sha256_bytes(raw)})
    return r.content,raw,coords,meta,label,trace,r


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True); ap.add_argument('--status-output',required=True); a=ap.parse_args()
    out=Path(a.output); st=Path(a.status_output); out.parent.mkdir(parents=True,exist_ok=True)
    artifact='JANUS-ECHO-COUSTEAU-EA-TPHASE-BLIND-CLUSTER-RUN-006-2026-08-21-v1.5'; started=utcnow()
    try:
        archive_raw,catalog_raw,coords,parse_meta,label,trace,response=acquire()

        # Frozen blind phase: no anchor argument enters clustering.
        blind=blind_cluster(coords)
        blind_hash=blind['freeze_sha256']

        # Only after blind freeze/hash do we reveal the frozen anchor.
        reveal=reveal_and_score(coords,blind,ANCHOR_LAT,ANCHOR_LON)
        nearest=[x['nearest_cluster']['anchor_to_center_km'] for x in reveal['configs'] if x.get('nearest_cluster')]
        p95=any(x.get('nearest_cluster') and x['nearest_cluster']['anchor_inside_cluster_p95_radius'] for x in reveal['configs'])
        maxr=any(x.get('nearest_cluster') and x['nearest_cluster']['anchor_inside_cluster_max_radius'] for x in reveal['configs'])
        verdict='ANCHOR_OVERLAPS_BLIND_CLUSTER_P95__TECTONIC_CONTROL_REQUIRED' if p95 else 'NO_BLIND_CLUSTER_P95_OVERLAP_WITH_FROZEN_ANCHOR'

        result={
          'artifact_id':artifact,'research_branch':'Janus-Echo-Кусто','started_at_utc':started,'completed_at_utc':utcnow(),
          'source':{
            'doi':DOI,'dataset':'EA_Hydroacoustics','data_set_uid':int(DATASET_UID),'file_uid':FILE_UID,
            'archive_filename':'MGDS_Download.tar','archive_sha256':sha256_bytes(archive_raw),'archive_bytes':len(archive_raw),
            'catalog_member_label':label,'catalog_member_sha256':sha256_bytes(catalog_raw),'catalog_member_bytes':len(catalog_raw),
            'raw_committed':False,'expected_event_count':EXPECTED_EVENT_COUNT,'parsed_valid_coordinate_count':int(len(coords)),
            'expected_count_exact_match':int(len(coords))==EXPECTED_EVENT_COUNT,'parse':parse_meta,
            'download':{'url':response.url,'content_type':response.headers.get('content-type'),'content_disposition':response.headers.get('content-disposition')},
            'acquisition_trace':trace,
          },
          'preregistration':{
            'anchor_hidden_during_clustering':True,'clustering_parameters_frozen_before_anchor_reveal':True,
            'dbscan_grid':DBSCAN_GRID,'null_domain':NULL_DOMAIN,'null_samples':NULL_SAMPLES,'null_seed':NULL_SEED,
            'parameter_change_from_runs_001_005':False,
          },
          'blind_phase':blind,'post_reveal':reveal,
          'summary':{
            'blind_freeze_sha256':blind_hash,'parsed_count':int(len(coords)),'expected_count_exact_match':int(len(coords))==EXPECTED_EVENT_COUNT,
            'nearest_catalog_event_to_anchor_km':reveal['nearest_event']['distance_km'],
            'nearest_blind_cluster_center_across_grid_km':round(min(nearest),3) if nearest else None,
            'anchor_inside_any_blind_cluster_p95_radius':p95,'anchor_inside_any_blind_cluster_max_radius':maxr,
            'verdict':verdict,'semantic_status':'UNCONFIRMED','tectonic_control_required':True,
          },
          'hard_rules':[
            'AUTHORITATIVE_DATASET_AND_FILE_UID_PREBOUND_BEFORE_CATALOG_READ','NESTED_ARCHIVE_UNPACK_IS_TRANSPORT_ONLY',
            'BLIND_CLUSTER_BEFORE_ANCHOR_REVEAL','NO_PARAMETER_RETUNING_AFTER_REVEAL',
            'MID_ATLANTIC_RIDGE_SEISMICITY_IS_MANDATORY_TECTONIC_CONTROL','RECTANGULAR_LOOK_ELSEWHERE_NULL_IS_DIAGNOSTIC_NOT_FORMAL',
            'DISTANCE_IS_NOT_CAUSATION','NO_RECENTERING','NO_UNDERWATER_PYRAMID_DETECTED_YET'],
          'status':'BLIND_CLUSTER_RUN_COMPLETE'}
        out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
        st.write_text(json.dumps({'artifact_id':artifact,'status':'SUCCESS','completed_at_utc':utcnow(),'parsed_count':int(len(coords)),
                                  'expected_count_exact_match':int(len(coords))==EXPECTED_EVENT_COUNT,'archive_sha256':sha256_bytes(archive_raw),
                                  'catalog_member_sha256':sha256_bytes(catalog_raw),'blind_freeze_sha256':blind_hash,'verdict':verdict,
                                  'result_path':str(out)},indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
        print(json.dumps(result['summary'],indent=2,ensure_ascii=False)); return 0
    except Exception as e:
        payload={'artifact_id':artifact,'status':'BLOCKED_DATA_ACQUISITION_OR_PARSE','started_at_utc':started,'completed_at_utc':utcnow(),
                 'error_type':type(e).__name__,'error':str(e),'scientific_interpretation':'ACQUISITION_OR_PARSE_BLOCKER_ONLY__NOT_A_NEGATIVE_CLUSTER_RESULT'}
        if 'trace' in locals(): payload['acquisition_trace']=trace
        st.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n',encoding='utf-8'); print(json.dumps(payload,indent=2)); return 2

if __name__=='__main__': raise SystemExit(main())
