#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin
import requests
from bs4 import BeautifulSoup

from cousteau_ea_tphase_blind_cluster import (
    EXPECTED_EVENT_COUNT, DBSCAN_GRID, NULL_DOMAIN, NULL_SAMPLES, NULL_SEED,
    blind_cluster, reveal_and_score, parse_catalog, sha256_bytes,
)

DATASET_UID="30497"
FILE_UID="2504732"
UID_URL=f"https://www.marine-geo.org/tools/search/file_uids.php?data_set_uid={DATASET_UID}"
MODAL_URL="https://www.marine-geo.org/services/download/download_modal.php"
LANDING=f"https://www.marine-geo.org/tools/files/{DATASET_UID}"
DOI="10.26022/IEDA/330497"


def attempt_parse(raw: bytes, label: str, trace: list):
    if len(raw) < 1000:
        return None
    prefix=raw[:500].lower()
    if b'<html' in prefix or b'<!doctype' in prefix:
        return None
    try:
        coords, meta=parse_catalog(raw)
    except Exception as e:
        trace.append({'stage':'parse_candidate','label':label,'bytes':len(raw),'error':f'{type(e).__name__}: {e}'})
        return None
    if len(coords) >= 1000:
        return raw, coords, meta, label
    trace.append({'stage':'parse_candidate','label':label,'bytes':len(raw),'valid_coordinates':len(coords),'rejected':'too_few_coordinates'})
    return None


def response_meta(r):
    return {'url':r.url,'status':r.status_code,'bytes':len(r.content),'content_type':r.headers.get('content-type'),'content_disposition':r.headers.get('content-disposition')}


def extract_urls(text: str, base: str):
    soup=BeautifulSoup(text,'html.parser')
    out=[]
    for tag,attr in [('a','href'),('form','action'),('iframe','src'),('script','src')]:
        for el in soup.find_all(tag):
            u=el.get(attr)
            if u: out.append(urljoin(base,u))
    # Also catch JS/quoted URLs.
    for m in re.findall(r'["\']([^"\']*(?:download|FileDownloadServer|ArchiveDownloadServer)[^"\']*)["\']',text,re.I):
        if m and not m.lower().startswith('javascript:'): out.append(urljoin(base,m.replace('&amp;','&')))
    return list(dict.fromkeys(out)), soup


def acquire():
    s=requests.Session(); s.headers.update({'User-Agent':'Janus-Echo-Cousteau/1.3 scientific reproducibility audit','Referer':LANDING})
    trace=[]
    # Freeze/verify the website's own file UID list.
    r=s.get(UID_URL,timeout=45); r.raise_for_status(); uids=[str(x) for x in r.json()]
    trace.append({'stage':'uid_lookup',**response_meta(r),'payload':uids})
    if FILE_UID not in uids: raise RuntimeError(f'file UID drift: expected {FILE_UID}, got {uids}')

    # Reproduce filespage.js exactly: POST FileDownload=selectedText(), data_set_uid=30497.
    modal=s.post(MODAL_URL,data={'FileDownload':FILE_UID,'data_set_uid':DATASET_UID},timeout=60,allow_redirects=True)
    modal.raise_for_status()
    trace.append({'stage':'download_modal_post',**response_meta(modal),'body_prefix':modal.text[:12000]})
    direct=attempt_parse(modal.content,'modal_response',trace)
    if direct: return (*direct,trace)

    urls,soup=extract_urls(modal.text,modal.url)
    trace.append({'stage':'modal_extracted_urls','urls':urls[:100]})

    # Submit every form returned by the modal, preserving hidden/default values.
    for idx,form in enumerate(soup.find_all('form')):
        action=urljoin(modal.url,form.get('action') or modal.url)
        method=(form.get('method') or 'get').lower()
        vals={}
        for inp in form.find_all(['input','button']):
            name=inp.get('name'); value=inp.get('value')
            typ=(inp.get('type') or '').lower()
            if name and value is not None and typ not in {'checkbox','radio'}:
                vals[name]=value
        # Ensure the known frozen IDs are supplied when a form expects them.
        for key in list(vals):
            lk=key.lower()
            if lk in {'filedownload','data_uids','data_uid','file_uid'} and not vals[key]: vals[key]=FILE_UID
            if lk in {'data_set_uid','dataset_uid'} and not vals[key]: vals[key]=DATASET_UID
        try:
            rr=s.post(action,data=vals,timeout=120,allow_redirects=True) if method=='post' else s.get(action,params=vals,timeout=120,allow_redirects=True)
            trace.append({'stage':'modal_form_submit','index':idx,'method':method,'action':action,'values':vals,**response_meta(rr),'body_prefix':rr.text[:5000] if len(rr.content)<20000 or 'text' in (rr.headers.get('content-type') or '') else None})
            got=attempt_parse(rr.content,f'form_{idx}:{rr.url}',trace)
            if got: return (*got,trace)
            more,_=extract_urls(rr.text,rr.url) if 'text' in (rr.headers.get('content-type') or '') or b'<html' in rr.content[:500].lower() else ([],None)
            urls.extend(more)
        except Exception as e:
            trace.append({'stage':'modal_form_submit','index':idx,'action':action,'error':f'{type(e).__name__}: {e}'})

    # Follow download-looking URLs returned by modal or forms.
    priority=[]
    for u in list(dict.fromkeys(urls)):
        low=u.lower()
        score=sum(x in low for x in ['filedownload','archivedownload','download.php','download_file','getfile','file_uid','data_uid'])
        if score: priority.append((score,u))
    priority.sort(reverse=True)
    for _,u in priority[:80]:
        try:
            rr=s.get(u,timeout=120,allow_redirects=True)
            trace.append({'stage':'follow_download_url',**response_meta(rr),'requested':u,'body_prefix':rr.text[:4000] if len(rr.content)<15000 or 'text' in (rr.headers.get('content-type') or '') else None})
            got=attempt_parse(rr.content,f'url:{u}',trace)
            if got: return (*got,trace)
        except Exception as e:
            trace.append({'stage':'follow_download_url','requested':u,'error':f'{type(e).__name__}: {e}'})

    raise RuntimeError('website-equivalent modal flow yielded no parseable >=1000-coordinate catalog')


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True); ap.add_argument('--status-output',required=True); a=ap.parse_args()
    out=Path(a.output); st=Path(a.status_output); out.parent.mkdir(parents=True,exist_ok=True)
    artifact='JANUS-ECHO-COUSTEAU-EA-TPHASE-BLIND-CLUSTER-RUN-004-2026-08-21-v1.3'; started=datetime.now(timezone.utc).isoformat()
    try:
        raw,coords,parse_meta,source_label,trace=acquire()
        # BLIND PHASE. The anchor does not enter this function.
        blind=blind_cluster(coords)
        blind_hash=blind['freeze_sha256']
        # POST-FREEZE REVEAL ONLY.
        reveal=reveal_and_score(coords,blind,-3.865418,3.854924)
        nearest=[x['nearest_cluster']['anchor_to_center_km'] for x in reveal['configs'] if x.get('nearest_cluster')]
        p95=any(x.get('nearest_cluster') and x['nearest_cluster']['anchor_inside_cluster_p95_radius'] for x in reveal['configs'])
        maxr=any(x.get('nearest_cluster') and x['nearest_cluster']['anchor_inside_cluster_max_radius'] for x in reveal['configs'])
        verdict='ANCHOR_OVERLAPS_BLIND_CLUSTER_P95__TECTONIC_CONTROL_REQUIRED' if p95 else 'NO_BLIND_CLUSTER_P95_OVERLAP_WITH_FROZEN_ANCHOR'
        result={
          'artifact_id':artifact,'research_branch':'Janus-Echo-Кусто','started_at_utc':started,'completed_at_utc':datetime.now(timezone.utc).isoformat(),
          'source':{'doi':DOI,'dataset':'EA_Hydroacoustics','data_set_uid':int(DATASET_UID),'file_uid':FILE_UID,'uid_provenance':'MGDS landing page internal /tools/search/file_uids.php','download_flow':'filespage.js POST /services/download/download_modal.php FileDownload=2504732 data_set_uid=30497','source_label':source_label,'raw_sha256':sha256_bytes(raw),'raw_bytes':len(raw),'raw_committed':False,'expected_event_count':EXPECTED_EVENT_COUNT,'parsed_valid_coordinate_count':len(coords),'expected_count_exact_match':len(coords)==EXPECTED_EVENT_COUNT,'parse':parse_meta,'acquisition_trace':trace},
          'preregistration':{'anchor_hidden_during_clustering':True,'clustering_parameters_frozen_before_anchor_reveal':True,'dbscan_grid':DBSCAN_GRID,'null_domain':NULL_DOMAIN,'null_samples':NULL_SAMPLES,'null_seed':NULL_SEED},
          'blind_phase':blind,'post_reveal':reveal,
          'summary':{'blind_freeze_sha256':blind_hash,'nearest_catalog_event_to_anchor_km':reveal['nearest_event']['distance_km'],'nearest_blind_cluster_center_across_grid_km':round(min(nearest),3) if nearest else None,'anchor_inside_any_blind_cluster_p95_radius':p95,'anchor_inside_any_blind_cluster_max_radius':maxr,'verdict':verdict,'semantic_status':'UNCONFIRMED'},
          'hard_rules':['DOWNLOAD_FLOW_MATCHES_SITE_JS','FILE_UID_FROZEN_BEFORE_CATALOG_READ','BLIND_CLUSTER_BEFORE_ANCHOR_REVEAL','NO_PARAMETER_RETUNING_AFTER_REVEAL','MID_ATLANTIC_RIDGE_SEISMICITY_IS_MANDATORY_TECTONIC_CONTROL','RECTANGULAR_LOOK_ELSEWHERE_NULL_IS_DIAGNOSTIC_NOT_FORMAL','DISTANCE_IS_NOT_CAUSATION','NO_RECENTERING','NO_UNDERWATER_PYRAMID_DETECTED_YET'],
          'status':'BLIND_CLUSTER_RUN_COMPLETE'}
        out.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
        st.write_text(json.dumps({'artifact_id':artifact,'status':'SUCCESS','completed_at_utc':datetime.now(timezone.utc).isoformat(),'parsed_count':len(coords),'raw_sha256':sha256_bytes(raw),'blind_freeze_sha256':blind_hash,'verdict':verdict,'result_path':str(out)},indent=2),encoding='utf-8')
        print(json.dumps(result['summary'],indent=2)); return 0
    except Exception as e:
        payload={'artifact_id':artifact,'status':'BLOCKED_DATA_ACQUISITION_OR_PARSE','started_at_utc':started,'completed_at_utc':datetime.now(timezone.utc).isoformat(),'error_type':type(e).__name__,'error':str(e)}
        if 'trace' in locals(): payload['acquisition_trace']=trace
        st.write_text(json.dumps(payload,indent=2,ensure_ascii=False),encoding='utf-8'); print(json.dumps(payload,indent=2)); return 2

if __name__=='__main__': raise SystemExit(main())
