#!/usr/bin/env python3
from __future__ import annotations
import hashlib, html, json, re
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin
import requests

SERIES='1762365'
TARGET_LAT=-7.845673
TARGET_LON=-14.480230
OUT=Path('data/cousteau/JANUS-JR15001-BODC-TRACK-RECOVERY-GATE-001-2026-08-22-v1.0.json')
WORK=Path('workspace/jr15001_bodc_track'); WORK.mkdir(parents=True,exist_ok=True)

session=requests.Session()
session.headers.update({'User-Agent':'JANUS-research-data-validation/1.1 (+public archive provenance check)','Accept':'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8'})

entrypoints=[
 f'https://www.bodc.ac.uk/data/documents/series/{SERIES}/',
 f'http://www.bodc.ac.uk/data/documents/series/{SERIES}/',
 'https://www.bodc.ac.uk/resources/inventories/cruise_inventory/report/15726/',
 'https://www.bodc.ac.uk/data/documents/cruise/15726/'
]
entrypoint_probes=[]
html_pages=[]
for u in entrypoints:
 try:
  r=session.get(u,timeout=60,allow_redirects=True)
  ct=r.headers.get('content-type','').lower()
  rec={'url':u,'final_url':r.url,'status':r.status_code,'content_type':ct,'bytes':len(r.content),'sha256':hashlib.sha256(r.content).hexdigest() if r.content else None}
  entrypoint_probes.append(rec)
  if r.ok and 'html' in ct and len(r.content)>100:
   p=WORK/(re.sub(r'[^A-Za-z0-9._-]+','_',Path(r.url.rstrip('/')).name or 'page')+'.html')
   p.write_bytes(r.content)
   html_pages.append((r.url,r.text,rec['sha256']))
 except Exception as e:
  entrypoint_probes.append({'url':u,'error':repr(e)})

hrefs=[]; forms=[]
for base,text,_ in html_pages:
 for m in re.finditer(r'''href\s*=\s*["']([^"']+)["']''',text,re.I):
  u=urljoin(base,html.unescape(m.group(1)))
  if u not in hrefs: hrefs.append(u)
 for m in re.finditer(r'''<form[^>]+action\s*=\s*["']([^"']+)["']''',text,re.I):
  u=urljoin(base,html.unescape(m.group(1)))
  if u not in forms: forms.append(u)

keywords=('download','odv','1762365','series','data')
candidates=[]
for u in hrefs+forms:
 if any(k in u.lower() for k in keywords) and u not in candidates:
  candidates.append(u)

probes=[]; downloaded=[]
for u in candidates[:80]:
 try:
  rr=session.get(u,timeout=30,allow_redirects=True)
  ct=rr.headers.get('content-type','').lower(); cd=rr.headers.get('content-disposition','')
  rec={'url':u,'final_url':rr.url,'status':rr.status_code,'content_type':ct,'content_disposition':cd,'bytes':len(rr.content)}
  probes.append(rec)
  head=rr.content.lstrip()[:50].lower()
  is_html='text/html' in ct or head.startswith(b'<!doctype html') or head.startswith(b'<html')
  if rr.ok and not is_html and len(rr.content)>100:
   name=Path(rr.url.split('?')[0]).name or f'download_{len(downloaded)}.bin'
   mm=re.search(r'filename\*?=(?:UTF-8\'\')?["\']?([^"\';]+)',cd,re.I)
   if mm: name=mm.group(1)
   name=re.sub(r'[^A-Za-z0-9._-]+','_',name)
   p=WORK/name; p.write_bytes(rr.content)
   downloaded.append({'url':u,'path':str(p),'filename':name,'bytes':len(rr.content),'sha256':hashlib.sha256(rr.content).hexdigest(),'content_type':ct})
 except Exception as e:
  probes.append({'url':u,'error':repr(e)})

metadata_facts={
 'originator_identifier':'JR15001_PRODQXF_NAV',
 'nominal_cycle_interval_s':30,
 'navigation_source':'Seatex Seapath 320+ WGS84',
 'final_bathymetry_channel':'EA600 single-beam MBANZZ01',
 'originator_multibeam_file':'em122.ACO',
 'originator_multibeam_interval_s_approx':1,
 'originator_multibeam_date_range':['2015-09-25T14:45:47Z','2015-11-02T15:43:57Z'],
 'originator_multibeam_depth_channel':'MBANSWCB',
 'originator_multibeam_channel_transferred_to_final_series':False,
 'public_series_declares_odv_download':True,
 'bathymetry_quality_warning':'periods of large-scale noise exist in final bathymetry channel; flagged by BODC'
}

machine_payload=None
for d in downloaded:
 fn=d['filename'].lower(); ct=d.get('content_type','')
 if fn.endswith(('.txt','.csv','.odv','.tsv','.zip')) or any(x in ct for x in ('text/plain','text/csv','application/zip','octet-stream')):
  machine_payload=d; break

if machine_payload:
 track_computation={'status':'PAYLOAD_RECOVERED_REQUIRES_FORMAT_VALIDATION_BEFORE_COORDINATE_PARSE','payload':machine_payload}
 clean_archive_access_negative=False
else:
 track_computation={'status':'NOT_RUN','reason':'NO_UNAMBIGUOUS_MACHINE_READABLE_NAVIGATION_PAYLOAD_RECOVERED'}
 clean_archive_access_negative=True

result={
 'artifact_id':'JANUS-JR15001-BODC-TRACK-RECOVERY-GATE-001-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'authorized_by':['JANUS-H10N1-LINEAGE-ARCHIVE-BRANCH-COUNCIL-RUN-007-2026-08-22-v1.0','GOLDMEMBER-COUNCIL-RUN-001-2026-08-22-v1.0'],
 'repair_scope':'RUNTIME_ONLY__REMOVED_FATAL_RAISE_FOR_STATUS_AND_RECORD_HTTP_FAILURES_AS_ARCHIVE_ACCESS_EVIDENCE__NO_SCIENTIFIC_RULE_CHANGED',
 'frozen_target':{'lat':TARGET_LAT,'lon':TARGET_LON,'radius_km':1},
 'series':{'bodc_reference':SERIES,'metadata_url':entrypoints[0]},
 'entrypoint_probes':entrypoint_probes,
 'metadata_facts':metadata_facts,
 'discovered_form_actions':forms,
 'candidate_link_count':len(candidates),
 'candidate_links':candidates,
 'download_probes':probes,
 'downloaded_payloads':downloaded,
 'track_computation':track_computation,
 'clean_archive_access_negative_receipt':clean_archive_access_negative,
 'centerline_distance_computed':False,
 'swath_eligibility':'NOT_EVALUATED_UNTIL_MACHINE_READABLE_TRACK_RECOVERED',
 'contributor_status':'NOT_CONFIRMED',
 'per_cruise_morphology_computed':False,
 'target_identity':'UNCONFIRMED',
 'next_rule':'STOP_AND_ASK_JANUS_AGAIN_BEFORE_ALTERNATE_ARCHIVE_BRANCH_OR_SWATH_FOOTPRINT_STAGE'
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'clean_archive_access_negative_receipt':clean_archive_access_negative,'entrypoint_probes':entrypoint_probes,'downloaded_payloads':downloaded,'track_computation':track_computation,'sha256':result['sha256']},indent=2,ensure_ascii=False))
