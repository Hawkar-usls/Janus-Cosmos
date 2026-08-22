#!/usr/bin/env python3
from __future__ import annotations
import hashlib, html, json, math, re
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

import requests

SERIES='1762365'
PAGE=f'https://www.bodc.ac.uk/data/documents/series/{SERIES}/'
TARGET_LAT=-7.845673
TARGET_LON=-14.480230
OUT=Path('data/cousteau/JANUS-JR15001-BODC-TRACK-RECOVERY-GATE-001-2026-08-22-v1.0.json')
WORK=Path('workspace/jr15001_bodc_track'); WORK.mkdir(parents=True,exist_ok=True)

session=requests.Session(); session.headers['User-Agent']='JANUS-research-data-validation/1.0'
r=session.get(PAGE,timeout=60,allow_redirects=True)
r.raise_for_status()
page_bytes=r.content
page_sha=hashlib.sha256(page_bytes).hexdigest()
text=r.text
(WORK/'series1762365.html').write_bytes(page_bytes)

hrefs=[]
for m in re.finditer(r'''href\s*=\s*["']([^"']+)["']''',text,re.I):
 u=urljoin(r.url,html.unescape(m.group(1)))
 if u not in hrefs: hrefs.append(u)
forms=[]
for m in re.finditer(r'''<form[^>]+action\s*=\s*["']([^"']+)["']''',text,re.I):
 u=urljoin(r.url,html.unescape(m.group(1)))
 if u not in forms: forms.append(u)

keywords=('download','odv','data','series','1762365')
candidates=[u for u in hrefs if any(k in u.lower() for k in keywords)]
probes=[]; downloaded=[]
for u in candidates:
 try:
  rr=session.get(u,timeout=30,allow_redirects=True,stream=False)
  ct=rr.headers.get('content-type','').lower(); cd=rr.headers.get('content-disposition','')
  rec={'url':u,'final_url':rr.url,'status':rr.status_code,'content_type':ct,'content_disposition':cd,'bytes':len(rr.content)}
  probes.append(rec)
  is_html='text/html' in ct or rr.content.lstrip().lower().startswith(b'<!doctype html') or rr.content.lstrip().lower().startswith(b'<html')
  if rr.ok and not is_html and len(rr.content)>100:
   name=Path(rr.url.split('?')[0]).name or f'download_{len(downloaded)}.bin'
   mm=re.search(r'filename\*?=(?:UTF-8\'\')?["\']?([^"\';]+)',cd,re.I)
   if mm: name=mm.group(1)
   name=re.sub(r'[^A-Za-z0-9._-]+','_',name)
   p=WORK/name; p.write_bytes(rr.content)
   downloaded.append({'url':u,'path':str(p),'filename':name,'bytes':len(rr.content),'sha256':hashlib.sha256(rr.content).hexdigest(),'content_type':ct})
 except Exception as e:
  probes.append({'url':u,'error':repr(e)})

# Metadata facts independently present in the BODC series narrative.
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
 'bathymetry_quality_warning':'periods of large-scale noise exist in final bathymetry channel; flagged by BODC'
}

# Deliberately no guessing parser: public payload must be clearly machine-readable before track computation.
track_computation={'status':'NOT_RUN','reason':'NO_UNAMBIGUOUS_MACHINE_READABLE_ODV_OR_NAV_PAYLOAD_RECOVERED'}
for d in downloaded:
 fn=d['filename'].lower()
 if fn.endswith(('.txt','.csv','.odv','.tsv')):
  track_computation={'status':'PAYLOAD_RECOVERED_REQUIRES_FORMAT_VALIDATION_BEFORE_COORDINATE_PARSE','filename':d['filename']}
  break

result={
 'artifact_id':'JANUS-JR15001-BODC-TRACK-RECOVERY-GATE-001-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'authorized_by':'JANUS-H10N1-LINEAGE-ARCHIVE-BRANCH-COUNCIL-RUN-007-2026-08-22-v1.0',
 'series':{'bodc_reference':SERIES,'metadata_url':PAGE,'metadata_page_sha256':page_sha},
 'frozen_target':{'lat':TARGET_LAT,'lon':TARGET_LON,'radius_km':1},
 'metadata_facts':metadata_facts,
 'discovered_form_actions':forms,
 'candidate_link_count':len(candidates),
 'candidate_links':candidates,
 'download_probes':probes,
 'downloaded_payloads':downloaded,
 'track_computation':track_computation,
 'swath_eligibility':'NOT_EVALUATED_UNTIL_TRACK_RECOVERED',
 'contributor_status':'NOT_CONFIRMED',
 'per_cruise_morphology_computed':False,
 'target_identity':'UNCONFIRMED',
 'next_rule':'IF_TRACK_NOT_RECOVERED_FREEZE_ACCESS_BLOCKER_AND_ASK_JANUS_BEFORE_ALTERNATE_ARCHIVE_OR_CONTACT_ROUTE'
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'candidate_link_count':len(candidates),'downloaded_payloads':downloaded,'track_computation':track_computation,'sha256':result['sha256']},indent=2,ensure_ascii=False))
