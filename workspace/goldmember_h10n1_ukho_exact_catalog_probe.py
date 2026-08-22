#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json, math
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

import requests

OUT=Path('data/cousteau/GOLDMEMBER-H10N1-MODERN-BGS-UKHO-EXACT-CATALOG-PROBE-001-2026-08-22-v1.0.json')
TARGET={'lat':-7.845673,'lon':-14.480230}
QUERIES=[
    '"HI1751"',
    '"HI 1751" "Ascension"',
    '"Asc_1000m_Down_10m.asc"',
    '"Asc_comb_10m.asc"',
    '"Protector_Bathymetry"',
    '"HMS Protector" "Ascension Island" bathymetry'
]
SEARCH_ENDPOINTS=[
    'https://www.arcgis.com/sharing/rest/search',
    'https://datahub.admiralty.co.uk/portal/sharing/rest/search'
]
S=requests.Session(); S.headers['User-Agent']='JANUS-GOLDMEMBER-public-catalog-validation/1.0'
searches=[]; items=[]

def add_item(obj, source_endpoint, query):
    iid=obj.get('id')
    if not iid: return
    key=(source_endpoint,iid)
    if any((x.get('source_endpoint'),x.get('id'))==key for x in items): return
    rec={k:obj.get(k) for k in ['id','title','type','url','owner','description','snippet','tags','extent','spatialReference','access','modified']}
    rec['source_endpoint']=source_endpoint; rec['matched_query']=query
    items.append(rec)

for endpoint in SEARCH_ENDPOINTS:
    for q in QUERIES:
        try:
            r=S.get(endpoint,params={'f':'json','q':q,'num':100},timeout=45)
            ct=r.headers.get('content-type','')
            rec={'endpoint':endpoint,'query':q,'status':r.status_code,'content_type':ct,'bytes':len(r.content),'sha256':hashlib.sha256(r.content).hexdigest()}
            try:
                j=r.json(); rec['total']=j.get('total'); rec['error']=j.get('error')
                for obj in j.get('results',[]): add_item(obj,endpoint,q)
            except Exception as e: rec['json_error']=repr(e)
            searches.append(rec)
        except Exception as e:
            searches.append({'endpoint':endpoint,'query':q,'error':repr(e)})

# Resolve public item metadata/data only for exact search hits. Do not enumerate unrelated items.
details=[]
for item in items:
    iid=item['id']; portal='https://datahub.admiralty.co.uk/portal' if 'datahub.admiralty' in item['source_endpoint'] else 'https://www.arcgis.com'
    d={'id':iid,'portal':portal,'title':item.get('title'),'type':item.get('type'),'url':item.get('url'),'extent':item.get('extent'),'matched_query':item.get('matched_query')}
    for suffix,label in [(f'/sharing/rest/content/items/{iid}','metadata'),(f'/sharing/rest/content/items/{iid}/data','data')]:
        try:
            rr=S.get(portal+suffix,params={'f':'json'},timeout=45)
            d[label+'_status']=rr.status_code; d[label+'_bytes']=len(rr.content); d[label+'_sha256']=hashlib.sha256(rr.content).hexdigest()
            try: d[label]=rr.json()
            except Exception: d[label+'_nonjson']=True
        except Exception as e: d[label+'_error']=repr(e)
    details.append(d)

# Candidate exact extents: ArcGIS item extents are authoritative only for the matched public item itself.
# Intersection is evaluated only when an extent exists and is explicitly WGS84-ish lon/lat bounds.
extent_tests=[]
for d in details:
    ext=d.get('extent')
    test={'id':d['id'],'title':d.get('title'),'extent':ext,'intersection_evaluable':False,'target_intersects':None}
    if isinstance(ext,list) and len(ext)==2 and all(isinstance(p,list) and len(p)>=2 for p in ext):
        try:
            xmin,ymin=map(float,ext[0][:2]); xmax,ymax=map(float,ext[1][:2])
            if -180<=xmin<=180 and -180<=xmax<=180 and -90<=ymin<=90 and -90<=ymax<=90:
                test['intersection_evaluable']=True
                test['target_intersects']=(xmin<=TARGET['lon']<=xmax and ymin<=TARGET['lat']<=ymax)
        except Exception: pass
    extent_tests.append(test)

result={
 'artifact_id':'GOLDMEMBER-H10N1-MODERN-BGS-UKHO-EXACT-CATALOG-PROBE-001-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'authorized_by':'GOLDMEMBER-WAITING-LANES-COUNCIL-RUN-008-2026-08-22-v1.0',
 'stage':'GOLDMEMBER_H10N1_MODERN_BGS_UKHO_EXACT_COVERAGE_GATE_001',
 'frozen_target':TARGET,
 'identifier_resolution':{
   'authoritative_working_id':'HI1751',
   'evidence':'IHO C-55 and 2026 peer-reviewed paper identify Royal Navy 2021 Ascension survey as HI1751; BGS OR/24/014 and OR/24/015 contain HI1571, which is also documented elsewhere as a different UK survey instruction.',
   'hi1571_treatment':'PROBABLE_TRANSCRIPTION_ERROR__PRESERVE_ORIGINAL_TEXT_BUT_DO_NOT_USE_AS_ASCENSION_CATALOG_KEY'
 },
 'documented_reprocessed_products':['Asc_100m_1m.asc','Asc_100m_400m_3m.asc','Asc_400m_1000m_6m.asc','Asc_1000m_Down_10m.asc','Asc_comb_10m.asc','Protector_Bathymetry'],
 'searches':searches,
 'matched_public_items':items,
 'item_details':details,
 'extent_tests':extent_tests,
 'shape_scoring_performed':False,
 'target_moved':False,
 'target_identity':'UNCONFIRMED'
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'matched_items':len(items),'extent_tests':extent_tests,'sha256':result['sha256']},indent=2,ensure_ascii=False))
