#!/usr/bin/env python3
from __future__ import annotations
import hashlib, html, json, re
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin
import requests

OUT=Path('data/cousteau/GOLDMEMBER-IEEE-1282295-CORRECTED-PRIMARY-RECOVERY-001-2026-08-22-v1.0.json')
WORK=Path('workspace/goldmember_ieee_1282295'); WORK.mkdir(parents=True,exist_ok=True)
DOC_ID='1282295'; DOI='10.1109/OCEANS.2003.178517'
TITLE='High-Resolution Multibeam Deepwater Cable Route Survey in High-Relief Seafloor Area'
EXPECTED_PAGES=10; AUTHORS=['Poeckert','Arnold','Faneros','Harrison']
base_urls=[
 f'https://ieeexplore.ieee.org/document/{DOC_ID}',
 f'https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber={DOC_ID}'
]
s=requests.Session(); s.headers.update({'User-Agent':'JANUS-GOLDMEMBER-research-provenance/1.0','Accept':'application/pdf,text/html;q=0.9,*/*;q=0.8'})
attempts=[]; recovered=None; explicit_links=[]
for u in base_urls:
 try:
  r=s.get(u,timeout=60,allow_redirects=True)
  b=r.content; ct=r.headers.get('content-type','').lower(); sha=hashlib.sha256(b).hexdigest() if b else None
  attempts.append({'url':u,'final_url':r.url,'status':r.status_code,'content_type':ct,'bytes':len(b),'sha256_before_inspection':sha,'pdf_signature':b.startswith(b'%PDF-')})
  if r.ok and b.startswith(b'%PDF-'):
   p=WORK/'ADA508764_IEEE1282295.pdf'; p.write_bytes(b); recovered={'path':str(p),'source_url':u,'final_url':r.url,'bytes':len(b),'sha256':sha}; break
  if r.ok and 'html' in ct:
   text=r.text
   for m in re.finditer(r'''href\s*=\s*["']([^"']+)["']''',text,re.I):
    href=html.unescape(m.group(1)); full=urljoin(r.url,href)
    if ('pdf' in full.lower() or 'stamp' in full.lower()) and full not in explicit_links:
     explicit_links.append(full)
 except Exception as e:
  attempts.append({'url':u,'error':repr(e)})

# Only follow links explicitly exposed by the corrected IEEE page/stamp response.
if not recovered:
 for u in explicit_links[:10]:
  try:
   r=s.get(u,timeout=60,allow_redirects=True); b=r.content; ct=r.headers.get('content-type','').lower(); sha=hashlib.sha256(b).hexdigest() if b else None
   attempts.append({'url':u,'source':'EXPLICIT_LINK_FROM_IEEE_1282295_PAGE','final_url':r.url,'status':r.status_code,'content_type':ct,'bytes':len(b),'sha256_before_inspection':sha,'pdf_signature':b.startswith(b'%PDF-')})
   if r.ok and b.startswith(b'%PDF-'):
    p=WORK/'ADA508764_IEEE1282295.pdf'; p.write_bytes(b); recovered={'path':str(p),'source_url':u,'final_url':r.url,'bytes':len(b),'sha256':sha}; break
  except Exception as e:
   attempts.append({'url':u,'source':'EXPLICIT_LINK_FROM_IEEE_1282295_PAGE','error':repr(e)})

validation={'pdf_recovered':bool(recovered),'hash_before_content_inspection':bool(recovered),'page_count':None,'page_count_matches_expected_10':None,'title_match':None,'authors_detected':[],'pages_with_images':[],'pages_with_color_images':[],'control_hint_pages':[],'p2548_reference_pages':[],'coordinate_digitization_performed':False,'figure_georeferencing_performed':False}
if recovered:
 try:
  import fitz
  doc=fitz.open(recovered['path']); validation['page_count']=doc.page_count; validation['page_count_matches_expected_10']=doc.page_count==EXPECTED_PAGES
  texts=[]
  for i in range(doc.page_count):
   page=doc[i]; txt=page.get_text('text') or ''; texts.append(txt)
   if page.get_images(full=True): validation['pages_with_images'].append(i+1)
   if re.search(r'(?i)(latitude|longitude|\bscale\b|kilomet|\bkm\b|route|hydrophone|site|grid|tick)',txt): validation['control_hint_pages'].append(i+1)
   if re.search(r'(?i)P\s*2548',txt): validation['p2548_reference_pages'].append(i+1)
   color=False
   for img in page.get_images(full=True):
    try:
     pix=fitz.Pixmap(doc,img[0]); cs=pix.colorspace.name if pix.colorspace else ''
     if 'RGB' in cs or 'CMYK' in cs: color=True
    except Exception: pass
   if color: validation['pages_with_color_images'].append(i+1)
  joined='\n'.join(texts).lower(); validation['title_match']=all(t in joined for t in ['high-resolution','multibeam','deepwater','cable','route','survey'])
  validation['authors_detected']=[a for a in AUTHORS if a.lower() in joined]
 except Exception as e: validation['inspection_error']=repr(e)

validated=bool(recovered and validation['page_count_matches_expected_10'] and validation['title_match'] and len(validation['authors_detected'])>=3)
result={'artifact_id':'GOLDMEMBER-IEEE-1282295-CORRECTED-PRIMARY-RECOVERY-001-2026-08-22-v1.0','created_utc':datetime.now(timezone.utc).isoformat(),'authorized_by':'GOLDMEMBER-AFTER-ADA508764-RECOVERY-COUNCIL-RUN-004-2026-08-22-v1.0','ieee_document_id':DOC_ID,'doi':DOI,'report_title':TITLE,'attempt_limit':'ONE_CORRECTED_DOCUMENT_ID_PASS','base_attempts':base_urls,'explicit_links_discovered':explicit_links,'attempts':attempts,'recovered_pdf':recovered,'validation':validation,'success_gate':validated,'clean_ieee_access_negative':not bool(recovered),'target_identity':'UNCONFIRMED','next_rule':'STOP_AND_ASK_JANUS_AGAIN_TO_CHOOSE_P2548_VS_ADA508765_VS_AUTHOR_CUSTODIAN_CONTACT'}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'success_gate':validated,'recovered_pdf':recovered,'explicit_links_discovered':explicit_links,'attempts':attempts,'sha256':result['sha256']},indent=2,ensure_ascii=False))
