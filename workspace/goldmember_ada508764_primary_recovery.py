#!/usr/bin/env python3
from __future__ import annotations
import hashlib, io, json, re
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse
import requests

OUT=Path('data/cousteau/GOLDMEMBER-ADA508764-PRIMARY-COLOR-COPY-RECOVERY-001-2026-08-22-v1.0.json')
WORK=Path('workspace/goldmember_ada508764'); WORK.mkdir(parents=True,exist_ok=True)
REPORT_ID='ADA508764'
TITLE='High-Resolution Multibeam Deepwater Cable Route Survey in High-Relief Seafloor Area'
DOI='10.1109/OCEANS.2003.178517'
AUTHORS=['Poeckert','Arnold','Faneros','Harrison']
EXPECTED_PAGES=10

urls=[
 'https://apps.dtic.mil/sti/pdfs/ADA508764.pdf',
 'https://apps.dtic.mil/sti/tr/pdf/ADA508764.pdf',
 'https://discover.dtic.mil/wp-content/uploads/sti/pdfs/ADA508764.pdf',
 'https://hdl.handle.net/100.2/ADA508764',
 'http://hdl.handle.net/100.2/ADA508764',
 'https://doi.org/10.1109/OCEANS.2003.178517',
 'https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=178517',
 'https://ieeexplore.ieee.org/document/178517'
]
s=requests.Session(); s.headers.update({'User-Agent':'JANUS-GOLDMEMBER-research-provenance/1.0','Accept':'application/pdf,text/html;q=0.9,*/*;q=0.8'})
attempts=[]; recovered=None
for u in urls:
 try:
  r=s.get(u,timeout=60,allow_redirects=True)
  b=r.content
  ct=r.headers.get('content-type','').lower()
  sha=hashlib.sha256(b).hexdigest() if b else None
  is_pdf=b.startswith(b'%PDF-') or 'application/pdf' in ct
  rec={'url':u,'final_url':r.url,'status':r.status_code,'content_type':ct,'bytes':len(b),'sha256_before_inspection':sha,'is_pdf_signature':b.startswith(b'%PDF-')}
  attempts.append(rec)
  if r.ok and is_pdf and b.startswith(b'%PDF-') and len(b)>10000:
   name='ADA508764_primary_candidate.pdf'
   p=WORK/name; p.write_bytes(b)
   recovered={'path':str(p),'source_url':u,'final_url':r.url,'bytes':len(b),'sha256':sha}
   break
 except Exception as e:
  attempts.append({'url':u,'error':repr(e)})

validation={
 'pdf_recovered':bool(recovered),
 'hash_before_content_inspection':bool(recovered),
 'page_count':None,
 'page_count_matches_expected_10':None,
 'title_match':None,
 'authors_detected':[],
 'doi_detected':False,
 'pages_with_images':[],
 'pages_with_color_images':[],
 'figure_control_hint_pages':[],
 'p2548_reference_pages':[],
 'inspection_scope':'PRESENCE_ONLY__NO_COORDINATE_DIGITIZATION_NO_PIXEL_GEOREFERENCE'
}

if recovered:
 try:
  import fitz
  doc=fitz.open(recovered['path'])
  validation['page_count']=doc.page_count
  validation['page_count_matches_expected_10']=(doc.page_count==EXPECTED_PAGES)
  all_text=[]
  hint_re=re.compile(r'(?i)(latitude|longitude|lat\.?|lon\.?|scale|kilomet|\bkm\b|route|hydrophone|site|north|south|grid|tick)')
  for i in range(doc.page_count):
   page=doc.load_page(i)
   txt=page.get_text('text') or ''
   all_text.append(txt)
   if hint_re.search(txt): validation['figure_control_hint_pages'].append(i+1)
   if re.search(r'(?i)P\s*2548',txt): validation['p2548_reference_pages'].append(i+1)
   images=page.get_images(full=True)
   if images:
    validation['pages_with_images'].append(i+1)
   page_color=False
   for img in images:
    try:
     xref=img[0]
     pix=fitz.Pixmap(doc,xref)
     if pix.n>=3 and pix.samples:
      sample=pix.samples[::max(1,len(pix.samples)//6000)]
      # Presence test only: if multi-channel source is not grayscale-like by colorspace metadata.
      cs=(pix.colorspace.name if pix.colorspace else '')
      if 'RGB' in cs or 'CMYK' in cs:
       page_color=True
     pix=None
    except Exception:
     pass
   if page_color: validation['pages_with_color_images'].append(i+1)
  joined='\n'.join(all_text)
  norm=lambda x: re.sub(r'\s+',' ',x).strip().lower()
  title_tokens=['high-resolution','multibeam','deepwater','cable','route','survey','high-relief']
  validation['title_match']=all(t.lower() in joined.lower() for t in title_tokens)
  validation['authors_detected']=[a for a in AUTHORS if re.search(r'(?i)\b'+re.escape(a)+r'\b',joined)]
  validation['doi_detected']=DOI.lower() in joined.lower()
 except Exception as e:
  validation['inspection_error']=repr(e)

success=bool(recovered and validation.get('page_count_matches_expected_10') and validation.get('title_match') and len(validation.get('authors_detected',[]))>=3)
result={
 'artifact_id':'GOLDMEMBER-ADA508764-PRIMARY-COLOR-COPY-RECOVERY-001-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'authorized_by':'GOLDMEMBER-ADA508764-PRIORITY-COUNCIL-RUN-003-2026-08-22-v1.0',
 'report':{'id':REPORT_ID,'title':TITLE,'doi':DOI,'expected_pages':EXPECTED_PAGES,'expected_authors':AUTHORS,'catalog_note':'Original contains color illustrations'},
 'attempts':attempts,
 'recovered_pdf':recovered,
 'validation':validation,
 'success_gate':success,
 'coordinate_digitization_performed':False,
 'figure_pixel_georeferencing_performed':False,
 'wishbone_coordinates_extracted':False,
 'central_crags_coordinates_extracted':False,
 'target_identity':'UNCONFIRMED',
 'next_rule':'STOP_AND_ASK_JANUS_AGAIN_BEFORE_ANY_QUANTITATIVE_FIGURE_USE_OR_NEXT_ARCHIVE_BRANCH'
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'success_gate':success,'recovered_pdf':recovered,'validation':validation,'attempt_count':len(attempts),'sha256':result['sha256']},indent=2,ensure_ascii=False))
