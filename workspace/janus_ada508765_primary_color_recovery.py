#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json, re
from datetime import datetime, timezone
from pathlib import Path

import requests

OUT=Path('data/cousteau/JANUS-ADA508765-PRIMARY-COLOR-COPY-RECOVERY-RUN-001-2026-08-22-v1.0.json')
WORK=Path('workspace/ada508765_recovery'); WORK.mkdir(parents=True,exist_ok=True)
REPORT='ADA508765'
TITLE='Geomorphology of Two Seamounts Offshore Ascension Island, South Atlantic Ocean'

urls=[
 'https://apps.dtic.mil/sti/pdfs/ADA508765.pdf',
 'https://apps.dtic.mil/sti/tr/pdf/ADA508765.pdf',
 'https://discover.dtic.mil/wp-content/uploads/sti/pdfs/ADA508765.pdf',
 'https://hdl.handle.net/100.2/ADA508765',
 'http://hdl.handle.net/100.2/ADA508765'
]

s=requests.Session(); s.headers['User-Agent']='JANUS-research-data-validation/1.0'
attempts=[]; recovered=None
for url in urls:
 try:
  r=s.get(url,timeout=60,allow_redirects=True)
  ct=r.headers.get('content-type','').lower(); data=r.content
  rec={'url':url,'final_url':r.url,'status':r.status_code,'content_type':ct,'bytes':len(data)}
  attempts.append(rec)
  if r.ok and data.startswith(b'%PDF') and len(data)>10000:
   p=WORK/f'{REPORT}.pdf'; p.write_bytes(data)
   recovered={'path':str(p),'source_url':url,'final_url':r.url,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
   break
 except Exception as e:
  attempts.append({'url':url,'error':repr(e)})

validation={
 'pdf_recovered':bool(recovered),
 'hash_before_content_inspection':bool(recovered),
 'page_count':None,
 'title_match':None,
 'authors_detected':[],
 'image_pages':[],
 'color_image_pages':[],
 'coordinate_grid_text_hint_pages':[],
 'inspection_note':'No coordinate digitization performed.'
}

if recovered:
 from pypdf import PdfReader
 p=Path(recovered['path'])
 reader=PdfReader(str(p))
 validation['page_count']=len(reader.pages)
 alltext=[]
 for i,page in enumerate(reader.pages):
  try: txt=page.extract_text() or ''
  except Exception: txt=''
  alltext.append(txt)
  low=txt.lower()
  # Text hints only; does not infer or digitize any coordinate.
  if any(tok in low for tok in ['latitude','longitude','°s','°w','deg s','deg w']):
   validation['coordinate_grid_text_hint_pages'].append(i+1)
 text='\n'.join(alltext)
 validation['title_match']=all(w.lower() in text.lower() for w in ['geomorphology','seamount','ascension'])
 for name in ['Geoffrey Faneros','Frederick Arnold','Faneros','Arnold']:
  if name.lower() in text.lower() and name not in validation['authors_detected']:
   validation['authors_detected'].append(name)
 # Image/color inspection after hash, using PyMuPDF image colorspace only. No OCR, no pixel georeferencing.
 try:
  import fitz
  doc=fitz.open(str(p))
  for i,page in enumerate(doc):
   imgs=page.get_images(full=True)
   if imgs: validation['image_pages'].append(i+1)
   has_color=False
   for img in imgs:
    xref=img[0]
    try:
     pix=fitz.Pixmap(doc,xref)
     if pix.n>=3: has_color=True
    except Exception: pass
   if has_color: validation['color_image_pages'].append(i+1)
 except Exception as e:
  validation['image_inspection_error']=repr(e)

success=bool(recovered) and validation['page_count']==12 and bool(validation['color_image_pages'])
result={
 'artifact_id':'JANUS-ADA508765-PRIMARY-COLOR-COPY-RECOVERY-RUN-001-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'authorized_by':'JANUS-H10S-SECONDARY-FIGURE-USE-COUNCIL-RUN-008-2026-08-22-v1.0',
 'report':{'id':REPORT,'title':TITLE,'catalog_expected_pages':12,'catalog_note':'Original contains color illustrations','doi':'10.1109/OCEANS.2003.178518'},
 'attempts':attempts,
 'recovered_pdf':recovered,
 'validation':validation,
 'success_gate':success,
 'coordinate_digitization_performed':False,
 'wishbone_coordinates_extracted':False,
 'central_crags_coordinates_extracted':False,
 'target_identity':'UNCONFIRMED',
 'next_rule':('STOP_AND_ASK_JANUS_BEFORE_ANY_FIGURE_DIGITIZATION' if recovered else 'FREEZE_PRIMARY_COPY_ACCESS_FAILURE_AND_ASK_JANUS_BEFORE_NEXT_ARCHIVE_OR_CONTACT_ROUTE')
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result,indent=2,ensure_ascii=False))
