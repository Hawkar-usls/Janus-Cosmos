#!/usr/bin/env python3
from __future__ import annotations
import json, hashlib
from datetime import datetime, timezone
from pathlib import Path
OUT=Path('data/cousteau/GOLDMEMBER-AFTER-ADA508764-RECOVERY-COUNCIL-RUN-004-2026-08-22-v1.0.json')
result={
 'artifact_id':'GOLDMEMBER-AFTER-ADA508764-RECOVERY-COUNCIL-RUN-004-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'question':'ADA508764 primary recovery returned no PDF. DTIC direct routes are access-failed, but DOI resolution revealed the exact IEEE document identifier 1282295; the previous IEEE stamp probes incorrectly used 178517 from the DOI suffix as an arnumber. Should GOLDMEMBER make a corrected IEEE-1282295 retrieval attempt before moving to P2548/ADA508765?',
 'evidence_snapshot':{
   'ada508764_pdf_recovered':False,
   'dtic_statuses':['403','403','404','handle_500'],
   'doi':'10.1109/OCEANS.2003.178517',
   'doi_resolved_ieee_document_id':'1282295',
   'previous_wrong_ieee_probe_document_id':'178517',
   'figure_inspection_performed':False,
   'coordinate_digitization_performed':False,
   'p2548_recovered':False,
   'ada508765_recovered':False,
   'target_identity':'UNCONFIRMED'
 },
 'council':{
   'HRain':'CORRECT_A_DISCOVERED_IDENTIFIER_BEFORE_ABANDONING_THE_PRESELECTED_PRIMARY_SOURCE',
   'iNaiHR':'DOCUMENT_ID_CORRECTION_HAS_ZERO_IDENTITY_WEIGHT',
   'DemiHead':'WRONG_ARNUMBER_FAILURE_DOES_NOT_COUNT_AS_IEEE_ACCESS_NEGATIVE',
   'Fast_CAT':'ONE_CORRECTED_IEEE_1282295_PASS_IS_ALLOWED_THEN_FREEZE_RESULT_NO_FURTHER_URL_GUESSING',
   'Aura':'ZERO_EVIDENCE_AUTHORITY',
   'Janus_Cosmos':'NO_CHANGE_TO_H10S_COORDINATES_OR_TARGETS',
   'Cousteau':'TRY_IEEE_METADATA_STAMP_AND_ANY_EXPLICIT_PUBLIC_PDF_LINK_FOR_1282295_ONLY',
   'Fundamentum':'VOLCANIC_BASELINE_REMAINS_PRIMARY',
   'AIFC':'HASH_BYTES_BEFORE_INSPECTION_AND_RECORD_REDIRECT_CHAIN',
   'Voice_of_Janus':'USE_THE_RIGHT_KEY_ONCE__IF_THE_DOOR_STAYS_CLOSED_MOVE_ON'
 },
 'janus_answer':{
   'AUTHORIZED_STAGE':'GOLDMEMBER_IEEE_1282295_CORRECTED_PRIMARY_RECOVERY_001',
   'ORDER':[
      'PROBE_IEEE_DOCUMENT_1282295_METADATA_AND_STAMP_ROUTES',
      'FOLLOW_ONLY_EXPLICIT_REDIRECTS_OR_PUBLIC_PDF_LINKS_DISCOVERED_FROM_THE_CORRECT_DOCUMENT_PAGE',
      'IF_PDF_RECOVERED_HASH_BEFORE_INSPECTION_AND_VALIDATE_TITLE_AUTHORS_PAGE_COUNT',
      'IF_NO_PDF_FREEZE_CLEAN_IEEE_ACCESS_NEGATIVE_FOR_1282295'
   ],
   'ATTEMPT_LIMIT':'ONE_CORRECTED_DOCUMENT_ID_PASS__NO_BRUTE_FORCE_OR_URL_ENUMERATION',
   'NOT_APPROVED':['NO_COORDINATE_DIGITIZATION','NO_FIGURE_GEOREFERENCE','NO_SHAPE_SCORING','NO_PYRAMID_OR_CRATER_LABEL'],
   'ON_COMPLETION':'STOP_AND_ASK_JANUS_AGAIN_TO_CHOOSE_P2548_VS_ADA508765_VS_AUTHOR_CUSTODIAN_CONTACT',
   'TARGET_IDENTITY':'UNCONFIRMED'
 },
 'hard_rules':['ASK_JANUS_BEFORE_NEXT_ARCHIVE_STAGE','CORRECT_IDENTIFIER_ONCE_ONLY','HASH_BEFORE_INSPECTION','NO_POSTHOC_URL_ENUMERATION','VOLCANIC_BASELINE_REMAINS_ACTIVE','NEGATIVE_RESULTS_REMAIN_NEGATIVE']
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result,indent=2,ensure_ascii=False))
