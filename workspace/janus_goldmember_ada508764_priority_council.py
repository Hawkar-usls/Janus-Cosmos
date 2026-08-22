#!/usr/bin/env python3
from __future__ import annotations
import json, hashlib
from datetime import datetime, timezone
from pathlib import Path
OUT=Path('data/cousteau/GOLDMEMBER-ADA508764-PRIORITY-COUNCIL-RUN-003-2026-08-22-v1.0.json')
result={
 'artifact_id':'GOLDMEMBER-ADA508764-PRIORITY-COUNCIL-RUN-003-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'question':'A new primary locator AD-A508764 has been identified: the 10-page color-illustrated OCEANS 2003 paper on the 2002 high-resolution multibeam cable-route survey. Should GOLDMEMBER recover this primary color copy before continuing ADA508765/P2548, and what may be extracted from it?',
 'evidence_snapshot':{
   'ada508764':{
      'title':'High-Resolution Multibeam Deepwater Cable Route Survey in High-Relief Seafloor Area',
      'pages':10,
      'doi':'10.1109/OCEANS.2003.178517',
      'dtic_id':'AD-A508764',
      'catalog_note':'Original contains color illustrations',
      'survey_scope':['multibeam','backscatter','sidescan','sub-bottom','cable routes','hydrophone sites'],
      'pdf_recovered':False
   },
   'ada508765_pdf_recovered':False,
   'p2548_recovered':False,
   'wishbone_georeferenced':False,
   'central_crags_georeferenced':False,
   'bodc_custodian_request_sent':True,
   'target_identity':'UNCONFIRMED'
 },
 'council':{
   'HRain':'PREFER_THE_PRIMARY_DOCUMENT_CLOSEST_TO_THE_MEASUREMENT_PROCESS',
   'iNaiHR':'USE_DOCUMENT_RELATIONSHIPS_ONLY_TO_FIND_PRIMARY_EVIDENCE',
   'DemiHead':'RECOVERY_DOES_NOT_AUTHORIZE_DIGITIZATION',
   'Fast_CAT':'ADA508764_PRECEDES_ADA508765_FOR_ROUTE_AND_SITE_GEOREFERENCE_BECAUSE_IT_DESCRIBES_THE_MEASUREMENT_SURVEY_DIRECTLY',
   'Aura':'ZERO_EVIDENCE_AUTHORITY',
   'Janus_Cosmos':'NO_BACKSOLVING_FROM_FROZEN_H10S_COORDINATES',
   'Cousteau':'RECOVER_COLOR_COPY_HASH_FIRST_THEN_INSPECT_FOR_ROUTE_SITE_COORDINATE_CONTROL_AND_SOURCE_REFERENCES',
   'Fundamentum':'VOLCANIC_BASELINE_REMAINS_PRIMARY',
   'AIFC':'RECORD_PDF_HASH_PAGE_COUNT_AUTHORS_DOI_AND_EVERY_FIGURE_PAGE_BEFORE_ANY_DERIVED_GEOMETRY',
   'Voice_of_Janus':'TAKE_THE_SURVEYORS_MAP_BEFORE_THE_INTERPRETERS_MAP'
 },
 'janus_answer':{
   'AUTHORIZED_STAGE':'GOLDMEMBER_ADA508764_PRIMARY_COLOR_COPY_RECOVERY_001',
   'PRIORITY':'ADA508764_FIRST__THEN_REASK_BEFORE_ADA508765_OR_P2548_NEXT_MOVE',
   'ORDER':[
      'TRY_OFFICIAL_DTIC_HANDLE_AND_PUBLIC_DTIC_DOWNLOAD_PATHS_FOR_ADA508764',
      'TRY_IEEE_OR_INSTITUTIONAL_OR_AUTHOR_PUBLIC_COPY_USING_EXACT_DOI_TITLE_AND_REPORT_ID',
      'HASH_ANY_RECOVERED_PDF_BEFORE_FIGURE_INSPECTION',
      'VERIFY_TITLE_AUTHORS_EXPECTED_10_PAGES_AND_COLOR_FIGURE_PRESENCE',
      'RECORD_WHETHER_FIGURES_HAVE_VISIBLE_LAT_LON_GRID_TICKS_SCALE_BARS_ROUTE_LABELS_OR_HYDROPHONE_SITE_MARKERS_WITHOUT_DIGITIZING_THEM',
      'RECORD_ANY_EXPLICIT_REFERENCE_TO_P2548_FIGURE_OR_SURVEY_PRODUCT_IDS'
   ],
   'NOT_APPROVED':[
      'NO_COORDINATE_DIGITIZATION',
      'NO_WISHBONE_OR_CRAG_POSITION_ESTIMATION',
      'NO_SCREENSHOT_BACKSOLVING',
      'NO_SHAPE_SCORING',
      'NO_PYRAMID_OR_CRATER_LABEL'
   ],
   'SUCCESS_GATE':'PRIMARY_COLOR_COPY_RECOVERED_AND_VALIDATED_OR_CLEAN_PRIMARY_ACCESS_FAILURE_FROZEN',
   'ON_SUCCESS':'STOP_AND_ASK_JANUS_AGAIN_BEFORE_ANY_QUANTITATIVE_USE_OR_NEXT_ARCHIVE_BRANCH',
   'TARGET_IDENTITY':'UNCONFIRMED'
 },
 'hard_rules':['ASK_JANUS_BEFORE_EACH_NEW_ARCHIVE_OR_PHYSICAL_STAGE','PRIMARY_SOURCE_FIRST','HASH_BEFORE_INSPECTION','NO_QUANTITATIVE_FIGURE_USE_WITHOUT_NEW_COUNCIL','VOLCANIC_BASELINE_REMAINS_ACTIVE','NEGATIVE_RESULTS_REMAIN_NEGATIVE']
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result,indent=2,ensure_ascii=False))
