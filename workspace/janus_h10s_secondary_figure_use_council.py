#!/usr/bin/env python3
from __future__ import annotations
import json, hashlib
from datetime import datetime, timezone
from pathlib import Path

DISC=Path('data/cousteau/JANUS-H10S-P2548-COLOR-FIGURE-ARCHIVE-DISCOVERY-001-2026-08-22-v1.0.json')
OUT=Path('data/cousteau/JANUS-H10S-SECONDARY-FIGURE-USE-COUNCIL-RUN-008-2026-08-22-v1.0.json')
d=json.loads(DISC.read_text(encoding='utf-8'))

result={
 'artifact_id':'JANUS-H10S-SECONDARY-FIGURE-USE-COUNCIL-RUN-008-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'question':'DTIC confirms a color-copy archive lead (AD-A508765), and Irving 2015 reproduces Southern Seamount Fig. 19. May secondary/reproduced figures be used for spatial extraction before the original color PDF/P2548 is recovered?',
 'evidence_snapshot':{
   'dtic_color_copy_catalogued':True,
   'dtic_pdf_bytes_recovered':False,
   'irving_fig19_secondary_reproduction_identified':True,
   'wishbone_exact_coordinates_known':False,
   'central_crags_exact_coordinates_known':False,
   'primary_2002_georeferenced_raster_recovered':False,
   'target_identity':'UNCONFIRMED'
 },
 'council':{
   'HRain':'USE_SECONDARY_FIGURES_TO_FIND_PRIMARY_SOURCES_NOT_TO_CREATE_NEW_PRIMARY_COORDINATES',
   'iNaiHR':'MAY_SUGGEST_ARCHIVE_CONNECTIONS_ONLY',
   'DemiHead':'NO_COORDINATE_DIGITIZATION_WITHOUT_VISIBLE_GEODETIC_CONTROL_AND_PRIMARY_OR_AUTHORITATIVE_SOURCE',
   'Fast_CAT':'SECONDARY_IMAGE_CAN_TEST_GROSS_LAYOUT_ONLY_IF_PREDECLARED__NOT_FINE_FEATURE_POSITION',
   'Aura':'ZERO_EVIDENCE_AUTHORITY',
   'Janus_Cosmos':'DO_NOT_BACKSOLVE_TARGET_COORDINATES_FROM_EXPECTED_H10S_POSITION',
   'Cousteau':'RECOVER_ADA508765_COLOR_PDF_OR_P2548_FIRST',
   'Fundamentum':'VOLCANIC_BASELINE_REMAINS_PRIMARY',
   'AIFC':'PRIMARY_SOURCE_PROVENANCE_REQUIRED_FOR_WISHBONE_OR_CRAG_POLYGON',
   'Voice_of_Janus':'THE_COPY_CAN_SHOW_WHAT_TO_LOOK_FOR__THE_ORIGINAL_MUST_TELL_YOU_WHERE_IT_IS'
 },
 'janus_answer':{
   'SECONDARY_FIGURE_COORDINATE_DIGITIZATION':'NOT_APPROVED',
   'SECONDARY_FIGURE_ALLOWED_USE':[ 
      'ARCHIVE_DISCOVERY_AND_FIGURE_IDENTITY_CONFIRMATION',
      'GROSS_NONQUANTITATIVE_LAYOUT_CONTEXT_WITHOUT_COORDINATES',
      'CROSSCHECK_THAT_A_LATER_PRIMARY_FIGURE_IS_THE_SAME_SCENE'
   ],
   'AUTHORIZED_STAGE':'ADA508765_PRIMARY_COLOR_COPY_RECOVERY_001',
   'ORDER':[ 
      'TRY_OFFICIAL_DTIC_HANDLE_AND_DTIC_PUBLIC_DOWNLOAD_PATHS_FOR_ADA508765',
      'TRY_INSTITUTIONAL_OR_AUTHOR_PUBLIC_COPIES_USING_EXACT_REPORT_ID_DOI_AND_TITLE',
      'HASH_ANY_RECOVERED_PDF_BEFORE_OPENING_FIGURES',
      'VERIFY_PAGE_COUNT_TITLE_AUTHORS_AND_COLOR_FIGURE_PRESENCE',
      'IF_VISIBLE_COORDINATE_GRID_OR_TICKS_EXIST__STOP_AND_ASK_JANUS_BEFORE_DIGITIZATION'
   ],
   'ON_FAILURE':'FREEZE_PRIMARY_COPY_ACCESS_FAILURE_AND_ASK_JANUS_BEFORE_CONTACTING_AUTHORS_OR_USING_SECONDARY_IMAGE_QUANTITATIVELY',
   'TARGET_IDENTITY':'UNCONFIRMED'
 },
 'hard_rules':['ASK_JANUS_BEFORE_NEXT_PHYSICAL_OR_ARCHIVE_STAGE','NO_SCREENSHOT_COORDINATE_ESTIMATION','NO_BACKSOLVING_FROM_H10S_COORDINATES','PRIMARY_SOURCE_REQUIRED_FOR_FINE_GEOREFERENCE','NEGATIVE_RESULTS_REMAIN_NEGATIVE']
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result,indent=2,ensure_ascii=False))
