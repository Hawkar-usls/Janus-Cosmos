#!/usr/bin/env python3
import json, hashlib
from datetime import datetime, timezone
from pathlib import Path
OUT=Path('data/cousteau/JANUS-HA10-DTIC-KEY-COUNCIL-RUN-001-2026-08-22-v1.0.json')
result={
 'artifact_id':'JANUS-HA10-DTIC-KEY-COUNCIL-RUN-001-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'question':'A DTIC archive identifier AD-A508765 has been identified for the 2003 Ascension seamount geomorphology paper, with catalog metadata stating that the original contains color illustrations. Should JANUS prioritize retrieval of the actual original report/figures, and what may be done if retrieved?',
 'evidence_snapshot':{
  'dtic_identifier':'AD-A508765',
  'dtic_handle':'http://hdl.handle.net/100.2/ADA508765',
  'title':'Geomorphology of Two Seamounts Offshore Ascension Island, South Atlantic Ocean',
  'doi':'10.1109/OCEANS.2003.178518',
  'catalog_says_original_contains_color_illustrations':True,
  'actual_report_bytes_recovered':False,
  'p2548_bytes_recovered':False,
  'wishbone_georeferenced':False,
  'target_identity':'UNCONFIRMED'
 },
 'council':{
  'HRain':'PRIORITIZE_ARCHIVE_KEY_CLOSEST_TO_ORIGINAL_SURVEY_GEOMETRY',
  'DemiHead':'CATALOG_METADATA_IS_NOT_FIGURE_EVIDENCE',
  'Fast_CAT':'RECOVER_BYTES_FIRST__FREEZE_HASH__THEN_DEFINE GEOREFERENCE_TEST BEFORE INSPECTION',
  'Cousteau':'COLOR_BATHYMETRY_OR_BACKSCATTER_FIGURE_MAY_SERVE_AS SECONDARY GEOREFERENCE ONLY IF SCALE_GRID_OR_CONTROL_POINTS ARE PRESENT',
  'Fundamentum':'NO_IDENTITY_PROMOTION_FROM FIGURE APPEARANCE',
  'AIFC':'REQUIRE FILE HASH PAGE NUMBER FIGURE NUMBER AND SOURCE CHAIN',
  'Voice_of_Janus':'RECOVER SOURCE__LOCATE FLOOR__THEN MEASURE'
 },
 'janus_answer':{
  'PRIORITY':'AD_A508765_ORIGINAL_RETRIEVAL_FIRST',
  'APPROVED_NEXT_ACTION':'RETRIEVE_ORIGINAL_REPORT_OR_COLOR_FIGURE_BYTES_FROM_DTIC_OR_A_TRACEABLE_MIRROR',
  'ON_RECOVERY_DO_FIRST':[ 
   'FREEZE_ORIGINAL_FILE_HASH_AND_SOURCE_URL',
   'IDENTIFY_PAGE_AND_FIGURE_NUMBERS',
   'CHECK_FOR_COORDINATE_GRID_SCALE_BAR_CONTOURS_STATION_MARKERS_OR_OTHER GEOREFERENCE CONTROL',
   'DO_NOT_DIGITIZE_OR_SCORE_SHAPE_YET'
  ],
  'IF_NO_GEOREFERENCE_PRESENT':'KEEP_AS_ILLUSTRATIVE_ONLY_AND_CONTINUE_P2548_NATIVE_SURVEY_RECOVERY',
  'IF_GEOREFERENCE_PRESENT':'STOP__FREEZE_SOURCE__ASK_JANUS_FOR_GEOREFERENCE_AND_EXTRACTION_GATE_BEFORE DIGITIZATION',
  'TARGET_IDENTITY':'UNCONFIRMED'
 },
 'hard_rules':['ASK_JANUS_BEFORE_DIGITIZATION','CATALOG_METADATA_IS_NOT_PHYSICAL_EVIDENCE','HASH_ORIGINAL_BYTES','NO_VISUAL_GUESSING','NO_POSTHOC_RETARGETING','NEGATIVE_RESULTS_REMAIN_NEGATIVE']
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result,indent=2,ensure_ascii=False))
