#!/usr/bin/env python3
import json, hashlib
from datetime import datetime, timezone
from pathlib import Path

OUT=Path('data/cousteau/JANUS-HA10-RECONSTRUCTION-NEXT-COUNCIL-RUN-001-2026-08-22-v1.0.json')

result={
  'artifact_id':'JANUS-HA10-RECONSTRUCTION-NEXT-COUNCIL-RUN-001-2026-08-22-v1.0',
  'created_utc':datetime.now(timezone.utc).isoformat(),
  'question':'Given the frozen HA10 reconstruction status, pending external replies, unrecovered P2548/BAS native grid/as-built survey, what is the next physical stage JANUS authorizes before any further interpretation?',
  'evidence_snapshot':{
    'status_json_frozen':True,
    'p2548_recovered':False,
    'wishbone_georeferenced':False,
    'as_built_anchor_seafloor_depths_recovered':False,
    'bas_native_grid_bytes_recovered':False,
    'bgs_ukho_h10n1_intersection_machine_verified':False,
    'external_specialist_replies_pending':True,
    'target_identity':'UNCONFIRMED'
  },
  'council':{
    'HRain':'PREFER_SOURCE_RECOVERY_OVER_NEW_HYPOTHESIS_GENERATION',
    'iNaiHR':'MAY_EXPAND_ARCHIVE_GRAPH_BUT_CANNOT_REWEIGHT_PHYSICAL_EVIDENCE',
    'DemiHead':'NO_NEW_IDENTITY_CLAIM_WHILE_PRIMARY_SURFACE_IS_MISSING',
    'Fast_CAT':'FREEZE_SEARCH_TERMS_AND_ARCHIVE_TARGETS_BEFORE_ARCHIVE_EXPANSION',
    'Cousteau':'RECOVER_PRIMARY_SURVEY_OR_NATIVE_GRID_BEFORE_FINE_GEOMORPHOLOGY',
    'Fundamentum':'KEEP_VOLCANIC_BASELINE_AND_RECORD_NULL_ARCHIVE_SEARCHES',
    'AIFC':'REQUIRE_SOURCE_URL_ARCHIVE_ID_DATE_AND_HASH_WHEN_BYTES_ARE_RECOVERED',
    'Voice_of_Janus':'LOCATE_FLOOR_THEN_SHAPE_THEN_ACOUSTICS'
  },
  'janus_answer':{
    'APPROVED_NEXT_STAGE':'PRIMARY_SOURCE_RECOVERY_SWEEP',
    'APPROVED_ACTIONS':[
      'SEARCH_PUBLIC_ARCHIVES_FOR_P2548_EXACT_TITLE_REPORT_NUMBER_AND_THALES_GEOSOLUTIONS_LINEAGE',
      'SEARCH_OCEANS_2003_AUTHOR_COPIES_AND_FIGURE_HOSTS_FOR_ORIGINAL_COLOR_MAPS',
      'SEARCH_HA10_FINAL_INSTALLATION_AS_BUILT_ROUTE_SURVEY_AND_ENGINEERING_REPORTS',
      'SEARCH_BAS_PDC_OR_RELATED_PUBLIC_ENDPOINTS_FOR_NATIVE_NETCDF_OR_ASCII_GRID_BYTES_WITHOUT_CHANGING_TARGET_WINDOWS',
      'SEARCH_BGS_MEDIN_OFFSHORE_GEOINDEX_FOR_ASCENSION_HA10_2002_SURVEY_RECORDS_AND_NATIVE_PRODUCT_IDS',
      'FREEZE_EACH_NULL_RESULT_AND_EACH_RECOVERED_SOURCE_BEFORE_ANY_SHAPE_SCORING'
    ],
    'NOT_APPROVED':[
      'DERIVE_WISHBONE_COORDINATES_FROM_TEXT_OR_SCREENSHOT',
      'USE_GEBCO_TO_CONFIRM_FINE_STRUCTURE',
      'MOVE_TARGET_WINDOWS_TO_MATCH_AVAILABLE_SURVEY_COVERAGE',
      'RUN_ACOUSTIC_IDENTITY_TESTS_BEFORE_PHYSICAL_LOCALIZATION'
    ],
    'ESCALATION_RULE':'IF_ONE_PRIMARY_GEOREFERENCED_SURFACE_OR_AS_BUILT_PRODUCT_IS_RECOVERED__STOP_ARCHIVE_EXPANSION__FREEZE_SOURCE__ASK_JANUS_FOR_EXTRACTION_AND_SCORING_GATE',
    'TARGET_IDENTITY':'UNCONFIRMED'
  },
  'hard_rules':[
    'ASK_JANUS_BEFORE_EACH_NEW_PHYSICAL_STAGE',
    'SOURCE_RECOVERY_PRECEDES_FINE_SCORING',
    'NO_POSTHOC_RETARGETING',
    'PRIMARY_VS_PROXY_LANES_MUST_NOT_MERGE',
    'NEGATIVE_RESULTS_REMAIN_NEGATIVE'
  ]
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result,indent=2,ensure_ascii=False))
