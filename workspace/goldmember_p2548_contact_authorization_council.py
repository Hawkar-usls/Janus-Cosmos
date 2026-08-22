#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json
from datetime import datetime, timezone
from pathlib import Path

OUT=Path('data/cousteau/GOLDMEMBER-P2548-CONTACT-AUTHORIZATION-COUNCIL-RUN-007-2026-08-22-v1.0.json')
result={
  'artifact_id':'GOLDMEMBER-P2548-CONTACT-AUTHORIZATION-COUNCIL-RUN-007-2026-08-22-v1.0',
  'created_utc':datetime.now(timezone.utc).isoformat(),
  'question':'Custodian lineage identification found Fugro as the documented corporate successor to Thales Geosolutions and CTBTO as the current HA10/IMS programme organization, while specific P2548 custody remains unconfirmed. May GOLDMEMBER contact one or both institutions, in what order, and with what exact scope?',
  'evidence_snapshot':{
    'fugro_documented_corporate_successor':True,
    'fugro_specific_p2548_custody_confirmed':False,
    'ctbto_current_ha10_program_responsibility':True,
    'ctbto_specific_p2548_custody_confirmed':False,
    'public_contact_routes_exist_for_both':True,
    'p2548_primary_bytes_recovered':False,
    'target_identity':'UNCONFIRMED'
  },
  'council':{
    'HRain':'ASK_FOR_THE_RECORD_NOT_FOR_VALIDATION_OF_A_STORY',
    'iNaiHR':'TWO_INDEPENDENT_CUSTODIAN_LANES_MAY_BE_QUERIED_WITHOUT_MERGING_THEIR_AUTHORITY',
    'DemiHead':'DO_NOT_CLAIM_EITHER_INSTITUTION_HOLDS_P2548__ASK_IF_THEY_DO_OR_CAN_ROUTE_US',
    'Fast_CAT':'PARALLEL_CONTACT_IS_ALLOWED_BECAUSE_THE_LANES_ARE_PREDECLARED_AND_NONCOMPETING',
    'Aura':'ZERO_EVIDENCE_AUTHORITY',
    'Janus_Cosmos':'NO_TARGET_COORDINATES_NEEDED_IN_THE_INITIAL_ARCHIVE_REQUEST',
    'Cousteau':'REQUEST_P2548_SCAN_OR_ACCESSION_AND_ASSOCIATED_2002_SURVEY_PRODUCT_IDENTIFIERS_ONLY',
    'Fundamentum':'NO_ARCHAEOLOGICAL_CLAIM_IN_CONTACT',
    'AIFC':'PRESERVE_SENT_MESSAGE_IDS_TIMESTAMPS_RECIPIENTS_AND_ANY_REPLY_AS_EXTERNAL_CALIBRATION',
    'Voice_of_Janus':'ASK_TWO_KEEPERS_THE_SAME_SIMPLE_QUESTION__DO_YOU_HAVE_THE_LEDGER_OR_KNOW_WHO_DOES'
  },
  'janus_answer':{
    'CONTACT_AUTHORIZED':True,
    'CONTACT_MODE':'PARALLEL_TWO_LANE_ARCHIVAL_PROVENANCE_REQUEST',
    'LANES':[
      {
        'institution':'Fugro',
        'reason':'documented corporate successor to Thales Geosolutions',
        'request':[ 
          'ASK_IF_FUGRO_ARCHIVES_HOLD_THALES_GEOSOLUTIONS_PACIFIC_REPORT_P2548_SURVEY_REPORT_FOR_HDAS_INSTALLATION_AT_ASCENSION_ISLAND_AUGUST_2002',
          'IF_YES_REQUEST_PUBLIC_OR_SHAREABLE_SCAN_PDF_OR_ARCHIVE_ACCESSION_REFERENCE',
          'ASK_FOR_IDENTIFIERS_OR_HOLDING_LOCATION_OF_ASSOCIATED_2002_MULTIBEAM_BACKSCATTER_SIDESCAN_SUBBOTTOM_PRODUCTS_IF_RETAINED',
          'IF_NO_ASK_WHETHER_ARCHIVE_CUSTODY_TRANSFERRED_AND_TO_WHOM'
        ]
      },
      {
        'institution':'CTBTO Preparatory Commission / IMS',
        'reason':'current HA10 IMS programme organization',
        'request':[
          'ASK_IF_TECHNICAL_OR_HISTORICAL_HA10_FILES_INCLUDE_OR_REFERENCE_THALES_GEOSOLUTIONS_REPORT_P2548_FROM_AUGUST_2002',
          'IF_YES_REQUEST_PUBLIC_OR_SHAREABLE_COPY_OR_ACCESSION_REFERENCE',
          'ASK_IF_ASSOCIATED_PREINSTALLATION_2002_SURVEY_PRODUCT_IDENTIFIERS_CAN_BE_SHARED',
          'IF_NOT_HELD_ASK_FOR_THE_APPROPRIATE_HISTORICAL_CUSTODIAN_ROUTE'
        ]
      }
    ],
    'CONTACT_STYLE':'HISTORICAL_ARCHIVE_AND_REPRODUCIBILITY_REQUEST_ONLY',
    'DO_NOT_INCLUDE':[
      'PYRAMID',
      'OSIRIS',
      'ANOMALY_CLAIM',
      '117_121HZ',
      'COSMOS_ALIGNMENT',
      'REQUEST_FOR_RESTRICTED_CURRENT_IMS_SECURITY_INFORMATION'
    ],
    'ON_SEND':'FREEZE_OUTBOUND_RECEIPTS_AND_CONTINUE_ONLY_ALREADY_AUTHORIZED_PUBLIC_ARCHIVE_WORK__DO_NOT_INFER_FROM_SILENCE',
    'ON_REPLY':'FREEZE_REPLY_AND_ASK_JANUS_BEFORE_OPENING_OR_QUANTITATIVELY_USING_ANY_NEW_PRIMARY_MATERIAL',
    'TARGET_IDENTITY':'UNCONFIRMED'
  },
  'hard_rules':[
    'ASK_JANUS_BEFORE_EACH_NEW_ARCHIVE_OR_PHYSICAL_STAGE',
    'NO_STORY_IN_ARCHIVE_REQUEST',
    'NO_RESTRICTED_SECURITY_DATA_REQUEST',
    'REPLY_IS_CALIBRATION_NOT_AUTHORITY_BY_ITSELF',
    'VOLCANIC_BASELINE_REMAINS_ACTIVE'
  ]
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result,indent=2,ensure_ascii=False))
