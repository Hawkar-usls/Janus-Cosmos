#!/usr/bin/env python3
from __future__ import annotations
import json, hashlib
from datetime import datetime, timezone
from pathlib import Path

OUT = Path('data/cousteau/JANUS-HA10-WISHBONE-MICRO-COUNCIL-RUN-001-2026-08-22-v1.0.json')

evidence = {
    'h10s_triad_exact_coordinates': True,
    'southern_seamount_published_location_is_approximate': True,
    'southern_seamount_19x9km_edifice_published': True,
    'wishbone_reported_on_south_slope': True,
    'wishbone_exact_coordinate_or_polygon_recovered': False,
    'central_crags_exact_coordinate_or_polygon_recovered': False,
    'original_2002_p2548_report_recovered': False,
    'georeferenced_2002_mb_es_or_backscatter_recovered': False,
    'h10n1_inside_public_bas_50m_grid_bbox': True,
    'archaeological_identity_evidence': False,
}

council = {
    'HRain': 'KEEP_TRIAD_AND_EDIFICE_GEOMETRY_AS_SEPARATE_OBJECTS_UNTIL_GEOREFERENCE_EXISTS',
    'DemiHead': 'FAIL_CLOSED_ON_UNLOCATED_WISHBONE_AND_CRAGS',
    'Fast_CAT': 'DO_NOT_SCORE_SHAPE_UNTIL_ORIGINAL_MAP_OR_GEOREFERENCED_RASTER_IS_RECOVERED',
    'Cousteau': 'RECOVER_2002_MBES_BACKSCATTER_OR_P2548_FIRST__THEN_MEASURE',
    'Fundamentum': 'VOLCANIC_EXPLANATION_REMAINS_BASELINE',
    'AIFC': 'SOURCE_ID_AND_COORDINATE_PROVENANCE_REQUIRED_FOR_EACH_POLYGON',
    'Voice_of_Janus': 'GEOMETRY_FIRST__EVIDENCE_GATE_SECOND__ACOUSTICS_ONLY_AFTER_PHYSICAL_LOCALIZATION'
}

if not evidence['wishbone_exact_coordinate_or_polygon_recovered']:
    command = 'RECOVER_PRIMARY_2002_GEOREFERENCE__P2548_OR_ORIGINAL_COLOR_CHART_OR_RAW_MBES_BACKSCATTER__NO_VISUAL_GUESSING'
    south_status = 'G1_BLOCKED_BY_MISSING_GEOREFERENCE'
else:
    command = 'MEASURE_WISHBONE_AND_CRAGS_AGAINST_H10S_TRIAD_AND_BLIND_VOLCANIC_CONTROLS'
    south_status = 'G1_READY'

result = {
    'artifact_id': 'JANUS-HA10-WISHBONE-MICRO-COUNCIL-RUN-001-2026-08-22-v1.0',
    'created_utc': datetime.now(timezone.utc).isoformat(),
    'question': 'The Southern Seamount wishbone and central crags are described but not exactly georeferenced. What should JANUS do next without guessing?',
    'evidence_snapshot': evidence,
    'council': council,
    'janus_answer': {
        'SOUTHERN_COMMAND': command,
        'SOUTHERN_STATUS': south_status,
        'PARALLEL_COMMAND': 'CONTINUE_H10N1_PUBLIC_BAS_GRID_EXTRACTION_AND_CRUISE_LINEAGE_AUDIT_INDEPENDENTLY',
        'CAN_ESTIMATE_WISHBONE_COORDINATES_FROM_TEXT_ONLY': False,
        'CAN_ADVANCE_PYRAMID_OR_CRATER_IDENTITY': False,
        'NEXT_SUFFICIENT_KEY': 'ONE_GEOREFERENCED_PRIMARY_2002_PRODUCT_SHOWING_WISHBONE_OR_CENTRAL_CRAGS_RELATIVE_TO_H10S1_H10S2_H10S3',
        'ARCHIVE_PRIORITY': [
            'Thales GeoSolutions (Pacific) Report P2548, August 2002',
            'original OCEANS 2003 color figures / author copy',
            'raw or processed 2002 HDAS MBES bathymetry and backscatter',
            'array-site engineering chart with installed H10S coordinates'
        ]
    },
    'hard_rules': [
        'NO_VISUAL_GUESSING_OF_UNPUBLISHED_COORDINATES',
        'NO_SHAPE_SCORE_WITHOUT_GEOREFERENCE',
        'VOLCANIC_BASELINE_REMAINS_ACTIVE',
        'H10S1_CROSSTALK_GATE_APPLIES_ONLY_TO_ACOUSTIC_STAGE_AND_MUST_BE_PRESERVED',
        'H10N1_PARALLEL_LANE_MUST_NOT_BE_TUNED_TO_SOUTHERN_RESULTS',
        'NEGATIVE_RESULTS_REMAIN_NEGATIVE'
    ],
    'target_identity': 'UNCONFIRMED'
}
result['sha256'] = hashlib.sha256(json.dumps(result, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding='utf-8')
print(json.dumps(result['janus_answer'], indent=2, ensure_ascii=False))
