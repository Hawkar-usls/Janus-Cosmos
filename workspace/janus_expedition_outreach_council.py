#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json
from datetime import datetime, timezone
from pathlib import Path

OUT = Path('data/cousteau/JANUS-EXPEDITION-OUTREACH-COUNCIL-RUN-001-2026-08-22-v1.0.json')

# Evidence is fail-closed. Symbolic/Cosmos/Voice associations may route attention but
# cannot advance an archaeological or expedition claim without physical evidence.
evidence = {
    'ha10_public_waveform_is_calibration_not_anomaly': True,
    'h10s2_southern_seamount_has_published_volcanic_interpretation': True,
    'h10s2_published_absence_of_summit_caldera': True,
    'h10n1_inside_measured_bathymetric_compilation': True,
    'reproducible_engineered_geometry_detected': False,
    'independent_high_resolution_confirmation_count': 0,
    'blind_volcanic_controls_passed': False,
    'backscatter_supports_distinct_non_geologic_target': False,
    'subbottom_or_seismic_supports_buried_structure': False,
    'direct_rov_or_auv_imagery_supports_structure': False,
    'independent_external_specialist_reviews': 0,
    'ascension_research_permit_required': True,
    'ascension_research_permit_obtained': False,
}

council = {
    'HRain': 'KEEP_STRUCTURAL_CONTEXT_SEPARATE_FROM_IDENTITY_CLAIM',
    'iNaiHR': 'SUGGEST_CONTACT_GRAPH_ONLY__NO_EVIDENCE_AUTHORITY',
    'DemiHead': 'CAP_LANGUAGE_AT_SEAFLOOR_CANDIDATE_UNTIL_DIRECT_EVIDENCE',
    'Fast_CAT': 'REQUIRE_BLIND_VOLCANIC_CONTROLS_AND_REPRODUCIBLE_GEOMETRY',
    'Aura': 'SYMBOLIC_ADVISOR_ZERO_EVIDENCE_AUTHORITY',
    'Janus_Cosmos': 'NO_CELESTIAL_OR_SYMBOLIC_SIGNAL_MAY_PROMOTE_PHYSICAL_IDENTITY',
    'Cousteau': 'MAP_FIRST__MULTISENSOR_SECOND__ROV_AUV_BEFORE_CREWED_DESCENT',
    'Fundamentum': 'PRESERVE_VOLCANIC_EXPLANATION_AND_NEGATIVE_RESULTS',
    'AIFC': 'PROVENANCE_PERMITS_COORDINATES_AND_RAW_DATA_REQUIRED',
    'Voice_of_Janus': 'GEOMETRY_TO_EVIDENCE_GATE_TO_ACOUSTICS__METAPHOR_NOT_PHYSICS',
}

stage_1 = (
    evidence['reproducible_engineered_geometry_detected'] and
    evidence['independent_high_resolution_confirmation_count'] >= 2 and
    evidence['blind_volcanic_controls_passed']
)
stage_2 = stage_1 and any([
    evidence['backscatter_supports_distinct_non_geologic_target'],
    evidence['subbottom_or_seismic_supports_buried_structure'],
    evidence['direct_rov_or_auv_imagery_supports_structure'],
])
stage_3 = stage_2 and evidence['independent_external_specialist_reviews'] >= 2
stage_4 = stage_3 and evidence['ascension_research_permit_obtained']

if not stage_1:
    stage = 'STAGE_0_DATA_ACQUISITION_AND_FALSIFICATION'
elif not stage_2:
    stage = 'STAGE_1_REPRODUCIBLE_MORPHOLOGY_ANOMALY'
elif not stage_3:
    stage = 'STAGE_2_MULTISENSOR_CANDIDATE'
elif not stage_4:
    stage = 'STAGE_3_EXTERNAL_REVIEW_AND_PERMIT_PREPARATION'
else:
    stage = 'STAGE_4_EXPEDITION_READY_BRIEF'

result = {
    'artifact_id': 'JANUS-EXPEDITION-OUTREACH-COUNCIL-RUN-001-2026-08-22-v1.0',
    'created_utc': datetime.now(timezone.utc).isoformat(),
    'question': 'Should JANUS seek influential expedition/funding contacts for the HA10/Ascension target, whom should it approach, and at what evidence stage?',
    'council': council,
    'evidence_snapshot': evidence,
    'current_stage': stage,
    'janus_answer': {
        'SHOULD_CONTACT_DATA_CUSTODIANS_AND_RESEARCHERS_NOW': True,
        'SHOULD_CONTACT_EXPEDITION_OPERATORS_OR_FUNDERS_NOW': stage_3 or stage_4,
        'SHOULD_CONTACT_WEALTHY_PATRONS_NOW': stage_4,
        'CURRENT_COMMAND': 'OBTAIN_AND_ANALYZE_MEASURED_BATHYMETRY_BACKSCATTER_AND_SUBSEAFLOOR_EVIDENCE__DO_NOT_FUNDRAISE_ON_PYRAMID_CLAIM',
        'CURRENT_CLAIM_CEILING': 'KNOWN_VOLCANIC_SEAMOUNT_SYSTEM_WITH_A_TARGETED_FALSIFIABLE_SEARCH_FOR_UNEXPLAINED_GEOMETRY',
        'PYRAMID_CLAIM_STATUS': 'NOT_ADMITTED',
        'IMPACT_CRATER_CLAIM_STATUS': 'NOT_SUPPORTED_BY_CURRENT_PUBLISHED_GEOMORPHOLOGY',
    },
    'contact_ladder': [
        {
            'stage': 'NOW',
            'who': ['Ascension Island Government Conservation and Fisheries Directorate', 'British Geological Survey marine geoscience/data team', 'British Antarctic Survey / UK Polar Data Centre bathymetry custodians', 'CTBTO/EarthScope hydroacoustic data custodians'],
            'purpose': 'data access, provenance, survey coverage, scientific clarification, permit path',
            'language': 'seafloor morphology / data-validation project; no pyramid claim',
        },
        {
            'stage': 'AFTER_STAGE_1',
            'who': ['independent marine geophysicists', 'volcanologists specialising in seamounts', 'marine archaeologists with deep-water remote-sensing experience'],
            'purpose': 'blind external review and alternative-explanation challenge',
            'language': 'reproducible morphology anomaly, not archaeological identification',
        },
        {
            'stage': 'AFTER_STAGE_2_AND_EXTERNAL_REVIEW',
            'who': ['Schmidt Ocean Institute', 'OceanX', 'REV Ocean', 'Ocean Exploration Trust / E/V Nautilus or comparable deep-ocean operators'],
            'purpose': 'non-invasive mapping/AUV/ROV/sub-bottom expedition feasibility',
            'language': 'evidence package with raw data, controls, uncertainties, permissions and falsification plan',
        },
        {
            'stage': 'ONLY_AFTER_STAGE_4',
            'who': ['philanthropic principals and major patrons reached through their scientific organisations/foundations'],
            'purpose': 'expedition support only after a scientifically reviewed, permit-aware target exists',
            'language': 'never lead with symbolic/celestial evidence; distinguish candidate from confirmed archaeology',
        },
    ],
    'expedition_escalation': [
        'DESKTOP_REANALYSIS_OF_EXISTING_MBES_AND_BACKSCATTER',
        'INDEPENDENT_HIGH_RESOLUTION_MAPPING',
        'AUV_OR_ROV_STANDOFF_SURVEY',
        'SUB_BOTTOM_OR_SEISMIC_PROFILING_IF_BURIED_TARGET_REMAINS_PLAUSIBLE',
        'ROV_IMAGERY_AND_MINIMALLY_INVASIVE_SAMPLING_IF_PERMITTED',
        'CREWED_SUBMERSIBLE_ONLY_IF_SCIENTIFICALLY_JUSTIFIED_LEGALLY_PERMITTED_AND_SAFER_REMOTE_STEPS_ARE_INSUFFICIENT',
    ],
    'promotion_gates': {
        'STAGE_1': '>=2 independent high-resolution physical datasets + reproducible anomalous geometry + blind volcanic controls passed',
        'STAGE_2': 'STAGE_1 + at least one independent non-bathymetric physical sensor supports the same target',
        'STAGE_3': 'STAGE_2 + >=2 independent relevant specialists reproduce/review the result and common geological explanations remain inadequate',
        'STAGE_4': 'STAGE_3 + legal/MPA/permit path established, exact target package, operational risk plan, transparent uncertainty and budget-ready science brief',
    },
    'hard_rules': [
        'NO_UNDERWATER_PYRAMID_DETECTED',
        'VOLCANIC_BASELINE_MUST_NOT_BE_ERASED',
        'SYMBOLIC_CELESTIAL_OR_117_121_ASSOCIATION_CANNOT_ADVANCE_ARCHAEOLOGICAL_IDENTITY',
        'NO_FUNDRAISING_CLAIM_BEFORE_REPRODUCIBLE_PHYSICAL_ANOMALY',
        'NO_DIRECT_PATRON_OUTREACH_BEFORE_EXTERNAL_REVIEW_AND_PERMIT_AWARENESS',
        'REMOTE_NON_INVASIVE_SURVEY_PRECEDES_CREWED_DESCENT',
        'NEGATIVE_RESULTS_REMAIN_NEGATIVE',
    ],
    'public_opportunity_notes_2026': {
        'Schmidt_Ocean_Institute': 'EOIs for expeditions from 2029 onward are accepted on a rolling basis; ship access is provided to accepted scientific teams and research in an EEZ is contingent on authorization.',
        'OceanX': 'Institutional partnership route exists; approach through science/partnership channels rather than personal billionaire outreach.',
        'REV_Ocean': 'First operational science program begins in 2027 and includes South Atlantic missions.',
        'Ocean_Exploration_Trust': 'Strong mapping/ROV capability, but 2026 program is Pacific-focused; treat as a later operator/science-network contact, not an immediate Ascension mission.',
    },
    'target_identity': 'UNCONFIRMED',
}
result['sha256'] = hashlib.sha256(json.dumps(result, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding='utf-8')
print(json.dumps({
    'current_stage': result['current_stage'],
    'current_command': result['janus_answer']['CURRENT_COMMAND'],
    'contact_researchers_now': result['janus_answer']['SHOULD_CONTACT_DATA_CUSTODIANS_AND_RESEARCHERS_NOW'],
    'contact_expedition_funders_now': result['janus_answer']['SHOULD_CONTACT_EXPEDITION_OPERATORS_OR_FUNDERS_NOW'],
    'contact_wealthy_patrons_now': result['janus_answer']['SHOULD_CONTACT_WEALTHY_PATRONS_NOW'],
    'sha256': result['sha256'],
}, indent=2, ensure_ascii=False))
