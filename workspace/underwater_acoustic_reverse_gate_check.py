#!/usr/bin/env python3
import json
from pathlib import Path

SRC = Path('data/love/JANUS-LOVE-EDEM-UNDERWATER-PYRAMID-REVERSE-ACOUSTIC-CLOSURE-GATE-2026-08-20-v1.0.json')
OUT = Path('data/love/JANUS-LOVE-EDEM-UNDERWATER-PYRAMID-REVERSE-ACOUSTIC-CLOSURE-RUN-001-RECEIPT.json')

def main():
    g = json.loads(SRC.read_text(encoding='utf-8'))
    findings = {x['finding']: x for x in g['reverse_audit_findings']}
    checks = {
        'point_frozen': g['frozen_point']['no_recenter'] is True,
        'reverse_starts_blind': g['reverse_operator'][0] == 'BLIND_SEARCH_OBSERVATION',
        'reverse_ends_forward_replay_closure': g['reverse_operator'][-2:] == ['FORWARD_REPLAY','HELDOUT_MATCH_OR_REJECT'],
        'unknown_scale_leak_repaired': findings['UNKNOWN_SCALE_LEAK']['status'] == 'REPAIRED',
        'fixed_frequency_circularity_repaired': findings['FIXED_FREQUENCY_CIRCULARITY']['status'] == 'REPAIRED',
        'inverse_nonuniqueness_preserved': findings['INVERSE_NONUNIQUENESS']['status'] == 'STRUCTURAL_LIMIT_PRESERVED',
        'buried_domain_branch_added': findings['BURIED_TARGET_DOMAIN_CHANGE']['status'] == 'BRANCH_ADDED',
        '119_not_blind_key': '117_121_hz_search_window' in g['blind_input_firewall']['forbidden'],
        '520_not_blind_key': '512_529_hz_search_window' in g['blind_input_firewall']['forbidden'],
        'target_label_not_allowed': 'label_pyramid' in g['blind_input_firewall']['forbidden'],
        'reference_controls_required': len(g['reference_library_required']['geometry_classes']) >= 5,
        'three_passes_required': 'at_least_3_independent_passes' in g['janus_backforth_closure']['advance_only_if'],
        'heldout_forward_replay_required': 'forward_replay_predicts_heldout_hydrophones_or_aspects' in g['janus_backforth_closure']['advance_only_if'],
        'no_underwater_pyramid_claim': 'NO_UNDERWATER_PYRAMID_DETECTED' in g['scientific_boundary']
    }
    status = 'PASS_REVERSE_OPERATOR_AUDIT__EMPIRICAL_IDENTIFIABILITY_BLOCKED_BY_REFERENCE_LIBRARY_AND_FIELD_DATA' if all(checks.values()) else 'FAIL_REVERSE_OPERATOR_AUDIT'
    receipt = {
        'artifact_id': 'JANUS-LOVE-EDEM-UNDERWATER-PYRAMID-REVERSE-ACOUSTIC-CLOSURE-RUN-001-2026-08-20-v1.0',
        'source': str(SRC),
        'status': status,
        'checks': checks,
        'reverse_conclusion': {
            'logical_reversibility': 'PASS',
            'unique_geometry_recovery_from_spectrum': 'NOT_ADMITTED',
            'blind_search_can_use_119_or_520_hz_prior': False,
            'scale_free_first_stage_required': True,
            'reference_library_required': True,
            'forward_replay_required': True,
            'current_field_detection_status': 'NO_TARGET_DETECTED'
        },
        'next_gate': 'BUILD_BLIND_REFERENCE_LIBRARY_AND_DEMONSTRATE_CLASS_STATE_SCALE_RECOVERY_ON_HELDOUT_CONTROLS'
    }
    OUT.write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(receipt, indent=2))
    if not all(checks.values()):
        raise SystemExit(1)

if __name__ == '__main__':
    main()
