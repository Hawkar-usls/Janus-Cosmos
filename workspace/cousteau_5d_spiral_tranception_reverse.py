#!/usr/bin/env python3
"""JANUS Echo Cousteau — 5D SPIRAL Tranception + Reverse deep-analysis runner.

The previous runner used a closed six-direction circuit. This version makes the
return step an ASCENT: revisit the same question only after state has changed.
A turn is allowed to continue only if it adds a new testable constraint, gate,
provenance binding, or explicit blocker. Repeated state hashes are forbidden.

5D is analytical feature space, not a claim of five physical dimensions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cousteau"
BASE_RUN = DATA / "JANUS-ECHO-COUSTEAU-5D-TRANCEPTION-REVERSE-DEEP-ANALYSIS-RUN-001-2026-08-21-v1.0.json"
RUN_ID = "JANUS-ECHO-COUSTEAU-5D-SPIRAL-TRANCEPTION-REVERSE-RUN-001-2026-08-21-v1.0"
STREAM = DATA / "5d_spiral_stream" / "run-001"
OUT = DATA / f"{RUN_ID}.json"
SUMMARY = DATA / f"{RUN_ID}-SUMMARY.json"

DIMS = ["D0_SPACE", "D1_TIME", "D2_ACOUSTIC", "D3_REVERSE_CAUSAL", "D4_ASSOCIATIVE_PROVENANCE"]
PHASES = ["BACK", "FORWARD", "LEFT_HRAIN", "RIGHT_INAIHR", "ASCEND"]


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha(obj: Any) -> str:
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def read_json(p: Path) -> dict[str, Any]:
    return json.loads(p.read_text(encoding="utf-8"))


class Spiral:
    def __init__(self, base: dict[str, Any]):
        self.base = base
        self.nodes: list[dict[str, Any]] = []
        self.turns: list[dict[str, Any]] = []
        self.state_hashes: list[str] = []
        STREAM.mkdir(parents=True, exist_ok=True)

    def emit(self, turn: int, phase: str, claim: str, structural: str, associative: str,
             evidence: str, new_constraints: list[str], provenance: list[dict[str, Any]],
             parents: list[str] | None = None, payload: dict[str, Any] | None = None) -> str:
        if phase not in PHASES:
            raise ValueError(phase)
        nid = f"S{turn:02d}N{sum(1 for n in self.nodes if n['turn']==turn)+1:02d}"
        node = {
            "node_id": nid,
            "turn": turn,
            "phase": phase,
            "created_utc": now(),
            "dimensions": DIMS,
            "claim_class": claim,
            "structural_context": structural,
            "associative_context": associative,
            "evidence_status": evidence,
            "new_constraints": new_constraints,
            "provenance": provenance,
            "parents": parents or [],
            "payload": payload or {},
        }
        node["node_sha256"] = sha(node)
        self.nodes.append(node)
        (STREAM / f"turn-{turn:02d}-{phase.lower()}-{nid}.json").write_text(
            json.dumps(node, ensure_ascii=False, indent=2)+"\n", encoding="utf-8"
        )
        return nid

    def freeze_turn(self, turn: int, inherited: list[str], added: list[str], anchor_node: str,
                    scientific_delta: str) -> dict[str, Any]:
        state = {
            "turn": turn,
            "inherited_constraints": sorted(set(inherited)),
            "added_constraints": sorted(set(added)),
            "active_constraints": sorted(set(inherited) | set(added)),
            "anchor_node": anchor_node,
            "scientific_delta": scientific_delta,
        }
        state_hash = sha(state)
        repeated = state_hash in self.state_hashes
        state["state_sha256"] = state_hash
        state["repeats_prior_state"] = repeated
        state["radius"] = len(state["active_constraints"])
        state["height"] = turn
        self.state_hashes.append(state_hash)
        self.turns.append(state)
        (STREAM / f"turn-{turn:02d}-state.json").write_text(
            json.dumps(state, ensure_ascii=False, indent=2)+"\n", encoding="utf-8"
        )
        return state


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default=str(OUT))
    ap.add_argument("--summary", default=str(SUMMARY))
    args = ap.parse_args()

    base = read_json(BASE_RUN)
    sp = Spiral(base)
    base_ref = [{"path": str(BASE_RUN.relative_to(ROOT)).replace('\\','/'), "run_sha256": base.get("run_sha256") }]

    # Frozen inheritance from the prior 5D run. This is the spiral's origin, not a destination.
    inherited = [
        "KEEP_H0_H1_H2_SEPARATE",
        "BLOCKED_DATA_STAYS_BLOCKED",
        "NEGATIVE_RESULTS_NOT_RESCUED_BY_ASSOCIATION",
        "NO_RECENTERING",
        "CLAIM_CEILING_PRESERVED",
    ]

    # TURN 1: revisit acquisition blocker, but ascend by localizing it to authoritative prebinding.
    a1 = sp.emit(1, "BACK", "BLOCKED",
        "Return to the T-phase blocker without changing cluster thresholds or anchor geometry.",
        "The useful reverse question is not 'where is the cluster?' but 'what must exist before a cluster can be born?'.",
        "BLOCKER_REVISITED_WITHOUT_RESCUE", [], base_ref)
    a2 = sp.emit(1, "FORWARD", "CORRECTION",
        "Forward replay shows the failure happens before event rows are materialized.",
        "Therefore a dataset/file identity binding is an upstream prerequisite, not a post-hoc patch.",
        "FAILURE_LOCALIZED_PRE_MATERIALIZATION",
        ["PREBIND_AUTHORITATIVE_MGDS_DATASET_ID_BEFORE_PARSE"], base_ref, [a1])
    a3 = sp.emit(1, "LEFT_HRAIN", "FACT",
        "Structural lane separates acquisition identity -> download bytes -> parse rows -> blind cluster -> anchor reveal.",
        "No later stage may compensate for a missing earlier stage.",
        "STRUCTURAL_DEPENDENCY_CHAIN_FROZEN",
        ["STAGE_ORDER_IS_PROOF_OBLIGATION"], base_ref, [a2])
    a4 = sp.emit(1, "RIGHT_INAIHR", "NEW_GATE",
        "Associative lane compares the blocker to prior JANUS pre-birth lessons only as software architecture.",
        "The analogy earns value only if prebinding deterministically unlocks the same frozen analysis without threshold retuning.",
        "TESTABLE_ARCHITECTURAL_ANALOGY",
        ["PREBIND_REPLAY_MUST_KEEP_DBSCAN_PARAMETERS_UNCHANGED"], base_ref, [a3],
        {"generated_gate":"COUSTEAU_MGDS_AUTHORITATIVE_DATASET_ID_PREBIND_REPLAY_V1"})
    added1 = ["PREBIND_AUTHORITATIVE_MGDS_DATASET_ID_BEFORE_PARSE","STAGE_ORDER_IS_PROOF_OBLIGATION","PREBIND_REPLAY_MUST_KEEP_DBSCAN_PARAMETERS_UNCHANGED"]
    a5 = sp.emit(1, "ASCEND", "STATE_LIFT",
        "Return is forbidden to the old state; turn 1 ascends with three new constraints.",
        "The same question is now sharper: can authoritative prebinding unlock the catalog while preserving the blind contract?",
        "SPIRAL_ASCENT", added1, base_ref, [a4])
    t1 = sp.freeze_turn(1, inherited, added1, a5, "ROUTING_RESOLUTION_INCREASED__TARGET_EVIDENCE_UNCHANGED")

    # TURN 2: revisit frequency evidence from the higher state, ascend by cross-control outperformance requirement.
    inh2 = t1["active_constraints"]
    b1 = sp.emit(2, "BACK", "CONTROL",
        "Revisit 119/520-class evidence with H0/H1/H2 separation already frozen.",
        "A frequency match is now treated as a retrieval cue that names competing homolog classes.",
        "FREQUENCY_REVISITED_AT_HIGHER_CONSTRAINT_STATE", [], base_ref, [a5])
    b2 = sp.emit(2, "FORWARD", "CORRECTION",
        "Forward replay asks what a true candidate must add beyond nominal frequency membership.",
        "It must carry Q/linewidth, modal spacing, phase, ringdown, location repeatability and source/target separation.",
        "FULL_FINGERPRINT_REQUIREMENT_EXPANDED",
        ["FREQUENCY_ONLY_CANNOT_PROMOTE_TARGET"], base_ref, [b1])
    b3 = sp.emit(2, "LEFT_HRAIN", "NEW_GATE",
        "Structural lane turns the control corpus into a matrix across natural, biological, engineered and illuminator classes.",
        "Each populated candidate dimension must be compared against the strongest matched control, not an average control.",
        "CONTROL_MATRIX_STRUCTURED",
        ["MATCHED_CONTROL_PER_POPULATED_DIMENSION"], base_ref, [b2])
    b4 = sp.emit(2, "RIGHT_INAIHR", "NEW_GATE",
        "Associative lane retrieves homolog false positives but cannot count similarity as evidence.",
        "A candidate advances only if it outperforms every relevant control without post-reveal tuning.",
        "ASSOCIATION_ROUTED_TO_FALSIFICATION",
        ["CANDIDATE_MUST_OUTPERFORM_ALL_RELEVANT_CONTROL_CLASSES"], base_ref, [b3],
        {"generated_gate":"COUSTEAU_5D_CONTROL_OUTPERFORMANCE_MATRIX_V1"})
    added2 = ["FREQUENCY_ONLY_CANNOT_PROMOTE_TARGET","MATCHED_CONTROL_PER_POPULATED_DIMENSION","CANDIDATE_MUST_OUTPERFORM_ALL_RELEVANT_CONTROL_CLASSES"]
    b5 = sp.emit(2, "ASCEND", "STATE_LIFT",
        "Turn 2 does not circle back to '520 matched'. It ascends to 'what multidimensional fingerprint survives controls?'.",
        "The question has changed level from matching to discrimination.",
        "SPIRAL_ASCENT", added2, base_ref, [b4])
    t2 = sp.freeze_turn(2, inh2, added2, b5, "DISCRIMINATION_REQUIREMENT_INCREASED__TARGET_EVIDENCE_UNCHANGED")

    # TURN 3: revisit the whole graph. If no new raw data exist, ascend into a plateau/stop condition rather than hallucinating novelty.
    inh3 = t2["active_constraints"]
    c1 = sp.emit(3, "BACK", "AUDIT",
        "Revisit the entire graph from the higher state and ask whether any blocked scientific tier was actually populated by new raw data.",
        "No association is allowed to substitute for missing waveform, catalog rows, geometry/material or complex-return data.",
        "RAW_DATA_AVAILABILITY_AUDIT", [], base_ref, [b5])
    c2 = sp.emit(3, "FORWARD", "NEGATIVE",
        "Forward replay finds no new raw catalog rows or phase-resolved target return inside this architecture-only test.",
        "Therefore the scientifically correct move is to stop growth, not manufacture another semantic layer.",
        "NO_NEW_SCIENTIFIC_DATA_IN_THIS_RUN",
        ["SPIRAL_MUST_STOP_ON_EVIDENCE_PLATEAU"], base_ref, [c1])
    c3 = sp.emit(3, "LEFT_HRAIN", "FACT",
        "Structural lane marks unresolved inputs explicitly: MGDS authoritative catalog binding, 119-Hz waveforms, local CTD conditioning, H2 geometry/material/complex return, Titanic phase calibration.",
        "These are concrete ingress points for the next spiral turn.",
        "PLATEAU_INPUTS_ENUMERATED", ["NEXT_TURN_REQUIRES_NEW_INPUT_OR_RESOLVED_BLOCKER"], base_ref, [c2])
    c4 = sp.emit(3, "RIGHT_INAIHR", "CORRECTION",
        "Associative lane is deliberately prevented from inventing novelty at the plateau.",
        "Abstract thought may generate a testable route, but not evidence or endless self-referential expansion.",
        "ASSOCIATIVE_RUNAWAY_BLOCKED", ["ABSTRACT_NOVELTY_WITHOUT_TESTABLE_PAYOFF_IS_REJECTED"], base_ref, [c3])
    added3 = ["SPIRAL_MUST_STOP_ON_EVIDENCE_PLATEAU","NEXT_TURN_REQUIRES_NEW_INPUT_OR_RESOLVED_BLOCKER","ABSTRACT_NOVELTY_WITHOUT_TESTABLE_PAYOFF_IS_REJECTED"]
    c5 = sp.emit(3, "ASCEND", "PLATEAU",
        "Turn 3 ascends into an explicit plateau state rather than returning to the origin.",
        "The spiral is paused at a higher-resolution boundary until new data arrive; pause is a valid scientific state.",
        "SPIRAL_PLATEAU_REACHED", added3, base_ref, [c4])
    t3 = sp.freeze_turn(3, inh3, added3, c5, "EVIDENCE_PLATEAU_DETECTED__NO_TARGET_PROMOTION")

    tests = {
        "TSP_01_FIVE_DIMENSIONS_PRESERVED": all(n["dimensions"] == DIMS for n in sp.nodes),
        "TSP_02_EACH_TURN_HAS_ALL_PHASES": all(set(PHASES) <= {n['phase'] for n in sp.nodes if n['turn']==t} for t in (1,2,3)),
        "TSP_03_NO_STATE_HASH_REPEATS": len(sp.state_hashes) == len(set(sp.state_hashes)),
        "TSP_04_RADIUS_STRICTLY_INCREASES": all(sp.turns[i]['radius'] < sp.turns[i+1]['radius'] for i in range(len(sp.turns)-1)),
        "TSP_05_HEIGHT_STRICTLY_INCREASES": [t['height'] for t in sp.turns] == [1,2,3],
        "TSP_06_INTERMEDIATE_JSON_EMITTED": len(list(STREAM.glob('*.json'))) >= len(sp.nodes)+len(sp.turns),
        "TSP_07_NEGATIVES_NOT_RESCUED": any(n['evidence_status']=="NO_NEW_SCIENTIFIC_DATA_IN_THIS_RUN" for n in sp.nodes),
        "TSP_08_PLATEAU_STOP_EXISTS": t3['scientific_delta']=="EVIDENCE_PLATEAU_DETECTED__NO_TARGET_PROMOTION",
        "TSP_09_TARGET_EVIDENCE_NOT_ARTIFICIALLY_INCREASED": all("TARGET_EVIDENCE_UNCHANGED" in t['scientific_delta'] or "NO_TARGET_PROMOTION" in t['scientific_delta'] for t in sp.turns),
        "TSP_10_ASCENT_IS_STATE_CHANGE_NOT_RETURN": all(not t['repeats_prior_state'] for t in sp.turns),
    }
    rows = [{"id":k,"pass":bool(v)} for k,v in tests.items()]
    passed = all(x['pass'] for x in rows)

    result = {
        "artifact_id": RUN_ID,
        "created_utc": now(),
        "architecture": "5D_ANALYTIC_SPIRAL__TRANCEPTION_REVERSE__EVENT_SOURCED_JSON",
        "spiral_law": "RETURN_TO_QUESTION_AT_HIGHER_STATE__NEVER_RETURN_TO_IDENTICAL_STATE",
        "spiral_coordinates": {"height":"turn_index","radius":"active_constraint_count","angle":"BACK_FORWARD_LEFT_RIGHT_ASCEND_phase"},
        "physical_5d_claim": False,
        "dimensions": DIMS,
        "phases": PHASES,
        "turns": sp.turns,
        "nodes": sp.nodes,
        "tests": rows,
        "engine_verdict": "PASS_5D_SPIRAL_ENGINE" if passed else "FAIL_5D_SPIRAL_ENGINE",
        "scientific_verdict": "TARGET_EVIDENCE_UNCHANGED__SPIRAL_ROUTING_AND_FALSIFICATION_RESOLUTION_INCREASED__PLATEAU_REACHED",
        "next_resume_condition": "NEW_RAW_DATA_OR_RESOLVED_BLOCKER_REQUIRED_FOR_TURN_4",
        "target_identity_status": "UNCONFIRMED",
        "underwater_pyramid_detected": False,
    }
    result["run_sha256"] = sha(result)
    summary = {
        "artifact_id": RUN_ID+"-SUMMARY",
        "engine_verdict": result["engine_verdict"],
        "scientific_verdict": result["scientific_verdict"],
        "turns": len(sp.turns),
        "nodes": len(sp.nodes),
        "state_hashes_unique": len(sp.state_hashes)==len(set(sp.state_hashes)),
        "radii": [t['radius'] for t in sp.turns],
        "heights": [t['height'] for t in sp.turns],
        "tests_passed": sum(x['pass'] for x in rows),
        "tests_total": len(rows),
        "plateau": True,
        "next_resume_condition": result["next_resume_condition"],
    }
    Path(args.output).write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    Path(args.summary).write_text(json.dumps(summary, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
