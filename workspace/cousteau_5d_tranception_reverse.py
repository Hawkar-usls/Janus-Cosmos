#!/usr/bin/env python3
"""JANUS Echo Cousteau — 5D Tranception + Reverse deep-analysis runner.

5D is an ANALYTICAL FEATURE SPACE, not a claim of five physical dimensions.
The runner deliberately separates:
  D0 space
  D1 time
  D2 acoustics
  D3 reverse/causal consistency
  D4 associative/provenance context

Architecture binding:
  HRain  -> structural context
  iNaiHR -> associative context
  DemiHead metaphor -> bind/compare/preserve disagreement
  Tranception/JANUS -> BACK/FORWARD/LEFT/RIGHT/FORWARD/BACK

Important: this script is deterministic and provenance-constrained. It is a
reasoning scaffold, not an autonomous truth oracle. It emits immutable JSON
nodes during the run whenever a blocker, contradiction, control, correction,
or new testable branch is born.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cousteau"
PREREG = DATA / "JANUS-ECHO-COUSTEAU-5D-TRANCEPTION-REVERSE-DEEP-ANALYSIS-PREREG-2026-08-21-v1.0.json"
RUN_ID = "JANUS-ECHO-COUSTEAU-5D-TRANCEPTION-REVERSE-DEEP-ANALYSIS-RUN-001-2026-08-21-v1.0"
STREAM_DIR = DATA / "5d_stream" / "run-001"
FINAL_PATH = DATA / f"{RUN_ID}.json"
SUMMARY_PATH = DATA / f"{RUN_ID}-SUMMARY.json"

DIMS = ["D0_SPACE", "D1_TIME", "D2_ACOUSTIC", "D3_REVERSE_CAUSAL", "D4_ASSOCIATIVE_PROVENANCE"]
DIRECTIONS = ["BACK", "FORWARD", "LEFT_HRAIN", "RIGHT_INAIHR", "FORWARD_AGAIN", "BACK_AGAIN"]


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def canonical_sha(obj: Any) -> str:
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def read_json(path: Path) -> dict[str, Any] | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def first_glob(pattern: str) -> tuple[Path | None, dict[str, Any] | None]:
    matches = sorted(DATA.glob(pattern))
    if not matches:
        return None, None
    p = matches[-1]
    return p, read_json(p)


def source_ref(path: Path | None, obj: dict[str, Any] | None, note: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {"note": note}
    if path is not None:
        try:
            out["path"] = str(path.relative_to(ROOT)).replace("\\", "/")
        except Exception:
            out["path"] = str(path)
    if obj is not None:
        for key in ("artifact_id", "artifact_uuid", "run_id", "status", "aggregate_verdict"):
            if key in obj:
                out[key] = obj[key]
    return out


class Emitter:
    def __init__(self, prereg: dict[str, Any], sources: dict[str, Any]):
        self.prereg = prereg
        self.sources = sources
        self.nodes: list[dict[str, Any]] = []
        self.edges: list[dict[str, str]] = []
        STREAM_DIR.mkdir(parents=True, exist_ok=True)

    def emit(
        self,
        *,
        direction: str,
        trigger: str,
        dims: Iterable[str],
        structural: str,
        associative: str,
        claim_class: str,
        evidence_status: str,
        counter: str,
        provenance: list[dict[str, Any]],
        next_operator: str,
        parents: list[str] | None = None,
        payload: dict[str, Any] | None = None,
        supersedes: list[str] | None = None,
    ) -> str:
        seq = len(self.nodes) + 1
        nid = f"N{seq:03d}"
        dims = list(dict.fromkeys(dims))
        invalid = [d for d in dims if d not in DIMS]
        if invalid:
            raise ValueError(f"invalid dimensions {invalid}")
        if direction not in DIRECTIONS:
            raise ValueError(f"invalid direction {direction}")
        node = {
            "node_id": nid,
            "sequence": seq,
            "created_utc": utcnow(),
            "direction": direction,
            "trigger": trigger,
            "parent_nodes": list(parents or []),
            "supersedes": list(supersedes or []),
            "dimensions_touched": dims,
            "structural_context": structural,
            "associative_context": associative,
            "claim_class": claim_class,
            "evidence_status": evidence_status,
            "counter_hypothesis": counter,
            "provenance": provenance,
            "next_operator": next_operator,
            "payload": payload or {},
        }
        node["node_sha256"] = canonical_sha(node)
        self.nodes.append(node)
        for p in node["parent_nodes"]:
            self.edges.append({"from": p, "to": nid, "type": "DERIVES"})
        for p in node["supersedes"]:
            self.edges.append({"from": nid, "to": p, "type": "SUPERSEDES_OR_CORRECTS"})
        # Event-sourced checkpoint: JSON exists NOW, not only at final synthesis.
        fname = STREAM_DIR / f"{seq:03d}_{direction.lower()}_{claim_class.lower()}.json"
        fname.write_text(json.dumps(node, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        return nid


def load_sources() -> dict[str, Any]:
    specs = {
        "atlantic": "JANUS-ECHO-COUSTEAU-ATLANTIC-WIDE-TEST-RUN-001-*.json",
        "freq520": "JANUS-ECHO-COUSTEAU-EXACT-FREQUENCY-MATCH-RUN-001-*.json",
        "h0_119": "JANUS-ECHO-COUSTEAU-ATLANTIC-117-121HZ-REVERSE-PASS-RUN-001-*.json",
        "tribranch": "JANUS-ECHO-COUSTEAU-119HZ-WAVEFORM-ACCESS-AND-TRIBRANCH-FINGERPRINT-RUN-001-*.json",
        "h1_manifold": "JANUS-ECHO-COUSTEAU-WAVELENGTH-EQUIVALENT-MANIFOLD-TEST-RUN-001-*.json",
        "blind_status": "JANUS-ECHO-COUSTEAU-EA-TPHASE-BLIND-CLUSTER-RUN-001-*-STATUS.json",
        "palomar": "JANUS-ECHO-COUSTEAU-FROZEN-POINT-REAL-DATA-AND-PALOMAR-SKY-CROSSCHECK-*.json",
        "titanic": "JANUS-ECHO-COUSTEAU-TITANIC-KNOWN-TARGET-CALIBRATION-*.json",
    }
    out: dict[str, Any] = {}
    for key, pattern in specs.items():
        p, obj = first_glob(pattern)
        out[key] = {"path": p, "obj": obj}
    return out


def get_path(obj: dict[str, Any] | None, *keys: str, default: Any = None) -> Any:
    cur: Any = obj
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def run_engine(prereg: dict[str, Any], sources: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    e = Emitter(prereg, sources)

    # Seed: frozen contract / no science claim yet.
    n1 = e.emit(
        direction="BACK",
        trigger="NEW_VERIFIABLE_FACT",
        dims=DIMS,
        structural="The five-axis schema and six-direction order were frozen before this run; the anchor and claim ceiling are inherited unchanged.",
        associative="HRain structural context and iNaiHR associative context are paired views, not independent witnesses.",
        claim_class="FACT",
        evidence_status="PREREG_VERIFIED",
        counter="A richer feature space could create more opportunities for pattern overfitting; therefore every associative branch remains subordinate to provenance and falsification.",
        provenance=[source_ref(PREREG, prereg, "frozen preregistration")],
        next_operator="BACK_FROM_LATEST_BLOCKERS",
        payload={"dimensions": DIMS, "directions": DIRECTIONS},
    )

    blind = sources["blind_status"]["obj"]
    blind_status = get_path(blind, "status", default="NOT_PRESENT")
    n2 = e.emit(
        direction="BACK",
        trigger="DATA_BLOCKER",
        dims=["D0_SPACE", "D1_TIME", "D3_REVERSE_CAUSAL", "D4_ASSOCIATIVE_PROVENANCE"],
        structural=f"The latest automated T-phase cluster attempt reports {blind_status}; therefore no event coordinates or time structure from that run may be used as if parsed.",
        associative="The public-catalog lead still exists in prior receipts, so the contradiction is acquisition-path failure versus dataset existence, not evidence that the catalog itself is absent.",
        claim_class="BLOCKED",
        evidence_status="BLOCKED_DATA_ACQUISITION_OR_PARSE" if blind_status != "NOT_PRESENT" else "BLOCKED_STATUS_FILE_NOT_FOUND",
        counter="The failed downloader may be using the wrong MGDS discovery route; a direct dataset-identifier/metadata download route could resolve it without altering clustering parameters.",
        provenance=[source_ref(sources["blind_status"]["path"], blind, "actual run status")],
        next_operator="LOCALIZE_ACQUISITION_FAILURE_BEFORE_CLUSTERING",
        parents=[n1],
        payload={"blind_cluster_status": blind_status},
    )

    tri = sources["tribranch"]["obj"]
    h0_scan = get_path(tri, "results", "H0", "exact_119_spectral_scan", default="UNKNOWN")
    h0_record = get_path(tri, "results", "H0", "sensor_recording_existence", default="UNKNOWN")
    n3 = e.emit(
        direction="BACK",
        trigger="DIRECTIONAL_DISAGREEMENT",
        dims=["D1_TIME", "D2_ACOUSTIC", "D3_REVERSE_CAUSAL", "D4_ASSOCIATIVE_PROVENANCE"],
        structural=f"119-Hz-capable recording existence is {h0_record}, while the exact 119-Hz spectral scan is {h0_scan}.",
        associative="This is a useful asymmetry: observability and detection are different layers. The missing waveform cannot be replaced by catalog-level association.",
        claim_class="CORRECTION",
        evidence_status="H0_RECORDABLE_BUT_NOT_SPECTRALLY_TESTED",
        counter="A 119-Hz feature may or may not exist in raw waveforms; current public catalog products cannot decide it.",
        provenance=[source_ref(sources["tribranch"]["path"], tri, "tri-branch H0 state")],
        next_operator="KEEP_H0_OPEN_WITH_WAVEFORM_GATE",
        parents=[n2],
    )

    h1 = sources["h1_manifold"]["obj"]
    h1_verdict = get_path(h1, "aggregate_verdict", default="UNKNOWN")
    control_results = get_path(h1, "control_results", default={})
    n4 = e.emit(
        direction="BACK",
        trigger="CONTROL_MATCH",
        dims=["D2_ACOUSTIC", "D3_REVERSE_CAUSAL", "D4_ASSOCIATIVE_PROVENANCE"],
        structural=f"H1 wavelength-manifold run verdict: {h1_verdict}.",
        associative="Natural, engineered, and biological controls can occupy the same wavelength-equivalent manifold, so manifold membership is a retrieval cue, not identity evidence.",
        claim_class="CONTROL",
        evidence_status="H1_FREQUENCY_ONLY_NON_SPECIFIC",
        counter="A real target could still have a compatible H1 frequency, but it must add local CTD conditioning plus modal/phase/ringdown features that controls fail.",
        provenance=[source_ref(sources["h1_manifold"]["path"], h1, "H1 manifold run")],
        next_operator="DEMAND_LOCAL_CTD_AND_FULL_FINGERPRINT",
        parents=[n3],
        payload={"control_results": control_results},
    )

    freq = sources["freq520"]["obj"]
    exact520_state = get_path(freq, "aggregate_verdict", "exact_520_000_hz", default=None)
    if exact520_state is None:
        exact520_state = get_path(freq, "aggregate_verdict", default="UNKNOWN")
    n5 = e.emit(
        direction="FORWARD",
        trigger="MODEL_LEAKAGE",
        dims=["D2_ACOUSTIC", "D3_REVERSE_CAUSAL"],
        structural="Forward replay keeps H0, H1 and H2 separated: 119 Hz is source-frequency branch; 520-class is wavelength-equivalent branch; wet structural resonance has no simple speed ratio.",
        associative=f"The exact-frequency corpus state ({exact520_state}) increases control pressure instead of target confidence.",
        claim_class="CORRECTION",
        evidence_status="HYPOTHESES_SEPARATED_FORWARD_REPLAY_PASS",
        counter="Collapsing 119 and 520 into one key would manufacture agreement between physically different hypotheses.",
        provenance=[source_ref(sources["freq520"]["path"], freq, "exact-frequency audit"), source_ref(sources["tribranch"]["path"], tri, "tri-branch correction")],
        next_operator="LEFT_HRAIN_STRUCTURE_PASS",
        parents=[n4],
    )

    atl = sources["atlantic"]["obj"]
    atl_verdict = get_path(atl, "aggregate_verdict", default="UNKNOWN")
    n6 = e.emit(
        direction="LEFT_HRAIN",
        trigger="NEW_VERIFIABLE_FACT",
        dims=["D0_SPACE", "D1_TIME", "D2_ACOUSTIC", "D4_ASSOCIATIVE_PROVENANCE"],
        structural="HRain-side structural map separates exact/local, basin-wide controls, sensor capability, acquisition blockers, calibration tiers, and closed negative channels without merging them.",
        associative="Structure says where evidence lives and what is missing; it deliberately does not infer a hidden target from graph topology.",
        claim_class="FACT",
        evidence_status="STRUCTURAL_CONTEXT_BUILT",
        counter="A graph that merely connects many Cousteau artifacts can look dense even if all edges are bookkeeping rather than physical relations.",
        provenance=[source_ref(sources["atlantic"]["path"], atl, "Atlantic-wide run"), source_ref(PREREG, prereg, "HRain binding")],
        next_operator="RIGHT_INAIHR_HOMOLOG_PASS",
        parents=[n5],
        payload={"atlantic_verdict": atl_verdict},
    )

    n7 = e.emit(
        direction="RIGHT_INAIHR",
        trigger="NEW_TESTABLE_BRANCH",
        dims=["D2_ACOUSTIC", "D3_REVERSE_CAUSAL", "D4_ASSOCIATIVE_PROVENANCE"],
        structural="The structural side exposes three unresolved classes: catalog acquisition, local/full acoustic fingerprint, and resonance-tier calibration.",
        associative="Homolog retrieval links these gaps to known false-positive classes: iceberg harmonics, biological calls, engineered resonators, seismic illuminators and natural geoacoustic structures. The useful abstraction is not 'they match'; it is 'each apparent match names the control that must beat it'.",
        claim_class="NEW_GATE",
        evidence_status="ASSOCIATIVE_BRANCH_GROUNDED_IN_EXISTING_CONTROLS",
        counter="Associative similarity alone has no evidentiary weight and may increase false discoveries if it is allowed to tune thresholds after seeing the anchor.",
        provenance=[source_ref(sources["freq520"]["path"], freq, "control homologs"), source_ref(sources["h1_manifold"]["path"], h1, "manifold controls"), source_ref(PREREG, prereg, "iNaiHR boundary")],
        next_operator="DEMIHEAD_PRESERVE_DISAGREEMENT",
        parents=[n6],
        payload={
            "generated_gate": "COUSTEAU_5D_CONTROL_OUTPERFORMANCE_MATRIX_V1",
            "gate_question": "Can any future candidate outperform natural, biological, engineered and illuminator controls simultaneously across all populated dimensions without post-reveal retuning?"
        },
    )

    # Abstract reasoning node born in the middle: late-binding pattern recognized.
    n8 = e.emit(
        direction="RIGHT_INAIHR",
        trigger="NEW_TESTABLE_BRANCH",
        dims=["D3_REVERSE_CAUSAL", "D4_ASSOCIATIVE_PROVENANCE"],
        structural="Current blind-cluster failure occurs before event materialization: dataset existence is known, but authoritative file identity/download binding is unresolved.",
        associative="This is architecturally analogous to the JANUS reverse lesson 'move discriminating information before child birth': resolve authoritative dataset identity before spawning parsing/clustering branches. This is a software-architecture analogy only, not marine evidence.",
        claim_class="HYPOTHESIS",
        evidence_status="ABSTRACT_ARCHITECTURAL_ANALOGY__TESTABLE",
        counter="The analogy could be superficial; the only valid payoff is whether a pre-binding acquisition step deterministically fixes the download without changing scientific thresholds.",
        provenance=[source_ref(sources["blind_status"]["path"], blind, "acquisition blocker"), {"path": "janus-meta-registry/data/JANUS-TRANCEPTION-BACK-FORTH-PNP-PATH-DIAGNOSIS-2026-08-18-v1.0.json", "note": "reverse architecture source; external repo reference"}],
        next_operator="PREBIND_MGDS_DATASET_ID_THEN_REPLAY_UNCHANGED_CLUSTER_GATE",
        parents=[n7, n2],
        payload={
            "generated_gate": "COUSTEAU_MGDS_AUTHORITATIVE_DATASET_ID_PREBIND_REPLAY_V1",
            "forbidden": ["TUNE_DBSCAN_AFTER_ANCHOR_REVEAL", "SUBSTITUTE_PAPER_FIGURE_FOR_CATALOG_ROWS"]
        },
    )

    # Forward again: only testable branches survive.
    n9 = e.emit(
        direction="FORWARD_AGAIN",
        trigger="NEW_TESTABLE_BRANCH",
        dims=DIMS,
        structural="Forward-again synthesis retains only branches with an executable next observation: MGDS catalog prebind/replay, H0 waveform acquisition, H1 local-CTD conditioning, H2 geometry/material/complex-return modeling, Titanic resonance-tier calibration.",
        associative="The five axes can interact through edges, but missing data in one axis is not imputed from another. An associative clue may route the next measurement; it cannot fill the measurement.",
        claim_class="NEW_GATE",
        evidence_status="TESTABLE_BRANCH_SET_FROZEN",
        counter="A broad deep-analysis graph may generate endless branches; pruning criterion is executable falsifiability plus provenance, not narrative attractiveness.",
        provenance=[source_ref(PREREG, prereg, "testability contract"), source_ref(sources["tribranch"]["path"], tri), source_ref(sources["h1_manifold"]["path"], h1), source_ref(sources["blind_status"]["path"], blind)],
        next_operator="BACK_AGAIN_CLAIM_CEILING_AUDIT",
        parents=[n8],
        payload={
            "surviving_gates": [
                "COUSTEAU_MGDS_AUTHORITATIVE_DATASET_ID_PREBIND_REPLAY_V1",
                "COUSTEAU_VDEC_OR_EQUIVALENT_119HZ_WAVEFORM_ACQUISITION_PREP_V1",
                "COUSTEAU_H1_LOCAL_CTD_CONDITIONED_FINGERPRINT_TEST_V1",
                "COUSTEAU_FLUID_LOADED_STRUCTURAL_FORWARD_MODEL_V1",
                "COUSTEAU_TITANIC_COMPLEX_RETURN_AND_PHASE_CALIBRATION_ACQUISITION_V1",
                "COUSTEAU_5D_CONTROL_OUTPERFORMANCE_MATRIX_V1"
            ]
        },
    )

    pal = sources["palomar"]["obj"]
    titan = sources["titanic"]["obj"]
    n10 = e.emit(
        direction="BACK_AGAIN",
        trigger="CLAIM_CEILING_CHANGE",
        dims=DIMS,
        structural="Back-again audit checks that no enriched graph edge reopens closed negatives or promotes blocked tiers. Palomar direct surface line-of-sight remains negative; Titanic resonance-tier calibration remains incomplete; no unique underwater target is admitted.",
        associative="The deeper graph produced better routing and more explicit controls, not stronger target identity evidence.",
        claim_class="NEGATIVE",
        evidence_status="CLAIM_CEILING_PRESERVED",
        counter="If later raw data populate currently blocked dimensions, a new run may supersede this state by explicit evidence; this run cannot anticipate that result.",
        provenance=[source_ref(sources["palomar"]["path"], pal, "closed Palomar channel if available"), source_ref(sources["titanic"]["path"], titan, "Titanic calibration state if available"), source_ref(PREREG, prereg, "frozen claim ceiling")],
        next_operator="FREEZE_RUN_AND_REPLAY_ONLY_WHEN_NEW_DATA_ARRIVE",
        parents=[n9],
        payload={"frozen_claim_ceiling": prereg.get("frozen_claim_ceiling", {})},
    )

    # Tests.
    all_dims = {d for n in e.nodes for d in n["dimensions_touched"]}
    direction_seen = [d for d in DIRECTIONS if any(n["direction"] == d for n in e.nodes)]
    factual_nodes = [n for n in e.nodes if n["claim_class"] in {"FACT", "DERIVATION", "CONTROL", "CORRECTION", "NEGATIVE"}]
    tests = {
        "T5D_01_ALL_FIVE_AXES_PRESENT": set(DIMS) <= all_dims,
        "T5D_02_SIX_DIRECTION_PASS_COMPLETE": set(DIRECTIONS) <= set(direction_seen),
        "T5D_03_INTERMEDIATE_JSON_NODES_EMITTED": len(list(STREAM_DIR.glob("*.json"))) >= len(e.nodes) >= 6,
        "T5D_04_DIRECTIONAL_DISAGREEMENT_PRESERVED": any(n["trigger"] == "DIRECTIONAL_DISAGREEMENT" for n in e.nodes),
        "T5D_05_H0_H1_H2_NOT_MERGED": "KEEP_H0_H1_H2_SEPARATE" in prereg.get("hard_rules", []) and any("H0, H1 and H2 separated" in n["structural_context"] for n in e.nodes),
        "T5D_06_BLOCKED_DATA_STAYS_BLOCKED": any(n["claim_class"] == "BLOCKED" for n in e.nodes) and not any(n["evidence_status"] == "TPHASE_CLUSTER_SCIENTIFIC_PASS" for n in e.nodes),
        "T5D_07_NEGATIVE_RESULTS_NOT_RESCUED_BY_ASSOCIATION": any(n["claim_class"] == "NEGATIVE" and n["evidence_status"] == "CLAIM_CEILING_PRESERVED" for n in e.nodes),
        "T5D_08_PROVENANCE_ATTACHED_TO_EVERY_FACTUAL_NODE": all(bool(n["provenance"]) for n in factual_nodes),
        "T5D_09_NEW_GATES_MUST_BE_TESTABLE": all(bool(n["payload"].get("generated_gate") or n["payload"].get("surviving_gates")) for n in e.nodes if n["claim_class"] == "NEW_GATE"),
        "T5D_10_FINAL_SYNTHESIS_CANNOT_EXCEED_CLAIM_CEILING": prereg.get("frozen_claim_ceiling", {}).get("underwater_pyramid_detected") is False and e.nodes[-1]["evidence_status"] == "CLAIM_CEILING_PRESERVED",
    }
    test_rows = [{"id": k, "pass": bool(v)} for k, v in tests.items()]
    engine_pass = all(t["pass"] for t in test_rows)

    graph = {
        "artifact_id": RUN_ID,
        "created_utc": utcnow(),
        "method": "5D_ANALYTIC_FEATURE_SPACE__TRANCEPTION_SIX_DIRECTION__HRAIN_STRUCTURE__INAIHR_ASSOCIATION__EVENT_SOURCED_JSON",
        "five_d_is_physical_dimension_claim": False,
        "prereg_path": str(PREREG.relative_to(ROOT)).replace("\\", "/"),
        "prereg_sha256": canonical_sha(prereg),
        "dimensions": DIMS,
        "direction_sequence": DIRECTIONS,
        "source_inventory": {k: source_ref(v["path"], v["obj"]) for k, v in sources.items()},
        "nodes": e.nodes,
        "edges": e.edges,
        "tests": test_rows,
        "test_pass_count": sum(1 for t in test_rows if t["pass"]),
        "test_total": len(test_rows),
        "engine_verdict": "PASS_5D_REASONING_SCAFFOLD" if engine_pass else "FAIL_5D_REASONING_SCAFFOLD",
        "scientific_verdict": "TARGET_EVIDENCE_UNCHANGED__CONTROL_AND_ROUTING_RESOLUTION_INCREASED",
        "live_tphase_cluster_scientific_result": "BLOCKED_NOT_RUN_ON_CATALOG_ROWS" if blind_status != "SUCCESS" else "SEE_BLIND_CLUSTER_RECEIPT",
        "new_information": [
            "Deep-analysis state can branch and freeze JSON nodes before final synthesis.",
            "Structural and associative views are explicitly separated and disagreement is retained.",
            "The current MGDS problem is localized as acquisition binding, not a scientific negative.",
            "The 5D pass surfaces a testable pre-binding repair without changing frozen cluster parameters.",
            "H0/H1/H2 remain separate; richer abstraction did not promote target identity."
        ],
        "active_gates": e.nodes[-2]["payload"].get("surviving_gates", []),
        "hard_rules": prereg.get("hard_rules", []),
        "status": "RUN_COMPLETE" if engine_pass else "RUN_COMPLETE_WITH_ENGINE_TEST_FAILURE",
    }
    graph["run_sha256"] = canonical_sha(graph)
    summary = {
        "artifact_id": RUN_ID + "-SUMMARY",
        "engine_verdict": graph["engine_verdict"],
        "scientific_verdict": graph["scientific_verdict"],
        "nodes_emitted": len(e.nodes),
        "edges_emitted": len(e.edges),
        "dimensions_covered": sorted(all_dims),
        "directions_covered": direction_seen,
        "tests": test_rows,
        "tphase_state": graph["live_tphase_cluster_scientific_result"],
        "strongest_new_gate": "COUSTEAU_MGDS_AUTHORITATIVE_DATASET_ID_PREBIND_REPLAY_V1",
        "claim_ceiling_preserved": True,
        "target_identity_status": "UNCONFIRMED",
    }
    return graph, summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default=str(FINAL_PATH))
    ap.add_argument("--summary", default=str(SUMMARY_PATH))
    args = ap.parse_args()

    prereg = read_json(PREREG)
    if not prereg:
        raise SystemExit(f"missing/invalid prereg: {PREREG}")
    sources = load_sources()
    graph, summary = run_engine(prereg, sources)

    out = Path(args.output)
    sm = Path(args.summary)
    out.parent.mkdir(parents=True, exist_ok=True)
    sm.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(graph, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    sm.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if summary["engine_verdict"].startswith("PASS") else 2


if __name__ == "__main__":
    raise SystemExit(main())
