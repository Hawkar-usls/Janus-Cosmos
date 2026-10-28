"""Adversarial zero-network tests of JANUS Cosmos Rowan source bridge."""
from __future__ import annotations

import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import rowan_cosmos_bridge as b


def encoded(obj):
    return b.render(obj)


def fixture():
    shadow = {
        "schema_id": "janus.rowan.bbc1.sequence_shadow_ridge.v1",
        "source": {"sha256": "a" * 64},
        "model": {"features": 44},
        "authority": {"deployment": "SHADOW_ONLY", "rowan_credits_used": 0},
    }
    topa = {
        "schema_id": "janus.topa.independent_shadow_replay.v1",
        "audit": "PASS_INDEPENDENT_WITHIN_ASSAY_REPLAY",
        "promotion_gate": "SHADOW_ONLY",
        "independent_remote_source_sha256": "a" * 64,
        "model_file_sha256": b.digest(encoded(shadow)),
        "heldout_variant_count": 259,
        "test_mae_independent": 1.0063159,
        "training_mean_baseline_mae_independent": 1.027715607,
        "model_family": "RIDGE_LINEAR_SEQUENCE_SHADOW",
        "external_lab_replication": "OPEN",
        "out_of_assay_generalization": "NOT_TESTED",
    }
    canonical = {
        "schema_id": "janus.rowan.observer_state.v1",
        "latest_durable_external_run": 12345,
        "last_completed_mutants": 8,
        "underlying_state_is_derived_from": "b" * 64,
    }
    csha = b.digest(encoded(canonical))
    aura = {
        "schema_id": "janus.aura.rowan_hypothesis_packet.v1",
        "source": {"state_sha256": csha, "upstream_run": 12345,
                   "completed_mutants": 8, "upstream_artifact_sha256": "b" * 64},
        "packet_digest_sha256": "c" * 64,
        "epistemic_status": "HYPOTHESIS_AND_FALSIFIER_PROPOSAL_ONLY",
        "source_has_unseen_outcome": False,
        "bounded_authority": {"may_submit_rowan_jobs": False},
        "next_discriminating_question_id": "EXTERNAL_HOLDOUT",
    }
    demi = {
        "schema_id": "janus.demiurge.rowan_council_proposal.v1",
        "scientific_authority": "NONE__PROPOSAL_NOT_EXPERIMENT",
        "gates": {"rowan_job_submission_new": False,
                  "production_model_weight_promotion": False},
        "sources": {"meta_registry": {
            "canonical_run": 12345, "sha256": csha,
            "source_artifact_sha256": "b" * 64},
            "aura_oracle": {"packet_sha256": b.digest(encoded(aura)),
                            "packet_digest": "c" * 64}},
        "next_independent_control": "external holdout",
    }
    coord = {
        "schema_id": "janus.meta_registry.rowan_council_state.v1",
        "meta_canonical_run": 12345,
        "source_binding_sha256": "d" * 64,
        "receipt_path": "data/rowan/coordination/receipts/ROWAN-COUNCIL-12345-dddddddddddddddd.json",
    }
    receipt = {
        "schema_id": "janus.meta_registry.rowan_council_provenance.v1",
        "source_binding_sha256": "d" * 64,
        "sources": {
            "meta_canonical_run": 12345,
            "meta_canonical_sha256": csha,
            "rowan_source_artifact_sha256": "b" * 64,
            "aura_packet_sha256": b.digest(encoded(aura)),
            "aura_packet_digest": "c" * 64,
            "demiurge_council_sha256": b.digest(encoded(demi)),
        },
        "authority": {
            "new_credit_spend": False,
            "new_rowan_jobs": False,
            "new_study_authorized": False,
            "hypotheses_are_empirical_proof": False,
        },
    }
    return {key: encoded(value) for key, value in
            {"canonical": canonical, "aura": aura, "demiurge": demi,
             "coordination": coord, "topa": topa,
             "shadow": shadow, "receipt": receipt}.items()}


class CosmosWitness(unittest.TestCase):
    def test_accept_exact_bound_independent_sources(self):
        packet = b.reconcile(fixture())
        self.assertEqual(packet["canonical_rowan_run"], 12345)
        self.assertEqual(packet["completed_mutants"], 8)
        self.assertEqual(packet["topa"]["test_variant_count"], 259)
        self.assertFalse(packet["astronomical_physical_evidence"])
        self.assertEqual(packet["rowan_provider_jobs_created"], 0)
        self.assertEqual(packet["rowan_credits_spent"], 0)
        self.assertEqual(packet["authority"], "OBSERVE_AND_ARCHIVE_ONLY")

    def test_canonical_drift_holds(self):
        raw = fixture()
        source = json.loads(raw["canonical"])
        source["latest_durable_external_run"] += 1
        raw["canonical"] = encoded(source)
        with self.assertRaisesRegex(b.GateHold, "AURA_CANONICAL_BINDING_DRIFT"):
            b.reconcile(raw)

    def test_aura_tamper_holds(self):
        raw = fixture()
        aura = json.loads(raw["aura"])
        aura["bounded_authority"]["may_submit_rowan_jobs"] = True
        raw["aura"] = encoded(aura)
        with self.assertRaises(b.GateHold):
            b.reconcile(raw)

    def test_demiurge_spend_authority_holds(self):
        raw = fixture()
        d = json.loads(raw["demiurge"])
        d["gates"]["rowan_job_submission_new"] = True
        raw["demiurge"] = encoded(d)
        with self.assertRaises(b.GateHold):
            b.reconcile(raw)

    def test_receipt_binding_drift_holds(self):
        raw = fixture()
        r = json.loads(raw["receipt"])
        r["source_binding_sha256"] = "f" * 64
        raw["receipt"] = encoded(r)
        with self.assertRaisesRegex(b.GateHold, "COORDINATION_RECEIPT_DRIFT"):
            b.reconcile(raw)

    def test_topa_promotion_leak_holds(self):
        raw = fixture()
        t = json.loads(raw["topa"])
        t["promotion_gate"] = "PRODUCTION"
        raw["topa"] = encoded(t)
        with self.assertRaisesRegex(b.GateHold, "TOPA_SHADOW_EVIDENCE_DRIFT"):
            b.reconcile(raw)

    def test_missing_source_holds(self):
        raw = fixture()
        del raw["receipt"]
        with self.assertRaisesRegex(b.GateHold, "INCOMPLETE_SOURCE_SET"):
            b.reconcile(raw)

    def test_idempotent_and_immutable(self):
        packet = b.reconcile(fixture())
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path, changed = b.persist(root, packet)
            self.assertTrue(changed)
            path2, changed = b.persist(root, packet)
            self.assertEqual(path, path2)
            self.assertFalse(changed)
            Path(path).write_text("tampered")
            with self.assertRaisesRegex(b.GateHold, "IMMUTABLE_RECEIPT_CHANGED"):
                b.persist(root, packet)

    def test_no_external_network_for_tests(self):
        self.assertTrue(set(b.SOURCES) >=
                        {"canonical", "aura", "demiurge", "coordination", "topa", "shadow"})
        self.assertNotIn("rowan_api_key", repr(b.SOURCES).lower())


if __name__ == "__main__":
    unittest.main()
