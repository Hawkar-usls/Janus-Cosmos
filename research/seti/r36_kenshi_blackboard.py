#!/usr/bin/env python3
"""Kenshi-inspired shared blackboard for JANUS distributed swarm.

The old Kenshi swarm idea is used here only as an engineering architecture:
each agent keeps a persistent local identity and local state while exchanging
append-only events through a bounded shared blackboard. State may synchronize;
identity may not collapse.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


FOUNDATION_LAWS = (
    "EVERY_AGENT_HAS_A_PERSISTENT_LOCAL_IDENTITY",
    "SYNC_STATE_NOT_IDENTITY",
    "SHARED_BLACKBOARD_IS_HANDOFF_NOT_HIVE_MIND",
    "EVENTS_ARE_APPEND_ONLY",
    "FACTS_MAY_CLUSTER__AGENTS_MUST_NOT",
    "FAILURE_BECOMES_EXPERIENCE_NOT_ERASURE",
    "STALE_AGENT_REMAINS_A_KNOWN_ANCESTOR",
    "LEASE_COORDINATES_WORK_NOT_AUTHORITY",
    "PEER_MESSAGE_IS_ROUTING_CONTEXT_NOT_EMPIRICAL_EVIDENCE",
)

EVENT_TYPES = {
    "HEARTBEAT",
    "STATE_DELTA",
    "OBSERVATION",
    "SOURCE_POINTER",
    "REQUEST_HELP",
    "HANDOFF",
    "CONFLICT",
    "WORK_CLAIM",
    "QUERY_EXECUTION",
    "FETCH_ACK",
    "CAUSAL_ACK",
    "GENOME_BIND",
    "HEALTH",
    "LEASE_STATE",
}


def _jsonable(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(v) for v in value]
    return repr(value)


def _hash(value: Any) -> str:
    raw = json.dumps(_jsonable(value), ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class BlackboardEvent:
    event_id: str
    sequence: int
    kind: str
    sender_id: str
    recipient_ids: List[str]
    payload: Any
    parent_event_ids: List[str]
    genome_node_ids: List[str]
    source_bound: bool
    fingerprint: str

    def to_dict(self) -> Dict[str, Any]:
        return _jsonable(asdict(self))


class KenshiSwarmBlackboard:
    """Append-only event bus with persistent local agent identities."""

    def __init__(self, active_window: int = 128):
        if active_window <= 0:
            raise ValueError("active_window must be positive")
        self.active_window = int(active_window)
        self.agents: Dict[str, Dict[str, Any]] = {}
        self.events: List[BlackboardEvent] = []
        self.archive_event_ids: List[str] = []
        self.active_event_ids: List[str] = []
        self.leases: Dict[str, Dict[str, Any]] = {}

    def register_agent(self, agent_id: str, role: str, *, lineage_id: Optional[str] = None, local_state: Any = None) -> Dict[str, Any]:
        agent_id = str(agent_id).strip()
        if not agent_id:
            raise ValueError("agent_id must be non-empty")
        if agent_id in self.agents:
            existing = self.agents[agent_id]
            if lineage_id and existing["lineage_id"] != str(lineage_id):
                raise ValueError("agent lineage cannot be replaced")
            return existing
        record = {
            "agent_id": agent_id,
            "lineage_id": str(lineage_id or f"JANUS_AGENT_LINEAGE::{agent_id}"),
            "role": str(role),
            "status": "ACTIVE",
            "local_state": _jsonable(local_state),
            "last_event_id": None,
            "event_ids": [],
        }
        self.agents[agent_id] = record
        return record

    def _require_agent(self, agent_id: str) -> Dict[str, Any]:
        if agent_id not in self.agents:
            raise KeyError(f"unknown agent: {agent_id}")
        return self.agents[agent_id]

    def append_event(
        self,
        sender_id: str,
        kind: str,
        payload: Any,
        *,
        recipient_ids: Optional[Iterable[str]] = None,
        parent_event_ids: Optional[Iterable[str]] = None,
        genome_node_ids: Optional[Iterable[str]] = None,
        source_bound: bool = False,
    ) -> BlackboardEvent:
        sender = self._require_agent(sender_id)
        kind = str(kind).upper()
        if kind not in EVENT_TYPES:
            raise ValueError(f"unsupported event type: {kind}")
        recipients: List[str] = []
        for value in recipient_ids or []:
            rid = str(value)
            self._require_agent(rid)
            if rid not in recipients:
                recipients.append(rid)
        parents = [str(x) for x in (parent_event_ids or [])]
        known = {event.event_id for event in self.events}
        for parent in parents:
            if parent not in known:
                raise ValueError(f"unknown parent event: {parent}")
        sequence = len(self.events)
        body = {
            "sequence": sequence,
            "kind": kind,
            "sender_id": sender_id,
            "recipient_ids": recipients,
            "payload": _jsonable(payload),
            "parent_event_ids": parents,
            "genome_node_ids": [str(x) for x in (genome_node_ids or [])],
            "source_bound": bool(source_bound),
        }
        fingerprint = _hash(body)
        event = BlackboardEvent(
            event_id=f"KSB::{sequence:08d}::{fingerprint[:16]}",
            fingerprint=fingerprint,
            **body,
        )
        self.events.append(event)
        self.active_event_ids.append(event.event_id)
        if len(self.active_event_ids) > self.active_window:
            self.archive_event_ids.append(self.active_event_ids.pop(0))
        sender["last_event_id"] = event.event_id
        sender["event_ids"].append(event.event_id)
        return event

    def sync_local_state(self, agent_id: str, delta: Dict[str, Any], *, parent_event_ids: Optional[Iterable[str]] = None) -> BlackboardEvent:
        agent = self._require_agent(agent_id)
        current = dict(agent.get("local_state") or {}) if isinstance(agent.get("local_state"), dict) else {}
        current.update(_jsonable(delta))
        agent["local_state"] = current
        return self.append_event(
            agent_id,
            "STATE_DELTA",
            {"delta": _jsonable(delta), "state_after": current},
            parent_event_ids=parent_event_ids,
        )

    def handoff(self, sender_id: str, recipient_ids: Iterable[str], payload: Any, *, parent_event_ids: Optional[Iterable[str]] = None) -> BlackboardEvent:
        return self.append_event(
            sender_id,
            "HANDOFF",
            payload,
            recipient_ids=recipient_ids,
            parent_event_ids=parent_event_ids,
        )

    def mark_stale(self, agent_id: str, reason: str) -> BlackboardEvent:
        agent = self._require_agent(agent_id)
        agent["status"] = "STALE_KNOWN_ANCESTOR"
        return self.append_event(agent_id, "HEALTH", {"status": agent["status"], "reason": str(reason)})

    def recover(self, agent_id: str, state: Any) -> BlackboardEvent:
        agent = self._require_agent(agent_id)
        agent["status"] = "RECOVERED"
        agent["local_state"] = _jsonable(state)
        return self.append_event(agent_id, "HEALTH", {"status": "RECOVERED", "state": _jsonable(state)})

    def claim_lease(self, agent_id: str, objective_id: str, *, now_turn: int, ttl_turns: int) -> BlackboardEvent:
        self._require_agent(agent_id)
        if ttl_turns <= 0:
            raise ValueError("ttl_turns must be positive")
        objective_id = str(objective_id)
        existing = self.leases.get(objective_id)
        if existing and int(existing["expires_after_turn"]) >= int(now_turn) and existing["owner_agent_id"] != agent_id:
            raise RuntimeError("objective already has an active lease")
        previous_event_id = existing.get("event_id") if existing else None
        event = self.append_event(
            agent_id,
            "LEASE_STATE",
            {
                "objective_id": objective_id,
                "owner_agent_id": agent_id,
                "state": "ACTIVE",
                "acquired_turn": int(now_turn),
                "expires_after_turn": int(now_turn) + int(ttl_turns),
                "supersedes_expired_lease": bool(existing and existing["owner_agent_id"] != agent_id),
            },
            parent_event_ids=[previous_event_id] if previous_event_id else [],
        )
        self.leases[objective_id] = {
            "objective_id": objective_id,
            "owner_agent_id": agent_id,
            "acquired_turn": int(now_turn),
            "expires_after_turn": int(now_turn) + int(ttl_turns),
            "event_id": event.event_id,
        }
        return event

    def inbox(self, agent_id: str) -> List[Dict[str, Any]]:
        self._require_agent(agent_id)
        return [
            event.to_dict()
            for event in self.events
            if agent_id in event.recipient_ids
        ]

    def bind_genome(self, agent_id: str, genome_node_id: str, *, parent_event_ids: Optional[Iterable[str]] = None) -> BlackboardEvent:
        return self.append_event(
            agent_id,
            "GENOME_BIND",
            {"genome_node_id": str(genome_node_id)},
            genome_node_ids=[str(genome_node_id)],
            parent_event_ids=parent_event_ids,
        )

    def validate(self) -> None:
        ids = [event.event_id for event in self.events]
        if len(ids) != len(set(ids)):
            raise ValueError("duplicate event identity")
        for expected, event in enumerate(self.events):
            if event.sequence != expected:
                raise ValueError("event sequence broken")
            body = event.to_dict()
            fingerprint = body.pop("fingerprint")
            body.pop("event_id")
            if _hash(body) != fingerprint:
                raise ValueError("event fingerprint mismatch")
        for agent_id, record in self.agents.items():
            if record["agent_id"] != agent_id:
                raise ValueError("agent identity mismatch")

    def to_dict(self) -> Dict[str, Any]:
        self.validate()
        return {
            "schema": "janus.swarm.kenshi_blackboard.v1",
            "model": "LOCAL_AGENTS_SHARED_APPEND_ONLY_BLACKBOARD",
            "laws": list(FOUNDATION_LAWS),
            "agents": _jsonable(self.agents),
            "events": [event.to_dict() for event in self.events],
            "active_event_ids": list(self.active_event_ids),
            "archive_event_ids": list(self.archive_event_ids),
            "leases": _jsonable(self.leases),
            "epistemic_boundary": {
                "blackboard_event_is_empirical_evidence": False,
                "same_source_multi_agent_is_independent_replication": False,
                "source_bound_flag_is_provenance_metadata_not_truth": True,
            },
        }

    def save(self, path: str | Path) -> None:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self.to_dict(), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def self_test() -> None:
    board = KenshiSwarmBlackboard(active_window=2)
    board.register_agent("A", "SCOUT")
    board.register_agent("B", "SCOUT")
    a = board.append_event("A", "OBSERVATION", {"fact": "x", "source": "s"}, source_bound=True)
    b = board.append_event("B", "OBSERVATION", {"fact": "x", "source": "s"}, source_bound=True)
    assert a.event_id != b.event_id
    handoff = board.handoff("A", ["B"], {"request": "verify x"}, parent_event_ids=[a.event_id])
    assert board.inbox("B")[-1]["event_id"] == handoff.event_id
    board.claim_lease("A", "OBJ", now_turn=1, ttl_turns=2)
    try:
        board.claim_lease("B", "OBJ", now_turn=2, ttl_turns=2)
        raise AssertionError("active lease conflict was not rejected")
    except RuntimeError:
        pass
    board.claim_lease("B", "OBJ", now_turn=4, ttl_turns=2)
    board.mark_stale("A", "heartbeat timeout")
    board.recover("A", {"status": "RECOVERED"})
    assert board.agents["A"]["lineage_id"] == "JANUS_AGENT_LINEAGE::A"
    assert len(board.archive_event_ids) > 0
    board.validate()
    print("JANUS_KENSHI_SWARM_BLACKBOARD_SELF_TEST=PASS")


if __name__ == "__main__":
    self_test()
