#!/usr/bin/env python3
"""Strict conformance adapter for OSIRIS v3 spiral compute.

The frozen parent OSIRIS v2/S𓂸ḥ/2 technical_forward keeps the full verified SAT
assignment inside its engine but intentionally exports only assignment_sha256 in
the public projection. The spiral contract requires a reusable SAT witness.

This adapter does not alter parent solving. It observes the already-created
engine object at the existing engine_projection verification boundary, copies a
SAT assignment only when assignment_verified is true, and passes that exact
witness to the spiral experience builder. Fundamentum remains untouched.
"""
from __future__ import annotations

from typing import Any

import osiris_spiral_runtime as base

_CAPTURED_ASSIGNMENT: dict[int, bool] | None = None
_original_projection = base.gate.v2.engine_projection
_original_make_experience = base.make_experience
_original_solve_spiral = base.solve_spiral


def capturing_projection(engine: dict[str, Any]) -> dict[str, Any]:
    global _CAPTURED_ASSIGNMENT
    if (
        engine.get("status") == "SAT"
        and engine.get("assignment_verified") is True
        and isinstance(engine.get("assignment"), dict)
    ):
        _CAPTURED_ASSIGNMENT = {int(k): bool(v) for k, v in engine["assignment"].items()}
    return _original_projection(engine)


def make_experience_with_parent_witness(formula, budget, solved, generation):
    if solved.get("status") == "SAT" and solved.get("assignment") is None and _CAPTURED_ASSIGNMENT is not None:
        enriched = dict(solved)
        enriched["assignment"] = dict(_CAPTURED_ASSIGNMENT)
        return _original_make_experience(formula, budget, enriched, generation)
    return _original_make_experience(formula, budget, solved, generation)


def solve_spiral_strict(formula, budget, state_path):
    global _CAPTURED_ASSIGNMENT
    _CAPTURED_ASSIGNMENT = None
    return _original_solve_spiral(formula, budget, state_path)


base.gate.v2.engine_projection = capturing_projection
base.make_experience = make_experience_with_parent_witness
base.solve_spiral = solve_spiral_strict


def main() -> None:
    base.main()


if __name__ == "__main__":
    main()
