#!/usr/bin/env python3
"""Compatibility repair for BRIDGE-R1 Palomar day materialization.

The frozen S0 release at commit 4005e200... contains the candidate/tile products,
but its historical data/plate_manifest.csv does not yet carry RA/Dec columns.
Use the current public plate manifest for angular centers while keeping the S0
candidate/tile release pinned and recording the fetched plate-manifest SHA-256
inside the generated ledger.
"""
from __future__ import annotations

import materialize_palomar_day_ledger as base

base.PLATE_MANIFEST_URL = "https://raw.githubusercontent.com/jannefi/poss1-plate-slice/main/data/plate_manifest.csv"

if __name__ == "__main__":
    raise SystemExit(base.main())
