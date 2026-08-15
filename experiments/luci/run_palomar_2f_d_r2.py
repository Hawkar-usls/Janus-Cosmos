#!/usr/bin/env python3
from __future__ import annotations

"""Second infrastructure repair for 2F-D archive inventory completeness.

R1 proved adaptive disjoint MJD partitioning works, but the NULL-MJD layer
contains >=5000 public imaging rows. R2 keeps the MJD logic and retrieves the
NULL-MJD layer with explicit TAP MAXREC, rejecting any VOTable OVERFLOW.
No scientific decision rule or source/morphology gate is changed.
"""

import io
import urllib.parse
import urllib.request

from astropy.table import Table

from experiments.luci import run_palomar_2f_d as base
from experiments.luci import run_palomar_2f_d_r1 as r1

NULL_MAXREC = 250000


def _tap_query_complete_null_mjd():
    adql = (
        f"SELECT {r1._select_columns()}{r1._base_where()}"
        "AND lbt.luci.mjd_obs IS NULL"
    )
    body = urllib.parse.urlencode({
        "REQUEST": "doQuery",
        "LANG": "ADQL",
        "QUERY": adql,
        "FORMAT": "votable",
        "MAXREC": str(NULL_MAXREC),
    }).encode("utf-8")
    req = urllib.request.Request(
        base.TAP_SYNC,
        data=body,
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "User-Agent": "JANUS-COSMOS-LUCI-JPFM-2F-D-R2/1.0",
        },
    )
    with urllib.request.urlopen(req, timeout=300) as resp:
        raw = resp.read()
    upper = raw.upper()
    if b'VALUE="OVERFLOW"' in upper or b"VALUE='OVERFLOW'" in upper:
        raise RuntimeError(
            f"NULL-MJD TAP query overflowed explicit MAXREC={NULL_MAXREC}; completeness not proven"
        )
    tab = Table.read(io.BytesIO(raw), format="votable")
    if len(tab) >= NULL_MAXREC:
        raise RuntimeError(
            f"NULL-MJD TAP query reached MAXREC={NULL_MAXREC} without explicit overflow marker; fail closed"
        )
    return tab


def query_full_luci_inventory_mjd_r2() -> list[dict]:
    by_name: dict[str, dict] = {}
    stack: list[tuple[float, float]] = [(r1.ROOT_MJD_MIN, r1.ROOT_MJD_MAX)]
    split_count = 0
    leaf_windows = 0

    while stack:
        lo, hi = stack.pop()
        tab = r1._query_window(lo, hi)
        if len(tab) >= r1.TOP:
            width = hi - lo
            if width <= r1.MIN_WINDOW_DAYS:
                raise RuntimeError(
                    f"LUCI MJD window remains saturated at safety floor: [{lo},{hi}) rows={len(tab)}"
                )
            split_count += 1
            if split_count > r1.MAX_SPLITS:
                raise RuntimeError("LUCI MJD inventory exceeded adaptive split safety bound")
            mid = 0.5 * (lo + hi)
            stack.append((mid, hi))
            stack.append((lo, mid))
            continue

        leaf_windows += 1
        for row in tab:
            fname = base._s(row["file_name"])
            if fname:
                by_name[fname] = r1._row_to_record(row)

    null_tab = _tap_query_complete_null_mjd()
    for row in null_tab:
        fname = base._s(row["file_name"])
        if fname:
            by_name[fname] = r1._row_to_record(row)

    if not by_name:
        raise RuntimeError("LUCI public imaging inventory unexpectedly empty")

    rows = sorted(by_name.values(), key=lambda x: x["file_name"])
    query_full_luci_inventory_mjd_r2.last_diagnostics = {
        "strategy": "adaptive_disjoint_mjd_windows_plus_complete_null_maxrec",
        "root_mjd_min": r1.ROOT_MJD_MIN,
        "root_mjd_max": r1.ROOT_MJD_MAX,
        "top_per_mjd_query": r1.TOP,
        "leaf_window_count": leaf_windows,
        "split_count": split_count,
        "null_mjd_maxrec": NULL_MAXREC,
        "null_mjd_row_count": len(null_tab),
        "unique_file_count": len(rows),
        "overflow_rejected": True,
    }
    return rows


query_full_luci_inventory_mjd_r2.last_diagnostics = {}
base.query_full_luci_inventory = query_full_luci_inventory_mjd_r2


if __name__ == "__main__":
    raise SystemExit(base.main())
