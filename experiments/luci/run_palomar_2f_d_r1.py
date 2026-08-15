#!/usr/bin/env python3
from __future__ import annotations

"""Infrastructure-only repair for 2F-D archive enumeration.

The scientific preregistration, full POSS-I population, crossmatch geometry,
exact-FITS chain, R1 source gate, matched controls, and claim ceiling are
unchanged. Only the public LUCI inventory enumeration strategy is replaced:
lexicographic ORDER BY pagination (rejected by the TAP service) becomes
adaptive disjoint MJD-window partitioning.
"""

import math

from experiments.luci import run_palomar_2f_d as base

TOP = 5000
ROOT_MJD_MIN = 0.0
ROOT_MJD_MAX = 100000.0
MIN_WINDOW_DAYS = 0.125
MAX_SPLITS = 256


def _select_columns() -> str:
    return (
        "lbt.luci.instrument,lbt.luci.telescope,lbt.luci.object,lbt.luci.filters,"
        "lbt.luci.gratname,lbt.luci.imagetyp,lbt.luci.file_name,lbt.luci.date_obs,"
        "lbt.luci.mjd_obs,lbt.luci.crval1,lbt.luci.crval2,lbt.luci.crpix1,lbt.luci.crpix2,"
        "lbt.luci.cd1_1,lbt.luci.cd1_2,lbt.luci.cd2_1,lbt.luci.cd2_2,"
        "lbt.luci.ctype1,lbt.luci.ctype2,lbt.luci.naxis1,lbt.luci.naxis2,"
        "lbt.luci.pixscale,lbt.lbt.file_url,lbt.lbt.policy"
    )


def _base_where() -> str:
    return (
        " FROM lbt.luci JOIN lbt.lbt ON lbt.luci.file_name=lbt.lbt.file_name "
        "WHERE lbt.lbt.policy='FREE' AND lbt.luci.imagetyp='SCIENCE' "
        "AND (lbt.luci.gratname='Mirror' OR lbt.luci.gratname='mirror') "
    )


def _row_to_record(row) -> dict:
    return {
        "file_name": base._s(row["file_name"]),
        "file_url": base._s(row["file_url"]),
        "instrument": base._s(row["instrument"]),
        "telescope": base._s(row["telescope"]),
        "target": base._s(row["object"]),
        "filters": base._s(row["filters"]),
        "date_obs": base._s(row["date_obs"]),
        "crval1": base._f(row["crval1"]),
        "crval2": base._f(row["crval2"]),
        "crpix1": base._f(row["crpix1"]),
        "crpix2": base._f(row["crpix2"]),
        "cd1_1": base._f(row["cd1_1"]),
        "cd1_2": base._f(row["cd1_2"]),
        "cd2_1": base._f(row["cd2_1"]),
        "cd2_2": base._f(row["cd2_2"]),
        "ctype1": base._s(row["ctype1"]),
        "ctype2": base._s(row["ctype2"]),
        "naxis1": base._f(row["naxis1"]),
        "naxis2": base._f(row["naxis2"]),
        "pixscale": base._f(row["pixscale"]),
    }


def _query_window(lo: float, hi: float):
    adql = (
        f"SELECT TOP {TOP} {_select_columns()}{_base_where()}"
        f"AND lbt.luci.mjd_obs >= {lo:.9f} AND lbt.luci.mjd_obs < {hi:.9f}"
    )
    return base.tap_query(adql, timeout=240)


def query_full_luci_inventory_mjd() -> list[dict]:
    """Enumerate all public imaging rows by recursively splitting saturated MJD windows.

    Every accepted leaf window contains fewer than TOP rows, so no result is
    truncated inside that disjoint time interval. NULL-MJD rows are queried
    separately and must also be below TOP or the run fails closed.
    """
    by_name: dict[str, dict] = {}
    stack: list[tuple[float, float]] = [(ROOT_MJD_MIN, ROOT_MJD_MAX)]
    split_count = 0
    leaf_windows = 0

    while stack:
        lo, hi = stack.pop()
        tab = _query_window(lo, hi)
        if len(tab) >= TOP:
            width = hi - lo
            if width <= MIN_WINDOW_DAYS:
                raise RuntimeError(
                    f"LUCI MJD window remains saturated at safety floor: [{lo},{hi}) rows={len(tab)}"
                )
            split_count += 1
            if split_count > MAX_SPLITS:
                raise RuntimeError("LUCI MJD inventory exceeded adaptive split safety bound")
            mid = 0.5 * (lo + hi)
            stack.append((mid, hi))
            stack.append((lo, mid))
            continue

        leaf_windows += 1
        for row in tab:
            fname = base._s(row["file_name"])
            if fname:
                by_name[fname] = _row_to_record(row)

    null_adql = (
        f"SELECT TOP {TOP} {_select_columns()}{_base_where()}"
        "AND lbt.luci.mjd_obs IS NULL"
    )
    null_tab = base.tap_query(null_adql, timeout=240)
    if len(null_tab) >= TOP:
        raise RuntimeError("NULL-MJD LUCI inventory is saturated; cannot prove completeness")
    for row in null_tab:
        fname = base._s(row["file_name"])
        if fname:
            by_name[fname] = _row_to_record(row)

    if not by_name:
        raise RuntimeError("LUCI public imaging inventory unexpectedly empty")

    rows = sorted(by_name.values(), key=lambda x: x["file_name"])
    # Add provenance only to runtime diagnostics, not to frozen CSV columns.
    query_full_luci_inventory_mjd.last_diagnostics = {
        "strategy": "adaptive_disjoint_mjd_windows",
        "root_mjd_min": ROOT_MJD_MIN,
        "root_mjd_max": ROOT_MJD_MAX,
        "top_per_query": TOP,
        "leaf_window_count": leaf_windows,
        "split_count": split_count,
        "null_mjd_row_count": len(null_tab),
        "unique_file_count": len(rows),
    }
    return rows


query_full_luci_inventory_mjd.last_diagnostics = {}
base.query_full_luci_inventory = query_full_luci_inventory_mjd


if __name__ == "__main__":
    raise SystemExit(base.main())
