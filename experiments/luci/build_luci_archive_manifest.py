from __future__ import annotations

import argparse
import io
import json
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

from astropy.table import Table

TAP_SYNC = "https://archive.lbto.org/tap/sync"


def tap_query(adql: str, *, timeout: int = 90) -> Table:
    payload = urllib.parse.urlencode({
        "REQUEST": "doQuery",
        "LANG": "ADQL",
        "QUERY": adql,
        "FORMAT": "votable",
    }).encode("utf-8")
    req = urllib.request.Request(
        TAP_SYNC,
        data=payload,
        headers={"User-Agent": "Janus-Cosmos-LUCI/1.0", "Content-Type": "application/x-www-form-urlencoded"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as response:
        raw = response.read()
    return Table.read(io.BytesIO(raw), format="votable")


def _s(value) -> str:
    if value is None:
        return ""
    try:
        if getattr(value, "mask", False):
            return ""
    except Exception:
        pass
    return str(value).strip()


def _find_column(columns: list[str], *candidates: str) -> str | None:
    by_lower = {c.lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in by_lower:
            return by_lower[candidate.lower()]
    return None


def discover_obscore() -> tuple[str, list[str]]:
    tables = tap_query("SELECT table_name, description FROM TAP_SCHEMA.tables")
    names = [_s(row["table_name"]) for row in tables]
    obscore = next((name for name in names if name.lower() == "ivoa.obscore"), None)
    if obscore is None:
        obscore = next((name for name in names if "obscore" in name.lower()), None)
    if obscore is None:
        raise RuntimeError(f"No ObsCore-like table exposed by LBT TAP. Tables: {names[:50]}")
    columns_table = tap_query(
        "SELECT column_name FROM TAP_SCHEMA.columns "
        f"WHERE table_name='{obscore.replace(chr(39), chr(39)*2)}'"
    )
    columns = [_s(row["column_name"]) for row in columns_table]
    return obscore, columns


def query_luci_rows(table_name: str, columns: list[str], *, limit: int) -> Table:
    instrument = _find_column(columns, "instrument_name", "instrument", "instrume")
    if instrument is None:
        raise RuntimeError("ObsCore table has no recognizable instrument column")
    where = f"({instrument} LIKE 'LUCI%' OR {instrument} LIKE 'LUCIFER%')"
    rights = _find_column(columns, "data_rights", "rights")
    if rights:
        where += f" AND ({rights}='public' OR {rights}='PUBLIC' OR {rights} IS NULL)"
    dataproduct = _find_column(columns, "dataproduct_type")
    if dataproduct:
        where += f" AND ({dataproduct}='image' OR {dataproduct}='IMAGE')"
    return tap_query(f"SELECT TOP {int(limit)} * FROM {table_name} WHERE {where}")


def build_manifest(rows: Table, *, table_name: str, max_targets: int = 5) -> dict:
    columns = list(rows.colnames)
    instrument_col = _find_column(columns, "instrument_name", "instrument", "instrume")
    access_col = _find_column(columns, "access_url", "access_reference", "url")
    target_col = _find_column(columns, "target_name", "object", "object_name", "obs_target")
    filter_col = _find_column(columns, "filter", "filter_name", "bandpass", "obs_bandpass")
    obsid_col = _find_column(columns, "obs_id", "observation_id", "filename", "file_name")
    em_min_col = _find_column(columns, "em_min")
    em_max_col = _find_column(columns, "em_max")

    missing = [name for name, col in (("instrument", instrument_col), ("access_url", access_col), ("target", target_col)) if col is None]
    if missing:
        raise RuntimeError(f"Cannot build downloadable LUCI manifest; missing columns: {missing}. Available: {columns}")

    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        target = _s(row[target_col]) or "UNKNOWN_TARGET"
        url = _s(row[access_col])
        if not url.startswith(("http://", "https://")):
            continue
        instrument = _s(row[instrument_col])
        if filter_col:
            filt = _s(row[filter_col]) or "UNKNOWN_FILTER"
        else:
            lo = _s(row[em_min_col]) if em_min_col else ""
            hi = _s(row[em_max_col]) if em_max_col else ""
            filt = f"EM_{lo}_{hi}" if (lo or hi) else "UNKNOWN_FILTER"
        grouped[target].append({
            "filter": filt,
            "band": filt,
            "instrument": instrument,
            "url": url,
            "archive_obs_id": _s(row[obsid_col]) if obsid_col else "",
        })

    targets = []
    for target in sorted(grouped):
        by_filter = {}
        for item in grouped[target]:
            by_filter.setdefault(item["filter"], item)
        items = list(by_filter.values())
        if len(items) < 2:
            continue
        targets.append({"target": target, "class": "archive_unspecified", "filters": items[:2]})
        if len(targets) >= max_targets:
            break

    if not targets:
        raise RuntimeError("LUCI rows were found, but no target had two downloadable distinct imaging bands in this query window")

    return {
        "schema": "janus.cosmos.luci.archive_manifest.v1",
        "source_archive": "LBT Archive / IA2 TAP",
        "tap_service": TAP_SYNC,
        "tap_table": table_name,
        "selection": "Public LUCI/LUCIFER imaging rows; first deterministic targets with >=2 distinct bands in TAP result.",
        "instrument_scope": ["LUCI1", "LUCI2", "legacy LUCIFER naming"],
        "targets": targets,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Discover public LUCI/LUCIFER imaging rows from the official LBT TAP archive")
    ap.add_argument("--limit", type=int, default=500)
    ap.add_argument("--max-targets", type=int, default=5)
    ap.add_argument("--output", default="data/runtime/luci/luci_archive_manifest.json")
    ap.add_argument("--metadata-only", action="store_true")
    args = ap.parse_args()

    table_name, columns = discover_obscore()
    report = {"tap": TAP_SYNC, "table": table_name, "column_count": len(columns), "columns": columns}
    print(json.dumps(report, indent=2))
    if args.metadata_only:
        return 0

    rows = query_luci_rows(table_name, columns, limit=args.limit)
    manifest = build_manifest(rows, table_name=table_name, max_targets=args.max_targets)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out} with {len(manifest['targets'])} targets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
