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


def discover_luci_table() -> tuple[str, list[str]]:
    tables = tap_query("SELECT table_name, description FROM TAP_SCHEMA.tables")
    names = [_s(row["table_name"]) for row in tables]
    table_name = next((name for name in names if name.lower() == "lbt.luci"), None)
    if table_name is None:
        table_name = next((name for name in names if "luci" in name.lower()), None)
    if table_name is None:
        raise RuntimeError(f"No LUCI table exposed by LBT TAP. Tables: {names[:50]}")
    columns_table = tap_query(
        "SELECT column_name FROM TAP_SCHEMA.columns "
        f"WHERE table_name='{table_name.replace(chr(39), chr(39)*2)}'"
    )
    columns = [_s(row["column_name"]) for row in columns_table]
    return table_name, columns


def _quote_adql(value: str) -> str:
    return value.replace("'", "''")


def query_luci_rows(*, limit: int, target: str | None = None) -> Table:
    where = [
        "lbt.luci.imagetyp='SCIENCE'",
        "(lbt.luci.gratname='Mirror' OR lbt.luci.gratname='mirror')",
        "lbt.lbt.policy='FREE'",
    ]
    if target:
        where.append(f"lbt.luci.object='{_quote_adql(target)}'")
    adql = (
        f"SELECT TOP {int(limit)} "
        "lbt.luci.instrument, lbt.luci.telescope, lbt.luci.object, "
        "lbt.luci.filters, lbt.luci.filter1, lbt.luci.filter2, "
        "lbt.luci.gratname, lbt.luci.imagetyp, lbt.luci.file_name, "
        "lbt.luci.date_obs, lbt.lbt.file_url, lbt.lbt.policy "
        "FROM lbt.luci JOIN lbt.lbt ON lbt.luci.file_name=lbt.lbt.file_name "
        "WHERE " + " AND ".join(where)
    )
    return tap_query(adql)


def _science_filter(row) -> str:
    f1 = _s(row["filter1"])
    f2 = _s(row["filter2"])
    excluded = {"", "clear", "blind", "open", "none"}
    if f2.lower() not in excluded:
        return f2
    if f1.lower() not in excluded:
        return f1
    combined = _s(row["filters"])
    return combined or "UNKNOWN_FILTER"


def build_manifest(rows: Table, *, target_requested: str | None, max_targets: int = 5) -> dict:
    grouped: dict[str, list[dict]] = defaultdict(list)
    rejected = 0
    for row in rows:
        target = _s(row["object"])
        url = _s(row["file_url"])
        policy = _s(row["policy"])
        grating = _s(row["gratname"])
        instrument = _s(row["instrument"])
        filt = _science_filter(row)
        if not target or target.upper() in {"NOTARGET", "UNKNOWN", "UNKNOWN_TARGET"}:
            rejected += 1
            continue
        if not url.startswith(("http://", "https://")) or policy.upper() != "FREE":
            rejected += 1
            continue
        if grating.lower() != "mirror":
            rejected += 1
            continue
        if instrument.upper() not in {"LUCI1", "LUCI2", "LUCI", "LUCIFER1", "LUCIFER2", "LUCIFER"}:
            rejected += 1
            continue
        grouped[target].append({
            "filter": filt,
            "band": filt,
            "instrument": instrument,
            "url": url,
            "archive_file_name": _s(row["file_name"]),
            "archive_date_obs": _s(row["date_obs"]),
            "archive_policy": policy,
            "archive_grating": grating,
            "archive_telescope": _s(row["telescope"]),
        })

    targets = []
    for target in sorted(grouped):
        # Sorting before de-duplication makes selection deterministic even if TAP
        # returns rows in a different physical order.
        items = sorted(
            grouped[target],
            key=lambda x: (x["filter"], x["archive_date_obs"], x["archive_file_name"]),
        )
        by_filter = {}
        for item in items:
            by_filter.setdefault(item["filter"], item)
        distinct = list(by_filter.values())
        if len(distinct) < 2:
            continue
        targets.append({"target": target, "class": "archive_unspecified", "filters": distinct})
        if len(targets) >= max_targets:
            break

    if not targets:
        scope = f" target={target_requested!r}" if target_requested else ""
        raise RuntimeError(
            "No LUCI imaging target with >=2 distinct downloadable FREE bands was found"
            f" in the queried rows.{scope}"
        )

    return {
        "schema": "janus.cosmos.luci.archive_manifest.v1",
        "source_archive": "LBT Archive / IA2 TAP",
        "tap_service": TAP_SYNC,
        "tap_tables": ["lbt.luci", "lbt.lbt"],
        "join_key": "file_name",
        "selection": (
            "FREE LUCI science rows with GRATNAME=Mirror only; deterministic first file per distinct filter; "
            "targets require at least two distinct imaging filters."
        ),
        "target_requested": target_requested,
        "instrument_scope": ["LUCI1", "LUCI2", "legacy LUCIFER naming"],
        "row_count_from_tap": len(rows),
        "row_count_rejected_pre_manifest": rejected,
        "targets": targets,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Build a public LUCI/LUCIFER imaging manifest from the official LBT TAP archive")
    ap.add_argument("--limit", type=int, default=5000)
    ap.add_argument("--max-targets", type=int, default=5)
    ap.add_argument("--target", default="", help="Optional exact LBT OBJECT name; useful for a deterministic infrastructure smoke test")
    ap.add_argument("--output", default="data/runtime/luci/luci_archive_manifest.json")
    ap.add_argument("--metadata-only", action="store_true")
    args = ap.parse_args()

    table_name, columns = discover_luci_table()
    report = {"tap": TAP_SYNC, "table": table_name, "column_count": len(columns), "columns": columns}
    print(json.dumps(report, indent=2))
    if args.metadata_only:
        return 0

    target = args.target.strip() or None
    rows = query_luci_rows(limit=args.limit, target=target)
    manifest = build_manifest(rows, target_requested=target, max_targets=args.max_targets)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(out),
        "target_count": len(manifest["targets"]),
        "targets": [
            {"target": t["target"], "filters": [x["filter"] for x in t["filters"]]}
            for t in manifest["targets"]
        ],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
