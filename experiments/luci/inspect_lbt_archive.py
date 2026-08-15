from __future__ import annotations

import json

from build_luci_archive_manifest import _s, tap_query


def rows_to_dicts(table):
    return [{name: _s(row[name]) for name in table.colnames} for row in table]


def main() -> int:
    report = {}
    for table_name in ("lbt.luci", "lbt.lbt"):
        cols = tap_query(
            "SELECT column_name, datatype, description FROM TAP_SCHEMA.columns "
            f"WHERE table_name='{table_name}'"
        )
        report[table_name] = rows_to_dicts(cols)

    keys = tap_query(
        "SELECT * FROM TAP_SCHEMA.keys "
        "WHERE from_table='lbt.luci' OR target_table='lbt.luci' "
        "OR from_table='lbt.lbt' OR target_table='lbt.lbt'"
    )
    report["keys"] = rows_to_dicts(keys)

    key_columns = tap_query("SELECT * FROM TAP_SCHEMA.key_columns")
    relevant_key_ids = {r.get("key_id", "") for r in report["keys"]}
    report["key_columns"] = [
        row for row in rows_to_dicts(key_columns)
        if row.get("key_id", "") in relevant_key_ids
    ]

    samples = tap_query(
        "SELECT TOP 5 instrument, telescope, object, filters, gratname, imagetyp, "
        "file_name, dataprod, date_obs FROM lbt.luci"
    )
    report["luci_samples"] = rows_to_dicts(samples)

    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
