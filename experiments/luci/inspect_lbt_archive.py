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

    report["gratname_values"] = rows_to_dicts(
        tap_query("SELECT DISTINCT gratname FROM lbt.luci")
    )
    report["optic_values"] = rows_to_dicts(
        tap_query("SELECT DISTINCT optic FROM lbt.luci")
    )
    report["science_mode_samples"] = rows_to_dicts(
        tap_query(
            "SELECT TOP 30 instrument, telescope, object, filters, filter1, filter2, "
            "gratname, optic, imagetyp, file_name, date_obs FROM lbt.luci "
            "WHERE imagetyp='SCIENCE'"
        )
    )
    report["joined_url_samples"] = rows_to_dicts(
        tap_query(
            "SELECT TOP 5 lbt.luci.file_name, lbt.lbt.file_url, lbt.lbt.policy "
            "FROM lbt.luci JOIN lbt.lbt ON lbt.luci.file_name=lbt.lbt.file_name"
        )
    )

    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
