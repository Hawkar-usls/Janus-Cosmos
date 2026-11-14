#!/usr/bin/env python3
"""Build a fail-closed, source-bound JANUS COSMOS FERMI-V7 browser dataset.

V7 consumes V4's previously extracted event-level diagnostics, NOT Table 7.
It audits the archival JSON against the pinned FITS manifest, but it does NOT
re-read FITS unless --raw-dir is explicitly provided. Internal row checks are
not independent replication and no LIV / source-emission constraint is fit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parents[1]
V4_MANIFEST = ROOT / "data/cosmos/COSMOS-ALI-V4-FERMI-RAW-FITS-MANIFEST-2026-10-18.json"
V4_EVENTS = ROOT / "data/cosmos/COSMOS-ALI-V4-FERMI-RAW-EVENT-DIAGNOSTICS-2026-10-18.json"
V6_GATES = ROOT / "data/cosmos/COSMOS-ALI-V6-JOINT-OBSERVATION-READINESS-2026-10-19.json"
EXPECTED = ("090510", "090902B", "130427A")
SCHEMA = "janus.cosmos.fermi.v7.event-evidence-browser.v1"


class EvidenceError(ValueError):
    pass


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise EvidenceError("JSON is not an object: " + str(path))
    return data


def assert_(condition: bool, reason: str) -> None:
    if not condition:
        raise EvidenceError(reason)


def number(value, key: str, *, minimum=None, maximum=None) -> float:
    assert_(not isinstance(value, bool) and isinstance(value, (int, float)), "invalid number: " + key)
    n = float(value)
    assert_(math.isfinite(n), "nonfinite " + key)
    if minimum is not None:
        assert_(n >= minimum, "below bound: " + key)
    if maximum is not None:
        assert_(n <= maximum, "above bound: " + key)
    return n


def trusted_fssc_url(url: str) -> bool:
    p = urlparse(url)
    return (
        p.scheme == "https"
        and p.netloc == "fermi.gsfc.nasa.gov"
        and not (p.params or p.query or p.fragment or p.username or p.password)
        and re.fullmatch(r"/FTP/fermi/data/lat/queries/L[A-Z0-9]+_(EV00|SC00)\.fits", p.path) is not None
    )


def canonical_event_ids(source: dict) -> list[tuple[int, int]]:
    return sorted((int(x["run_id"]), int(x["event_id"])) for x in source["events"])


def audited_snapshot(manifest: dict, raw: dict, v6: dict, raw_dir: Path | None = None) -> dict:
    assert_(manifest.get("schema") == "janus.cosmos.fermi-lat-fits-manifest.v4", "unexpected FITS manifest schema")
    assert_(raw.get("schema") == "janus.cosmos.fermi-lat-extended-raw-ingest.v4", "unexpected V4 event schema")
    assert_(v6.get("schema") == "janus.cosmos.fermi-lat-gbm-cross-instrument-readiness.v6", "unexpected V6 gate schema")
    assert_(raw.get("gates", {}).get("LIV_OBSERVED") == "NOT_ESTABLISHED", "historical physical ceiling changed")
    assert_(raw.get("gates", {}).get("LAT_EVENT_ASSOCIATION") == "NOT_ESTIMATED", "association status changed")
    assert_(raw.get("gates", {}).get("IRF_FORWARD_FOLDING") == "NOT_IMPLEMENTED", "instrument gate changed")
    assert_(v6.get("scientific_gates", {}).get("OBSERVED_MULTIBURST_LIV_LRT", "").startswith("BLOCKED"), "V6 science gate changed")
    calc_manifest_digest = hashlib.sha256(
        json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    assert_(calc_manifest_digest == raw.get("manifest_sha256"), "V4 source manifest link mismatch")
    configs = manifest.get("bursts", [])
    sources = raw.get("sources", [])
    assert_(len(configs) == len(sources) == len(EXPECTED), "unexpected number of source bursts")
    assert_(tuple(x.get("grb") for x in configs) == EXPECTED, "manifest cohort is not frozen")
    assert_(tuple(x.get("grb") for x in sources) == EXPECTED, "V4 events cohort is not frozen")
    assert_(raw.get("counts", {}).get("bursts") == len(EXPECTED), "aggregate burst count mismatch")
    assert_(manifest.get("claim_ceiling") == "NO_OBSERVATIONAL_LIV_INFERENCE_WITHOUT_IRF_GBM_SOURCE_LAG_MODEL",
            "source claim ceiling missing")
    assert_(len({s["query_id"] for s in configs}) == len(configs), "duplicate query ID")
    assert_(all(s["selection"] == manifest["common_selection"] for s in configs), "manifest selection changed")
    assert_(all(s["selection"] == manifest["common_selection"] for s in sources), "V4 selection changed")
    assert_(all(s.get("trigger_time_not_emission_time") is True for s in sources),
            "trigger misrepresented as emission time")
    total_raw, total_selected = 0, 0
    rows = []
    view_sources = []
    for config, source in zip(configs, sources):
        grb = config["grb"]
        assert_(config["query_id"] == source["query_id"], "query ID mismatch: " + grb)
        assert_(config["class"] == source["class"], "GRB class mismatch: " + grb)
        assert_(str(config["redshift"]) == str(source["redshift"]), "redshift mismatch: " + grb)
        trigger = float(config["trigger_met_estimated"])
        assert_(math.isfinite(trigger) and trigger > 0 and abs(trigger - float(source["trigger_met_estimated"])) < 1e-8,
                "estimated GBM trigger changed: " + grb)
        prov = source.get("fits_provenance", {})
        for role in ("event", "spacecraft"):
            m = config[role]
            p = prov.get(role, {})
            assert_(trusted_fssc_url(m["url"]), "untrusted FITS source URL")
            assert_(m["url"] == p.get("source_url") and m["sha256"] == p.get("sha256") and m["bytes"] == p.get("bytes"),
                    "archived FITS provenance mismatch: " + grb + "/" + role)
            assert_(re.fullmatch(r"[a-f0-9]{64}", m["sha256"]) is not None, "invalid FITS sha256")
            assert_(m["filename"] == grb + ("_EV00.fits" if role == "event" else "_SC00.fits"),
                    "FITS local path mismatch")
        counts = source["counts"]
        actual = source["events"]
        assert_(len(actual) == counts["after_diagnostic_selection"], "selected count mismatch: " + grb)
        count_rejected = counts["rejected_by_first_failed_cut"]
        assert_(set(count_rejected) == {"not_gti","bad_spacecraft_or_gap","class","energy","zenith","theta","roi"},
                "unexpected rejection accounting categories")
        assert_(all(isinstance(x,int) and x>=0 for x in count_rejected.values()), "invalid rejection count")
        assert_(len(actual) + sum(count_rejected.values()) == counts["events_in_archive_roi"],
                "ROI rejection accounting mismatch: " + grb)
        ordered = [(number(e["met_s"], "met", minimum=0), int(e["run_id"]), int(e["event_id"])) for e in actual]
        assert_(ordered == sorted(ordered), "event order not monotonically sorted: " + grb)
        assert_(len(set(canonical_event_ids(source))) == len(actual), "duplicate event identity: " + grb)
        for event in actual:
            energy = number(event["energy_MeV"], "energy", minimum=1000, maximum=300000)
            seconds = number(event["after_trigger_s"], "trigger-relative time")
            met = number(event["met_s"], "MET", minimum=0)
            assert_(abs((met-trigger)-seconds) <= 0.000003, "event time/trigger inconsistency: " + grb)
            angular = number(event["angular_sep_deg"], "angular ROI", minimum=0,
                             maximum=config["selection"]["max_roi_deg"]+0.000001)
            theta = number(event["theta_deg"], "incidence", minimum=0,
                           maximum=config["selection"]["max_theta_deg"]+0.000001)
            zenith = number(event["zenith_deg"], "zenith", minimum=0,
                            maximum=config["selection"]["max_zenith_deg"]+0.000001)
            mask = event["event_class_mask"]
            assert_(isinstance(mask,int) and (mask & 128) != 0, "event lost SOURCE bit")
            event_type = event["event_type_mask"]
            assert_(isinstance(event_type,int) and event_type>0, "event type invalid")
            run, eid = event["run_id"], event["event_id"]
            assert_(type(run) is int and type(eid) is int and run>=0 and eid>=0, "invalid run/event identity")
            rows.append({"grb":grb,"run_id":run,"event_id":eid,
                         "met_s":met,"after_trigger_s":seconds,"energy_MeV":energy,
                         "angular_sep_deg":angular,"theta_deg":theta,"zenith_deg":zenith,
                         "event_class_mask":mask,"event_type_mask":event_type})
        total_raw += counts["events_in_archive_roi"]
        total_selected += len(actual)
        view_sources.append({
            "grb":grb,"class":source["class"],"redshift":source["redshift"],
            "query_id":source["query_id"],"trigger_met_estimated":config["trigger_met_estimated"],
            "counts":{"events_in_query_roi":counts["events_in_archive_roi"],
                      "diagnostic_selected":len(actual)},
            "selection":config["selection"],
            "event_file":{"url":config["event"]["url"],"sha256":config["event"]["sha256"],
                          "bytes":config["event"]["bytes"],"filename":config["event"]["filename"]},
            "spacecraft_file":{"url":config["spacecraft"]["url"],"sha256":config["spacecraft"]["sha256"],
                               "bytes":config["spacecraft"]["bytes"],"filename":config["spacecraft"]["filename"]},
        })
    assert_(total_raw == raw["counts"]["raw_events_in_query_rois"], "aggregate input events mismatch")
    assert_(total_selected == raw["counts"]["diagnostic_source_region_candidates"], "aggregate selected mismatch")
    assert_(total_selected == 175 and total_raw == 6520, "unexpected historical V4 cohort")
    # Frozen event identity ordering is insensitive to reordering entire input sources and rows.
    # This is an INTERNAL integrity control, not an independent scientific experiment.
    canon = json.dumps(sorted(rows,key=lambda x:(x["grb"],x["run_id"],x["event_id"])),
                       sort_keys=True,separators=(",",":")).encode("utf-8")
    digest = hashlib.sha256(canon).hexdigest()

    raw_status = "RAW_FITS_NOT_RECHECKED_THIS_BUILD"
    if raw_dir is not None:
        raw_dir = raw_dir.resolve()
        for config in configs:
            for role in ("event", "spacecraft"):
                m = config[role]
                candidate = raw_dir / m["filename"]
                assert_(candidate.is_file() and candidate.stat().st_size == m["bytes"] and
                        sha256_file(candidate) == m["sha256"], "missing / mismatched pinned raw FITS: " + m["filename"])
        # This is an internal replay of the pre-existing V4 parser. Not a second independent parser.
        sys.path.insert(0, str(ROOT / "workspace"))
        from janus_cosmos_fermi_raw_v4 import run as replay_v4
        replay = replay_v4(manifest, raw_dir)
        assert_(replay["counts"] == raw["counts"], "raw FITS replay aggregate mismatch")
        for computed, saved in zip(replay["sources"], sources):
            assert_(computed["grb"] == saved["grb"] and
                    computed["counts"] == saved["counts"] and
                    computed["events"] == saved["events"], "raw FITS replay changed event records")
        raw_status = "PINNED_FITS_BYTES_AND_INTERNAL_V4_REPLAY_MATCHED_THIS_BUILD"

    return {
        "schema":SCHEMA,
        "dataset_label":"FERMI-V7 SOURCE-BOUND EVENT INSPECTOR / DERIVED FROM V4",
        "data_class":"DERIVED_ARCHIVAL_EVENT_DIAGNOSTICS_NOT_NEW_INDEPENDENT_OBSERVATIONS",
        "source": {
            "v4_manifest_path":str(V4_MANIFEST.relative_to(ROOT)),
            "v4_event_diagnostics_path":str(V4_EVENTS.relative_to(ROOT)),
            "v6_readiness_path":str(V6_GATES.relative_to(ROOT)),
            "v4_manifest_file_sha256":sha256_file(V4_MANIFEST),
            "v4_event_diagnostics_file_sha256":sha256_file(V4_EVENTS),
            "v6_readiness_file_sha256":sha256_file(V6_GATES),
            "fits_manifest_canonical_sha256":calc_manifest_digest,
            "event_rows_canonical_sha256":digest,
            "nasa_fssc_data_help":"https://fermi.gsfc.nasa.gov/ssc/LATDataQuery_help.html",
            "nasa_event_columns":"https://fermi.gsfc.nasa.gov/ssc/data/access/lat/lat_data_columns.html"
        },
        "counts":{"grbs":len(EXPECTED),"archive_query_roi_events":total_raw,
                  "diagnostic_selected":total_selected,
                  "note":"175 diagnostic candidates are NOT 175 confirmed GRB-associated photons"},
        "gates":{
            "archived_manifest_link":"PASS_RECOMPUTED_FROM_JSON",
            "event_identity_and_rejection_accounting":"PASS_RECOMPUTED_FROM_JSON",
            "source_order_and_event_id_audit":"PASS_RECOMPUTED_FROM_JSON",
            "source_byte_recheck_this_build":raw_status,
            "independent_external_replication":"NOT_DONE",
            "full_LAT_PSF_EDISP_background_forward_fold":"NOT_DONE",
            "GBM_source_emission_offset_calibration":"NOT_DONE",
            "multi_burst_LIV_inference":"BLOCKED",
            "observational_LIV_p_value":None,
            "lattice_spacing_fit":None,
            "retrocausality_detection":"NOT_ESTABLISHED",
        },
        "scientific_warning":
          "Event TIME/MET is arrival at LAT; time relative to estimated GBM trigger is NOT emission time. "
          "Diagnostic SOURCE-class event selection is NOT event-by-event source association. "
          "No effective-area/PSF/EDISP, diffuse background and source-lag model is applied, "
          "so no astrophysical LIV constraint, p-value or lattice-spacing estimate is authorized.",
        "sources":view_sources,
        "events":rows,
    }


def main() -> int:
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--out",required=True,type=Path)
    cli.add_argument("--raw-dir",type=Path,default=None,
                     help="Optional exact six-file V4 FITS directory for SHA256+legacy parser replay")
    args = cli.parse_args()
    try:
        data = audited_snapshot(load(V4_MANIFEST),load(V4_EVENTS),load(V6_GATES),args.raw_dir)
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(data,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    except (OSError,ValueError,KeyError,TypeError,ImportError) as e:
        print("FERMI-V7 BLOCKED: "+str(e),file=sys.stderr)
        return 2
    print("FERMI-V7 JSON AUDIT PASS: %s events / %s query ROI rows; raw=%s" %
          (data["counts"]["diagnostic_selected"],data["counts"]["archive_query_roi_events"],
           data["gates"]["source_byte_recheck_this_build"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
