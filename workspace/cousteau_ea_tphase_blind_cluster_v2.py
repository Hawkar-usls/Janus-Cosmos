#!/usr/bin/env python3
from __future__ import annotations

import argparse, hashlib, io, json, re, tarfile, zipfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

import requests
from bs4 import BeautifulSoup

from cousteau_ea_tphase_blind_cluster import (
    EXPECTED_EVENT_COUNT, DBSCAN_GRID, NULL_DOMAIN, NULL_SAMPLES, NULL_SEED,
    blind_cluster, reveal_and_score, parse_catalog, sha256_bytes,
)

DOI = "10.26022/IEDA/330497"
DATASET = "EA_Hydroacoustics"
DATASET_UID = "30497"
LANDING = "https://www.marine-geo.org/tools/files/30497"
SEARCH = "https://www.marine-geo.org/services/search/datasets"
FILESERVER = "https://www.marine-geo.org/services/FileServer"
FDS = "https://www.marine-geo.org/services/FileDownloadServer"
ADS = "https://www.marine-geo.org/services/ArchiveDownloadServer"


def brief(r):
    return {"url": r.url, "status": r.status_code, "bytes": len(r.content),
            "content_type": r.headers.get("content-type"),
            "prefix": r.text[:700] if "text" in (r.headers.get("content-type") or "") or len(r.content)<5000 else None}


def try_parse_catalog_bytes(raw: bytes, label: str, attempts: list):
    if len(raw) < 1000: return None
    low = raw[:1000].lower()
    if b"<html" in low or b"<!doctype" in low: return None
    try:
        coords, meta = parse_catalog(raw)
    except Exception as e:
        attempts.append({"candidate": label, "parse_error": f"{type(e).__name__}: {e}", "bytes": len(raw)})
        return None
    n = len(coords)
    # A real event catalog must be substantial; exact 6843 is preferred but not required for acquisition.
    if n < 1000:
        attempts.append({"candidate": label, "rejected": "too_few_valid_coordinates", "valid_coordinates": n})
        return None
    return raw, coords, meta, label


def try_archive(raw: bytes, label: str, attempts: list):
    members = []
    try:
        if raw[:2] == b"PK":
            with zipfile.ZipFile(io.BytesIO(raw)) as z:
                for name in z.namelist():
                    if name.endswith("/"): continue
                    b = z.read(name)
                    members.append({"name": name, "bytes": len(b)})
                    got = try_parse_catalog_bytes(b, label+"::"+name, attempts)
                    if got: return got, members
        else:
            with tarfile.open(fileobj=io.BytesIO(raw), mode="r:*") as t:
                for m in t.getmembers():
                    if not m.isfile(): continue
                    f = t.extractfile(m)
                    if not f: continue
                    b = f.read()
                    members.append({"name": m.name, "bytes": len(b)})
                    got = try_parse_catalog_bytes(b, label+"::"+m.name, attempts)
                    if got: return got, members
    except Exception as e:
        attempts.append({"archive": label, "archive_error": f"{type(e).__name__}: {e}", "bytes": len(raw)})
    return None, members


def candidate_urls_from_text(text: str, base: str):
    out = []
    for u in re.findall(r'(?:https?://|/)[^"\'<>\s]+', text):
        u = u.replace("&amp;", "&").rstrip(")]}>,;")
        low = u.lower()
        if "filedownloadserver" in low or "archivedownloadserver" in low or low.endswith((".txt",".dat",".csv",".tsv",".zip",".tar",".tgz",".gz")):
            out.append(urljoin(base, u))
    return list(dict.fromkeys(out))


def acquire(session: requests.Session):
    attempts = []
    candidates = []

    # 1) Search API with several documented/legacy-friendly query spellings.
    query_variants = [
        {"query": DATASET}, {"q": DATASET}, {"term": DATASET}, {"id": DATASET},
        {"doi": DOI}, {"query": DOI}, {"data_set_uid": DATASET_UID},
    ]
    for params in query_variants:
        try:
            r = session.get(SEARCH, params=params, timeout=45)
            attempts.append({"stage": "search_api", "params": params, **brief(r)})
            candidates += candidate_urls_from_text(r.text, r.url)
            # Recursively mine obvious UID-ish values from JSON/XML/text.
            for m in re.finditer(r'(?:data[_-]?uid|file[_-]?uid|uid)["\'\s:=]+([A-Za-z0-9._-]+)', r.text, re.I):
                candidates.append(FDS + "?data_uid=" + m.group(1))
        except Exception as e:
            attempts.append({"stage": "search_api", "params": params, "error": str(e)})

    # 2) Non-spatial FileServer query: catalogs may have no file boundary geometry.
    for params in [
        {"format":"full_info", "data_type":"Earthquake:Catalog:Microseismicity"},
        {"format":"full_info", "data_type":"Seismic:Passive"},
        {"format":"full_info", "data_set":DATASET},
        {"format":"full_info", "dataset":DATASET},
    ]:
        try:
            r = session.get(FILESERVER, params=params, timeout=60)
            attempts.append({"stage":"fileserver_nonspatial", "params":params, **brief(r)})
            candidates += candidate_urls_from_text(r.text, r.url)
            for m in re.finditer(r'(?:data[_-]?uid|file[_-]?uid|uid)["\'\s:=]+([A-Za-z0-9._-]+)', r.text, re.I):
                candidates.append(FDS + "?data_uid=" + m.group(1))
        except Exception as e:
            attempts.append({"stage":"fileserver_nonspatial", "params":params, "error":str(e)})

    # 3) Landing page: inspect hrefs, scripts, inputs and download-like forms.
    try:
        r = session.get(LANDING, timeout=60)
        attempts.append({"stage":"landing", **brief(r)})
        candidates += candidate_urls_from_text(r.text, r.url)
        soup = BeautifulSoup(r.text, "html.parser")
        for tag in soup.find_all(["a","script","link"]):
            u = tag.get("href") or tag.get("src")
            if u and ("download" in u.lower() or "file" in u.lower()):
                candidates.append(urljoin(r.url, u))
        form_summaries = []
        for form in soup.find_all("form"):
            action = urljoin(r.url, form.get("action") or r.url)
            method = (form.get("method") or "get").lower()
            vals = {}
            for inp in form.find_all("input"):
                name, val = inp.get("name"), inp.get("value")
                if name and val is not None: vals[name] = val
            form_summaries.append({"action":action,"method":method,"inputs":vals})
            if "download" in action.lower() or any("file" in k.lower() or "uid" in k.lower() for k in vals):
                try:
                    rr = session.post(action, data=vals, timeout=60, allow_redirects=True) if method=="post" else session.get(action, params=vals, timeout=60, allow_redirects=True)
                    attempts.append({"stage":"landing_form", "action":action, "method":method, "inputs":vals, **brief(rr)})
                    got = try_parse_catalog_bytes(rr.content, "landing_form:"+rr.url, attempts)
                    if got: return got, attempts
                    arch, members = try_archive(rr.content, "landing_form_archive:"+rr.url, attempts)
                    if arch: return arch, attempts
                    candidates += candidate_urls_from_text(rr.text, rr.url)
                except Exception as e:
                    attempts.append({"stage":"landing_form", "action":action, "error":str(e)})
        attempts.append({"stage":"landing_forms_inventory", "forms":form_summaries[:20]})
    except Exception as e:
        attempts.append({"stage":"landing", "error":str(e)})

    # 4) Common dataset/archive parameter spellings, including dataset page numeric UID.
    direct = [
        ADS+"?data_set_uid="+DATASET_UID,
        ADS+"?dataset_uid="+DATASET_UID,
        ADS+"?data_set="+DATASET,
        ADS+"?dataset="+DATASET,
        FDS+"/metadata?data_set_uid="+DATASET_UID,
        FDS+"/metadata?dataset_uid="+DATASET_UID,
    ]
    candidates += direct
    candidates = list(dict.fromkeys(candidates))

    # 5) Try every discovered/direct candidate. Accept only a payload that parses as >=1000 coords.
    for u in candidates[:120]:
        try:
            r = session.get(u, timeout=90, allow_redirects=True)
            attempts.append({"stage":"candidate", **brief(r)})
            got = try_parse_catalog_bytes(r.content, r.url, attempts)
            if got: return got, attempts
            arch, members = try_archive(r.content, r.url, attempts)
            if members: attempts.append({"stage":"archive_members", "url":r.url, "members":members[:100]})
            if arch: return arch, attempts
        except Exception as e:
            attempts.append({"stage":"candidate", "url":u, "error":f"{type(e).__name__}: {e}"})

    raise RuntimeError("MGDS acquisition v2 exhausted Search API, non-spatial FileServer, landing forms, and archive candidates")


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--output",required=True); ap.add_argument("--status-output",required=True); a=ap.parse_args()
    out=Path(a.output); status=Path(a.status_output); out.parent.mkdir(parents=True,exist_ok=True)
    base={"artifact_id":"JANUS-ECHO-COUSTEAU-EA-TPHASE-BLIND-CLUSTER-RUN-002-2026-08-21-v1.1","started_at_utc":datetime.now(timezone.utc).isoformat()}
    try:
        s=requests.Session(); s.headers.update({"User-Agent":"Janus-Echo-Cousteau/1.1 scientific reproducibility audit"})
        acquired, attempts = acquire(s)
        raw, coords, parse_meta, source_label = acquired
        blind=blind_cluster(coords)
        # ANCHOR REVEAL IS DELIBERATELY AFTER blind phase is finalized and hashed.
        reveal=reveal_and_score(coords, blind, -3.865418, 3.854924)
        nd=[x["nearest_cluster"]["anchor_to_center_km"] for x in reveal["configs"] if x.get("nearest_cluster")]
        anyp=any(x.get("nearest_cluster") and x["nearest_cluster"]["anchor_inside_cluster_p95_radius"] for x in reveal["configs"])
        nearest_event=reveal["nearest_event"]["distance_km"]
        verdict="ANCHOR_OVERLAPS_BLIND_CLUSTER_P95__REQUIRES_TECTONIC_CONTROL" if anyp else "NO_BLIND_CLUSTER_P95_OVERLAP_WITH_FROZEN_ANCHOR"
        result={
          "artifact_id":base["artifact_id"],"research_branch":"Janus-Echo-Кусто","completed_at_utc":datetime.now(timezone.utc).isoformat(),
          "source":{"doi":DOI,"dataset":DATASET,"dataset_page":LANDING,"source_label":source_label,"raw_sha256":sha256_bytes(raw),"raw_bytes":len(raw),"raw_committed":False,"license":"CC BY-NC-SA 3.0","expected_event_count":EXPECTED_EVENT_COUNT,"parsed_count":len(coords),"expected_count_exact_match":len(coords)==EXPECTED_EVENT_COUNT,"parse":parse_meta,"acquisition_attempts":attempts},
          "preregistration":{"anchor_hidden_during_clustering":True,"clustering_parameters_frozen_before_anchor_reveal":True,"dbscan_grid":DBSCAN_GRID,"null_domain":NULL_DOMAIN,"null_samples":NULL_SAMPLES,"null_seed":NULL_SEED},
          "blind_phase":blind,"post_reveal":reveal,
          "summary":{"nearest_catalog_event_to_anchor_km":nearest_event,"nearest_blind_cluster_center_across_grid_km":round(min(nd),3) if nd else None,"anchor_inside_any_blind_cluster_p95_radius":anyp,"verdict":verdict,"semantic_status":"UNCONFIRMED"},
          "hard_rules":["BLIND_CLUSTER_BEFORE_ANCHOR_REVEAL","NO_PARAMETER_RETUNING_AFTER_REVEAL","MID_ATLANTIC_RIDGE_SEISMICITY_IS_MANDATORY_CONTROL","DISTANCE_IS_NOT_CAUSATION","NO_RECENTERING","NO_UNDERWATER_PYRAMID_DETECTED_YET"],
          "status":"BLIND_CLUSTER_RUN_COMPLETE"}
        out.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
        status.write_text(json.dumps({**base,"status":"SUCCESS","completed_at_utc":datetime.now(timezone.utc).isoformat(),"parsed_count":len(coords),"source_sha256":sha256_bytes(raw),"verdict":verdict,"result_path":str(out)},indent=2),encoding="utf-8")
        print(json.dumps(result["summary"],indent=2)); return 0
    except Exception as e:
        # attempts may contain enough diagnostic detail to make the blocker reproducible.
        payload={**base,"status":"BLOCKED_DATA_ACQUISITION_OR_PARSE","completed_at_utc":datetime.now(timezone.utc).isoformat(),"error_type":type(e).__name__,"error":str(e)}
        if "attempts" in locals(): payload["acquisition_attempts"]=attempts
        status.write_text(json.dumps(payload,indent=2,ensure_ascii=False),encoding="utf-8")
        print(json.dumps(payload,indent=2),flush=True); return 2

if __name__=="__main__": raise SystemExit(main())
