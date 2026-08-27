#!/usr/bin/env python3
"""JANUS Archiv: TOPA-inspired append-only provenance corpus + local FTS search.

This is intentionally a discovery/provenance sidecar. Search rank, snippets and
indexed metadata never become scientific evidence. Primary catalog/archive values
must remain attached to their source/query provenance in experiment receipts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sqlite3
from pathlib import Path
from typing import Any, Dict, Iterable, List

SCHEMA = "janus.cosmos.archiv.record.v1"
TOKEN_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.+\-/]{1,}")


def canonical_json(obj: Any) -> str:
    return json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def seal_record(record: Dict[str, Any]) -> Dict[str, Any]:
    out=dict(record)
    out.setdefault("schema",SCHEMA)
    out.setdefault("scientific_authority","DISCOVERY_OR_PROVENANCE_METADATA_NOT_TRUTH")
    clean={k:v for k,v in out.items() if k!="record_sha256"}
    out["record_sha256"]=sha256_text(canonical_json(clean))
    return out


def append_jsonl(path: Path, records: Iterable[Dict[str, Any]]) -> int:
    path.parent.mkdir(parents=True,exist_ok=True); n=0
    with path.open("a",encoding="utf-8") as f:
        for record in records:
            f.write(canonical_json(seal_record(record))+"\n"); n+=1
    return n


def read_jsonl(path: Path) -> Iterable[Dict[str, Any]]:
    if not path.exists(): return
    with path.open("r",encoding="utf-8") as f:
        for line in f:
            if line.strip(): yield json.loads(line)


def meaningful_tokens(q: str) -> List[str]:
    out=[]
    for token in TOKEN_RE.findall(q):
        t=token.lower()
        if len(t)>=2 and t not in out: out.append(t)
    return out[:16]


def build_index(corpus: Path, db: Path) -> Dict[str,Any]:
    records=list(read_jsonl(corpus)); db.parent.mkdir(parents=True,exist_ok=True)
    if db.exists(): db.unlink()
    conn=sqlite3.connect(db)
    try:
        conn.execute("CREATE TABLE records (id INTEGER PRIMARY KEY, payload TEXT NOT NULL)")
        conn.execute("CREATE VIRTUAL TABLE records_fts USING fts5(title, source, locator, query, notes, payload UNINDEXED)")
        for rec in records:
            payload=canonical_json(rec)
            conn.execute("INSERT INTO records(payload) VALUES (?)",(payload,))
            conn.execute("INSERT INTO records_fts(title,source,locator,query,notes,payload) VALUES (?,?,?,?,?,?)",(
                str(rec.get("title") or ""),str(rec.get("source") or ""),str(rec.get("locator") or ""),
                str(rec.get("query") or ""),str(rec.get("notes") or ""),payload))
        conn.commit()
    finally: conn.close()
    return {"engine":"sqlite_fts5","documents":len(records),"index_path":str(db),"rank_is_truth":False}


def search(db: Path, query: str, limit: int=20) -> List[Dict[str,Any]]:
    tokens=meaningful_tokens(query)
    match=" OR ".join(f'"{t}"' for t in tokens) if tokens else f'"{query.strip()}"'
    conn=sqlite3.connect(db); conn.row_factory=sqlite3.Row
    try:
        rows=conn.execute("SELECT payload,bm25(records_fts,3.0,1.5,1.2,1.0,0.5,0.0) score FROM records_fts WHERE records_fts MATCH ? ORDER BY score LIMIT ?",(match,limit)).fetchall()
        out=[]
        for rank,row in enumerate(rows,1):
            rec=json.loads(row["payload"]); rec["local_search"]={"engine":"sqlite_fts5","rank":rank,"bm25_native":float(row["score"]),"rank_is_truth":False}; out.append(rec)
        return out
    finally: conn.close()


def self_test() -> None:
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        root=Path(td); corpus=root/"c.jsonl"; db=root/"i.sqlite"
        append_jsonl(corpus,[{"title":"AllWISE epoch","source":"IRSA","locator":"MJD 55400","query":"AllWISE RA_pm"}])
        meta=build_index(corpus,db); hits=search(db,"AllWISE epoch")
        assert meta["documents"]==1 and hits and hits[0]["scientific_authority"].endswith("NOT_TRUTH")
    print("JANUS_ARCHIV_SEARCH_SELF_TEST=PASS")


def main() -> int:
    p=argparse.ArgumentParser(); p.add_argument("--self-test",action="store_true"); a=p.parse_args()
    if a.self_test: self_test(); return 0
    p.error("use as a library or extend with explicit corpus/index paths"); return 2


if __name__=="__main__": raise SystemExit(main())
