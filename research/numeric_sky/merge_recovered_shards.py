#!/usr/bin/env python3
"""Merge immutable R0 scan and targeted R0R corrections into eight reproducible source shards.
Source artifacts remain untouched; only explicit original FAILED rows may be replaced.
"""
import contextlib
import gzip
import hashlib
import json
from pathlib import Path

ROOT=Path("runtime/numeric_sky")
ORIGINAL=ROOT/"original"
REPAIRED=ROOT/"repaired"
MERGED=ROOT/"collected"
MERGED.mkdir(parents=True,exist_ok=True)

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def save(path,v):path.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
def main():
    source={}
    coord={}
    original_inv=[]
    original_sha={}
    for i in range(8):
        folder=ORIGINAL/f"numsky-shard-{i}"
        invent=load(folder/"inventory.json");co=load(folder/"explicit_icrs_candidates.json")
        assert invent["shard"]==i and invent["shards"]==8
        assert all(r["repo"] not in source for r in invent["repos"])
        for r in invent["repos"]:source[r["repo"]]=(i,r)
        original_inv.append(invent);coord[i]=co["candidates"]
        original_sha[i]=sha(folder/"numbers.jsonl.gz")
    assert len(source)==48
    failed={r for r,(i,e) in source.items() if e["status"]=="FAILED"}
    updated={}
    repair_coords={i:[] for i in range(8)}
    repair_sha={}
    for i in range(4):
        folder=REPAIRED/f"numsky-repaired-{i}"
        invent=load(folder/"inventory.json")
        co=load(folder/"explicit_icrs_candidates.json")
        assert invent["shard"]==i and invent["shards"]==4 and invent["source_run"]==44833839206
        assert len(co["candidates"])==sum(x["candidate_pairs"] for x in invent["repos"])
        for r in invent["repos"]:
            repo=r["repo"]
            assert repo in failed and repo not in updated
            shard,old=source[repo]
            assert old["head_sha"]==r["head_sha"],"REPAIRED_REVISION_MISMATCH"
            updated[repo]=r
        for c in co["candidates"]:
            assert c["repo"] in updated and c["repo"] in source
            repair_coords[source[c["repo"]][0]].append(c)
        repair_sha[i]=sha(folder/"numbers.jsonl.gz")
    assert set(updated)==failed, "ALL_FAILED_REPOSITORIES_MUST_HAVE_REPAIR_RECORD"
    with contextlib.ExitStack() as context:
        sinks={}
        for i in range(8):
            d=MERGED/f"numsky-shard-{i}";d.mkdir(parents=True,exist_ok=True)
            sinks[i]=context.enter_context(gzip.open(d/"numbers.jsonl.gz","wt",encoding="utf-8",compresslevel=6))
            with gzip.open(ORIGINAL/f"numsky-shard-{i}"/"numbers.jsonl.gz","rt",encoding="utf-8") as f:
                for line in f:
                    # The original FAILED repositories have no numeric records.
                    sinks[i].write(line)
        for j in range(4):
            with gzip.open(REPAIRED/f"numsky-repaired-{j}"/"numbers.jsonl.gz","rt",encoding="utf-8") as f:
                for line in f:
                    obj=json.loads(line)
                    assert obj["repo"] in updated and updated[obj["repo"]]["status"]!="FAILED"
                    src_idx,old=source[obj["repo"]]
                    assert obj["revision"]==old["head_sha"]
                    sinks[src_idx].write(line)
    totals={"source_run":44833839206,"recovery_run":44834866914,"source_shas":original_sha,"repair_shas":repair_sha,
            "input_repos":len(source),"source_failed":len(failed),
            "recovery_fixed":sum(x["status"]!="FAILED" for x in updated.values()),"recovery_remaining_failed":sum(x["status"]=="FAILED" for x in updated.values())}
    full=[]
    for i in range(8):
        new_recs=[]
        for old in original_inv[i]["repos"]:
            rep=updated.get(old["repo"],old)
            if old["repo"] in updated and rep["status"]!="FAILED":
                rep={**rep,"recovery_origin":"R0R_FROM_ORIGINAL_FAILED","initial_status":"FAILED"}
            new_recs.append(rep);full.append(rep)
        d=MERGED/f"numsky-shard-{i}"
        save(d/"inventory.json",{"shard":i,"shards":8,"repos":new_recs,"private_excluded":original_inv[i]["private_excluded"]})
        save(d/"explicit_icrs_candidates.json",{"shard":i,"candidates":coord[i]+repair_coords[i]})
        save(d/"audit_summary.json",{"shard":i,"repositories":len(new_recs),
            "complete":sum(r["status"]=="COMPLETE_SELECTED_CODE" for r in new_recs),
            "partial":sum(r["status"]=="PARTIAL" for r in new_recs),
            "failed":sum(r["status"]=="FAILED" for r in new_recs),
            "selected_files":sum(r["selected_files"] for r in new_recs),
            "numbers":sum(r["numbers"] for r in new_recs),
            "coordinate_candidates":len(coord[i]+repair_coords[i])})
    totals.update({"full_complete":sum(r["status"]=="COMPLETE_SELECTED_CODE" for r in full),
        "full_partial":sum(r["status"]=="PARTIAL" for r in full),
        "full_failed":sum(r["status"]=="FAILED" for r in full),
        "numbers":sum(r["numbers"] for r in full),
        "labelled_coordinates":sum(len(load(MERGED/f"numsky-shard-{i}"/"explicit_icrs_candidates.json")["candidates"]) for i in range(8))})
    save(ROOT/"numeric_source_merge_receipt.json",totals)
    print("NUMSKY_R0R_MERGED",json.dumps({k:v for k,v in totals.items() if k not in ("source_shas","repair_shas")}),flush=True)
if __name__=="__main__":main()
