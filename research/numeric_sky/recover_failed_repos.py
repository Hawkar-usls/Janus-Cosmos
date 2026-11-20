#!/usr/bin/env python3
"""R0R: rescan only failed original-code repositories, selected by immutable R0 receipts."""
import argparse,gzip,hashlib,json,pathlib,sys
import scan_numeric_sources as base
P=base.PROTO
ROOT=pathlib.Path("runtime/numeric_sky")
SOURCE=ROOT/"collected"
def read_failure_list():
    r=[]
    for i in range(8):
        f=SOURCE/f"numsky-shard-{i}"/"inventory.json"
        assert f.is_file(),str(f)
        j=json.loads(f.read_text())
        assert j["shard"]==i and j["shards"]==8
        r.extend(j["repos"])
    expected=P["provenance"]["public_repositories"]
    assert len(r)==len(expected)==48
    assert sorted(x["repo"] for x in r)==sorted(expected)
    failures=[x for x in r if x["status"]=="FAILED"]
    assert len(failures)>=1
    return r,failures
def main():
    arg=argparse.ArgumentParser()
    arg.add_argument("--shard",required=True,type=int)
    arg.add_argument("--shards",type=int,default=4)
    a=arg.parse_args()
    assert a.shards==4 and a.shard in range(a.shards)
    originals,failures=read_failure_list()
    path=ROOT/f"repair-{a.shard:02d}"
    path.mkdir(parents=True,exist_ok=True)
    targets=[v for i,v in enumerate(failures) if i%a.shards==a.shard]
    results=[];candidates=[]
    with gzip.open(path/"numbers.jsonl.gz","wb",compresslevel=6) as out:
        for v in targets:
            info,cand=base.scan_repository(v["repo"],path,out)
            assert info["repo"]==v["repo"]
            assert info["head_sha"]==v["head_sha"],"FAILED_SOURCE_REVISION_DRIFT:"+v["repo"]
            results.append(info);candidates.extend(cand)
            print("NUMSKY_RECOVERY_REPO",json.dumps({"repo":info["repo"],"status":info["status"],"files":info["selected_files"],
                 "numbers":info["numbers"],"candidates":info["candidate_pairs"],"errors":info["errors"][:3]}),flush=True)
    (path/"inventory.json").write_text(json.dumps({"shard":a.shard,"shards":4,"source_run":44833839206,
        "source_failures":len(failures),"source_selected":len(targets),"repos":results},indent=2))
    (path/"explicit_icrs_candidates.json").write_text(json.dumps({"shard":a.shard,"candidates":candidates},indent=2))
    summary={"shard":a.shard,"sources":len(targets),"recovered":sum(x["status"]!="FAILED" for x in results),
         "remaining_failed":sum(x["status"]=="FAILED" for x in results),
         "files":sum(x["selected_files"] for x in results),"numbers":sum(x["numbers"] for x in results),
         "coordinates":len(candidates)}
    (path/"audit_summary.json").write_text(json.dumps(summary,indent=2))
    print("NUMSKY_RECOVERY_SHARD",json.dumps(summary),flush=True)
if __name__=="__main__":main()
