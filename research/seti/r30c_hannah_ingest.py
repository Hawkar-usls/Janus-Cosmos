#!/usr/bin/env python3
"""Incremental real-observation ingestion; never learns origin from detector outputs."""
import base64,hashlib,json,pathlib,time
import r30c_hannah_memory as h
import r26c_reference_alignment as r
import r9_expanded_search as eng
P=pathlib.Path("research/seti/r30c_hannah_ingest_protocol.json")
S=json.loads(P.read_text());OUT=pathlib.Path("runtime/r30c_ingest");OUT.mkdir(parents=True,exist_ok=True)
h.OUT=OUT;r.OUT=OUT;r.BLC_OUT=OUT;eng.OUT=OUT
def run():
    t0=time.perf_counter();checks=h.selftest();h.emit("HANNAH_INGEST_SELFTEST",checks)
    p=S["parent"];n=S["new_observations"]
    r.BLC_P["artifacts"]["parent"]={"id":p["id"],"digest":p["digest"]}
    z,idx,zp=r.archive("parent");raw=z.read(idx["library.json"]);z.close();zp.unlink()
    assert hashlib.sha256(raw).hexdigest()==p["library_sha256"]
    prev,_=h.append(json.loads(raw),[]);assert len(prev["entries"])==p["expected_entries"]
    r.BLC_P["artifacts"]["fresh_manifest"]={"id":n["id"],"digest":n["digest"]}
    z,idx,zp=r.archive("fresh_manifest");m=r.get(z,idx,"manifest.json");z.close();zp.unlink()
    q=m["selected"][n["selected_index"]];assert q["split"]=="component_fit"
    records=[];obs=[];seen=set();new_scans=set()
    for f in q["files"]:
        scan=f["name"];assert scan not in new_scans;new_scans.add(scan)
        for t,raw_tile in r.tiles(f):
            assert t["slice_sha256"] not in seen;seen.add(t["slice_sha256"])
            ids=[]
            for score,path in h.search_tile(raw_tile):
                meta=h.meta(t,scan,{"id":n["id"],"digest":n["digest"]})
                meta.update({"source_artifact_role":"FROZEN_METADATA_MANIFEST_NOT_RAW_TILE","archive_url":f["url"],"target":q["target"],"day":q["day"],"role":f["role"]})
                rec=h.make_record(raw_tile,path,meta);records.append(rec);ids.append(rec["id"])
            obs.append({"scan":scan,"target":q["target"],"day":q["day"],"role":f["role"],"tile_sha256":t["slice_sha256"],"candidate_ids":ids})
    assert len(obs)==16 and len(new_scans)==4
    assert not seen&{x["provenance"]["source_tile_sha256"] for x in prev["entries"]}
    after,delta=h.append(prev,records);rep,again=h.append(after,records)
    assert after["entries"][:len(prev["entries"])]==prev["entries"] and rep==after and again["added"]==0
    before_review=h.review(prev);review=h.review(after)
    assert review["verified_rfi"]==before_review["verified_rfi"]==1
    h.save("library.json",after);h.save("observations.json",obs);h.save("review.json",review)
    lib=OUT/"library.json";payload=lib.read_bytes();enc=base64.b64encode(payload).decode("ascii");chunks=[enc[i:i+3000] for i in range(0,len(enc),3000)]
    # Public measured fingerprints only. No environment or secrets in export.
    for i,chunk in enumerate(chunks):print("HANNAH_LEDGER_CHUNK",i,chunk,flush=True)
    result={"status":"PASS_INCREMENTAL_MEMORY_DEVELOPMENT" if delta["added"]>0 else "HOLD_NO_NEW_CANDIDATES",
       "parent":p,"source_manifest":n,"target":q["target"],"day":q["day"],"new_observations":len(obs),"new_scans":len(new_scans),"delta":delta,"replay":again,
       "parent_review_counts":{k:before_review[k] for k in ("verified_rfi","rfi_candidates","unknown_without_recurrence")},
       "review_counts":{k:review[k] for k in ("verified_rfi","rfi_candidates","unknown_without_recurrence")},"recurrence_pairs":len(review["recurrence_pairs"]),
       "library_sha256":hashlib.sha256(payload).hexdigest(),"export_chunks":len(chunks),"selftests":checks,"seconds":time.perf_counter()-t0,"production_adoption":False}
    h.save("result.json",result);h.emit("HANNAH_INGEST_SUMMARY",result)
if __name__=="__main__":run()
