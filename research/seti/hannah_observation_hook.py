"""Measured-only observation hook; stage data until the frozen exam is sealed."""
import json,hashlib,os
import numpy as np
import r30c_hannah_memory as h
class ObservationAccumulator:
    def __init__(self,parent,manifest_artifact):
        self.parent,_=h.append(parent,[])
        self.manifest_artifact=manifest_artifact
        self.seen=set()
        self.old_hashes={x["provenance"]["source_tile_sha256"] for x in self.parent["entries"]}
        self.records=[];self.observations=[]
    def ingest(self,q,f,t,raw):
        ident=t["slice_sha256"];assert h.ishash(ident)
        assert hashlib.sha256(np.ascontiguousarray(raw).tobytes()).hexdigest()==ident
        if ident in self.seen or ident in self.old_hashes:raise ValueError("Repeated source tile; not a new observation")
        self.seen.add(ident)
        meta=h.meta(t,f["name"],self.manifest_artifact)
        meta.update({"source_artifact_role":"FROZEN_METADATA_MANIFEST_NOT_RAW_TILE","archive_url":f["url"],"target":q["target"],"day":q["day"],"role":f["role"]})
        ids=[]
        for score,path in h.search_tile(raw):
            rec=h.make_record(raw,path,meta)
            assert rec["label"]=="UNKNOWN"
            self.records.append(rec);ids.append(rec["id"])
        self.observations.append({"target":q["target"],"day":q["day"],"scan":f["name"],"role":f["role"],"tile_sha256":ident,"candidate_ids":ids})
    def seal(self):
        new,delta=h.append(self.parent,self.records)
        assert new["entries"][:len(self.parent["entries"])]==self.parent["entries"]
        again,rep=h.append(new,self.records)
        assert again==new and rep["added"]==0
        assert sum(x["label"]=="VERIFIED_RFI" for x in new["entries"])==sum(x["label"]=="VERIFIED_RFI" for x in self.parent["entries"])
        return new,delta,rep
def manifest_artifact():
    i=os.environ.get("HANNAH_MANIFEST_ARTIFACT_ID","")
    d=os.environ.get("HANNAH_MANIFEST_ARTIFACT_SHA256","").removeprefix("sha256:")
    if not i.isdecimal() or not h.ishash(d):raise ValueError("Missing frozen manifest artifact authority")
    return {"id":int(i),"digest":d}
