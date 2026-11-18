#!/usr/bin/env python3
"""Recover exact sealed R32 G4 bytes without repository write credentials."""
import base64,hashlib,json,pathlib
import r30c_hannah_memory as mem
import r26c_reference_alignment as r
OUT=pathlib.Path("runtime/r33_export");OUT.mkdir(parents=True,exist_ok=True)
PARENT=pathlib.Path("research/seti/memory/hannah_library_v0003.json")
G4_SHA="17acaca5bf50a5872b51647984b49dadb46b3c2a3aa4536dd150ae70673d7bb3"
G3_SHA="ec382441e707ebf131e3cb84b7037a9729f31ef3163412e99eca4462dcf0d807"
ART={"id":14411063806,"digest":"3cc29b108a54ca2c41eabcb6fde54a3c39e9d9064e09f2f74ea1d662774db0eb"}
def sha(b):return hashlib.sha256(b).hexdigest()
def restore():
    r.BLC_OUT=OUT;r.BLC_P["artifacts"]["hannah_g4"]=ART
    z,n,p=r.archive("hannah_g4")
    raw=z.read(n["library.json"]);result=r.get(z,n,"hook_result.json")
    z.close();p.unlink()
    assert sha(raw)==G4_SHA==result["library_sha256"]
    parent=PARENT.read_bytes();assert sha(parent)==G3_SHA==result["parent_sha256"]
    old,_=mem.append(json.loads(parent),[]);new,_=mem.append(json.loads(raw),[])
    assert len(old["entries"])==109 and len(new["entries"])==169 and new["entries"][:109]==old["entries"]
    assert len({e["id"] for e in new["entries"]})==169
    assert all(e["label"]=="UNKNOWN" for e in new["entries"][109:])
    replay,counts=mem.append(new,new["entries"][109:]);assert replay==new and counts["added"]==0
    rev=mem.review(new);assert rev["verified_rfi"]==1 and rev["rfi_candidates"]==14 and rev["unknown_without_recurrence"]==154
    (OUT/"library.json").write_bytes(raw)
    return new,raw
def main():
    new,raw=restore();enc=base64.b64encode(raw).decode();chunks=[enc[i:i+3000] for i in range(0,len(enc),3000)]
    for i,c in enumerate(chunks):print("HANNAH_G4_CHUNK",i,c,flush=True)
    print("HANNAH_G4_EXPORT",json.dumps({"status":"PASS_EXACT_G4_RECOVERY","entries":len(new["entries"]),"sha256":sha(raw),"bytes":len(raw),"chunks":len(chunks),"parent_preserved":109,"new_entries":60,"verified":1}),flush=True)
if __name__=="__main__":main()
