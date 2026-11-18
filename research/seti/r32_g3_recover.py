#!/usr/bin/env python3
"""Recover exact sealed R31 G3 bytes without repository write credentials."""
import base64,hashlib,json,pathlib
import r30c_hannah_memory as mem
import r26c_reference_alignment as r
OUT=pathlib.Path("runtime/r32_export");OUT.mkdir(parents=True,exist_ok=True)
PARENT=pathlib.Path("research/seti/memory/hannah_library_v0002.json")
G3_SHA="ec382441e707ebf131e3cb84b7037a9729f31ef3163412e99eca4462dcf0d807"
G2_SHA="95f90afbb07e1bca4d73fb28d1da64e778f59d019782587ef1ec7a868ddebc1e"
ART={"id":14411051155,"digest":"0f1d3d785f566622773150eefc5d262a9613dd90b438b94ffcf378e61cc9c1db"}
def sha(b):return hashlib.sha256(b).hexdigest()
def restore():
    r.BLC_OUT=OUT;r.BLC_P["artifacts"]["hannah_g3"]=ART
    z,n,p=r.archive("hannah_g3")
    raw=z.read(n["library.json"]);result=r.get(z,n,"hook_result.json")
    z.close();p.unlink()
    assert sha(raw)==G3_SHA==result["library_sha256"]
    parent=PARENT.read_bytes();assert sha(parent)==G2_SHA==result["parent_sha256"]
    old,_=mem.append(json.loads(parent),[]);new,_=mem.append(json.loads(raw),[])
    assert len(old["entries"])==66 and len(new["entries"])==109 and new["entries"][:66]==old["entries"]
    assert len({e["id"] for e in new["entries"]})==109
    assert all(e["label"]=="UNKNOWN" for e in new["entries"][66:])
    replay,counts=mem.append(new,new["entries"][66:]);assert replay==new and counts["added"]==0
    rev=mem.review(new);assert rev["verified_rfi"]==1 and rev["rfi_candidates"]==12 and rev["unknown_without_recurrence"]==96
    (OUT/"library.json").write_bytes(raw)
    return new,raw
def main():
    new,raw=restore();enc=base64.b64encode(raw).decode();chunks=[enc[i:i+3000] for i in range(0,len(enc),3000)]
    for i,c in enumerate(chunks):print("HANNAH_G3_CHUNK",i,c,flush=True)
    print("HANNAH_G3_EXPORT",json.dumps({"status":"PASS_EXACT_G3_RECOVERY","entries":len(new["entries"]),"sha256":sha(raw),"bytes":len(raw),"chunks":len(chunks),"parent_preserved":66,"new_entries":43,"verified":1}),flush=True)
if __name__=="__main__":main()
