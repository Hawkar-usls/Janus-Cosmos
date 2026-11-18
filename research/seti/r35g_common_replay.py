#!/usr/bin/env python3
"""Replay exact frozen R35B manifest; export common calibration, evaluation and shared-reference inputs."""
import hashlib,io,json,os,pathlib,platform,resource,time,zipfile
import numpy as np
import requests
import r35b_grand_benchmark as b
P=pathlib.Path("research/seti/r35g_common_native_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r35g");OUT.mkdir(parents=True,exist_ok=True)
IN=OUT/"inputs";IN.mkdir(exist_ok=True)
def sha(x):return hashlib.sha256(x).hexdigest()
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def frozen():
    a=S["source"]["manifest"];u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['id']}/zip"
    r=requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=60);r.raise_for_status()
    raw=r.content;assert len(raw)<16*1024**2 and sha(raw)==a["digest"]
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        n=[i for i in z.infolist() if pathlib.PurePosixPath(i.filename).name=="manifest.json"];assert len(n)==1 and n[0].file_size<8*1024**2
        data=z.read(n[0])
    m=json.loads(data);assert m["protocol_sha256"]==sha(b.P.read_bytes())
    (OUT/"manifest.json").write_bytes(data)
    return m
def encode(d,label,name,L,marg,th):
    rows=d["rows"];arr={"V":d["Y"].astype(np.float32),"y":b.labels(d),"truth":np.stack([r[7] for r in rows]),"active":np.stack([r[8] for r in rows]),"amp":np.asarray([r[4] for r in rows]),"model":np.asarray([r[2] for r in rows],dtype="U32"),"tile_sha":np.asarray([r[6] for r in rows],dtype="U128")}
    scores={};per={}
    for k,a in b.allarms(d,L,marg).items():
        per[k],rec=b.metrics(d,a,th[k]);scores[k]=rec
        arr[k+"_alert"]=rec["alert"];arr[k+"_localized"]=rec["localized_positive"];arr[k+"_score"]=rec["scores"];arr[k+"_paths"]=rec["paths"]
    path=IN/(label+"_"+name+".npz")
    np.savez_compressed(path,**arr)
    return {"name":path.name,"sha256":sha(path.read_bytes()),"array_sha256":sha(arr["V"].tobytes()),"examples":len(rows),"metrics":per}
def main():
    m=frozen();selected=m["selected"];assert len(selected)==8 and len({q["day"] for q in selected})==8
    assert [q["split"] for q in selected]==["component_fit"]*2+["null_calibration"]*2+["evaluation"]*4
    f=S["freshness_audit"];overlap=[{"target":q["target"],"day":q["day"],"split":q["split"]} for q in selected if q["target"] in f["known_r32_targets"] or q["day"] in f["known_r32_days"]]
    emit("R35G_FRESHNESS",{"overlaps":overlap,"status":"HOLD_R32_OVERLAP" if overlap else "NO_R32_GROUP_OVERLAP","targets":[{k:q[k] for k in ("target","day","split")} for q in selected]})
    init=b.taylor.init_tree();checks={**b.h.selftest(),**b.selftest()};seen=set();t0=time.perf_counter()
    def gather(q):
        rows,R,ch,hs=b.ref.gather(q);assert not seen&hs;seen.update(hs)
        return b.panel(rows,R),ch
    fits=[gather(q)[0] for q in selected[:2]]
    L=b.h.fit([(d["rows"],d["pack"]) for d in fits])
    marg={}
    for k in b.B["arms"][:-1]:
        v=b.controls(fits,k,L);t=b.cut(v);marg[k]=[t,max(.1,t-float(np.quantile(v,.8)))]
    null=[gather(q)[0] for q in selected[2:4]]
    th={k:b.cut(np.concatenate([b.allarms(d,L,marg)[k][0][b.labels(d)==0] for d in null])) for k in b.B["arms"]}
    lock={"component_fit":L,"margins":marg,"thresholds":th,"manifest_sha256":sha((OUT/"manifest.json").read_bytes())}
    save("lock.json",lock);lockhash=sha((OUT/"lock.json").read_bytes());emit("R35G_LOCK",{"sha256":lockhash,"thresholds":th})
    records=[encode(d,"cal",q["target"],L,marg,th) for q,d in zip(selected[2:4],null)]
    assert all(sum(t["metrics"][k]["control_alerts"] for t in records)<=12 for k in b.B["arms"])
    for q in selected[4:]:
        d,ch=gather(q);shared=b.panel([v[0] for v in ch],np.stack([v[1] for v in ch]))
        records.append(encode(d,"eval",q["target"],L,marg,th));records.append(encode(shared,"shared",q["target"],L,marg,th))
        save("partial.json",records);emit("R35G_INPUT",{"target":q["target"],"record":records[-2],"shared":records[-1]})
    assert sha((OUT/"lock.json").read_bytes())==lockhash and len(records)==10
    save("receipt.json",{"status":"PASS_COMMON_REPLAY","source_manifest":S["source"]["manifest"],"lock_sha256":lockhash,"records":records,"r32_overlap":overlap,"freshness_status":"HOLD_R32_OVERLAP" if overlap else "NO_R32_GROUP_OVERLAP_NOT_EXTERNAL_BLIND","checks":checks,"compile_seconds":init,"total_seconds":time.perf_counter()-t0,"cpu":platform.processor(),"peak_rss_kib_process":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"production_adoption":False})
    emit("R35G_SUMMARY",{"files":len(records),"lock_sha256":lockhash,"overlap":overlap,"checks":checks,"production_adoption":False})
if __name__=="__main__":main()
