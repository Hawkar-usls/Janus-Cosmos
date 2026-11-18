#!/usr/bin/env python3
"""Fresh paired HANNAH advisory HOLD exam, stage all measured tiles independently."""
import hashlib,itertools,json,pathlib,sys,time,math
import numpy as np
import r30_false_alarm_shield as shield
import r26c_reference_alignment as ref
import r30c_hannah_memory as mem
import hannah_observation_hook as hook
P=pathlib.Path("research/seti/r31_hannah_paired_protocol.json");S=json.loads(P.read_text());H=S["hannah"]
OUT=pathlib.Path("runtime/r31");OUT.mkdir(parents=True,exist_ok=True)
shield.P=P;shield.S=S;shield.C=S["shield"];shield.OUT=OUT
shield.h.S=S;shield.rb.S=S;ref.S=S;shield.od.S=S;ref.OUT=OUT;shield.eng.OUT=OUT
shield.hold.S={**shield.hold.S,"seed":S["seed"],"exclude":sorted(set(shield.hold.S["exclude"]+S["exclude_extra"]))}
shield.hold.P=P;shield.hold.OUT=OUT;mem.OUT=OUT
def save(n,v):mem.save(n,v)
def emit(k,v):mem.emit(k,v)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def gather(q,acc):
    fs=q["files"];assert [f["role"] for f in fs]==["ON1","OFF1","ON2","OFF2"]
    rng=np.random.default_rng(int.from_bytes(hashlib.sha256((S["seed"]+"|"+q["target"]).encode()).digest()[:8],"big"))
    plan=list(itertools.product(S["amplitudes"],S["injection_models"],S["gap_levels"]))*4;rng.shuffle(plan)
    rows=[];refs=[];challenges=[];sources=set();k=0
    for i in (0,2):
        ti=ref.tiles(fs[i]);ri=ref.tiles(fs[i+1]);assert len(ti)==len(ri)==4
        for (a,X),(b,R) in zip(ti,ri):
            acc.ingest(q,fs[i],a,X);acc.ingest(q,fs[i+1],b,R)
            assert a["center_mhz"]==b["center_mhz"]; X,R,lo=ref.common(a,b,X,R);assert X.shape[0]==R.shape[0]==16
            sources.update([a["slice_sha256"],b["slice_sha256"]])
            # Hard patches picked from ON background before injection, not used to estimate labels.
            z=ref.norm(X.astype(float));rank=np.argsort(-np.median(z,axis=0),kind="stable");hard=[]
            for j in rank:
                if 64<=j<X.shape[1]-64 and all(abs(j-v)>=128 for v in hard):hard.append(int(j))
                if len(hard)==4:break
            assert len(hard)==4
            for trial in range(16):
                start=hard[trial-12]-64 if trial>=12 else int(rng.integers(0,X.shape[1]-128+1))
                x=X[:,start:start+128].astype(float);r=R[:,start:start+128].astype(float)
                amp,model,gap=plan[k];k+=1;t=np.arange(16);u=t-7.5
                if model=="stationary":d=np.zeros(16,int)
                elif model=="straight":d=np.rint(t*rng.uniform(-1.8,1.8)).astype(int)
                elif model=="off_grid_quadratic":d=np.rint(rng.uniform(-.9,.9)*u+rng.uniform(-.12,.12)*u*u).astype(int)
                else:d=np.rint(rng.uniform(-.5,.5)*u+3.4*np.sin(2*np.pi*t/16+rng.uniform(0,6.28))).astype(int)
                d-=d[0];p=int(rng.integers(8-int(d.min()),120-int(d.max())))+d
                assert p.min()>=8 and p.max()<120 and np.max(abs(np.diff(p)))<=3
                active=np.ones(16,bool)
                if gap:active[rng.choice(16,int(16*gap),replace=False)]=False
                scale=1.4826*np.median(abs(x-np.median(x,axis=1)[:,None]),axis=1)
                assert np.all(scale>0)
                v=x.copy();v[t[active],p[active]]+=amp*scale[active]
                for label,xx in ((0,x),(1,v)):
                    rows.append((xx,label,model,gap,amp,trial>=12,a["slice_sha256"],p.copy(),active.copy(),{"channel0_mhz":a["channel0_mhz"]+(lo+start)*a["df_hz"]/1e6,"df_hz":a["df_hz"],"ts":a["ts"]}));refs.append(r)
                if model=="stationary":
                    # Same artificial stationary source in reference: explicit leakage stress, not ordinary calibration.
                    rr=r.copy();ss=1.4826*np.median(abs(r-np.median(r,axis=1)[:,None]),axis=1)
                    rr[t[active],p[active]]+=amp*ss[active]
                    challenges.append((rows[-1],rr))
    assert len(rows)==256 and k==128 and len(challenges)==32
    return rows,np.stack(refs),challenges,sources

def load():
    p=pathlib.Path(H["parent_path"]);assert sha(p)==H["parent_sha256"]
    ledger,_=mem.append(json.loads(p.read_text()),[]);assert len(ledger["entries"])==H["parent_entries"]
    return ledger
def evidence(row,path,entries):
    # Query fingerprint from raw patch along selected hybrid path. No labels read.
    m=row[9]
    try:f=mem.features(row[0],np.asarray(path),m["df_hz"],m["ts"])
    except ValueError:return {"status":"NO_FEATURE_SUPPORT","matches":[]}
    freq=m["channel0_mhz"]*1e6+int(path[0])*m["df_hz"]
    q={"features":f,"frequency_hz":freq,"provenance":{"source_scan":"QUERY_ONLY","source_tile_sha256":row[6]}}
    found=[]
    for e in entries:
        if e["provenance"]["source_tile_sha256"]==row[6]:raise ValueError("Exam memory leakage")
        if abs(freq-e["frequency_hz"])>H["match"]["frequency_hz"]:continue
        if abs(f["drift_hz_s"]-e["features"]["drift_hz_s"])>H["match"]["drift_hz_s"]:continue
        d=mem.compare(q,e)
        if d["grid_compatible"] and d["shape_similarity"]>=H["match"]["shape_similarity"]:found.append(e["id"])
    return {"status":"COMPARED","matches":sorted(found)}
def panel(rows,R):
    pre,_=shield.rb.panel(rows,R)
    return pre
def exam(rows,R,L,entries,groups):
    pre=panel(rows,R);q=pre[0];z,p,_,_=shield.h.fuse(q,L,"registered_hybrid")
    alarm=z>=L["thresholds"]["registered_hybrid"];y=shield.h.labels(rows)
    loc=shield.rb.loc(p,rows);weak=np.array([v[4]<=1.25 for v in rows]);baseline=alarm&loc&(y==1)
    matches=[evidence(row,path,entries) if a else {"status":"BELOW_BASELINE","matches":[]} for row,path,a in zip(rows,p,alarm)]
    out={};rec={}
    for name in H["arms"]:
        ids=groups[name];hold=np.array([bool(set(d["matches"])&ids) for d in matches]);keep=alarm&~hold;hit=keep&loc&(y==1)
        out[name]={"localized":int(hit.sum()),"weak_localized":int((hit&weak).sum()),"control_alerts":int((keep&(y==0)).sum()),"controls":int((y==0).sum()),"positives":int((y==1).sum()),
           "lost_localized":int((baseline&~hit).sum()),"lost_weak":int((baseline&~hit&weak).sum()),"held_controls":int((alarm&hold&(y==0)).sum()),"held_positives":int((alarm&hold&(y==1)).sum())}
        rec[name]={"alert":keep.tolist(),"hold":(alarm&hold).tolist()}
    return out,{"y":y.tolist(),"weak":weak.tolist(),"localized":loc.tolist(),"score":z.tolist(),"path":p.tolist(),"source_tile":[v[6] for v in rows],"evidence":matches,"arms":rec}
def main():
    t0=time.perf_counter();checks=mem.selftest();shield.a.init_tree();checks.update(shield.h.selftest())
    parent=load();entries=parent["entries"];parent_ids={e["id"] for e in entries};review=mem.review(parent)
    groups={"baseline_full":set(),"verified_only_hold":{e["id"] for e in entries if e["label"]=="VERIFIED_RFI"},"recurrent_only_hold":set(review["derived_candidate_ids"]),"all_memory_hold":parent_ids}
    m=json.loads((OUT/"manifest.json").read_text());assert m["protocol_sha256"]==sha(P)
    selected=m["selected"];assert len(selected)==6 and len({q["day"] for q in selected})==6
    assert [q["split"] for q in selected]==["component_fit"]*2+["evaluation"]*4
    assert not set(q["day"] for q in selected)&set(S["shield"]["exclude_days"])
    acc=hook.ObservationAccumulator(parent,hook.manifest_artifact());seen=set()
    def fetch(q):
        rows,R,ch,hs=gather(q,acc);assert len(hs)==16 and not seen&hs;seen.update(hs)
        emit("R31_ACQUIRED",{"target":q["target"],"day":q["day"],"split":q["split"],"tiles":len(hs)})
        return rows,R,ch
    cal=[]
    for q in selected[:2]:
        rows,R,ch=fetch(q);cal.append((rows,panel(rows,R)))
    L=shield.h.fit(cal);save("lock.json",{"baseline":L,"memory_sha256":H["parent_sha256"],"manifest_sha256":sha(OUT/"manifest.json"),"policy":H["match"]});lock=sha(OUT/"lock.json");emit("R31_LOCK",{"sha256":lock,"memory_entries":len(entries),"recurrence":len(groups["recurrent_only_hold"])})
    out={};rec={};leak={}
    for q in selected[2:]:
        rows,R,ch=fetch(q);out[q["target"]],rec[q["target"]]=exam(rows,R,L,entries,groups)
        leak[q["target"]],_=exam([v[0] for v in ch],np.stack([v[1] for v in ch]),L,entries,groups)
        save("partial.json",out);emit("R31_TARGET",{"target":q["target"],"metrics":out[q["target"]]})
    assert sha(OUT/"lock.json")==lock and sha(pathlib.Path(H["parent_path"]))==H["parent_sha256"]
    fields=list(next(iter(out.values()))["baseline_full"]);total={k:{f:sum(v[k][f] for v in out.values()) for f in fields} for k in H["arms"]}
    lt={k:{f:sum(v[k][f] for v in leak.values()) for f in fields} for k in H["arms"]}
    assert all(v["controls"]==512 and v["positives"]==512 for v in total.values())
    assert all(v["positives"]==128 for v in lt.values())
    b=total["baseline_full"];v=total["recurrent_only_hold"];gate={"less_controls":v["control_alerts"]<b["control_alerts"],"no_lost_localized":v["lost_localized"]==0,"no_lost_weak":v["lost_weak"]==0,"no_lost_leak":lt["recurrent_only_hold"]["lost_localized"]==0}
    result={"status":"FRESH_ARCHIVE_PAIRED_DEVELOPMENT","gate":"PASS_DEVELOPMENT" if all(gate.values()) else "HOLD","gate_checks":gate,"modes":total,"leak_challenge":lt,"per_target":out,"targets":[{k:q[k] for k in ("target","day","split")} for q in selected],"lock_sha256":lock,"memory_sha256":H["parent_sha256"],"protocol_sha256":sha(P),"selftests":checks,"scope":S["scope"],"production_adoption":False}
    save("scores.json",rec);save("result.json",result);emit("R31_SUMMARY",{k:v for k,v in result.items() if k!="per_target"})
    # Only now make staged new observation memory available as a separate output.
    assert len(acc.observations)==96;new,delta,replay=acc.seal();save("library.json",new);save("observations.json",acc.observations)
    nextreview=mem.review(new);save("review.json",nextreview)
    summary={"status":"PASS_HOOK_DEVELOPMENT","delta":delta,"replay":replay,"observations":len(acc.observations),"library_sha256":sha(OUT/"library.json"),"parent_sha256":H["parent_sha256"],"verified_rfi":nextreview["verified_rfi"],"rfi_candidates":nextreview["rfi_candidates"],"unknown_without_recurrence":nextreview["unknown_without_recurrence"],"total_seconds":time.perf_counter()-t0,"production_adoption":False}
    save("hook_result.json",summary);emit("R31_HOOK_SUMMARY",summary)
if __name__=="__main__":
    if sys.argv[1]=="freeze":shield.freeze()
    elif sys.argv[1]=="scan":main()
    else:raise ValueError("phase")
