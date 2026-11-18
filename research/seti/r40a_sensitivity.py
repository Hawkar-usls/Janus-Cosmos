#!/usr/bin/env python3
"""Development only: robust local contrast and temporal consistency at archived candidate paths."""
import hashlib,io,json,os,pathlib,time,zipfile
import numpy as np
import requests
P=pathlib.Path("research/seti/r40a_sensitivity_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r40a");OUT.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def restore():
    a=S["source"];u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['id']}/zip"
    r=requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=120)
    r.raise_for_status();b=r.content;assert len(b)<128*1024**2 and sha(b)==a["sha256"]
    z=zipfile.ZipFile(io.BytesIO(b));items=[i for i in z.infolist() if not i.is_dir()]
    assert sum(i.file_size for i in items)<512*1024**2
    n={pathlib.PurePosixPath(i.filename).name:i.filename for i in items};assert len(n)==len(items)
    d=json.loads(z.read(n["receipt.json"]));lk=z.read(n["lock.json"])
    assert d["status"]=="PASS_COMMON_REPLAY" and sha(lk)==d["lock_sha256"]
    assert len(d["records"])==10
    return z,n,d,json.loads(lk)
def load(z,n,r):
    b=z.read(n[r["name"]]);assert sha(b)==r["sha256"]
    with np.load(io.BytesIO(b),allow_pickle=False) as q:
        d={k:q[k].copy() for k in q.files}
    assert sha(d["V"].tobytes())==r["array_sha256"]
    return d
def local(V,p):
    V=np.asarray(V,dtype=np.float64);p=np.asarray(p,dtype=int);N,T,W=V.shape
    assert (T,W)==(16,128) and p.shape==(N,T) and np.all((p>=0)&(p<W))
    i=np.arange(N)[:,None];t=np.arange(T)[None,:]
    signal=V[i,t,p];coords=p[:,:,None]+np.asarray(S["side_offsets"])
    valid=(coords>=0)&(coords<W);assert (valid.sum(2)>=3).all()
    side=V[i[:,:,None],t[:,:,None],np.clip(coords,0,W-1)]
    side=np.where(valid,side,np.nan);bg=np.nanmedian(side,axis=2)
    scale=np.maximum(S["side_scale_floor"],1.4826*np.nanmedian(np.abs(side-bg[:,:,None]),axis=2))
    v=(signal-bg)/scale;assert np.isfinite(v).all()
    a=v.sum(1)/np.sqrt(T)
    # Fixed winsorized four-block consistency, not a veto; allow gaps and retain signed evidence.
    cl=np.clip(v,-S["clip"],S["clip"]).reshape(N,4,4).sum(2)/2
    b=.5*cl.sum(1)/2+.5*np.sort(cl,axis=1)[:,1:3].mean(1)*2
    return a,b
def scores(q):
    V=q["V"];paths=np.stack([q[k+"_paths"] for k in S["candidate_paths"]],axis=1)
    assert np.issubdtype(paths.dtype,np.integer)
    z=[];u=[]
    for j in range(paths.shape[1]):
        a,b=local(V,paths[:,j]);z.append(a);u.append(b)
    ix=np.arange(len(V));a=np.stack(z,1);b=np.stack(u,1);j=a.argmax(1);k=b.argmax(1)
    s=q["hybrid_r28_score"];p=q["hybrid_r28_paths"]
    return {"r28_frozen":(s,p),"local_contrast":(a[ix,j],paths[ix,j]),"block_balanced":(b[ix,k],paths[ix,k])}
def cut(z):
    x=np.sort(np.asarray(z,dtype=float))[::-1];assert len(x)==256 and np.isfinite(x).all()
    # explicit empirical budget, no floating point boundary ambiguity in float32
    t=float(np.nextafter(x[int(len(x)*S["cal_control_fraction"])],np.inf))
    assert int((x>=t).sum())<=12;return t
def combine(z,L):
    k=["r28_frozen","block_balanced"];a=np.stack([(np.asarray(z[v][0],dtype=float)-L[v]["center"])/L[v]["scale"] for v in k],axis=1)
    j=a.argmax(1);i=np.arange(len(j));p=np.stack([z[v][1] for v in k],axis=1)
    return a[i,j],p[i,j]
def met(q,s,p,th):
    y=q["y"];amp=q["amp"];act=q["active"].astype(bool);truth=q["truth"]
    pred=np.asarray(s,dtype=float)>=th
    loc=((abs(p-truth)<=2)&act).sum(1)>=.75*act.sum(1)
    hit=pred&loc&(y==1);weak=np.isin(amp,S["weak_amplitudes"])
    return {"localized":int(hit.sum()),"weak":int((hit&weak).sum()),"control_alerts":int((pred&(y==0)).sum()),"positives":int((y==1).sum()),"controls":int((y==0).sum())},hit,pred
def selftest():
    rng=np.random.default_rng(401);x=rng.normal(size=(4,16,128));p=np.full((4,16),64);p[0]=0;p[1]=127
    a,b=local(x,p);xx=x.copy();xx[np.arange(4)[:,None],np.arange(16),p]+=5
    aa,bb=local(xx,p);assert (aa>a).all() and (bb>b).all()
    # A strong one-row impulse must be penalized relative to distributed evidence.
    x=np.zeros((2,16,128));p=np.full((2,16),64);x[0,0,64]=16;x[1,:,64]=1
    a,b=local(x,p);assert np.isclose(a[0],a[1]) and b[1]>b[0]
    return {"edge_no_wrap":"PASS","distributed_signal":"PASS","one_row_impulse":"PASS"}
def main():
    begin=time.perf_counter();checks=selftest();z,n,rec,orig=restore()
    cal=[r for r in rec["records"] if r["name"].startswith("cal_")];assert len(cal)==2
    parts=[]
    for r in cal:
        q=load(z,n,r);parts.append((q,scores(q)))
    L={}
    for k in S["arms"][:-1]:
        x=np.concatenate([np.asarray(v[k][0],dtype=float)[q["y"]==0] for q,v in parts])
        t=cut(x);L[k]={"threshold":t,"center":t,"scale":max(.1,t-float(np.quantile(x,.8))),"controls":len(x),"alerts":int((x>=t).sum())}
    # Original R35G operates at its own archived lock; do not relabel baseline.
    old=float(orig["thresholds"]["hybrid_r28"])
    assert abs(old-L["r28_frozen"]["threshold"])<1e-6
    L["r28_frozen"]["threshold"]=old
    x=np.concatenate([combine(v,L)[0][q["y"]==0] for q,v in parts]);t=cut(x)
    L["r28_plus_balanced"]={"threshold":t,"alerts":int((x>=t).sum()),"controls":len(x)}
    save("lock.json",{"thresholds":L,"source_lock":rec["lock_sha256"],"protocol":sha(P.read_bytes())})
    lh=sha((OUT/"lock.json").read_bytes());emit("R40A_LOCK",L)
    rows=[]
    for r in rec["records"]:
        if r["name"].startswith("cal_"):continue
        assert sha((OUT/"lock.json").read_bytes())==lh
        q=load(z,n,r);st=time.perf_counter();v=scores(q);v["r28_plus_balanced"]=combine(v,L);elapsed=time.perf_counter()-st
        per={};hits={};alerts={};arr={"y":q["y"],"truth":q["truth"],"active":q["active"],"amp":q["amp"]}
        for k,(s,p) in v.items():
            per[k],hits[k],alerts[k]=met(q,s,p,L[k]["threshold"])
            arr[k+"_score"]=s;arr[k+"_path"]=p;arr[k+"_alert"]=alerts[k];arr[k+"_localized"]=hits[k]
        assert np.array_equal(hits["r28_frozen"],q["hybrid_r28_localized"]) and np.array_equal(alerts["r28_frozen"],q["hybrid_r28_alert"])
        b=hits["r28_frozen"];a=alerts["r28_frozen"];y=q["y"];weak=np.isin(q["amp"],S["weak_amplitudes"])
        for k in per:
            h=hits[k];p=alerts[k];per[k].update({"added":int((h&~b).sum()),"lost":int((b&~h).sum()),"weak_added":int((h&~b&weak).sum()),"weak_lost":int((b&~h&weak).sum()),"controls_added":int((p&~a&(y==0)).sum()),"controls_removed":int((a&~p&(y==0)).sum())})
        np.savez_compressed(OUT/r["name"],**arr)
        row={"record":r["name"],"metrics":per,"score_wall_seconds":elapsed};rows.append(row);save("partial.json",rows);emit("R40A_GROUP",row)
    def total(prefix):
        rr=[v["metrics"] for v in rows if v["record"].startswith(prefix)];assert len(rr)==4
        return {k:{m:sum(x[k][m] for x in rr) for m in rr[0][k]} for k in S["arms"]}
    e=total("eval_");shared=total("shared_");b=e["r28_frozen"];bs=shared["r28_frozen"]
    assert b["positives"]==b["controls"]==512 and bs["positives"]==128
    g={}
    for k in S["arms"][1:]:
        v=e[k];w=shared[k]
        g[k]={"weak_gain":v["weak"]-b["weak"],"control_reduction":b["control_alerts"]-v["control_alerts"],"total_gain":v["localized"]-b["localized"],"shared_gain":w["localized"]-bs["localized"],"development_dual_gate":v["weak"]>b["weak"] and v["control_alerts"]<b["control_alerts"] and v["localized"]>=b["localized"] and w["localized"]>=bs["localized"]}
    d={"status":"DEVELOPMENT_ONLY","checks":checks,"eval":e,"shared":shared,"gates":g,"rows":rows,"lock_sha256":lh,"wall_seconds":time.perf_counter()-begin,"scope":S["limits"],"production_adoption":False,"discovery_relaunch_allowed":False}
    save("result.json",d);emit("R40A_SUMMARY",{k:v for k,v in d.items() if k!="rows"})
if __name__=="__main__":main()
