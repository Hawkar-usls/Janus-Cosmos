#!/usr/bin/env python3
"""R26C: physical-grid checked ON/OFF registration with independent references."""
import hashlib,json,pathlib,sys,time,statistics,itertools
import numpy as np
import r24_odontocrystal as od
import r17c_frozen_holdout as hold
import r9_expanded_search as eng
import os,io,zipfile,requests
P=pathlib.Path("research/seti/r26c_reference_protocol.json")
S=json.loads(P.read_text());OUT=pathlib.Path("runtime/r26c");OUT.mkdir(parents=True,exist_ok=True)
hold.S={**hold.S,"seed":S["seed"],"exclude":sorted(set(hold.S["exclude"]+S["exclude_extra"]))}
hold.OUT=OUT;hold.P=P;eng.OUT=OUT
def emit(k,d):print(k,json.dumps(d),flush=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def norm(x):
    m=np.median(x,axis=-1,keepdims=True);s=1.4826*np.median(abs(x-m),axis=-1,keepdims=True)
    assert np.all(s>0) and np.isfinite(s).all()
    return (x-m)/s
def register(X,R):
    a=np.median(X,axis=1);b=np.median(R,axis=1);A=np.maximum(a-1,0);B=np.maximum(b-1,0)
    cor=[];prof=[]
    for shift in S["registration"]["shift_bins"]:
        v=np.roll(B,shift,axis=1);p=np.roll(b,shift,axis=1)
        if shift>0:v[:,:shift]=0;p[:,:shift]=0
        if shift<0:v[:,shift:]=0;p[:,shift:]=0
        cor.append((A*v).sum(1)/(np.sqrt((A*A).sum(1)*(v*v).sum(1))+1e-9));prof.append(p)
    C=np.stack(cor,axis=1);j=C.argmax(1);conf=C.max(1)
    clean=(a.max(1)<1.5)&(b.max(1)<1.5);g=np.where(clean,0,np.where(conf>=.8,1,2))
    profile=np.stack(prof,axis=1)[np.arange(len(X)),j]
    res=X-np.maximum(profile,0)[:,None,:];res[g!=1]=X[g!=1]
    return res,g,np.asarray(S["registration"]["shift_bins"])[j],conf
def common(a,b,X,R):
    # Constrain to exact overlapping integer bins; no resampling of measured power.
    df=a["df_hz"];assert np.isclose(df,b["df_hz"],rtol=1e-12,atol=1e-12)
    off=(a["channel0_mhz"]-b["channel0_mhz"])*1e6/df;k=int(round(off));assert abs(off-k)<.001
    lo=max(0,-k);hi=min(X.shape[1],R.shape[1]-k)
    assert hi-lo>=128 and np.isclose(a["ts"],b["ts"],rtol=1e-9)
    return X[:,lo:hi],R[:,lo+k:hi+k],lo
def tiles(f):
    out=[]
    for t in eng.get_tiles(f):
        x=np.load(OUT/t["artifact"],allow_pickle=False)
        assert list(x.shape)==t["shape"] and hashlib.sha256(x.tobytes()).hexdigest()==t["slice_sha256"]
        out.append((t,x))
    return out
def gather(q):
    fs=q["files"];assert [f["role"] for f in fs]==["ON1","OFF1","ON2","OFF2"]
    rng=np.random.default_rng(int.from_bytes(hashlib.sha256((S["seed"]+"|"+q["target"]).encode()).digest()[:8],"big"))
    plan=list(itertools.product(S["amplitudes"],S["injection_models"],S["gap_levels"]))*4;rng.shuffle(plan)
    rows=[];refs=[];challenges=[];sources=set();k=0
    for i in (0,2):
        ti=tiles(fs[i]);ri=tiles(fs[i+1]);assert len(ti)==len(ri)==4
        for (a,X),(b,R) in zip(ti,ri):
            assert a["center_mhz"]==b["center_mhz"]; X,R,lo=common(a,b,X,R);assert X.shape[0]==R.shape[0]==16
            sources.update([a["slice_sha256"],b["slice_sha256"]])
            # Hard patches picked from ON background before injection, not used to estimate labels.
            z=norm(X.astype(float));rank=np.argsort(-np.median(z,axis=0),kind="stable");hard=[]
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
                    rows.append((xx,label,model,gap,amp,trial>=12,a["slice_sha256"],p.copy(),active.copy()));refs.append(r)
                if model=="stationary":
                    # Same artificial stationary source in reference: explicit leakage stress, not ordinary calibration.
                    rr=r.copy();ss=1.4826*np.median(abs(r-np.median(r,axis=1)[:,None]),axis=1)
                    rr[t[active],p[active]]+=amp*ss[active]
                    challenges.append((rows[-1],rr))
    assert len(rows)==256 and k==128 and len(challenges)==32
    return rows,np.stack(refs),challenges,sources
def calc(X,R):
    X=norm(X);R=norm(R);V,g,sh,conf=register(X,R);unaligned=X-np.maximum(np.median(R,axis=1),0)[:,None,:]
    a={}
    for name,w in (("raw",X),("unaligned",unaligned),("registered",V)):
        t=time.perf_counter();s,p,m=od.sliced(w,0,None);a[name]=(s,p,m,time.perf_counter()-t)
    return a,g,sh,conf
def cut(x):
    z=np.sort(np.asarray(x,dtype=float))[::-1];assert len(z)>0
    return float(np.nextafter(z[int(len(z)*.05)],np.inf))
def lock_fit(groups):
    y=np.concatenate([np.asarray([r[1] for r in rows]) for rows,R,_,_ in groups]);sc={};gg=[]
    for rows,R,_,_ in groups:
        a,g,_,_=calc(np.stack([r[0] for r in rows]),R);gg.append(g)
        for k,v in a.items():sc.setdefault(k,[]).append(v[0])
    g=np.concatenate(gg);counts=[int(np.sum((y==0)&(g==i))) for i in range(3)];pool={k:cut(np.concatenate(v)[y==0]) for k,v in sc.items()}
    strata={k:[cut(np.concatenate(v)[(y==0)&(g==i)]) if counts[i]>=30 else pool[k] for i in range(3)] for k,v in sc.items()}
    return {"pooled":pool,"strata":strata,"counts":counts,"sparse_fallback":[i for i,n in enumerate(counts) if n<30]}
def score(rows,R,L):
    t=time.perf_counter();a,g,sh,conf=calc(np.stack([r[0] for r in rows]),R);out={};dec={}
    for k in S["modes"]:
        base,rule=k.rsplit("_",1);s,p,m,w=a[base];th=L["pooled"][base] if rule=="pooled" else np.asarray(L["strata"][base])[g]
        pred=s>=th;out[k]={**od.measure(pred,p,m,rows),"dp_seconds":w}
        y=np.asarray([r[1] for r in rows]);truth=np.stack([r[7] for r in rows]);active=np.stack([r[8] for r in rows])
        loc=((abs(p-truth)<=2)&active).sum(1)>=.75*active.sum(1);dec[k]=(pred,loc)
        out[k]["regimes"]={str(i):{"localized":int(np.sum(pred&loc&(y==1)&(g==i))),"positives":int(np.sum((y==1)&(g==i))),"control_alerts":int(np.sum(pred&(y==0)&(g==i))),"controls":int(np.sum((y==0)&(g==i)))} for i in range(3)}
    for k,(p,l) in dec.items():
        b,bl=dec["raw_strata"];y=np.asarray([r[1] for r in rows]);q=p&l&(y==1);ref=b&bl&(y==1)
        out[k]["paired_localized"]={"added":int(np.sum(q&~ref)),"lost":int(np.sum(ref&~q))}
    return out,{"regimes":np.bincount(g,minlength=3).tolist(),"shift_counts":{str(i):int(np.sum((g==1)&(sh==i))) for i in range(-3,4)},"seconds_full_panel":time.perf_counter()-t}
def selftest():
    rng=np.random.default_rng(260301);x=rng.normal(0,.1,(1,16,128));r=rng.normal(0,.1,(1,16,128))
    x[:,:,30]+=6;x[:,:,60]+=4;r[:,:,32]+=6;r[:,:,62]+=4
    v,g,sh,c=register(x,r);assert g[0]==1 and sh[0]==-2 and np.max(abs(np.median(v,axis=1)))<.2
    x2=x.copy();r2=r.copy();r2[:,:,32]-=6;r2[:,:,62]-=4;r2[:,:,90]+=8
    v,g,_,_=register(x2,r2);assert g[0]==2 and np.array_equal(v,x2)
    t={"df_hz":2.,"channel0_mhz":1000.,"ts":1.};b={**t,"channel0_mhz":1000.000004}
    a,rr,lo=common(t,b,np.zeros((16,140)),np.ones((16,140)));assert lo==2 and a.shape==rr.shape==(16,138)
    return {"shift_recovery":"PASS","unmatched_fallback":"PASS","physical_grid_overlap":"PASS"}
BLC_P=json.loads(pathlib.Path("research/seti/r15d_r16d_r17d_protocol.json").read_text())
BLC_OUT=OUT/"blc"
def archive(key):
    a=BLC_P["artifacts"][key];path=BLC_OUT/(key+".zip");h=hashlib.sha256();n=0
    with requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['id']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},stream=True,timeout=(15,120)) as r:
        r.raise_for_status()
        with path.open("wb") as f:
            for b in r.iter_content(1024*1024):
                n+=len(b);assert n<=500*1024**2;h.update(b);f.write(b)
    assert h.hexdigest()==a["digest"]
    z=zipfile.ZipFile(path);assert sum(i.file_size for i in z.infolist())<=1200*1024**2
    names={}
    for i in z.infolist():
        if i.is_dir():continue
        b=pathlib.PurePosixPath(i.filename).name
        assert b not in names;names[b]=i.filename
    return z,names,path
def get(z,n,k):return json.loads(z.read(n[k]))
def tile(z,n,t,k=None):
    x=np.load(io.BytesIO(z.read(n[k or t["artifact"]])),allow_pickle=False)
    assert list(x.shape)==t["shape"] and hashlib.sha256(x.tobytes()).hexdigest()==t["slice_sha256"]
    return x

def blc_search(y,t):
    s=BLC_P["r15d"]
    nr,nc=y.shape;df=t["df_hz"];ts=t["ts"];best=[]
    for width in sorted(set(min(nr,w) for w in s["windows"])):
        ans=None
        for start in sorted(set(list(range(0,nr-width+1,max(1,width//2)))+[nr-width])):
            tt=np.arange(width)*ts
            for d in s["drift_grid_hz_s"]:
                for a in s["drift_change_grid_hz_s"]:
                    # d is initial rate of segment; a is change over full scan.
                    shifts=np.rint((d*tt+0.5*a*tt**2/max((nr-1)*ts,1))/df).astype(int)
                    lo=max(0,-int(shifts.min()));hi=min(nc,nc-int(shifts.max()))
                    if hi<=lo:continue
                    vals=sum(y[start+i,lo+k:hi+k] for i,k in enumerate(shifts))
                    med=np.median(vals);sd=1.4826*np.median(np.abs(vals-med))
                    if sd<=0:continue
                    scores=(vals-med)/sd
                    freq=t["channel0_mhz"]+np.arange(lo,hi)*df/1e6
                    ix=np.flatnonzero(np.abs(freq-s["center_mhz"])*1e6<=s["half_hz"])
                    if not len(ix):continue
                    j=int(ix[np.argmax(scores[ix])]);v=float(scores[j])
                    record={"score":v,"frequency_mhz":float(freq[j]),"drift_hz_s":d,"drift_change_hz_s":a,"start":start,"rows":width}
                    if ans is None or v>ans["score"]:ans=record
        best.append(ans)
    return best

def blc_replay():
    BLC_OUT.mkdir(exist_ok=True)
    z,n,p=archive("blc");m,mn,mp=archive("blc_manifest")
    manifest=get(m,mn,"manifest.json");result=get(z,n,"result.json")
    rr=[r for r in result["rows"] if "tile" in r];assert len(rr)==15
    hs={f["name"]:f["header"] for f in manifest["selected"]}
    data={r["name"]:tile(z,n,r["tile"],r["artifact"]) for r in rr}
    results=[]
    for row in rr:
        off=[r for r in rr if r["role"].startswith("OFF") and r["name"]!=row["name"]]
        off.sort(key=lambda r:abs(hs[r["name"]]["tstart"]-hs[row["name"]]["tstart"]))
        ref=off[0];x,r,lo=common(row["tile"],ref["tile"],data[row["name"]],data[ref["name"]])
        tt={**row["tile"],"channel0_mhz":row["tile"]["channel0_mhz"]+lo*row["tile"]["df_hz"]/1e6}
        X=eng.normalized(x)[None].astype(float);R=eng.normalized(r)[None].astype(float);v,g,sh,c=register(X,R)
        raw=blc_search(X[0],tt);clean=blc_search(v[0],tt)
        results.append({"name":row["name"],"role":row["role"],"reference":ref["name"],"regime":int(g[0]),"shift":int(sh[0]),"similarity":float(c[0]),"raw":raw,"registered":clean})
    anchor=next(r for r in results if r["name"]==BLC_P["r15d"]["anchor"])
    def match(r):
        return bool(any(x is not None and x["score"]>=10 and abs(x["frequency_mhz"]-982.002472)*1e6<=10 and abs(x["drift_hz_s"]-.028)<=.012 for x in r))
    d={"rows":results,"anchor":anchor,"anchor_raw_recovered":match(anchor["raw"]),"anchor_registered_recovered":match(anchor["registered"]),"scope":"Targeted replay of known BLC1, terrestrial RFI per published investigation; not blind or source classification. Any failure is retained. OFF scans are not certified signal-free."}
    save("blc_result.json",d);emit("R26C_BLC",{k:v for k,v in d.items() if k!="rows"})
    z.close();m.close();p.unlink();mp.unlink()
def main():
    checks=selftest();m=json.loads((OUT/"manifest.json").read_text());assert m["protocol_sha256"]==hashlib.sha256(P.read_bytes()).hexdigest() and len(m["selected"])==6
    assert len({q["target"] for q in m["selected"]})==6 and all(q["target"] not in hold.S["exclude"] for q in m["selected"])
    cal=[gather(q) for q in m["selected"][:2]];L=lock_fit(cal);L["targets"]=[q["target"] for q in m["selected"][:2]];save("lock.json",L);emit("R26C_LOCK",L)
    seen=set()
    for _,_,_,h in cal:assert not seen&h;seen|=h
    out={};diag={};leak={}
    for q in m["selected"][2:]:
        t=time.perf_counter();rows,R,ch,h=gather(q);load=time.perf_counter()-t;assert not seen&h;seen|=h
        out[q["target"]],diag[q["target"]]=score(rows,R,L);diag[q["target"]]["load_seconds"]=load
        leak[q["target"]]=score([x[0] for x in ch],np.stack([x[1] for x in ch]),L)[0]
        save("partial.json",{"metrics":out,"diagnostics":diag});emit("R26C_TARGET",{"target":q["target"],"metrics":out[q["target"]]})
    keys=("detected","localized","positives","control_alerts","controls","dp_seconds")
    total={k:{a:sum(x[k][a] for x in out.values()) for a in keys} for k in S["modes"]}
    for k in total:
        for f,vals in (("amplitudes",S["amplitudes"]),("gaps",S["gap_levels"]),("models",S["injection_models"]),("regimes",[0,1,2])):
            total[k][f]={str(v):{a:sum(x[k][f][str(v)][a] for x in out.values()) for a in (("localized","positives","control_alerts","controls") if f=="regimes" else keys[:5])} for v in vals}
        total[k]["paired_localized"]={a:sum(x[k]["paired_localized"][a] for x in out.values()) for a in ("added","lost")}
    d={"status":"FRESH_ARCHIVE_SYNTHETIC_INJECTIONS","samples":1024,"targets":list(out),"selftests":checks,"lock":L,"modes":total,"diagnostics":diag,"leak_challenge":{k:{a:sum(v[k][a] for v in leak.values()) for a in ("detected","localized","positives")} for k in total},"per_target":out,"scope":"Actual measured ON/OFF backgrounds, synthetic labels; no externally blind labels. OFF independently acquired, not physically independent of shared receiver. Registration can subtract sky signal; raw retained, no automatic veto. Source unknown.","production_adoption":False}
    save("result.json",d);emit("R26C_SUMMARY",{k:v for k,v in d.items() if k!="per_target"})
if __name__=="__main__":
    if sys.argv[1]=="freeze":hold.freeze()
    elif sys.argv[1]=="scan":main()
    elif sys.argv[1]=="blc":blc_replay()
    elif sys.argv[1]=="selftest":emit("R26C_SELFTEST",selftest())
    else:raise ValueError("phase")
