#!/usr/bin/env python3
"""R40N authentic raw ON/OFF rival trajectories and physically mapped repeat witness."""
import sys
MODE=sys.argv[1];assert MODE in ("replay","controls")
import json,pathlib,time,itertools,hashlib
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import ExtraTreesClassifier
import r40l_joint_shield as l
k=l.k;h=l.h;c=l.c;t=l.t;f=l.f;b=l.b;ref=b.ref
P=pathlib.Path("research/seti/r40n_raw_rivals_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r40n")/MODE;O.mkdir(parents=True,exist_ok=True)
l.O=O;k.O=O
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def archive(key):return c.artifact(S["sources"][key],32)
def read(z,n,k):return json.loads(z.read(n[k]))
def gather(q,tag):
    w=O/'raw'/tag;w.mkdir(parents=True,exist_ok=True);ref.OUT=w;b.sh.eng.OUT=w
    pre=[ref.tiles(f) for f in q['files']]
    meta=[];cm=[]
    fs=q["files"];assert [f["role"] for f in fs]==["ON1","OFF1","ON2","OFF2"]
    rng=np.random.default_rng(int.from_bytes(hashlib.sha256((ref.S["seed"]+"|"+q["target"]).encode()).digest()[:8],"big"))
    plan=list(itertools.product(ref.S["amplitudes"],ref.S["injection_models"],ref.S["gap_levels"]))*4;rng.shuffle(plan)
    rows=[];refs=[];challenges=[];sources=set();k=0
    for i in (0,2):
        ti=pre[i];ri=pre[i+1];assert len(ti)==len(ri)==4
        for wi,((a,X),(bb,R)) in enumerate(zip(ti,ri)):
            assert a["center_mhz"]==bb["center_mhz"]; X,R,lo=ref.common(a,bb,X,R);assert X.shape[0]==R.shape[0]==16
            sources.update([a["slice_sha256"],bb["slice_sha256"]])
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
                    rows.append((xx,label,model,gap,amp,trial>=12,a["slice_sha256"],p.copy(),active.copy()));refs.append(r);meta.append((i,wi,a['channel0_mhz']+(lo+start)*a['df_hz']/1e6,a['df_hz']))
                if model=="stationary":
                    # Same artificial stationary source in reference: explicit leakage stress, not ordinary calibration.
                    rr=r.copy();ss=1.4826*np.median(abs(r-np.median(r,axis=1)[:,None]),axis=1)
                    rr[t[active],p[active]]+=amp*ss[active]
                    challenges.append((rows[-1],rr));cm.append(meta[-1])
    assert len(rows)==256 and k==128 and len(challenges)==32
    return rows,np.stack(refs),challenges,sources,meta,cm,pre

def set_config(seed,windows,half,budget):
    b.S["seed"]=seed;ref.S=b.S
    b.sh.eng.SPEC={**b.sh.eng.SPEC,"window_centers_mhz":windows,"window_half_hz":half,"limits":{**b.sh.eng.SPEC["limits"],"spectrum_budget_per_file_mib":budget}}
def normalized_contrast(V):
    idx=np.arange(5,123);u=np.zeros_like(V);side=V[:,:,idx[:,None]+np.array([-5,-4,-3,3,4,5])]
    u[:,:,idx]=np.clip(V[:,:,idx]-np.median(side,axis=3),-4,4)
    return u
def rivals(X,D,p):
    # Enumerate every nonwrapping R35B path, excluding the selected candidate basin.
    U=normalized_contrast(X);V=normalized_contrast(D);n=len(X);row=np.arange(16);ix=np.arange(n)
    base=np.minimum(U[ix[:,None],row[None,:],p].mean(1),V[ix[:,None],row[None,:],p].mean(1))
    best=np.full((n,3),-np.inf);bestp=np.zeros((n,16),int)
    for d in b.D:
        lo=5-int(d.min());hi=123-int(d.max())
        if hi<=lo:continue
        starts=np.arange(lo,hi);pp=starts[:,None]+d
        distinct=(abs(pp[None,:,:]-p[:,None,:])>2).sum(2)>=5
        a=sum(U[:,j,lo+v:hi+v] for j,v in enumerate(d))/16
        b0=sum(V[:,j,lo+v:hi+v] for j,v in enumerate(d))/16
        for j,score in enumerate((np.minimum(a,b0),a,b0)):
            s=np.where(distinct,score,-np.inf);pos=s.argmax(1);v=s[ix,pos];m=v>best[:,j]
            best[m,j]=v[m]
            if j==0:bestp[m]=pp[pos[m]]
    assert np.isfinite(best).all() and ((abs(bestp-p)>2).sum(1)>=5).all()
    witness=np.minimum(U[ix[:,None],row[None,:],bestp].mean(1),V[ix[:,None],row[None,:],bestp].mean(1))
    assert np.allclose(witness,best[:,0],atol=1e-12,rtol=1e-12)
    return np.column_stack([best,base[:,None]-best]),bestp
def shape(V,p):
    i=np.arange(len(V))[:,None];j=np.arange(16)[None,:];v=V[i,j,p]
    def band(ds):
        idx=p[:,:,None]+np.array(ds);valid=(idx>=0)&(idx<128);x=V[i[:,:,None],j[:,:,None],np.clip(idx,0,127)]
        cnt=valid.sum((1,2));assert (cnt>0).all()
        return (np.where(valid,np.clip(x,-6,6),0).sum((1,2))/cnt)
    narrow=band([-2,-1,1,2]);near=band([-5,-4,-3,3,4,5]);wide=band(list(range(-12,-7))+list(range(8,13)))
    return np.stack([np.clip(v,-6,6).mean(1)-narrow,narrow-near,near-wide,wide,(V>3).mean((1,2))],1)
def repeat(meta,pre,p):
    out=np.zeros((len(meta),5));freq=[]
    for j,(i,wi,f0,df) in enumerate(meta):
        assert (p[j]==p[j,0]).all()
        fq=f0+p[j,0]*df/1e6;freq.append(fq)
        oth=2-i;rr=[]
        for slot in (oth,oth+1):
            a,V=pre[slot][wi];vdf=a["df_hz"]
            if not np.isclose(df,vdf,rtol=1e-12,atol=1e-12):break
            off=(f0-a["channel0_mhz"])*1e6/vdf;idx=int(round(off))
            if abs(off-idx)>=.001 or idx<0 or idx+128>V.shape[1]:break
            Z=ref.norm(V[:,idx:idx+128].astype(float));path=int(p[j,0]);near=Z[:,path+np.array([-5,-4,-3,3,4,5])]
            x=np.clip(Z[:,path]-np.median(near,axis=1),-4,4);rr.append((x.mean(),(x>0).mean()))
        if len(rr)==2:out[j]=[1,rr[0][0],rr[1][0],rr[0][1],rr[1][1]]
    return out,np.asarray(freq)
def encode(rows,R,meta,pre,old=None):
    X=ref.norm(np.stack([v[0] for v in rows]));D=ref.norm(R);Y,_,_,_=ref.register(X,D)
    q={"V":Y.astype(np.float32),"y":np.array([v[1] for v in rows]),"truth":np.stack([v[7] for v in rows]),"active":np.stack([v[8] for v in rows]),"amp":np.array([v[4] for v in rows])}
    s,p,sup,reg,F=l.evidence(rows,R,q)
    correct=l.correct(q,p);weak=np.isin(q["amp"],t.S["weak_amplitudes"])
    if old is not None:
        for key,v in {"s":s,"p":p,"support":sup,"regime":reg,"F":F,"y":q["y"],"correct":correct,"weak":weak}.items():
            assert np.array_equal(v,old[key]),"R40L parity "+key
    rv,rp=rivals(X,D,p);rep,freq=repeat(meta,pre,p)
    a=np.column_stack([F,rv,shape(X,p),shape(D,p)])
    assert np.isfinite(a).all() and np.isfinite(rep).all()
    return {**(old or {}),"s":s,"p":p,"support":sup,"regime":reg,"F":F,"y":q["y"],"correct":correct,"weak":weak,"raw_features":a,"repeat_features":rep,"rival_path":rp,"frequency_mhz":freq,"source_tile_sha":np.array([v[6] for v in rows])}
def selftest():
    rng=np.random.default_rng(4019);X=rng.normal(0,.1,(2,16,128));D=X.copy();p=np.full((2,16),64)
    X[:,:,64]+=4;D[:,:,64]+=4
    alt=30+np.rint(np.arange(16)/2).astype(int);X[1,np.arange(16),alt]+=5;D[1,np.arange(16),alt]+=5
    f0,rp=rivals(X,D,p);assert f0[0,0]<1 and (abs(rp[1]-alt)<=1).sum()>=12
    a={"df_hz":-2.8,"channel0_mhz":3000.};v=np.random.default_rng(7).normal(size=(16,256));v[:,74]+=10;pre=[[(a,v)]]*4
    rr,fq=repeat([(0,0,3000.-10*2.8/1e6,-2.8)],pre,np.full((1,16),64));assert rr[0,0]==1 and rr[0,1]>3 and abs(fq[0]-(3000.-74*2.8/1e6))<1e-9
    emit("R40N_SELFTEST",{"rival_exclusion":"PASS","alternative_drift":"PASS","signed_frequency_mapping":"PASS"})
def manifests():
    hz,hn=archive("historical_manifest");hm=read(hz,hn,"manifest.json")
    zz,nn=archive("r40h_freeze");m=read(zz,nn,"manifest.json");fr=read(zz,nn,"freeze.json");sel=read(zz,nn,"selection_receipt.json")
    assert c.sha(zz.read(nn["manifest.json"]))==sel["manifest_sha256"]
    return hm,m,fr,sel
def replay():
    st=time.perf_counter();selftest();hm,m,fr,sel=manifests()
    lz,ln=archive("r40l");raw=lz.read(ln["lock.json"]);assert c.sha(raw)==S["sources"]["r40l"]["lock_sha256"];L=json.loads(raw)
    kz,kn=archive("r40k");kl=read(kz,kn,"lock.json");kr=read(kz,kn,"result.json");th=L["parent_joint_min_threshold"]
    train=[];cal=[];seen=set()
    def group(g,historical):
        name=g["target"];rows,R,ch,hs,meta,cm,pre=gather(g,("hist_" if historical else "")+name)
        assert len(hs)==16 and not seen&hs;seen.update(hs)
        split=g.get("split")
        if historical and split=="null_calibration":assert sorted(hs)==next(v["raw_tile_hashes"] for v in kl["calibration"] if v["target"]==name)
        if not historical:assert sorted(hs)==next(v["raw_tile_hashes"] for v in kr["rows"] if v["target"]==name)
        pairs=[("cal" if historical and split=="null_calibration" else "eval",rows,R,meta)]
        if pairs[0][0]!="cal":pairs.append(("shared",[v[0] for v in ch],np.stack([v[1] for v in ch]),cm))
        for kind,rr,ref0,mm in pairs:
            key=("historical_" if historical else "")+kind+"_"+name+".npz";old=f.load(lz,ln,key)
            d=encode(rr,ref0,mm,pre,old)
            assert np.array_equal((~d["parent_alert"])&(d["s"]>=th),d["joint_rescue"])
            np.savez_compressed(O/key,**d);emit("R40N_RAW_GROUP",{"record":key,"rows":len(rr),"prior_exact":True,"recurrence_coverage":int(d["repeat_features"][:,0].sum()),"tiles":len(hs)})
            yield key,d
    for g in hm["selected"]:
        if g["split"]=="component_fit":continue
        for key,d in group(g,True):(cal if g["split"]=="null_calibration" else train).append((key,d))
    assert len(train)==8 and len(cal)==2
    def use(d):return ~d["parent_alert"]
    def null(d):return use(d)&(d["y"]==0)
    def pos(d):return d["joint_rescue"]&d["correct"]
    y=np.concatenate([d["correct"][use(d)] for _,d in train]);assert y.any() and (~y).any()
    models={};cuts={};serial={}
    for feature in ("raw","repeat"):
        def fv(d,feature=feature):return d["raw_features"] if feature=="raw" else np.column_stack([d["raw_features"],d["repeat_features"]])
        X=np.concatenate([fv(d)[use(d)] for _,d in train])
        scaler=StandardScaler().fit(X)
        for kind in ("logistic","tree"):
            model=LogisticRegression(C=.1,class_weight="balanced",max_iter=1000,random_state=401).fit(scaler.transform(X),y) if kind=="logistic" else ExtraTreesClassifier(n_estimators=128,max_depth=4,min_samples_leaf=8,class_weight="balanced",random_state=401,n_jobs=1).fit(X,y)
            def score(d,model=model,scaler=scaler,fv=fv,kind=kind):
                return model.decision_function(scaler.transform(fv(d))) if kind=="logistic" else model.predict_proba(fv(d))[:,1]
            name=feature+"_"+kind;models[name]=score
            gp=np.concatenate([score(d)[pos(d)] for _,d in train]);gc=np.concatenate([score(d)[null(d)] for _,d in cal]);assert len(gp)>0 and len(gc)==246
            cuts[name]={"preserve":float(np.nextafter(gp.min(),-np.inf)),"null025":l.cut(gc,.025)}
            serial[name]={"mean":scaler.mean_.tolist(),"scale":scaler.scale_.tolist(),"coef":model.coef_.tolist(),"intercept":model.intercept_.tolist()} if kind=="logistic" else {"trees":[{"children_left":e.tree_.children_left.tolist(),"children_right":e.tree_.children_right.tolist(),"feature":e.tree_.feature.tolist(),"threshold":e.tree_.threshold.tolist(),"value":e.tree_.value.tolist()} for e in model.estimators_]}
    lock={"protocol_sha256":c.sha(P.read_bytes()),"source_lock_sha256":c.sha(raw),"cuts":cuts,"models":serial,"train_rows":len(y),"train_positive":int(y.sum()),"calibration_nulls":246,"bank_paths":len(b.D)}
    save("lock.json",lock);lh=c.sha((O/"lock.json").read_bytes());emit("R40N_LOCK",{k:v for k,v in lock.items() if k!="models"})
    set_config(h.S["seed"],h.S["windows_mhz"],h.S["half_hz"],h.S["max_mib_per_scan"]);results=[]
    for g in m["selected"]:
        for key,d in group(g,False):
            a=d["parent_alert"];base=d["parent_hit"];res=d["joint_rescue"];old=d["joint_min_hit"];y=d["y"];weak=d["weak"];pre={"joint_min":np.ones(len(y),bool),"r40g":np.zeros(len(y),bool),"learned05":d["learned05_alert"]&~a}
            for name,score in models.items():
                z=score(d)
                for cut,v in cuts[name].items():pre[name+"_"+cut]=z>=v
            out={}
            for name,guard in pre.items():
                rescue=res&guard;aa=a|rescue;hh=base|(rescue&d["correct"])
                out[name]={"localized":int(hh.sum()),"weak":int((hh&weak).sum()),"control_alerts":int((aa&(y==0)).sum()),"lost_parent":int((base&~hh).sum()),"lost_joint":int((old&~hh).sum()),"lost_joint_weak":int((old&~hh&weak).sum()),"positives":int((y==1).sum()),"controls":int((y==0).sum())}
                d[name+"_alert"]=aa;d[name+"_hit"]=hh
            np.savez_compressed(O/key,**d);row={"record":key,"metrics":out};results.append(row);save("partial.json",results);emit("R40N_GROUP",row)
        assert c.sha((O/"lock.json").read_bytes())==lh
    def total(kind):
        vv=[r["metrics"] for r in results if r["record"].startswith(kind+"_")];assert len(vv)==5
        return {a:{k:sum(x[a][k] for x in vv) for k in vv[0][a]} for a in vv[0]}
    e=total("eval");sh=total("shared");assert [e["joint_min"][k] for k in ("localized","weak","control_alerts")]==[215,43,43] and [sh["joint_min"][k] for k in ("localized","weak")]==[111,36]
    gates={a:{"strict":e[a]["localized"]>=215 and e[a]["weak"]>=43 and e[a]["control_alerts"]<=35 and sh[a]["localized"]>=111 and sh[a]["weak"]>=36 and e[a]["lost_joint"]==sh[a]["lost_joint"]==0,"minimum":e[a]["localized"]>=214 and e[a]["weak"]>=41 and e[a]["control_alerts"]<=35 and sh[a]["weak"]>=21 and e[a]["lost_parent"]==sh[a]["lost_parent"]==0} for a in e}
    d={"status":"DEVELOPMENT_REPLAY_ONLY","eval":e,"shared":sh,"gates":gates,"rows":results,"lock_sha256":lh,"seconds":time.perf_counter()-st,"production_adoption":False};save("result.json",d);emit("R40N_SUMMARY",{k:v for k,v in d.items() if k!="rows"})
def controls():
    st=time.perf_counter();selftest();hm,hm2,fr,sel=manifests();rz,rn=archive("r37_manifest");rm=read(rz,rn,"manifest.json")
    C=S["new_controls"];prior=hm["selected"]+hm2["selected"]+rm["selected"]["Parkes"]+rm["selected"].get("GBT",[])
    cde=json.loads(pathlib.Path("research/seti/r40cde_protocol.json").read_text())
    names=set(sel["excluded_names"])|{g["target"] for g in prior}|set(cde["targets"]);days=set(sel["excluded_days"])|{int(g["day"]) for g in prior}|set(cde["days"])
    b.sh.S={**b.S,"seed":C["seed"]};b.sh.P=P;b.sh.C={**b.sh.C,"fit_targets":0,"null_targets":C["targets"],"eval_targets":0,"max_queries":C["max_queries"],"exclude_days":sorted(set(b.sh.C["exclude_days"])|days)}
    b.sh.hold.S={**b.sh.hold.S,"exclude":sorted(set(b.sh.hold.S["exclude"])|names)};b.sh.OUT=O;b.sh.freeze()
    m=json.loads((O/"manifest.json").read_text());assert len(m["selected"])==C["targets"] and not {g["target"] for g in m["selected"]}&names and not {int(g["day"]) for g in m["selected"]}&days
    def coords(g):
        hd=g["files"][0]["header"];return str(hd["src_raj"]),str(hd["src_dej"])
    known={coords(g) for g in prior if g.get("files") and "header" in g["files"][0]};new=[coords(g) for g in m["selected"]];assert len(set(new))==len(new) and not set(new)&known
    mh=c.sha((O/"manifest.json").read_bytes());emit("R40N_CONTROL_MANIFEST",{"manifest_sha256":mh,"targets":[[g["target"],g["day"]] for g in m["selected"]],"coordinate_exact_match_excluded":True})
    set_config(C["seed"],C["windows_mhz"],C["half_hz"],C["max_mib_per_scan"]);L=h.readlocks();b.taylor.init_tree()
    lz,ln=archive("r40l");raw=lz.read(ln["lock.json"]);assert c.sha(raw)==S["sources"]["r40l"]["lock_sha256"];old=json.loads(raw)
    kz,kn=archive("r40k");kr=read(kz,kn,"result.json");kl=read(kz,kn,"lock.json")
    seen=set(fr["excluded_raw_hashes"]);seen.update(v for r in kr["rows"] for v in r["raw_tile_hashes"]);seen.update(v for r in kl["calibration"] for v in r["raw_tile_hashes"])
    rowsout=[];allreg=[];allparent=[];alltarget=[];allscore=[]
    for g in m["selected"]:
        name=g["target"];rr,R,ch,hs,meta,cm,pre=gather(g,name);assert len(hs)==16 and not seen&hs;seen.update(hs)
        ix=np.array([v[1]==0 for v in rr]);assert ix.sum()==128
        rows=[v for v,x in zip(rr,ix) if x];mm=[v for v,x in zip(meta,ix) if x];R=R[ix]
        d=encode(rows,R,mm,pre);pan=b.panel(rows,R);arms=b.base_arms(pan,L["r40c"]["original"]["component_fit"])
        q={"V":pan["Y"].astype(np.float32),"y":np.zeros(128,int),"truth":np.stack([v[7] for v in rows]),"active":np.stack([v[8] for v in rows]),"amp":np.array([v[4] for v in rows])}
        for a in t.S["paths"]:q[a+"_score"],q[a+"_paths"]=arms[a][:2]
        hh=f.historical(q,L["r40c"]);_,_,_,a=h.replay(hh,L);a=a.astype(bool);reg=d["regime"]>old["reference_cut"];d.update({"parent_alert":a,"parent_hit":np.zeros(128,bool),"joint_rescue":(~a)&(d["s"]>=old["parent_joint_min_threshold"]),"regime_high":reg})
        np.savez_compressed(O/(name+".npz"),**d)
        allreg.extend(reg.tolist());allparent.extend(a.tolist());alltarget.extend([name]*len(a));allscore.extend(d["s"].tolist())
        row={"target":name,"day":g["day"],"controls":128,"parent_alerts":int(a.sum()),"joint_extra_alerts":int(d["joint_rescue"].sum()),"joint_max":float(d["s"].max()),"repeat_coverage":int(d["repeat_features"][:,0].sum()),"raw_tile_hashes":sorted(hs)}
        rowsout.append(row);save("partial.json",rowsout);emit("R40N_CONTROL_GROUP",row);assert c.sha((O/"manifest.json").read_bytes())==mh
    rg=np.asarray(allreg);pa=np.asarray(allparent);ta=np.asarray(alltarget);sc=np.asarray(allscore);assert len(rg)==1024
    stats={str(i):{"controls":int(((rg==i)&~pa).sum()),"targets":len(set(ta[(rg==i)&~pa].tolist()))} for i in (0,1)}
    qualified=all(v["controls"]>=128 and v["targets"]>=3 for v in stats.values())
    d={"status":"PASS_KUSTO_STRATUM_COUNTS_ONLY" if qualified else "HOLD_SPARSE_KUSTO_STRATUM","manifest_sha256":mh,"controls":len(rg),"targets":len(rowsout),"parent_control_alerts":int(pa.sum()),"strata":stats,"reference_cut":old["reference_cut"],"stratum_maxima":{str(i):float(sc[(rg==i)&~pa].max()) if ((rg==i)&~pa).any() else None for i in (0,1)},"groups":rowsout,"formal_fap":False,"alias_completeness_proven":False,"seconds":time.perf_counter()-st,"production_adoption":False}
    save("result.json",d);emit("R40N_CONTROL_SUMMARY",{k:v for k,v in d.items() if k!="groups"})
if __name__=="__main__":
    if MODE=="replay":replay()
    else:controls()
