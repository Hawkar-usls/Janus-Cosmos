"""R42A3B full diagnostic score contract, frozen R41 Minimal routing."""
import pathlib,json,hashlib,io,zipfile,os
import numpy as np
import requests
import r42a_common_export as src
base=src.base;a=src.a;n=src.n;l=base.l
P=pathlib.Path("research/seti/r42a3_score_contract_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r42a3b");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(k,v):(O/k).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def transform(m):
    m=np.asarray(m,float);assert np.isfinite(m).all()
    return m/(1+abs(m))
def unpack():
    q=S["source"];u="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/"+str(q["artifact"])
    h={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"}
    def get(u):
        r=requests.get(u,headers=h,timeout=120);r.raise_for_status();return r
    x=get(u).json();assert x["workflow_run"]["id"]==q["run"] and x["digest"]=="sha256:"+q["digest"]
    b=get(u+"/zip").content;assert len(b)<128*1024**2 and sha(b)==q["digest"]
    z=zipfile.ZipFile(io.BytesIO(b));f=[x for x in z.infolist() if not x.is_dir()]
    assert sum(x.file_size for x in f)<512*1024**2
    names={pathlib.PurePosixPath(x.filename).name:x for x in f};assert len(names)==len(f)
    r=json.loads(z.read(names["receipt.json"]));assert len(r["records"])==22 and sum(x["rows"] for x in r["records"])==3168
    return r,lambda k:z.read(names[k])
def group(i):
    assert 1<=i<11
    r,rd=unpack()
    gs=base.plan()["groups"];g=gs[i];num=f"{i+1:02d}"
    assert len([x for x in r["records"] if x["group"]==i+1 and x["fit"]==g["fit"]])==2
    seed={"n6":a.v.C["seed"],"n7":base.w.C["seed"],"n9":base.old9.C["seed"]}[g["source"]]
    n.set_config(seed,a.C["windows_mhz"],a.C["half_hz"],a.C["max_mib_per_scan"]);base.b.taylor.init_tree()
    m,L,old=base.locks("minimal",base.plan())
    rr,R,ch,hs,meta,cm,pre=n.gather(g["group"],"r42a3b_"+num)
    assert sorted(hs)==g["raw_tile_hashes"]
    records=[]
    for kind,rows,refs,mm in (("eval",rr,R,meta),("shared",[x[0] for x in ch],np.stack([x[1] for x in ch]),cm)):
        cap={};orig=l.parent
        def observe(q,locks):
            hh=l.f.historical(q,locks["r40c"]);ss,z,_,alert=l.h.replay(hh,locks)
            fused=hh["r40c_fused_score"]-locks["r40c"]["thresholds"]["r40c_fused"]
            screen=ss-locks["r40f"]["threshold"];tree=z-locks["r40g"]["thresholds"]["r40g_tree"]
            parent=np.minimum(transform(fused),np.maximum(transform(screen),transform(tree)))
            assert np.array_equal(parent>=0,alert)
            cap.update(fused=fused,screen=screen,tree=tree,parent_margin=parent)
            aa,pp,hit=orig(q,locks)
            assert np.array_equal(aa,alert) and np.array_equal(pp,hh["r40c_fused_path"])
            return aa,pp,hit
        l.parent=observe
        try:d=base.infer("minimal",rows,refs,mm,pre,L,old,m,None,kind)
        finally:l.parent=orig
        oldarr=a.read_npz(base.C/("dev_"+num+"_"+kind+".npz"));a.verify_rows(oldarr,rows);base.check(d,oldarr,m)
        name=num+"-"+kind+".npz";rec=next(x for x in r["records"] if x["name"]==name)
        raw=rd(name);assert sha(raw)==rec["sha256"]
        with np.load(io.BytesIO(raw),allow_pickle=False) as q:
            X=n.ref.norm(np.stack([v[0] for v in rows]));Y,_,_,_=n.ref.register(X,n.ref.norm(refs));V=Y.astype(np.float32)
            assert sha(V.tobytes())==rec["array_sha256"] and np.array_equal(q["V"],V)
            assert np.array_equal(q["minimal_alert"],d["alert"]) and np.array_equal(q["minimal_hit"],d["hit"])
            for key in ("y","truth","active","amp","weak"):assert np.array_equal(q[key],d[key])
        pa=d["parent_alert"];idx=d["idx"];joint=d["joint_score"]-old["parent_joint_min_threshold"]
        prune=np.zeros(len(pa),float);prune[idx]=d["score"]-m["cuts"]["preserve"]
        supported=np.zeros(len(pa),bool);supported[idx]=d["use"]
        # A missing prune stage is explicitly ineligible; the original fail-open condition applies.
        parent=np.minimum(cap["parent_margin"],np.where(supported,transform(prune),1.))
        score=np.where(pa,parent,transform(joint))
        path=np.where(pa[:,None],d["parent_path"],d["p"])
        assert np.isfinite(score).all() and np.array_equal(score>=0,d["alert"])
        assert np.array_equal((score>=0)&l.correct({"truth":d["truth"],"active":d["active"],"y":d["y"]},path),d["hit"])
        np.savez_compressed(O/(num+"-"+kind+"-scores.npz"),score=score,path=path,alert=d["alert"],parent_route=pa,prune_supported=supported,prune_evaluated=pa,prune_margin=prune,joint_margin=joint,source_tile_sha=d["source_tile_sha"],**cap)
        v={"kind":kind,"rows":len(score),"finite_scores":int(np.isfinite(score).sum()),"original_alerts":int(d["alert"].sum()),"original_hits":int(d["hit"].sum()),"parent_rows":int(pa.sum()),"prune_supported":int(supported.sum()),"decision_changes":0,"registered_sha256":rec["array_sha256"],"score_file_sha256":sha((O/(num+"-"+kind+"-scores.npz")).read_bytes())}
        records.append(v);print("R42A3B_KIND",json.dumps(v),flush=True)

    result={"group":i+1,"fit":g["fit"],"rows":sum(x["rows"] for x in records),"records":records}
    assert result["rows"]==288
    save("group-"+num+".json",result);print("R42A3B_GROUP",json.dumps(result),flush=True)

B=pathlib.Path("research/seti/r42a3b_full_score_protocol.json");Q=json.loads(B.read_text())
def prepare():
    assert sha(P.read_bytes())==Q["pilot"]["protocol_sha256"]
    assert sha(pathlib.Path("research/seti/r42a3_score_contract.py").read_bytes())==Q["pilot"]["code_sha256"]
    src.prepare()
    q=Q["pilot"];u="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/"+str(q["artifact"])
    h={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"}
    def get(u):
        r=requests.get(u,headers=h,timeout=120);r.raise_for_status();return r
    a=get(u).json();assert a["workflow_run"]["id"]==q["run"] and a["digest"]=="sha256:"+q["digest"]
    b=get(u+"/zip").content;assert len(b)<32*1024**2 and sha(b)==q["digest"]
    z=zipfile.ZipFile(io.BytesIO(b));files=[x for x in z.infolist() if not x.is_dir()]
    names={pathlib.PurePosixPath(x.filename).name:x.filename for x in files};assert len(names)==len(files) and sum(x.file_size for x in files)<64*1024**2
    raw=z.read(names["result.json"]);res=json.loads(raw)
    assert res["group"]==1 and res["fit"] and res["rows"]==288 and res["protocol_sha256"]==q["protocol_sha256"] and res["code_sha256"]==q["code_sha256"]
    for r in res["records"]:
        name="01-"+r["kind"]+"-scores.npz";raw=z.read(names[name]);assert sha(raw)==r["score_file_sha256"];(O/name).write_bytes(raw)
    save("group-01.json",{"group":1,"fit":True,"rows":288,"records":res["records"],"pilot_reused":True})
    save("preparation.json",{"groups":11,"pilot":q,"protocol_sha256":sha(B.read_bytes())})
    print("R42A3B_PREPARED",flush=True)

def summarize():
    r,rd=unpack();groups=[json.loads((O/f"group-{i:02d}.json").read_text()) for i in range(1,12)]
    assert sum(g["rows"] for g in groups)==3168 and sum(g["fit"] for g in groups)==6 and groups[0]["pilot_reused"]
    hold=[];fit=[]
    for g in groups:
        assert len(g["records"])==2
        for v in g["records"]:
            name=f'{g["group"]:02d}-'+v["kind"]+".npz";rec=next(x for x in r["records"] if x["name"]==name)
            assert rec["fit"]==g["fit"] and rec["array_sha256"]==v["registered_sha256"]
            raw=rd(name);assert sha(raw)==rec["sha256"]
            b=(O/(name[:-4]+"-scores.npz")).read_bytes();assert sha(b)==v["score_file_sha256"]
            with np.load(io.BytesIO(raw),allow_pickle=False) as q,np.load(io.BytesIO(b),allow_pickle=False) as z:
                s=z["score"];p=z["path"];y=q["y"].astype(bool);w=q["weak"].astype(bool)
                assert len(s)==len(y)==v["rows"] and np.isfinite(s).all() and np.array_equal(s>=0,z["alert"]) and np.array_equal(z["alert"],q["minimal_alert"])
                loc=((abs(p-q["truth"])<=2)&q["active"]).sum(1)>=.75*q["active"].sum(1)
                assert np.array_equal((s>=0)&loc&y,q["minimal_hit"])
                if "tile_sha" in q:assert np.array_equal(q["tile_sha"],z["source_tile_sha"])
                (fit if g["fit"] else hold).append((s,loc,y,w,np.full(len(y),v["kind"]=="shared",bool),np.full(len(y),g["group"])))
    def merge(a):return tuple(np.concatenate(x) for x in zip(*a))
    s,loc,y,w,shared,group=merge(hold)
    fs,fl,fy,fw,fsh,fg=merge(fit)
    assert len(s)==1440 and y.sum()==800 and (~y).sum()==640 and len(fs)==1728 and (~fy).sum()==768
    def metrics(t,x=s,l=loc,yy=y,ww=w,sh=shared):
        a=x>=t; hit=a&l&yy
        return {"localized":int(hit.sum()),"weak":int((hit&ww).sum()),"controls":int((a&~yy).sum()),"ordinary":int((hit&~sh).sum()),"shared":int((hit&sh).sum())}
    locked=metrics(0);assert all(locked[k]==v for k,v in Q["expected_audit"].items())
    thresholds=np.unique(s)[::-1]
    points=[{"localized":0,"weak":0,"controls":0,"ordinary":0,"shared":0,"fpr":0.,"tpr":0.,"fnr":1.}]
    for t in thresholds:
        d=metrics(t);points.append(d|{"fpr":d["controls"]/640,"tpr":d["localized"]/800,"fnr":1-d["localized"]/800})
    caps=[{"cap":c,"best":max((x for x in points if x["fpr"]<=c+1e-12),key=lambda x:(x["localized"],x["weak"],-x["controls"]))} for c in Q["caps"]]
    # Fit-only threshold, independent of audit score ordering. Not adopted.
    budget=int((~fy).sum()*35/640);cal=np.sort(fs[~fy])[::-1];t=float(np.nextafter(cal[budget],np.inf))
    calmet=metrics(t,fs,fl,fy,fw,fsh);assert calmet["controls"]<=budget
    auditmet=metrics(t)
    per_group=[]
    for i in sorted(set(group)):
        k=group==i;per_group.append({"group":int(i),"positives":int(y[k].sum()),"negatives":int((~y[k]).sum()),"locked":metrics(0,s[k],loc[k],y[k],w[k],shared[k]),"fit_only_diagnostic":metrics(t,s[k],loc[k],y[k],w[k],shared[k])})
    result={"status":"PASS_FULL_HISTORICAL_FIXED_ROUTING_SCORE_EXPORT__HOLD_FULL_NOISE_AND_REROUTING","groups":11,"rows":3168,"finite_scores":3168,"pilot_reused":True,"decision_changes":0,"locked":locked,"caps":caps,"fit_only_diagnostic":{"null_budget":budget,"fit":calmet,"audit":auditmet,"not_adopted":True},"per_group":per_group,"full_pipeline_roc":False,"production_adoption":False,"external_validation":False}
    save("result.json",result);save("curve.json",{"points":points,"definition":"fixed routing Minimal diagnostic, not full JANUS+NOISE ROC"})
    print("R42A3B_SUMMARY",json.dumps(result),flush=True)
