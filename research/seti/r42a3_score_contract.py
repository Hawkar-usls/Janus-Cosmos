"""R42A3A diagnostic score contract pilot, frozen R41 Minimal routing."""
import pathlib,json,hashlib,io,zipfile,os
import numpy as np
import requests
import r42a_common_export as src
base=src.base;a=src.a;n=src.n;l=base.l
P=pathlib.Path("research/seti/r42a3_score_contract_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r42a3");O.mkdir(parents=True,exist_ok=True)
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
def run():
    src.prepare();r,rd=unpack()
    gs=base.plan()["groups"];i=next(i for i,g in enumerate(gs) if g["fit"]);g=gs[i];num=f"{i+1:02d}"
    assert len([x for x in r["records"] if x["group"]==i+1 and x["fit"]])==2
    seed={"n6":a.v.C["seed"],"n7":base.w.C["seed"],"n9":base.old9.C["seed"]}[g["source"]]
    n.set_config(seed,a.C["windows_mhz"],a.C["half_hz"],a.C["max_mib_per_scan"]);base.b.taylor.init_tree()
    m,L,old=base.locks("minimal",base.plan())
    rr,R,ch,hs,meta,cm,pre=n.gather(g["group"],"r42a3_pilot_"+num)
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
        records.append(v);print("R42A3A_KIND",json.dumps(v),flush=True)
    result={"status":"PASS_COMPONENT_SCORE_EXPORT_FIXED_ROUTING__HOLD_FULL_NOISE_SCORE","group":i+1,"fit":True,"rows":sum(x["rows"] for x in records),"records":records,"protocol_sha256":sha(P.read_bytes()),"code_sha256":sha(pathlib.Path(__file__).read_bytes()),"production_adoption":False,"external_validation":False,"full_pipeline_roc":False}
    assert result["rows"]==288
    save("result.json",result);print("R42A3A_SUMMARY",json.dumps(result),flush=True)
if __name__=="__main__":run()
