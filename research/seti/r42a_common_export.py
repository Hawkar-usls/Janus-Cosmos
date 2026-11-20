"""Exact R41 common-input export; R41S1F historical NOISE masks are archived, not retrained."""
import sys,pathlib,json,time,hashlib
import numpy as np
import r41s1f_unclipped_joint as f
base=f.base;a=f.a;n=f.n
O=pathlib.Path("runtime/r42a");O.mkdir(parents=True,exist_ok=True)
IN=O/"inputs";IN.mkdir(exist_ok=True)
P=pathlib.Path("research/seti/r42a_competitor_protocol.json");L=json.loads(P.read_text())
def sha(b):return hashlib.sha256(b).hexdigest()
def save(name,d):(O/name).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def read(name):return json.loads((O/name).read_text())
def prepare():
    base.prepare()
    d,rd=f.download(L["source"])
    assert d["baseline"]=={"controls":53,"localized":318,"weak":69,"ordinary":205,"shared":113}
    assert d["models"]["raw_noise"]["hypothetical"]["controls"]==51 and d["models"]["raw_noise"]["hypothetical"]["localized"]==318
    save("noise-model.json",d["models"]["raw_noise"])
    for i in range(1,12):
        name=f"group-{i:02d}.json";(O/name).write_bytes(rd(name))
    save("preparation.json",{"protocol_hash":sha(P.read_bytes()),"source":L["source"],"groups":11,"production_adoption":False})
def group(i):
    assert read("preparation.json")["protocol_hash"]==sha(P.read_bytes())
    g=base.plan()["groups"][i];num=f"{i+1:02d}";r=read("group-"+num+".json")
    assert g["fit"]==r["fit_viewed"] and g["group"]["target"]==r["target"]
    seed={"n6":a.v.C["seed"],"n7":base.w.C["seed"],"n9":base.old9.C["seed"]}[g["source"]]
    n.set_config(seed,a.C["windows_mhz"],a.C["half_hz"],a.C["max_mib_per_scan"]);base.b.taylor.init_tree()
    m,LL,old0=base.locks("minimal",base.plan())
    t=time.perf_counter();rr,R,ch,hs,meta,cm,pre=n.gather(g["group"],"r42a_"+num)
    acq=time.perf_counter()-t;assert sorted(hs)==g["raw_tile_hashes"]
    noise=read("noise-model.json")
    keep={(p["group"],p["kind"],p["row"]):p["would_keep"] for lane in ("parent","rescue") for p in noise[lane]["predictions"]}
    records=[]
    for kind,rows,refs,mm in (("eval",rr,R,meta),("shared",[x[0] for x in ch],np.stack([x[1] for x in ch]),cm)):
        t=time.perf_counter();d=base.infer("minimal",rows,refs,mm,pre,LL,old0,m,None,kind);janus_wall=time.perf_counter()-t
        oldarr=a.read_npz(base.C/("dev_"+num+"_"+kind+".npz"));a.verify_rows(oldarr,rows);base.check(d,oldarr,m)
        archived=[v for v in r["records"] if v["kind"]==kind]
        assert [v["row"] for v in archived]==np.flatnonzero(d["alert"]).tolist()
        for v in archived:
            j=v["row"];assert v["label_evaluation_only"]==int(d["y"][j]) and v["localized_evaluation_only"]==bool(d["hit"][j]) and v["weak"]==bool(d["weak"][j])
        t=time.perf_counter()
        X=n.ref.norm(np.stack([v[0] for v in rows]));Y,_,_,_=n.ref.register(X,n.ref.norm(refs));V=Y.astype(np.float32)
        assert V.shape==(len(rows),16,128) and np.isfinite(V).all()
        prep_wall=time.perf_counter()-t
        mask=np.array([keep.get((i+1,kind,j),True) for j in range(len(rows))],bool)
        if g["fit"]:assert mask.all()
        arr={"V":V,"y":d["y"],"truth":d["truth"],"active":d["active"],"amp":d["amp"],"weak":d["weak"],"minimal_alert":d["alert"],"minimal_hit":d["hit"],"noise_alert":d["alert"]&mask,"noise_hit":d["hit"]&mask,"tile_sha":d["source_tile_sha"]}
        name=f"{num}-{kind}.npz";path=IN/name;np.savez_compressed(path,**arr)
        records.append({"name":name,"sha256":sha(path.read_bytes()),"array_sha256":sha(V.tobytes()),"group":i+1,"target":r["target"],"fit":bool(g["fit"]),"kind":kind,"rows":len(rows),"janus_wall_s":janus_wall,"registration_wall_s":prep_wall})
    save(f"export-{num}.json",{"group":i+1,"acquisition_wall_s":acq,"records":records});print("R42A_GROUP",num,flush=True)
def summary():
    groups=[read(f"export-{i:02d}.json") for i in range(1,12)]
    records=[r for g in groups for r in g["records"]];assert len(records)==22 and sum(r["rows"] for r in records)==3168
    def count(fit):
        out={}
        for method in ("minimal","noise"):
            d={k:0 for k in ("ordinary","shared","weak","controls","positives","negatives")}
            for r in records:
                if r["fit"]!=fit:continue
                with np.load(IN/r["name"],allow_pickle=False) as q:
                    h=q[method+"_hit"];d["ordinary" if r["kind"]=="eval" else "shared"]+=int(h.sum());d["weak"]+=int((h&q["weak"]).sum())
                    d["controls"]+=int((q[method+"_alert"]&(q["y"]==0)).sum());d["positives"]+=int((q["y"]==1).sum());d["negatives"]+=int((q["y"]==0).sum())
            d["localized"]=d["ordinary"]+d["shared"];out[method]=d
        return out
    fit=count(True);audit=count(False)
    assert audit["minimal"]["localized"]==318 and audit["minimal"]["weak"]==69 and audit["minimal"]["controls"]==53
    assert audit["noise"]["localized"]==318 and audit["noise"]["weak"]==69 and audit["noise"]["controls"]==51
    d={"status":"PASS_HISTORICAL_COMMON_INPUT","records":records,"fit":fit["minimal"],"audit":audit,"acquisition_wall_s":sum(g["acquisition_wall_s"] for g in groups),"historical_only":True,"production_adoption":False}
    save("receipt.json",d);print("R42A_BASE",json.dumps({k:v for k,v in d.items() if k!="records"}),flush=True)
if __name__=="__main__":
    cmd=sys.argv[1]
    if cmd=="prepare":prepare()
    elif cmd=="group":group(int(sys.argv[2]))
    elif cmd=="summary":summary()
