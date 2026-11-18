#!/usr/bin/env python3
"""R40H: numeric-only frozen inference and prospective archive transfer."""
import sys
MODE=sys.argv[1]
import json,pathlib,time,hashlib
import numpy as np
import r40g_recovery as g
f=g.f;c=g.c;t=g.t
P=pathlib.Path("research/seti/r40h_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r40h")/MODE;O.mkdir(parents=True,exist_ok=True)
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def readlocks():
    L={}
    for k,a in S["locks"].items():
        z,n=c.artifact(a,32);raw=z.read(n[a["file"]]);assert c.sha(raw)==a["lock_sha256"]
        L[k]=json.loads(raw);(O/(k+".json")).write_bytes(raw)
    assert L["r40f"]["parent_model_sha256"]==S["locks"]["r40c"]["lock_sha256"] and L["r40g"]["parent_screen_sha256"]==S["locks"]["r40f"]["lock_sha256"]
    return L
def tree(X,L):
    X=np.asarray(X,dtype=np.float32);acc=np.zeros(len(X),float)
    for e in L["trees"]:
        left=np.asarray(e["children_left"],int);right=np.asarray(e["children_right"],int);fe=np.asarray(e["feature"],int);th=np.asarray(e["threshold"],float);v=np.asarray(e["value"],float)
        ix=np.zeros(len(X),int)
        for step in range(6):
            m=left[ix]>=0
            if not m.any():break
            i=ix[m];ix[m]=np.where(X[np.flatnonzero(m),fe[i]]<=th[i],left[i],right[i])
        assert (left[ix]<0).all() and v.shape[1:]==(1,2)
        leaf=v[ix,0,:];assert (leaf.sum(1)>0).all()
        acc+=leaf[:,1]/leaf.sum(1)
    return acc/len(L["trees"])
def replay(q,L):
    ss=g.screen(q,L["r40f"]);old=ss>=L["r40f"]["threshold"]
    z=tree(g.features(q),L["r40g"])
    alert=q["r40c_fused_alert"].astype(bool)
    return ss,z,alert&old,alert&(old|(z>=L["r40g"]["thresholds"]["r40g_tree"]))
def freeze():
    L=readlocks();z,n=c.artifact(S["locks"]["r40g"],32)
    checks=[];seen=set()
    for target,a in f.S["artifacts"].items():
        zz,nn=c.artifact(a,128);res=json.loads(zz.read(nn["result.json"]));assert res["model_sha256"]==S["locks"]["r40c"]["lock_sha256"]
        seen.update(res["source_tile_hashes"])
        for kind in ("eval","shared"):
            q=f.load(zz,nn,kind+".npz");prior=f.load(z,n,target+"_"+kind+".npz");ss,v,old,new=replay(q,L)
            assert np.array_equal(prior["original_path"],q["r40c_fused_path"]) and np.array_equal(prior["original_alert"],q["r40c_fused_alert"])
            assert np.allclose(ss,prior["screen_score"],atol=1e-12,rtol=0) and np.allclose(v,prior["r40g_tree"],atol=1e-12,rtol=0)
            assert np.array_equal(old,prior["r40f_screen_priority_alert"]) and np.array_equal(new,prior["r40g_tree_priority_alert"])
            checks.append({"target":target,"kind":kind,"rows":len(v),"max_probability_delta":float(np.max(abs(v-prior["r40g_tree"])))})
    assert sum(x["rows"] for x in checks)==1440
    hz,hn,rec,old=t.parent.restore();ms=json.loads(hz.read(hn["manifest.json"]))
    for r in rec["records"]:
        q=t.parent.load(hz,hn,r);seen.update(str(v) for v in q["tile_sha"])
    save("freeze.json",{"status":"PASS_NUMERIC_SERIALIZATION","protocol_sha256":c.sha(P.read_bytes()),"locks":S["locks"],"rows":1440,"exact_decisions":True,"checks":checks,"excluded_raw_hashes":sorted(seen)})
    emit("R40H_FREEZE",{"rows":1440,"max_probability_delta":max(x["max_probability_delta"] for x in checks),"exact_decisions":True})
    # Select by metadata, after model is locked, before accessing any new power.
    import r35b_grand_benchmark as b
    rz,rn=c.artifact(S["r37_manifest"],32);rm=json.loads(rz.read(rn["manifest.json"]))
    prev=ms["selected"]+rm["selected"]["Parkes"]+rm["selected"].get("GBT",[])
    audit=json.loads(pathlib.Path("research/seti/r35g_common_native_protocol.json").read_text())["freshness_audit"]
    names={q["target"] for q in prev}|set(audit["known_r32_targets"])
    days={int(q["day"]) for q in prev}|set(audit["known_r32_days"])
    b.S["seed"]=S["seed"];b.sh.S=b.S
    b.sh.C={**b.sh.C,"fit_targets":0,"null_targets":0,"eval_targets":S["targets"],"max_queries":S["max_queries"],"exclude_days":sorted(set(b.sh.C["exclude_days"])|days)}
    b.sh.hold.S={**b.sh.hold.S,"exclude":sorted(set(b.sh.hold.S["exclude"])|names)}
    b.sh.OUT=O
    b.sh.freeze()
    m=json.loads((O/"manifest.json").read_text())
    assert len(m["selected"])==S["targets"] and not {q["target"] for q in m["selected"]}&names and not {int(q["day"]) for q in m["selected"]}&days
    # Pairwise near-identical coordinates are not independent target names.
    def coords(q):
        hd=q["files"][0]["header"];a=str(hd["src_raj"]);d=str(hd["src_dej"]);return (a,d)
    known={coords(q) for q in prev if q.get("files") and "header" in q["files"][0]}
    new=[coords(q) for q in m["selected"]];assert len(set(new))==len(new) and not set(new)&known
    save("selection_receipt.json",{"status":"PASS_METADATA_SELECTION","manifest_sha256":c.sha((O/"manifest.json").read_bytes()),"excluded_names":sorted(names),"excluded_days":sorted(days),"targets":[{"target":q["target"],"day":q["day"]} for q in m["selected"]],"coordinate_exact_match_excluded":True,"alias_completeness_proven":False})
    emit("R40H_SELECTION",{"targets":[[q["target"],q["day"]] for q in m["selected"]],"manifest_sha256":c.sha((O/"manifest.json").read_bytes())})
def evaluate():
    import os
    A={"id":int(os.environ["FREEZE_ID"]),"sha256":os.environ["FREEZE_DIGEST"].removeprefix("sha256:")}
    zz,nn=c.artifact(A,32);receipt=json.loads(zz.read(nn["freeze.json"]));sel=json.loads(zz.read(nn["selection_receipt.json"]));raw=zz.read(nn["manifest.json"]);assert c.sha(raw)==sel["manifest_sha256"] and receipt["protocol_sha256"]==c.sha(P.read_bytes())
    L={}
    for k,a in S["locks"].items():
        bb=zz.read(nn[k+".json"]);assert c.sha(bb)==a["lock_sha256"];L[k]=json.loads(bb)
    m=json.loads(raw);assert len(m["selected"])==S["targets"]
    import r35b_grand_benchmark as b
    b.S["seed"]=S["seed"];b.ref.S=b.S;b.ref.OUT=O;b.sh.eng.OUT=O
    b.sh.eng.SPEC={**b.sh.eng.SPEC,"window_centers_mhz":S["windows_mhz"],"window_half_hz":S["half_hz"],"limits":{**b.sh.eng.SPEC["limits"],"spectrum_budget_per_file_mib":S["max_mib_per_scan"]}}
    b.taylor.init_tree();seen=set(receipt["excluded_raw_hashes"]);records=[];rows=[]
    (O/"manifest.json").write_bytes(raw)
    def encode(rows,R):
        d=b.panel(rows,R);arms=b.base_arms(d,L["r40c"]["original"]["component_fit"])
        v={"V":d["Y"].astype(np.float32),"y":b.labels(d),"truth":np.stack([r[7] for r in rows]),"active":np.stack([r[8] for r in rows]),"amp":np.asarray([r[4] for r in rows]),"tile_sha":np.asarray([r[6] for r in rows])}
        for k in t.S["paths"]:v[k+"_score"],v[k+"_paths"]=arms[k][:2]
        return v
    for group in m["selected"]:
        start=time.perf_counter();name=group["target"];assert name.isalnum()
        rr,R,ch,hs=b.ref.gather(group);assert not seen&hs;seen.update(hs)
        target=[]
        for kind,q in (("eval",encode(rr,R)),("shared",encode([v[0] for v in ch],np.stack([v[1] for v in ch])))):
            h=f.historical(q,L["r40c"]);ss,v,old,new=replay(h,L)
            p=h["r40c_fused_path"];y=h["y"];weak=np.isin(h["amp"],t.S["weak_amplitudes"]);hit=h["r40c_fused_localized"].astype(bool);base=h["r28_frozen_localized"].astype(bool)
            a={"r28":(h["r28_frozen_alert"],base),"r40c":(h["r40c_fused_alert"],hit),"r40f":(old,old&hit),"r40g":(new,new&hit)}
            met={}
            for k,(alert,hits) in a.items():
                met[k]={"localized":int(hits.sum()),"weak":int((hits&weak).sum()),"control_alerts":int((alert&(y==0)).sum()),"lost_vs_r40c":int((hit&~hits).sum()),"lost_vs_r28":int((base&~hits).sum()),"weak_lost_vs_r40c":int((hit&~hits&weak).sum()),"positives":int((y==1).sum()),"controls":int((y==0).sum())}
            fn=kind+"_"+name+".npz";np.savez_compressed(O/fn,V=h["V"],y=y,truth=h["truth"],active=h["active"],amp=h["amp"],r28_hit=base,r40g_hit=new&hit,r40g_alert=new,r40g_path=p,r40g_probability=v,r40f_score=ss)
            records.append({"name":fn,"sha256":c.sha((O/fn).read_bytes()),"array_sha256":c.sha(h["V"].tobytes())})
            target.append({"kind":kind,"metrics":met})
        row={"target":name,"day":group["day"],"metrics":target,"wall_seconds":time.perf_counter()-start,"raw_tile_hashes":sorted(hs)};rows.append(row);save("partial.json",rows);emit("R40H_TARGET",row)
    def tot(kind):
        rr=[v["metrics"] for row in rows for v in row["metrics"] if v["kind"]==kind]
        return {k:{j:sum(x[k][j] for x in rr) for j in rr[0][k]} for k in rr[0]}
    e=tot("eval");sh=tot("shared");x=e["r40g"];b=e["r28"]
    assert b["positives"]==b["controls"]==640 and sh["r28"]["positives"]==160
    gate=x["weak"]>b["weak"] and x["control_alerts"]<b["control_alerts"] and x["localized"]>=b["localized"] and sh["r40g"]["localized"]>=sh["r28"]["localized"] and x["lost_vs_r40c"]==0 and sh["r40g"]["lost_vs_r40c"]==0
    d={"status":"PROSPECTIVE_ARCHIVE_INJECTION_ONLY","eval":e,"shared":sh,"gate":gate,"records":records,"rows":rows,"freeze":A,"manifest_sha256":sel["manifest_sha256"],"production_adoption":False,"discovery_relaunch_allowed":False}
    save("receipt.json",d);emit("R40H_SUMMARY",{k:v for k,v in d.items() if k not in ("records","rows")})
if __name__=="__main__":
    if MODE=="freeze":freeze()
    elif MODE=="evaluate":evaluate()
    else:raise ValueError(MODE)
