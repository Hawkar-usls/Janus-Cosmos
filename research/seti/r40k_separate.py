#!/usr/bin/env python3
"""R40K independently normalized ON/OFF stationary rescue, frozen control threshold."""
import sys
sys.argv[1:]=["separate"]
import json,pathlib,time
import numpy as np
import r40h_frozen as h
import r35b_grand_benchmark as b
c=h.c;t=h.t;f=h.f
P=pathlib.Path("research/seti/r40k_separate_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r40k");O.mkdir(parents=True,exist_ok=True)
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def candidate(X,R):
    def contrast(V):
        idx=np.arange(5,123);s=V[:,:,idx];side=V[:,:,idx[:,None]+np.array([-5,-4,-3,3,4,5])]
        return np.clip(s-np.median(side,axis=3),-4,4).mean(1)
    a=contrast(X);d=contrast(R)
    z={"on_only":a,"joint_min":np.minimum(a,d),"joint_mean":np.where((a>0)&(d>0),(a+d)/2,-4.)}
    out={}
    for k,v in z.items():
        j=v.argmax(1);out[k]=(v[np.arange(len(v)),j],np.repeat((j+5)[:,None],16,axis=1))
    return out
def regen(rows,R,q):
    X=b.ref.norm(np.stack([v[0] for v in rows]));Y,_,_,_=b.ref.register(X,b.ref.norm(R))
    assert np.array_equal(Y.astype(np.float32),q["V"]), "registered common array mismatch"
    for k,v in {"y":np.array([r[1] for r in rows]),"truth":np.stack([r[7] for r in rows]),"active":np.stack([r[8] for r in rows]),"amp":np.array([r[4] for r in rows])}.items():assert np.array_equal(v,q[k]),k
    if "tile_sha" in q:assert np.array_equal(np.array([r[6] for r in rows]),q["tile_sha"])
    return candidate(X,b.ref.norm(R))
def gather(group,tag):
    w=O/"raw"/tag;w.mkdir(parents=True,exist_ok=True);b.ref.OUT=w;b.sh.eng.OUT=w
    return b.ref.gather(group)
def main():
    st=time.perf_counter();L=h.readlocks();z,n,rec,_=t.parent.restore();m=json.loads(z.read(n["manifest.json"]))
    null={k:[] for k in S["arms"]};count=0;calcount=0;calrows=[]
    for group in m["selected"]:
        name=group["target"]
        if group["split"]!="null_calibration":continue
        rr=next(x for x in rec["records"] if x["name"]=="cal_"+name+".npz");q=t.parent.load(z,n,rr)
        rows,R,ch,hs=gather(group,"cal_"+name);ca=regen(rows,R,q)
        hh=f.historical(q,L["r40c"]);_,_,_,a=h.replay(hh,L);mask=q["y"]==0
        count+=int((a&mask).sum());calcount+=int(mask.sum())
        for k,(s,p) in ca.items():null[k].extend(s[~a&mask].tolist())
        calrows.append({"target":name,"raw_tile_hashes":sorted(hs),"registered_exact":True,"controls":int(mask.sum())})
        emit("R40K_CAL_GROUP",calrows[-1])
    assert count==10 and calcount==256
    th={k:float(np.nextafter(max(v),np.inf)) for k,v in null.items()}
    lock={"thresholds":th,"calibration_parent_controls":count,"calibration_controls":calcount,"additional_calibration_controls":0,"calibration":calrows,"protocol_sha256":c.sha(P.read_bytes())}
    save("lock.json",lock);lh=c.sha((O/"lock.json").read_bytes());emit("R40K_LOCK",lock)
    zz,nn=c.artifact(S["sources"]["r40h"],192);rec=json.loads(zz.read(nn["receipt.json"]));m=json.loads(zz.read(nn["manifest.json"]))
    assert len(m["selected"])==5 and rec["manifest_sha256"]==c.sha(zz.read(nn["manifest.json"]))
    b.S["seed"]=h.S["seed"];b.ref.S=b.S
    b.sh.eng.SPEC={**b.sh.eng.SPEC,"window_centers_mhz":h.S["windows_mhz"],"window_half_hz":h.S["half_hz"],"limits":{**b.sh.eng.SPEC["limits"],"spectrum_budget_per_file_mib":h.S["max_mib_per_scan"]}}
    results=[]
    for group in m["selected"]:
        name=group["target"];rows,R,ch,hs=gather(group,name)
        expected=next(v for v in rec["rows"] if v["target"]==name)
        assert sorted(hs)==expected["raw_tile_hashes"]
        for kind,rr,ref in (("eval",rows,R),("shared",[v[0] for v in ch],np.stack([v[1] for v in ch]))):
            key=kind+"_"+name+".npz";ar=next(v for v in rec["records"] if v["name"]==key)
            assert c.sha(zz.read(nn[key]))==ar["sha256"];q=f.load(zz,nn,key)
            assert c.sha(q["V"].tobytes())==ar["array_sha256"]
            ca=regen(rr,ref,q);a=q["r40g_alert"].astype(bool);hit=q["r40g_hit"].astype(bool);y=q["y"];weak=np.isin(q["amp"],t.S["weak_amplitudes"]);act=q["active"].astype(bool)
            out={};arr={"parent_alert":a,"parent_hit":hit,"parent_path":q["r40g_path"]}
            def met(aa,hh):
                return {"localized":int(hh.sum()),"weak":int((hh&weak).sum()),"control_alerts":int((aa&(y==0)).sum()),"parent_lost":int((hit&~hh).sum()),"added_localized":int((hh&~hit).sum()),"positives":int((y==1).sum()),"controls":int((y==0).sum())}
            out["r40g"]=met(a,hit)
            for k,(s,p) in ca.items():
                rescue=(~a)&(s>=th[k]);correct=(((abs(p-q["truth"])<=2)&act).sum(1)>=.75*act.sum(1))&(y==1)
                aa=a|rescue;hh=hit|(rescue&correct)
                out[k]=met(aa,hh);assert out[k]["parent_lost"]==0
                arr.update({k+"_score":s,k+"_path":p,k+"_rescue":rescue,k+"_localized":hh})
            np.savez_compressed(O/(key),**arr)
            row={"target":name,"kind":kind,"registered_exact":True,"raw_tile_hashes":sorted(hs),"metrics":out};results.append(row);save("partial.json",results);emit("R40K_GROUP",row)
        assert c.sha((O/"lock.json").read_bytes())==lh
    def total(kind):
        rr=[v["metrics"] for v in results if v["kind"]==kind];assert len(rr)==5
        return {k:{j:sum(v[k][j] for v in rr) for j in rr[0][k]} for k in rr[0]}
    e=total("eval");sh=total("shared")
    assert [e["r40g"][k] for k in ("localized","weak","control_alerts")]==[213,41,31] and sh["r40g"]["weak"]==18
    gates={k:e[k]["localized"]>=214 and e[k]["weak"]>=41 and sh[k]["weak"]>=21 and e[k]["control_alerts"]<=35 and e[k]["parent_lost"]==sh[k]["parent_lost"]==0 for k in S["arms"]}
    d={"status":"DEVELOPMENT_REPLAY_ONLY","eval":e,"shared":sh,"gates":gates,"rows":results,"lock_sha256":lh,"seconds":time.perf_counter()-st,"production_adoption":False}
    save("result.json",d);emit("R40K_SUMMARY",{k:v for k,v in d.items() if k!="rows"})
if __name__=="__main__":main()
