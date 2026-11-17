#!/usr/bin/env python3
"""Controlled family evidence, not a classifier of extraterrestrial origin."""
import pathlib,json,hashlib,sys,time,collections
import numpy as np
import r26c_reference_alignment as r
import r24_odontocrystal as od
P=json.loads(pathlib.Path("research/seti/r26d_family_protocol.json").read_text())
OUT=pathlib.Path("runtime/r26d");OUT.mkdir(parents=True,exist_ok=True)
def emit(k,d):print(k,json.dumps(d),flush=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def restore():
    r.BLC_OUT=OUT;r.BLC_P["artifacts"]["family_source"]=P["source"];r.BLC_P["artifacts"]["family_manifest"]=P["manifest"]
    z,n,p=r.archive("family_source");m,mn,mp=r.archive("family_manifest")
    manifest=r.get(m,mn,"manifest.json")
    assert z.read(n["manifest.json"])==m.read(mn["manifest.json"])
    files={f["name"]:f for f in P["files"]};assert len(files)==24
    data={}
    for q in manifest["selected"]:
        assert [f["role"] for f in q["files"]]==["ON1","OFF1","ON2","OFF2"]
        data[q["target"]]=[]
        for i in (0,2):
            a,b=[files[q["files"][j]["name"]] for j in (i,i+1)]
            assert len(a["tiles"])==len(b["tiles"])==4
            for ta,tb in zip(a["tiles"],b["tiles"]):
                assert ta["center_mhz"]==tb["center_mhz"]
                x,y,lo=r.common(ta,tb,r.tile(z,n,ta),r.tile(z,n,tb))
                t={**ta,"channel0_mhz":ta["channel0_mhz"]+lo*ta["df_hz"]/1e6}
                data[q["target"]].append((x,y,t,a["name"],b["name"]))
    z.close();m.close();p.unlink();mp.unlink()
    assert set(data)==set(P["calibration"]+P["evaluation"]) and all(len(x)==8 for x in data.values())
    return data
def path(rng):
    t=np.arange(16);d=np.rint(rng.uniform(-.65,.65)*t+4*np.sin(2*np.pi*t/16+rng.uniform(0,6.28))).astype(int);d-=d[0]
    assert np.max(abs(np.diff(d)))<=3 and np.max(abs(d))<=20
    return d
def inject(x,q,a):
    # Power-domain addition in robust row-scale units, prior to normalization.
    sc=1.4826*np.median(abs(x-np.median(x,axis=1)[:,None]),axis=1)
    assert np.all(sc>0) and np.all((q>=0)&(q<x.shape[1]))
    x[np.arange(16),q]+=a*sc
def sample(data,target):
    rng=np.random.default_rng(int.from_bytes(hashlib.sha256((str(P["seed"])+target).encode()).digest()[:8],"big"))
    X=[];R=[];labels=[];truth=[]
    for a,b,t,on,off in data:
        for _ in range(P["patches_per_on_tile"]):
            start=int(rng.integers(0,a.shape[1]-384+1));x=a[:,start:start+384].astype(float);y=b[:,start:start+384].astype(float)
            d=path(rng);spacing=int(rng.choice(P["spacing_bins"]));amp=float(rng.choice(P["sibling_amplitudes"]))
            for label in P["classes"]:
                u=x.copy();v=y.copy();inject(u,48+d,P["anchor_amplitude"])
                if label in ("family","mirror","on_only_comb","unrelated"):
                    for j in (1,2):
                        dd=path(rng) if label=="unrelated" else -d if label=="mirror" and j==1 else d
                        inject(u,48+j*spacing+dd,amp)
                        if label in ("family","mirror"):inject(v,48+j*spacing+dd,amp)
                if label=="repeat_only":inject(v,48+d,amp)
                X.append(u);R.append(v);labels.append(label);truth.append(48+d)
    return r.norm(np.stack(X)),r.norm(np.stack(R)),np.array(labels),np.stack(truth)
def matched(V,center,d):
    # Whole-path shift of at most two bins; mirrored morphology tested explicitly.
    ans=np.full(len(V),-np.inf);t=np.arange(16)[None,:];i=np.arange(len(V))[:,None]
    for sign in (-1,1):
        for delta in (-2,-1,0,1,2):
            q=center[:,None]+sign*d+delta
            ok=np.all((q>=0)&(q<V.shape[2]),axis=1)
            z=V[i,t,np.clip(q,0,V.shape[2]-1)].sum(1)/4
            ans=np.maximum(ans,np.where(ok,z,-np.inf))
    return ans
def unrestricted(V,cent):
    # Same start-frequency region, no shared morphology constraint.
    rows=np.arange(len(V))[:,None,None];t=np.arange(16)[None,:,None]
    cols=cent[:,None]+np.arange(-24,25)[None,:]
    assert np.all((cols>=0)&(cols<V.shape[2]))
    z=V[rows,t,cols[:,None,:]]
    return od.sliced(z,0,None)[0]
def scores(X,R):
    s,p,_=od.sliced(X[:,:,16:80],0,None);p=p+16;cent=p[:,0];d=p-cent[:,None]
    v={k:np.full(len(X),-np.inf) for k in P["modes"]}
    for spacing in P["spacing_bins"]:
        a=cent+spacing;b=cent+2*spacing
        freq=np.minimum(unrestricted(X,a),unrestricted(X,b))
        ma=matched(X,a,d);mb=matched(X,b,d);m=np.minimum(ma,mb)
        off=np.maximum(matched(R,a,d),matched(R,b,d))
        for k,z in (("frequency_pair",freq),("morphology_pair",m),("morphology_off",np.minimum(m,off))):v[k]=np.maximum(v[k],z)
    assert all(np.isfinite(x).all() for x in v.values())
    return v,p,s
def cut(v):
    a=np.sort(v)[::-1];return float(np.nextafter(a[int(len(a)*P["control_fraction"])],np.inf))
def test():
    rng=np.random.default_rng(123);x=rng.normal(0,.1,(4,16,384));v=x.copy();d=np.rint(4*np.sin(np.arange(16)*np.pi/8)).astype(int)
    for i in range(4):inject(x[i],48+d,80)
    for j in (1,2):
        inject(x[0],48+80*j+d,80);inject(v[0],48+80*j+d,80)
        inject(x[1],48+80*j-d,80);inject(v[1],48+80*j-d,80)
        inject(x[2],48+80*j+d,80)
    a,p,s=scores(r.norm(x),r.norm(v))
    assert np.all(abs(p-(48+d))<=1)
    assert a["morphology_off"][0]>10 and a["morphology_off"][1]>10
    assert a["morphology_off"][2]<10 and a["morphology_pair"][2]>10 and a["morphology_pair"][3]<10
    return "PASS"
def catalog(data):
    # Bounded fixed linear drift bank, top two distinct lines per measured scan/tile.
    hits=[]
    for target,items in data.items():
        for a,b,t,on,off in items:
            for X,name,role in ((a,on,"ON"),(b,off,"OFF")):
                y=r.norm(X.astype(float));w=y.shape[1];best=[]
                for shift in range(-32,33):
                    d=np.rint(np.arange(16)*shift/15).astype(int);lo=max(0,-int(d.min()));hi=min(w,w-int(d.max()))
                    z=sum(y[i,lo+k:hi+k] for i,k in enumerate(d))/4
                    for jj in np.argpartition(z,-2)[-2:]:
                        best.append((float(z[jj]),lo+int(jj),shift))
                keep=[]
                for sc,j,sh in sorted(best,reverse=True):
                    if sc<10:break
                    if any(abs(j-k)<8 for _,k,_ in keep):continue
                    keep.append((sc,j,sh))
                    if len(keep)==2:break
                for sc,j,sh in keep:
                    hits.append({"target":target,"scan":name,"role":role,"center":t["center_mhz"],"frequency_hz":t["channel0_mhz"]*1e6+j*t["df_hz"],"drift_hz_s":sh*t["df_hz"]/(15*t["ts"]),"score":sc})
    pairs=[]
    for i,a in enumerate(hits):
        for b in hits[i+1:]:
            if a["scan"]==b["scan"] or a["center"]!=b["center"]:continue
            if abs(a["frequency_hz"]-b["frequency_hz"])<=8.4 and abs(a["drift_hz_s"]-b["drift_hz_s"])<=.011:
                pairs.append({"a":a,"b":b})
    return {"hits":hits,"recurrence_pairs":pairs,"scope":"Top two per tile above fixed score, same-frequency candidates at scan starts only; no whole-band coverage, labels, or origin inference."}
def main():
    assert test()=="PASS";data=restore();cal=[];ev=[]
    for target in P["calibration"]+P["evaluation"]:
        X,R,l,t=sample(data[target],target);start=time.perf_counter();s,p,a=scores(X,R)
        z={"target":target,"label":l,"scores":s,"anchor":np.mean(abs(p-t)<=2,axis=1)>=.75,"seconds":time.perf_counter()-start}
        (cal if target in P["calibration"] else ev).append(z)
    th={}
    for k in P["modes"]:
        v=np.concatenate([a["scores"][k][np.isin(a["label"],P["calibration_negative"])] for a in cal])
        th[k]=cut(v)
    save("lock.json",{"thresholds":th,"calibration_examples":sum(len(a["label"]) for a in cal)})
    out={}
    for k in P["modes"]:
        out[k]={}
        for label in P["classes"]:
            n=anchor=alert=conditional=0
            for a in ev:
                m=a["label"]==label;q=a["scores"][k]>=th[k];h=a["anchor"]
                n+=int(m.sum());anchor+=int((m&h).sum());alert+=int((m&q).sum());conditional+=int((m&q&h).sum())
            out[k][label]={"examples":n,"anchor_localized":anchor,"family_alerts":alert,"alerts_on_localized_anchor":conditional}
    cat=catalog(data);save("catalog.json",cat)
    d={"status":"DEVELOPMENT_REPLAY","selftest":"PASS","thresholds":th,"modes":out,"eval_examples":sum(len(a["label"]) for a in ev),"calibration_examples":sum(len(a["label"]) for a in cal),"evaluation_seconds":sum(a["seconds"] for a in ev),"real_catalog":{"hits":len(cat["hits"]),"pairs":len(cat["recurrence_pairs"]),"cross_target_pairs":sum(x["a"]["target"]!=x["b"]["target"] for x in cat["recurrence_pairs"]),"on_off_pairs":sum(x["a"]["role"]!=x["b"]["role"] for x in cat["recurrence_pairs"])},"scope":P["scope"],"production_adoption":False}
    save("result.json",d);emit("R26D_SUMMARY",d)
if __name__=="__main__":
    if len(sys.argv)>1 and sys.argv[1]=="selftest":emit("R26D_SELFTEST",test())
    else:main()
