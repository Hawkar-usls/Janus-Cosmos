#!/usr/bin/env python3
import json,pathlib,time,numpy as np
import r27a_benchmark as a
import r26c_reference_alignment as r
import r24_odontocrystal as od
P=json.loads(pathlib.Path("research/seti/r27a2_equal_protocol.json").read_text());OUT=pathlib.Path("runtime/r27a2");OUT.mkdir(parents=True,exist_ok=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def restore():
    r.BLC_OUT=OUT;r.BLC_P["artifacts"]["equal_source"]=P["source"];r.BLC_P["artifacts"]["equal_manifest"]=P["manifest"]
    z,n,p=r.archive("equal_source");m,mn,mp=r.archive("equal_manifest");manifest=r.get(m,mn,"manifest.json")
    assert z.read(n["manifest.json"])==m.read(mn["manifest.json"])
    by={f["name"]:f for f in P["files"]};assert len(by)==24
    def tiles(f):
        return [(t,r.tile(z,n,t)) for t in by[f["name"]]["tiles"]]
    r.tiles=tiles
    return z,m,p,mp,manifest
def compute(q):
    rows,R,_,h=r.gather(q);X=r.norm(np.stack([x[0] for x in rows]));Y=r.register(X,r.norm(R))[0];res={}
    for prefix,V in (("raw",X),("reg",Y)):
        for k,fn in (("exact",lambda v:od.sliced(v,0,None)[:2]),("taylor",a.tree),("linear",a.straight)):
            start=time.perf_counter();s,p=fn(V);res[prefix+"_"+k]=(s,p,time.perf_counter()-start)
    y=np.array([x[1] for x in rows]);truth=np.stack([x[7] for x in rows]);active=np.stack([x[8] for x in rows]);amp=np.array([x[4] for x in rows]);data={}
    for k,(s,p,t) in res.items():
        loc=((np.abs(p-truth)<=2)&active).sum(1)>=.75*active.sum(1)
        data[k]={"scores":s,"localized":loc,"time":t}
    return {"target":q["target"],"y":y,"amp":amp,"arms":data}
def main():
    a.init_tree();a.selftest();z,m,p,mp,manifest=restore()
    cal=[compute(q) for q in manifest["selected"][:2]]
    modes=P["modes"];th={k:r.cut(np.concatenate([x["arms"][k]["scores"][x["y"]==0] for x in cal])) for k in modes};save("lock.json",th)
    ev=[compute(q) for q in manifest["selected"][2:]];d={};y=np.concatenate([x["y"] for x in ev]);amp=np.concatenate([x["amp"] for x in ev])
    for k in modes:
        s=np.concatenate([x["arms"][k]["scores"] for x in ev]);l=np.concatenate([x["arms"][k]["localized"] for x in ev]);ok=s>=th[k]
        # Ranking diagnostic: evaluation-control threshold, explicitly retrospective.
        roc={}
        for n in (5,10,17,25):
            v=np.sort(s[y==0])[::-1];t=np.nextafter(v[n],np.inf);roc[str(n)]={"localized":int(np.sum((s>=t)&l&(y==1))),"controls":int(np.sum((s>=t)&(y==0)))}
        d[k]={"localized":int(np.sum(ok&l&(y==1))),"controls":int(np.sum(ok&(y==0))),"seconds":sum(x["arms"][k]["time"] for x in ev),"weak_075_125":int(np.sum(ok&l&(y==1)&(amp<=1.25))),"retrospective_curve":roc}
    rows=[{"target":x["target"],"y":x["y"].tolist(),"amp":x["amp"].tolist(),"arms":{k:{"scores":v["scores"].tolist(),"localized":v["localized"].tolist()} for k,v in x["arms"].items()}} for x in cal+ev];save("scores.json",rows)
    result={"status":"DEVELOPMENT_EQUAL_INPUT_REPLAY","modes":d,"thresholds":th,"scope":P["scope"],"eval_examples":len(y),"production_adoption":False};save("result.json",result);print("R27A2_SUMMARY",json.dumps(result),flush=True);z.close();m.close();p.unlink();mp.unlink()
if __name__=="__main__":main()
