#!/usr/bin/env python3
"""Separate DC-touching hypotheses, without removing measurements or stationary sources."""
import json,pathlib,time,numpy as np,h5py,hdf5plugin
import r26f_voyager as voy
import r9_expanded_search as eng
OUT=pathlib.Path("runtime/r27c");OUT.mkdir(parents=True,exist_ok=True)
def bank(V,forbidden):
    n,w=V.shape;best=None;clean=None;t=np.arange(n)
    for sh in range(-32,33):
        ds=np.rint(t*sh/(n-1)).astype(int);lo=max(0,-int(ds.min()));hi=min(w,w-int(ds.max()))
        a=sum(V[i,lo+k:hi+k] for i,k in enumerate(ds))/4
        ids=np.arange(lo,hi);touch=np.any(np.isin(ids[:,None]+ds,forbidden),axis=1)
        for mode,mask in (("raw",np.ones(len(ids),bool)),("separated",~touch)):
            j=np.flatnonzero(mask)
            if not len(j):continue
            k=int(j[a[j].argmax()]);rec={"score":float(a[k]),"start_bin":int(ids[k]),"shift":sh,"touches_center":bool(touch[k])}
            if mode=="raw" and (best is None or rec["score"]>best["score"]):best=rec
            if mode=="separated" and (clean is None or rec["score"]>clean["score"]):clean=rec
    return best,clean
def main():
    p,load=voy.download()
    with h5py.File(p,"r") as f:
        d=f["data"];a=d.attrs;x=np.asarray(d[:,0,:],np.float32);f0=float(a["fch1"]);df=float(a["foff"])*1e6;ts=float(a["tsamp"])
    assert x.shape==(16,1048576);V=eng.normalized(x);n,w=V.shape;dc=w//2;forbidden=np.arange(dc-1,dc+2)
    # Full band and full R26F drift range; separate, not delete, DC-touching tracks.
    best={"raw":None,"separated":None,"center_touching":None};t0=time.perf_counter()
    for sh in range(-68,69):
        z,lo,d=eng.scores(V,sh);ids=np.arange(lo,lo+len(z));touch=np.zeros(len(ids),bool)
        for k in np.unique(d):
            for f in forbidden:
                j=f-k-lo
                if 0<=j<len(z):touch[j]=True
        for key,mask in (("raw",np.ones(len(z),bool)),("separated",~touch),("center_touching",touch)):
            if not mask.any():continue
            j=int(np.argmax(np.where(mask,z,-np.inf)));k=lo+j;rec={"score":float(z[j]),"frequency_mhz":f0+k*df/1e6,"drift_hz_s":sh*df/(15*ts),"start_bin":k,"touches_center":bool(touch[j])}
            if best[key] is None or rec["score"]>best[key]["score"]:best[key]=rec
    # Synthetic control: stationary lines at/away center, and crossing center.
    rng=np.random.default_rng(2703);tests=[]
    for name,track in (("stationary_center",np.full(16,64)),("stationary_off_center",np.full(16,40)),("crossing_center",np.arange(56,72))):
        v=rng.normal(0,.1,(16,128));v[np.arange(16),track]+=8
        raw,clean=bank(v,np.arange(63,66));tests.append({"case":name,"raw":raw,"separated":clean,"raw_recovered":abs(raw["start_bin"]-track[0])<=1 and abs(raw["shift"]-(track[-1]-track[0]))<=1})
    assert all(q["raw_recovered"] for q in tests)
    q=best["separated"];ok=abs(q["frequency_mhz"]-voy.P["expected"]["frequency_mhz"])*1e6<=100 and abs(q["drift_hz_s"]-voy.P["expected"]["drift_hz_s"])<=.06
    d={"status":"SEPARATE_RANKING_ONLY","center_bin":dc,"center_frequency_mhz":f0+dc*df/1e6,"best":best,"voyager_first_after_separation":ok,"synthetic_controls":tests,"seconds":time.perf_counter()-t0,"source_sha256":voy.P["sha256"],"scope":"Fixed central 3-bin mask on candidate trajectories, raw ranking retained. Instrumental origin of center spike is not independently verified. Crossing or stationary physical signals may touch center; no automatic physical veto.","production_adoption":False}
    (OUT/"result.json").write_text(json.dumps(d,indent=2)+"\n");print("R27C_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
