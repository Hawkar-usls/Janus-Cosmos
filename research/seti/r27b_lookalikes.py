#!/usr/bin/env python3
import hashlib,json,pathlib,re,sys
from urllib.parse import unquote
import numpy as np
import r4_chunk_replay as base
import r9_expanded_search as eng
import r15c_r16c_spectral as sig
import r26c_reference_alignment as r
P=json.loads(pathlib.Path("research/seti/r27b_lookalike_protocol.json").read_text())
OUT=pathlib.Path("runtime/r27b");OUT.mkdir(parents=True,exist_ok=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def main():
    r.BLC_OUT=OUT;z,n,p=r.archive("blc_manifest");m=r.get(z,n,"manifest.json");z.close();p.unlink()
    anchor=next(f for f in m["selected"] if f["name"]==r.BLC_P["r15d"]["anchor"])
    controls=sorted([f for f in m["selected"] if f["role"]=="CONTROL"],key=lambda f:abs(f["header"]["tstart"]-anchor["header"]["tstart"]))
    picks=[anchor]+controls[:2];times={f["name"].split("_",1)[1]:f["role"] for f in picks};records={}
    for target in ("ProxCen","1421-490"):
        raw,url=base.fetch(eng.API,{"target":target,"telescope":"Parkes","file-type":"HDF5"},limit=32*base.MIB)
        obj=json.loads(raw);rs=obj if isinstance(obj,list) else next((obj[k] for k in ("data","files","results") if isinstance(obj.get(k),list)),[])
        for v in rs:
            u=v.get("url","")
            if not isinstance(u,str) or not eng.allowed(u):continue
            name=unquote(base.basename(u))
            if not re.match(r"blc\d+_",name) or name.split("_",1)[1] not in times:continue
            records[name]={"name":name,"url":u,"role":times[name.split("_",1)[1]],"bytes":int(v["size"])}
    assert len(records)<=P["max_files"]
    centers=[P["anchor_mhz"]+k*P["oscillator_mhz"] for k in P["multipliers"]];out=[];coverage=[]
    for f in sorted(records.values(),key=lambda f:f["name"]):
        try:
            with base.RangeFile(f["url"],8*base.MIB) as rr:
                h=sig.parse_header(rr);assert h["size"]==f["bytes"] and int(h["tstart"])==58602
                lo,hi=sorted([h["fch1"],h["fch1"]+(h["nchans"]-1)*h["foff"]])
                coverage.append({"name":f["name"],"role":f["role"],"lo":lo,"hi":hi,"rows":h["rows"]})
                for fc in centers:
                    if not lo+P["half_hz"]/1e6<fc<hi-P["half_hz"]/1e6:continue
                    x,t=sig.fil_window(rr,h,fc,P["half_hz"]);y=eng.normalized(x).astype(float);nr,nc=y.shape;best=None
                    for dr in P["drifts"]:
                        ds=np.rint(np.arange(nr)*t["ts"]*dr/t["df_hz"]).astype(int);l=max(0,-int(ds.min()));u=min(nc,nc-int(ds.max()))
                        if u<=l:continue
                        a=sum(y[i,l+k:u+k] for i,k in enumerate(ds));med=np.median(a);s=1.4826*np.median(abs(a-med))
                        if s<=0:continue
                        a=(a-med)/s;freq=t["channel0_mhz"]+np.arange(l,u)*t["df_hz"]/1e6;ok=np.flatnonzero(abs(freq-fc)*1e6<=P["search_half_hz"])
                        if not len(ok):continue
                        j=int(ok[a[ok].argmax()]);q={"score":float(a[j]),"frequency_mhz":float(freq[j]),"offset_hz":float((freq[j]-fc)*1e6),"drift_hz_s":dr}
                        if best is None or q["score"]>best["score"]:best=q
                    a={"name":f["name"],"role":f["role"],"center":fc,"tile":t,"best":best,"strong":bool(best and best["score"]>=P["minimum_score"])}
                    out.append(a);np.save(OUT/(hashlib.sha256((f["name"]+str(fc)).encode()).hexdigest()[:24]+".npy"),x,allow_pickle=False)
                    print("R27B_WINDOW",json.dumps(a),flush=True)
        except Exception as e:coverage.append({"name":f["name"],"error":str(e)})
    d={"status":"TARGETED_PUBLISHED_COMPANION_REPLAY","expected_centers":centers,"files_indexed":len(records),"coverage":coverage,"results":out,"strong_on_centers":sorted({v["center"] for v in out if v["role"]=="ON" and v["strong"]}),"scope":P["scope"],"production_adoption":False}
    save("result.json",d);print("R27B_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
