#!/usr/bin/env python3
import collections,hashlib,io,json,math,os,pathlib,re,sys,zipfile
from urllib.parse import unquote
import numpy as np
import requests,h5py,hdf5plugin
import r15_r16_archive_pilots as old
import r9_expanded_search as eng
import r4_chunk_replay as base
S=json.loads(pathlib.Path("research/seti/r15b_r16b_r17b_protocol.json").read_text())
def emit(k,x):print(k,json.dumps(x),flush=True)
def save(p,x):p.write_text(json.dumps(x,indent=2)+"\n")
def freeze15(out):
    p=S["r15b"];records={};queries=[];examples=[]
    for target in p["targets"]:
        for typ in p["types"]:
            try:
                raw,url=base.fetch(eng.API,{"target":target,"telescope":"Parkes","file-type":typ},limit=32*base.MIB)
                d=json.loads(raw);rs=d if isinstance(d,list) else next((d[k] for k in ("data","files","results") if isinstance(d.get(k),list)),[])
                (out/("index-"+target+"-"+typ+".json")).write_bytes(raw)
                q={"target":target,"type":typ,"rows":len(rs),"sha256":hashlib.sha256(raw).hexdigest(),"date_matches":0};queries.append(q)
                for v in rs:
                    if not isinstance(v,dict):continue
                    u=v.get("url","")
                    if not isinstance(u,str) or not eng.allowed(u):continue
                    b=unquote(base.basename(u))
                    # Match explicit UTC day, numeric MJD, or exact MJD filename token.
                    txt=json.dumps(v,ensure_ascii=False)
                    date=bool(re.search(r"(?<!\d)58602(?:[._-]|\b)|2019-04-29|20190429",txt))
                    if not date:continue
                    q["date_matches"]+=1
                    label=str(v.get("target","")).upper().replace("_","").replace(" ","")
                    role="ON" if "PROX" in label else "CALIBRATOR" if "1421" in label or "1934" in label else None
                    if role is None:continue
                    records[u]={"url":u,"name":b,"role":role,"metadata":v,"type":typ}
                examples.extend([{"target":target,"type":typ,"keys":list(v),"name":base.basename(v.get("url","")),"mjd":v.get("mjd"),"target_label":v.get("target")} for v in rs[:2] if isinstance(v,dict)])
            except Exception as e:queries.append({"target":target,"type":typ,"error":str(e)})
    selected=[];attempts=[]
    for role in ("ON","CALIBRATOR"):
        # No filename resolution assumptions. Try actual HDF5 names in order.
        pool=sorted((v for v in records.values() if v["role"]==role and v["name"].lower().endswith((".h5",".hdf5"))),key=lambda v:v["name"])
        for f in pool[:80]:
            try:
                h=old.header(f["url"])
                if int(h["mjd"])!=58602 or not old.eligible(h,982.002571,2500):raise ValueError("header date/coverage/resolution")
                selected.append({**f,"header":h})
                if sum(v["role"]==role for v in selected)>=p["limit_per_role"]:break
            except Exception as e:attempts.append({"name":f["name"],"error":str(e)})
    m={"queries":queries,"examples":examples,"exact_records":list(records.values()),"attempts":attempts,"selected":selected}
    save(out/"manifest.json",m)
    emit("R15B_MANIFEST",{"queries":queries,"examples":examples,"exact_names":[{"name":v["name"],"role":v["role"],"type":v["type"]} for v in records.values()][:300],"attempts":attempts,"selected":[{"name":v["name"],"role":v["role"],"url":v["url"],"header":v["header"]} for v in selected]})
def scan15(out):
    m=json.loads((out/"manifest.json").read_text());rows=[]
    for f in m["selected"]:
        try:
            x,t=old.slice_read(f["url"],f["header"],982.002571,2500,out)
            h=old.blc_top(x,t,.15)
            r={"name":f["name"],"role":f["role"],"tile":t,"best":h,"strong":bool(h and h["score"]>=12)}
        except Exception as e:r={"name":f["name"],"role":f["role"],"error":str(e)}
        rows.append(r);emit("R15B_SCAN",r)
    s={"indexed_matching_products":len(m["exact_records"]),"selected":len(m["selected"]),"read":sum("tile" in r for r in rows),"strong_on":sum(r.get("strong",False) and r["role"]=="ON" for r in rows),"strong_calibrator":sum(r.get("strong",False) and r["role"]=="CALIBRATOR" for r in rows),"status":"TARGETED_REPLAY_ONLY" if any("tile" in r for r in rows) else "HOLD_NO_SPECTRAL_READ","confirmed_technosignatures":0}
    save(out/"result.json",{"summary":s,"rows":rows});emit("R15B_SUMMARY",s)
def archive(a):
    u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['id']}/zip"
    buf=bytearray()
    with requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},stream=True,timeout=(15,60)) as r:
        r.raise_for_status()
        for b in r.iter_content(1048576):
            buf.extend(b);assert len(buf)<=160*1024**2
    assert hashlib.sha256(buf).hexdigest()==a["digest"]
    z=zipfile.ZipFile(io.BytesIO(buf));assert sum(i.file_size for i in z.infolist())<=240*1024**2
    names={}
    for i in z.infolist():
        b=pathlib.PurePosixPath(i.filename).name
        if b in ("manifest.json","result.json","r9-result.json") or re.fullmatch("[0-9a-f]{24}\\.npy",b):
            assert b not in names;names[b]=i.filename
    return z,names
def wide(y,t,f,dr):
    df=t["df_hz"];n,m=y.shape;ts=t["ts"];c=(f-t["channel0_mhz"])*1e6/df;s0=round(dr*(n-1)*ts/df)
    dsmax=math.ceil(.2*(n-1)*ts/abs(df));rad=math.ceil(100/abs(df))
    best=None
    for s in range(s0-dsmax,s0+dsmax+1):
        r=eng.scores(y,s)
        if r is None:continue
        v,lo,ds=r
        a=max(lo,round(c)-rad);b=min(lo+len(v),round(c)+rad+1)
        if b<=a:continue
        for q in np.argsort(v[a-lo:b-lo])[-5:]:
            k=a+int(q);pers=float(np.mean(y[np.arange(n),k+ds]>3))
            if pers<.5:continue
            h={"score":float(v[k-lo]),"persistence":pers,"frequency_mhz":t["channel0_mhz"]+k*df/1e6,"drift_hz_s":s*df/((n-1)*ts)}
            if best is None or h["score"]>best["score"]:best=h
    return best
def run16(out):
    z,names=archive(S["r16b"]["artifact"]);d=json.loads(z.read(names["result.json"]));m=json.loads(z.read(names["manifest.json"]))
    candidates=[x for x in d["forward"] if x["status"]=="PRELIMINARY_DIRECTIONAL_HOLD"];assert len(candidates)==1
    h=candidates[0];F=h["frequency_mhz"];D=h["drift_hz_s"];anchor=m["selected"][0]
    # The forward candidate may come from ON1 or ON2, infer anchor from recorded bin/score by exact replay.
    matches=[]
    for i in (0,2):
        r=d["rows"][i]
        if any(a["bin"]==h["bin"] and a["shift"]==h["shift"] and abs(a["score"]-h["score"])<1e-5 for a in r["hits"]):matches.append(i)
    assert len(matches)==1
    anchor=m["selected"][matches[0]];T=anchor["header"]["mjd"]
    frozen={"candidate":h,"anchor_role":anchor["role"],"anchor_mjd":T,"source_artifact":S["r16b"]["artifact"],"status":"FROZEN_BEFORE_REPLAY"}
    save(out/"frozen_candidate.json",frozen);emit("R16B_CANDIDATE",frozen)
    rows=[]
    for r,o in zip(d["rows"],m["selected"]):
        t=r["tile"];x=np.load(io.BytesIO(z.read(names[t["artifact"]])),allow_pickle=False)
        assert list(x.shape)==t["shape"] and hashlib.sha256(x.tobytes()).hexdigest()==t["slice_sha256"]
        y=eng.normalized(x);dt=(o["header"]["mjd"]-T)*86400
        models=[("fixed",F,D),("zero",F,0),("continuous",F+D*dt/1e6,D)]
        tests={}
        for k,f,dr in models:
            a=eng.matched(y,t,f,dr);w=wide(y,t,f,dr)
            tests[k]={"predicted_mhz":f,"match":a,"strong":eng.present(a),"exploratory_wide":w,"wide_strong":eng.present(w)}
        pointing={}
        try:
            with base.RangeFile(o["url"],8*base.MIB) as rem:
                with h5py.File(rem,"r") as hf:
                    aa=hf["data"].attrs
                    pointing={k:str(aa.get(k,"")) for k in ("source_name","src_raj","src_dej","tstart")}
                    assert float(aa["tstart"])==o["header"]["mjd"]
        except Exception as e:pointing={"error":str(e)}
        row={"role":o["role"],"name":o["name"],"pointing":pointing,"dt_seconds":dt,"tests":tests};rows.append(row);emit("R16B_SCAN",row)
    z.close()
    off=[r for r in rows if r["role"].startswith("OFF")];on=[r for r in rows if r["role"].startswith("ON")]
    strong=lambda r:any(t["strong"] for t in r["tests"].values())
    widehit=lambda r:any(t["wide_strong"] for t in r["tests"].values())
    s={"scans":len(rows),"anchor":h,"on_strong":sum(strong(r) for r in on),"off_strong":sum(strong(r) for r in off),"off_wide_strong":sum(widehit(r) for r in off),"status":"RFI_CONCERN_REVIEW" if any(strong(r) for r in off) else "HOLD_INDEPENDENT_CONFIRMATION","pointing_headers":[r["pointing"] for r in rows],"confirmed_technosignatures":0,"independently_confirmed":0}
    save(out/"result.json",{"summary":s,"rows":rows});emit("R16B_SUMMARY",s)
if __name__=="__main__":
    lane=sys.argv[1];out=pathlib.Path("runtime/"+lane);out.mkdir(parents=True,exist_ok=True)
    if lane=="r15b": {"freeze":freeze15,"scan":scan15}[sys.argv[2]](out)
    else:run16(out)
