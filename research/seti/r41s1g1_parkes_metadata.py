"""Metadata-only audit of zero-pointing Parkes cadences; never infer actual pointing from target name."""
import hashlib,json,pathlib,os,struct,math,urllib.parse
import numpy as np,h5py,requests
import r40n10c_archive_adapter as old
eng=old.eng;prev=old.prev
P=pathlib.Path("research/seti/r41s1g1_parkes_metadata_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r41s1g1");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(k,d):(O/k).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def read(k):return json.loads((O/k).read_text())
def simple(v):
    if isinstance(v,bytes):return v.decode("utf-8",errors="replace")
    if isinstance(v,np.ndarray):return [simple(x) for x in v.tolist()] if v.ndim else simple(v.item())
    if isinstance(v,np.generic):return simple(v.item())
    if isinstance(v,float) and not math.isfinite(v):return str(v)
    if isinstance(v,(str,int,float,bool,type(None))):return v
    return str(v)
def attrs(url):
    with eng.base.RangeFile(url,8*eng.base.MIB) as r:
        with h5py.File(r,"r") as h:
            out={}
            for k,v in (("root",h),("data",h["data"])):
                out[k]={a:simple(b) for a,b in v.attrs.items()}
            return out
def filindex(target):
    d={}
    for q in (target,target+"_S",target+"_R"):
        try:
            raw,url=eng.base.fetch(eng.API,{"target":q,"telescope":"Parkes","file-type":"filterbank"},limit=16*eng.base.MIB)
            (O/("filindex-"+sha(q.encode())[:16]+".json")).write_bytes(raw)
            j=json.loads(raw);rows=j if isinstance(j,list) else next((j[k] for k in ("data","files","results") if isinstance(j.get(k),list)),[])
            for x in rows:
                u=x.get("url","") if isinstance(x,dict) else ""
                if isinstance(u,str) and eng.allowed(u) and urllib.parse.urlparse(u).path.lower().endswith(".fil"):
                    name=urllib.parse.unquote(urllib.parse.urlparse(u).path.rsplit("/",1)[-1]);d.setdefault(name,{})[u]=x
        except Exception as e:print("R41S1G1_INDEX_HOLD",q,type(e).__name__,str(e)[:160],flush=True)
    return d
def sigproc(url):
    # Strict SIGPROC keyword types, bounded header reads. Unknown keyword => HOLD.
    ints=set("machine_id telescope_id data_type barycentric pulsarcentric nbits nsamples nchans nifs nbeams ibeam sumifs signed".split())
    doubles=set("az_start za_start src_raj src_dej tstart tsamp fch1 foff refdm period fchannel".split())
    strings={"source_name","rawdatafile"}
    with eng.base.RangeFile(url,S["max_header_bytes"]) as r:
        def item(n):
            assert r.tell()+n<=S["max_header_bytes"];v=r.read(n);assert len(v)==n;return v
        def word():
            n=struct.unpack("<i",item(4))[0];assert 0<n<=4096;return item(n).decode("ascii")
        assert word()=="HEADER_START";out={}
        for _ in range(256):
            k=word()
            if k=="HEADER_END":return out
            if k in ("FREQUENCY_START","FREQUENCY_END"):continue
            if k in ints:v=struct.unpack("<i",item(4))[0]
            elif k in doubles:v=struct.unpack("<d",item(8))[0];assert math.isfinite(v)
            elif k in strings:v=word()
            else:raise ValueError("unknown SIGPROC key "+k)
            if k=="fchannel":out.setdefault(k,[]).append(v)
            else:assert k not in out;out[k]=v
        raise ValueError("missing header end")
def coordinate(h):
    def packed(v,ra):
        v=float(v);a=int(abs(v)//10000);b=int(abs(v)//100)%100;c=abs(v)%100
        assert 0<=b<60 and 0<=c<60 and (a<24 if ra else a<=90)
        return (a+b/60+c/3600)*(15 if ra else (-1 if v<0 else 1))
    pp=(packed(h["src_raj"],True),packed(h["src_dej"],False));assert pp!=(0,0)
    return pp
def prepare():
    z,n=old.c.artifact(S["source"],128)
    m=json.loads(z.read(n["same_receiver_manifest.json"]));ctx=json.loads(z.read(n["exclusions.json"]))
    assert m["diagnosis"]["matched"]==80 and m["diagnosis"]["exact_index_snapshot_matches"]==80
    plan=[]
    for item in m["attempts"]:
        for e in item["records"]:
            if e.get("matched") and e.get("error",{}).get("frames",[]) and "r!=0 or d!=0" in e["error"]["frames"][-1].get("code",""):
                plan.append({"target":item["target"],"ordinal":e["ordinal"],"day":e["day"],"files":e["files"],"index_sha256":item["index_sha256"],"previous_headers":e["headers"]});break
    plan=sorted(plan,key=lambda a:sha(("R41S1G1|"+a["target"]).encode()))[:S["max_cadences"]]
    assert plan
    save("plan.json",{"plan":plan,"exclusions":ctx,"protocol_sha256":sha(P.read_bytes()),"original_rejections":80})
    print("R41S1G1_PLAN",len(plan),flush=True)
def group(i):
    lock=read("plan.json");assert lock["protocol_sha256"]==sha(P.read_bytes())
    if i>=len(lock["plan"]):return
    q=lock["plan"][i];ctx=lock["exclusions"];target=q["target"];out={"target":target,"day":q["day"],"records":[],"status":"HOLD_NO_INDEPENDENT_POINTING"}
    try:
        rr,errs=prev.query(target);assert [x["index_sha256"] for x in rr]==q["index_sha256"]
        fs=old.qual.entries(rr,"Parkes",target,set(ctx["days"]),set(ctx["files"]),set(ctx["urls"]))
        g=old.qual.cadence(fs)[q["ordinal"]];assert [f["name"] for f in g]==q["files"]
        fi=filindex(target)
        for f,prior in zip(g,q["previous_headers"]):
            x={"name":f["name"],"url":f["url"],"catalog":simple(f["metadata"]),"original_header":prior["header"]}
            try:
                hd=eng.hdr(f);assert hd==prior["header"]
                x["hdf5_attrs"]=attrs(f["url"])
                matches=fi.get(f["name"][:-3]+".fil",{})
                x["exact_fil_catalog"]=matches
                if len(matches)==1:
                    u=next(iter(matches));h=sigproc(u);x["fil_header"]=h
                    assert str(h["source_name"]).strip().upper()==hd["source"].upper()
                    assert abs(h["tstart"]-hd["mjd"])*86400<=3 and abs(h["foff"]-hd["df"])<=1e-12 and abs(h["tsamp"]-hd["ts"])<=1e-7 and abs(h["fch1"]-hd["f0"])<=1e-8
                    x["fil_nonzero_pointing_candidate_deg"]=coordinate(h)
                    x["candidate_only_not_independently_corroborated"]=True
            except Exception as e:x["hold_reason"]=type(e).__name__+": "+str(e)[:300]
            out["records"].append(x)
        out["status"]="METADATA_DIAGNOSTIC_ONLY"
    except Exception as e:out["error"]=type(e).__name__+": "+str(e)[:300]
    save(f"group-{i:02d}.json",out)
    print("R41S1G1_GROUP",i,target,len(out["records"]),sum("fil_nonzero_pointing_candidate_deg" in r for r in out["records"]),flush=True)
def summary():
    plan=read("plan.json")["plan"];g=[read(f"group-{i:02d}.json") for i in range(len(plan))]
    r=[x for a in g for x in a["records"]]
    d={"status":"HOLD_INDEPENDENT_POINTING_CORROBORATION","cadences":len(g),"files":len(r),"hdf5_attrs_recovered":sum("hdf5_attrs" in x for x in r),"exact_fil_match":sum(len(x.get("exact_fil_catalog",{}))==1 for x in r),"candidate_nonzero_fil_pointings":sum("fil_nonzero_pointing_candidate_deg" in x for x in r),"qualified_fresh_groups":0,"no_power_read":True,"production_adoption":False,"frozen_model_lock":S["model_lock"]}
    save("result.json",d);print("R41S1G1_SUMMARY",json.dumps(d),flush=True)
