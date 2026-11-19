#!/usr/bin/env python3
"""R40N2: independent archive acquisition; never pool incompatible instruments."""
import sys
sys.argv[1:]=["controls"]
import json,pathlib,time,hashlib,re,io,math
import numpy as np
import r40n_raw_rivals as n
import r37b_source_header as qual
import r37_source_scout as scout
c=n.c;b=n.b;h=n.h;l=n.l;f=n.f;t=n.t
P=pathlib.Path("research/seti/r40n2_instrument_controls_protocol.json");S=json.loads(P.read_text());C=S["selection"]
O=pathlib.Path("runtime/r40n2");O.mkdir(parents=True,exist_ok=True);n.O=O;n.k.O=O;n.l.O=O
eng=b.sh.eng
def save(k,v):(O/k).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def emit(k,v):print(k,json.dumps(v,allow_nan=False),flush=True)
def archive(k):return c.artifact(S["sources"][k],64)
def index_query(tel,target,query):
    raw,u=scout.body("https://seti.berkeley.edu/opendata/api/query-files",{"target":query,"telescope":tel,"file-type":"HDF5"})
    d=json.loads(raw);rows=d if isinstance(d,list) else next((d[k] for k in ("data","files","results") if isinstance(d.get(k),list)),[])
    out=[]
    for row in rows:
        if not isinstance(row,dict):continue
        u0=row.get("url","")
        if not isinstance(u0,str) or not eng.allowed(u0):continue
        name=scout.urllib.parse.unquote(scout.urllib.parse.urlparse(u0).path.rsplit("/",1)[-1])
        try:mjd=float(row.get("mjd"));day=int(mjd)
        except (ValueError,TypeError,OverflowError):continue
        if not name.lower().endswith(".h5"):continue
        out.append({"name":name,"url":u0,"day":day,"mjd":mjd,"metadata":row})
    return {"telescope":tel,"target":target,"query":query,"index_sha256":c.sha(raw),"files":out,"rows":len(rows),"url":u}
def select(old):
    names,days,ns,urls,coords=qual.history()
    hm,m,fr,sel=n.manifests();zz,nn=n.archive("r37_manifest");rm=n.read(zz,nn,"manifest.json")
    cde=json.loads(pathlib.Path("research/seti/r40cde_protocol.json").read_text())
    names.update(sel["excluded_names"]);names.update(cde["targets"]);days.update(sel["excluded_days"]);days.update(cde["days"])
    prior=hm["selected"]+m["selected"]+rm["selected"]["Parkes"]+rm["selected"].get("GBT",[])+old["selected"]
    # Identical or ambiguous historical declinations use R37B conservative sky exclusion.
    for g in prior:
        names.add(g["target"]);days.add(int(g["day"]))
        for v in g["files"]:ns.add(v["name"]);urls.add(v["url"])
        hd=g["files"][0]["header"];mt=g["files"][0].get("metadata",{})
        coords.append({"target":g["target"],**qual.history_pointing(hd,mt)})
    raw,u=scout.body(scout.S["index"]);strings=set(scout.walk(json.loads(raw)))
    (O/"target-index.json").write_bytes(raw)
    aliases={x[:-2].upper() for x in strings if re.fullmatch(r"[A-Za-z0-9_-]{2,48}_S",x,re.I) and not re.fullmatch(r"HIP\d+_S",x,re.I)}
    gb={x[:-4].upper() for x in strings if re.fullmatch(r"[A-Za-z0-9_-]{2,48}_OFF",x,re.I)}
    published=set(scout.S["gbt_targets"])
    order=lambda pool:sorted(pool-names,key=lambda x:(x not in published,hashlib.sha256((C["seed"]+"|"+x).encode()).hexdigest()))
    queues={"GBT":order(gb|published)[:C["max_gbt_targets"]],"Parkes":order(aliases)[:C["max_parkes_aliases"]]}
    # Override only metadata windows, not cadence or pointing checks.
    qual.S["parkes"]={**qual.S["parkes"],"windows_mhz":C["parkes_windows_mhz"],"half_hz":C["half_hz"]}
    manifest={"protocol_sha256":c.sha(P.read_bytes()),"target_index_sha256":c.sha(raw),"old_manifest_sha256":c.sha(json.dumps(old,sort_keys=True).encode()),"candidate_counts":{k:len(v) for k,v in queues.items()},"candidate_order":queues,"selected":[],"attempts":[]}
    save("manifest.json",manifest)
    for tel in ("GBT","Parkes"):
        for target in queues[tel]:
            if len(manifest["selected"])>=C["groups"]:break
            if target in names:continue
            a={"telescope":tel,"target":target,"errors":[]};manifest["attempts"].append(a)
            try:
                rr=[index_query(tel,target,q) for q in ((target,target+"_OFF") if tel=="GBT" else (target,))]
                fs=qual.entries(rr,tel,target,days,ns,urls);cand=qual.cadence(fs);a["files"]=len(fs);a["cadences"]=len(cand);a["index_sha256"]=[v["index_sha256"] for v in rr]
                for group in cand[:C["header_attempts"]]:
                    try:
                        hs,pp=qual.qualify(tel,target,group,coords)
                        day=int(group[0]["day"]);assert day not in days and target not in names
                        g={"telescope":tel,"target":target,"day":day,"pointing":pp[0],"split":"null_calibration","files":[{**v,"header":hd,"role":role} for v,hd,role in zip(group,hs,("ON1","OFF1","ON2","OFF2"))]}
                        manifest["selected"].append(g);names.add(target);days.add(day);coords.append({"target":target,"pointing":pp[0],"status":"FULL_HEADER_POINTING"});ns.update(v["name"] for v in group);urls.update(v["url"] for v in group)
                        emit("R40N2_SELECTED",{"telescope":tel,"target":target,"day":day,"count":len(manifest["selected"])});break
                    except Exception as e:a["errors"].append({"day":group[0]["day"],"error":type(e).__name__+": "+str(e)[:200]})
            except Exception as e:a["errors"].append({"error":type(e).__name__+": "+str(e)[:200]})
            save("manifest.json",manifest)
    assert len({g["day"] for g in manifest["selected"]})==len({g["target"] for g in manifest["selected"]})==len(manifest["selected"])
    save("manifest.json",manifest);emit("R40N2_MANIFEST",{"targets":[[g["telescope"],g["target"],g["day"]] for g in manifest["selected"]],"queries":len(manifest["attempts"]),"candidate_counts":manifest["candidate_counts"],"sha256":c.sha((O/"manifest.json").read_bytes())})
    return manifest,fr
def selftest():
    a=np.array([1.,2.,2.,4.]);assert l.cut(a,.25)>2
    assert qual.sep((0,0),(0,.2))>.1
    assert qual.sep((359.99,0),(.01,0))<.1
    emit("R40N2_SELFTEST",{"tie_safe_threshold":"PASS","spherical_separation_wrap":"PASS"})
def main():
    start=time.perf_counter();selftest()
    oz,on=archive("r40n1");old=n.read(oz,on,"manifest.json");receipt=n.read(oz,on,"controls.json");assert receipt["groups"]==4 and receipt["controls"]==512
    assert c.sha(oz.read(on["manifest.json"]))==receipt["manifest_sha256"]
    m,fr=select(old);mh=c.sha((O/"manifest.json").read_bytes())
    if not m["selected"]:
        d={"status":"HOLD_NO_NEW_HEADER_QUALIFIED_GROUPS","groups":0,"new_controls":0,"total_controls":512,"manifest_sha256":mh,"production_adoption":False};save("result.json",d);emit("R40N2_SUMMARY",d);return
    L=h.readlocks();b.taylor.init_tree()
    lz,ln=n.archive("r40l");raw=lz.read(ln["lock.json"]);assert c.sha(raw)==n.S["sources"]["r40l"]["lock_sha256"];lock=json.loads(raw)
    rz,rn=archive("r40n");raw=rz.read(rn["lock.json"]);assert c.sha(raw)==S["sources"]["r40n"]["lock_sha256"];rlock=json.loads(raw)
    kz,kn=n.archive("r40k");kr=n.read(kz,kn,"result.json");kl=n.read(kz,kn,"lock.json")
    seen=set(fr["excluded_raw_hashes"]);seen.update(v for r in kr["rows"] for v in r["raw_tile_hashes"]);seen.update(v for r in kl["calibration"] for v in r["raw_tile_hashes"]);seen.update(v for r in receipt["records"] for v in r["raw_tile_hashes"])
    for name in rn:
        if name.endswith(".npz"):
            with np.load(io.BytesIO(rz.read(rn[name])),allow_pickle=False) as q:seen.update(q["source_tile_sha"].tolist())
    controls={"Parkes":[],"GBT":[]};records=[]
    def raw_score(d):
        v=rlock["models"]["raw_logistic"];X=d["raw_features"];return ((X-np.asarray(v["mean"]))/np.asarray(v["scale"]))@np.asarray(v["coef"])[0]+v["intercept"][0]
    for g in old["selected"]:
        name=g["target"];d=f.load(oz,on,"control_"+name+".npz");assert len(d["y"])==128 and not d["y"].any()
        controls["Parkes"].append((name,int(g["day"]),d))
    for g in m["selected"]:
        tel=g["telescope"];name=g["target"];wins=C["windows_mhz"] if tel=="GBT" else C["parkes_windows_mhz"]
        n.set_config(C["seed"]+"|"+tel,wins,C["half_hz"],C["max_mib_per_scan"])
        rr,R,ch,hs,meta,cm,pre=n.gather(g,tel+"_"+name);assert len(hs)==16 and not seen&hs;seen.update(hs)
        ix=np.array([v[1]==0 for v in rr]);assert ix.sum()==128
        rows=[v for v,x in zip(rr,ix) if x];mm=[v for v,x in zip(meta,ix) if x];R=R[ix];d=n.encode(rows,R,mm,pre)
        pan=b.panel(rows,R);arms=b.base_arms(pan,L["r40c"]["original"]["component_fit"])
        q={"V":pan["Y"].astype(np.float32),"y":np.zeros(128,int),"truth":np.stack([v[7] for v in rows]),"active":np.stack([v[8] for v in rows]),"amp":np.array([v[4] for v in rows]),"tile_sha":np.array([v[6] for v in rows])}
        for a in t.S["paths"]:q[a+"_score"],q[a+"_paths"]=arms[a][:2]
        hh=f.historical(q,L["r40c"]);_,_,_,a=h.replay(hh,L);a=a.astype(bool)
        d.update(parent_alert=a,parent_hit=np.zeros(128,bool),joint_rescue=(~a)&(d["s"]>=lock["parent_joint_min_threshold"]),raw_logistic_score=raw_score(d))
        np.savez_compressed(O/("control_"+tel+"_"+name+".npz"),**d);controls[tel].append((name,int(g["day"]),d))
        row={"telescope":tel,"target":name,"day":g["day"],"controls":128,"parent_alerts":int(a.sum()),"joint_extra_alerts":int(d["joint_rescue"].sum()),"raw_tile_hashes":sorted(hs),"native_df_hz":float(g["files"][0]["header"]["df"]*1e6),"native_ts_s":g["files"][0]["header"]["ts"],"windows_mhz":wins}
        records.append(row);save("partial.json",records);emit("R40N2_CONTROL",row);assert c.sha((O/"manifest.json").read_bytes())==mh
    instruments={}
    for tel,pp in controls.items():
        if not pp:continue
        vv=np.concatenate([d["regime"] for _,_,d in pp]);pa=np.concatenate([d["parent_alert"] for _,_,d in pp]);day=np.repeat([v[1] for v in pp],128)
        joint=np.concatenate([d["s"] for _,_,d in pp]);learn=np.concatenate([raw_score(d) for _,_,d in pp]);reject=~pa
        # Parkes existing reference split is not changed; GBT has independent development reference median.
        rc=float(lock["reference_cut"] if tel=="Parkes" else np.median(vv));reg=vv>rc
        stats={str(i):{"controls":int((reject&(reg==i)).sum()),"days":len(set(day[reject&(reg==i)].tolist()))} for i in (0,1)}
        qualified=len(pp)>=4 and all(v["controls"]>=128 and v["days"]>=3 for v in stats.values())
        thresholds={}
        for key,score in (("joint",joint),("raw_logistic",learn)):
            pooled=l.cut(score[reject],.025);th=[l.cut(score[reject&(reg==i)],.025) for i in (0,1)] if qualified else [pooled,pooled]
            thresholds[key]={"pooled":pooled,"conditional":th,"pooled_calibration_alerts":int((reject&(score>=pooled)).sum()),"conditional_calibration_alerts":int((reject&(score>=np.asarray(th)[reg.astype(int)])).sum()),"day_maxima":{str(j):float(score[(day==j)&reject].max()) for j in sorted(set(day.tolist())) if ((day==j)&reject).any()}}
        v={"groups":len(pp),"controls":len(vv),"reference_cut":rc,"strata":stats,"fallback":not qualified,"qualified_counts_only":qualified,"parent_alerts_diagnostic":int(pa.sum()),"joint_extra_alerts_diagnostic":int(sum(d["joint_rescue"].sum() for _,_,d in pp)),"thresholds_development_only":thresholds,"formal_fap":False,"positive_sensitivity_tested":False}
        instruments[tel]=v;emit("R40N2_INSTRUMENT",{"telescope":tel,**v})
    d={"status":"PASS_ACQUISITION_COUNTS_ONLY" if len(records)==4 else "HOLD_PARTIAL_ACQUISITION","manifest_sha256":mh,"groups":len(records),"new_controls":128*len(records),"total_controls":512+128*len(records),"total_groups":4+len(records),"instruments":instruments,"records":records,"cross_instrument_pooling":False,"physical_null_verified":False,"alias_completeness_proven":False,"detector_replay_run":False,"production_adoption":False,"seconds":time.perf_counter()-start}
    save("result.json",d);emit("R40N2_SUMMARY",{k:v for k,v in d.items() if k!="records"})
if __name__=="__main__":main()
