#!/usr/bin/env python3
"""Frozen R40N10 transfer; metadata-only qualification, two nonpooled receiver lanes."""
import sys
ARGS=sys.argv[1:]
import collections,concurrent.futures,hashlib,io,json,math,pathlib,re,time,traceback
import numpy as np
import r40n10_artmoney_turbo as a
n=a.n;c=a.c;h=a.h;b=a.b;qual=a.qual;scout=a.scout;eng=a.eng
P=pathlib.Path("research/seti/r40n10b_expansion_protocol.json");S=json.loads(P.read_text());C=S["selection"]
O=pathlib.Path("runtime/r40n10b");O.mkdir(parents=True,exist_ok=True);CACHE=O/"cache";CACHE.mkdir(exist_ok=True);PROG=O/"progress";PROG.mkdir(exist_ok=True)
a.P=P;a.S=S;a.C=C;a.O=O;a.CACHE=CACHE;a.PROG=PROG
n.O=O;n.k.O=O;n.l.O=O;eng.OUT=O
W=C["windows_mhz"];HALF=C["half_hz"]/1e6+.001
ROLES=("ON1","OFF1","ON2","OFF2")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
a.emit=lambda k,d:emit(k.replace("R40N10_","R40N10B_"),d)
def save(k,d):a.write(O/k,d)
def load(k):return json.loads((O/k).read_text())
def cp(k,d):a.write(PROG/(k+".json"),d);emit("R40N10B_PROGRESS",d)
def require(ok,why):
    if not ok:raise ValueError(why)
def msg(e):return type(e).__name__+": "+str(e)[:240]
def key(x):return re.sub(r"\s+","",str(x)).upper()
def contains(h,w):return h["bounds"][0]<w-HALF and h["bounds"][1]>w+HALF
def query(target):
    rows=[];errors=[]
    for q in (target,target+"_S",target+"_R"):
        try:rows.append(a.index_query("Parkes",target,q))
        except Exception as e:errors.append({"query":q,"error":msg(e)})
    return rows,errors
def old_source():
    z,nn=c.artifact(S["source"],128)
    for name in ("lock","manifest"):
        require(c.sha(z.read(nn[name+".json"]))==S["source"][name+"_sha256"],"N10 "+name+" digest")
    plan=n.read(z,nn,"plan.json");receipt=n.read(z,nn,"freeze_receipt.json")
    require(receipt["plan_sha256"]==c.sha(z.read(nn["plan.json"])) and receipt["lock_sha256"]==S["source"]["lock_sha256"],"N10 plan provenance")
    require(len(plan["excluded_raw_hashes"])==864,"prior exclusion hash count")
    return z,nn,plan
def prepare():
    z,nn,prior=old_source();m=n.read(z,nn,"manifest.json")
    (O/"lock.json").write_bytes(z.read(nn["lock.json"]))
    save("parent_manifest.json",m)
    save("old_rejection.json",{"queries":len(m["attempts"]),"targets_with_cadence":sum(x.get("cadences",0)>0 for x in m["attempts"]),"cadences":sum(x.get("cadences",0) for x in m["attempts"]),"reasons":dict(collections.Counter(e["error"] for x in m["attempts"] for e in x["errors"])),"zero_eligible_files":sum(x.get("files")==0 for x in m["attempts"])})
    names,days,ns,urls,coords=qual.history()
    hm,mh,fr,sel=n.manifests();zz,nn0=n.archive("r37_manifest");rm=n.read(zz,nn0,"manifest.json")
    cde=json.loads(pathlib.Path("research/seti/r40cde_protocol.json").read_text())
    names.update(sel["excluded_names"]);names.update(cde["targets"]);days.update(sel["excluded_days"]);days.update(cde["days"])
    prev=hm["selected"]+mh["selected"]+rm["selected"]["Parkes"]+rm["selected"].get("GBT",[])
    for k in ("n1","n2","n6","n7","n8","n9"):
        zz,nn0=a.checked_archive(k);prev+=n.read(zz,nn0,"manifest.json")["selected"]
    for g in prev:
        names.add(g["target"]);days.add(int(g["day"]))
        for v in g["files"]:ns.add(v["name"]);urls.add(v["url"])
        f=g["files"][0];coords.append({"target":g["target"],**qual.history_pointing(f["header"],f.get("metadata",{}))})
    names.update(S["published_uwl"]["excluded"]);names.update(["3I","3IATLAS","3I-ATLAS","ATLAS","3I_ATLAS","A11PL3Z"])
    raw,_=scout.body(scout.S["index"]);idx=sorted(set(scout.walk(json.loads(raw))));(O/"target-index.json").write_bytes(raw)
    reference=prior["development"][0]["group"]["files"][0]["header"]
    context={"names":sorted(names),"days":sorted(days),"files":sorted(ns),"urls":sorted(urls),"coords":coords,"excluded_raw_hashes":prior["excluded_raw_hashes"],"reference":{"df":reference["df"],"ts":reference["ts"]},"index_sha256":c.sha(raw)}
    save("exclusions.json",context)
    h.O=CACHE;L=h.readlocks();a.write(CACHE/"parent_locks.json",L)
    for k in ("n7","n8","n9"):
        zz,nn0=a.checked_archive(k);(CACHE/(k+"_lock.json")).write_bytes(zz.read(nn0["lock.json"]))
    lz,ln=n.archive("r40l");raw=lz.read(ln["lock.json"]);require(c.sha(raw)==n.S["sources"]["r40l"]["lock_sha256"],"joint lock digest");(CACHE/"joint_lock.json").write_bytes(raw)
    require(c.sha((O/"lock.json").read_bytes())==S["source"]["lock_sha256"],"frozen N10 lock")
    cp("00-prepare",{"stage":"PREPARE","excluded_names":len(names),"excluded_days":len(days),"excluded_hashes":len(prior["excluded_raw_hashes"]),"original_rejection":load("old_rejection.json"),"model_sha256":S["source"]["lock_sha256"]})
def exclude(context,selected):
    names=set(context["names"]);days=set(context["days"]);ns=set(context["files"]);urls=set(context["urls"]);coords=list(context["coords"])
    for g in selected:
        names.add(g["target"]);days.add(g["day"]);coords.append({"target":g["target"],"pointing":g["pointing"]})
        for v in g["files"]:
            for f in v.get("subbands",[v]):ns.add(f["name"]);urls.add(f["url"])
    return names,days,ns,urls,coords
def cadence_check(pp,coords):
    require(qual.sep(pp[0],pp[2])<.001,"ON pointing mismatch")
    require(qual.sep(pp[0],pp[1])>=.01 and qual.sep(pp[0],pp[3])>=.01,"ON/OFF not separated")
    require(pp[1]==pp[3],"OFF pointing mismatch")
    require(all(qual.history_distance(pp[0],x)>=qual.A["pointing_separation_degrees"] for x in coords),"Prior or duplicate sky pointing")
def headers_same(target,group,coords,reference,cache):
    hs=[]
    for f in group:
        if f["url"] not in cache:cache[f["url"]]=eng.hdr(f)
        hs.append(cache[f["url"]])
    for hd,v in zip(hs,group):
        require(len(hd["shape"])==3 and hd["shape"][0]==16 and hd["shape"][1]>=1,"not 16 native integrations")
        require(0<abs(hd["df"])*1e6<=10 and 0<hd["ts"]<=60,"spectral/time resolution invalid")
        require(abs((hd["mjd"]-v["day"])*86400-v["seconds"])<=3,"filename/header time mismatch")
        require(all(contains(hd,w) for w in W),"missing original four frequency windows")
        require(key(hd["source"])==target+"_"+v["position"],"source header mismatch")
        require(abs(hd["df"]-reference["df"])<=1e-12 and abs(hd["ts"]-reference["ts"])<=1e-7,"different native 10cm resolution/time")
    require(len({(x["shape"][0],round(x["ts"],7),round(x["df"],12)) for x in hs})==1,"within cadence sampling mismatch")
    df=hs[0]["df"];require(all(abs((x["f0"]-hs[0]["f0"])/df-round((x["f0"]-hs[0]["f0"])/df))<.001 for x in hs),"native grid mismatch")
    pp=[qual.pointing(hd,v["metadata"]) for hd,v in zip(hs,group)]
    cadence_check(pp,coords)
    return hs,pp
def uwl_scan_candidates(rows,target,days,ns,urls):
    out={}
    for q in rows:
        for f in q["files"]:
            if f["day"] in days or f["name"] in ns or f["url"] in urls or not eng.allowed(f["url"]):continue
            m=f["metadata"];role=key(m.get("target",""))
            if role not in (target+"_S",target+"_R"):continue
            # Native high-resolution products only. Header resolution checked later.
            if not re.search(r"(?:fine|rawspec\.0000|seti)",f["name"],re.I):continue
            try:fc=float(m.get("center_freq"))
            except (TypeError,ValueError):continue
            if not math.isfinite(fc) or fc<=0:continue
            if fc>1e6:fc/=1e6
            k=(f["day"],round(float(f["mjd"]),7),role)
            out.setdefault(k,{})[f["url"]]={**f,"center_meta":fc,"position":role[-1],"kind":"ON" if role[-1]=="S" else "OFF"}
    scans=[]
    for (day,mjd,role),ff in sorted(out.items()):
        bywin=[]
        for w in W:
            cc=sorted((v for v in ff.values() if abs(v["center_meta"]-w)<=80),key=lambda v:(abs(v["center_meta"]-w),v["name"],v["url"]))
            if not cc:break
            bywin.append(cc[:3])
        if len(bywin)==4:scans.append({"day":day,"mjd":mjd,"kind":"ON" if role.endswith("_S") else "OFF","position":role[-1],"candidates":bywin})
    return scans
def uwl_cadences(scans):
    out=[];by=collections.defaultdict(list)
    for s in scans:by[s["day"]].append(s)
    for day,ss in sorted(by.items()):
        ss.sort(key=lambda q:q["mjd"])
        for i in range(len(ss)-3):
            v=ss[i:i+4]
            if [x["kind"] for x in v]!=["ON","OFF","ON","OFF"]:continue
            if all(0<(y["mjd"]-x["mjd"])*86400<=qual.A["max_gap_seconds"]+5 for x,y in zip(v,v[1:])):out.append(v)
    return out
def qualify_uwl(target,cad,coords,cache):
    files=[];pp=[]
    for s,role in zip(cad,ROLES):
        parts=[];p0=[];times=[]
        for w,candidates in zip(W,s["candidates"]):
            errs=[];part=None
            for f in candidates:
                try:
                    if f["url"] not in cache:cache[f["url"]]=eng.hdr(f)
                    hd=cache[f["url"]]
                    require(key(hd["source"])==target+"_"+s["position"],"UWL source mismatch")
                    require(len(hd["shape"])==3 and 16<=hd["shape"][0]<=256 and hd["shape"][1]==1,"UWL expected >=16 native rows of one intensity plane")
                    require(0<abs(hd["df"])*1e6<=10 and 0<hd["ts"]<=60,"UWL spectral/time resolution")
                    require(contains(hd,w),"UWL native window not covered")
                    require(hd["bounds"][1]-hd["bounds"][0]<=150,"not UWL subband")
                    require(abs((hd["mjd"]-s["mjd"])*86400)<=5,"UWL metadata/header timestamp mismatch")
                    pt=qual.pointing(hd,f["metadata"])
                    part={**f,"header":hd,"role":role,"window_mhz":w};p0.append(pt);times.append(hd["mjd"]);break
                except Exception as e:errs.append(msg(e))
            require(part is not None,"UWL window "+str(w)+": "+" | ".join(errs)[:300]);parts.append(part)
        require(len({x["url"] for x in parts})==4,"four distinct UWL window files")
        require(max(times)-min(times)<1e-8,"UWL subbands not simultaneous")
        require(all(qual.sep(p0[0],x)<.001 for x in p0),"UWL subbands pointings differ")
        require(len({(x["header"]["shape"][0],round(x["header"]["ts"],7),round(x["header"]["df"],12)) for x in parts})==1,"UWL subband integration mismatch")
        files.append({**parts[0],"subbands":parts,"role":role});pp.append(p0[0])
    require(len({x["url"] for f in files for x in f["subbands"]})==16,"UWL source file uniqueness")
    for j in range(4):
        hh=[f["subbands"][j]["header"] for f in files]
        require(len({(round(x["ts"],7),round(x["df"],12)) for x in hh})==1,"UWL cadence resolution mismatch")
        df=hh[0]["df"];require(all(abs((x["f0"]-hh[0]["f0"])/df-round((x["f0"]-hh[0]["f0"])/df))<.001 for x in hh),"UWL cadence native grid mismatch")
    hh=[f["header"] for f in files]
    require(all(0<(y["mjd"]-x["mjd"])*86400<=qual.A["max_gap_seconds"] for x,y in zip(hh,hh[1:])),"UWL native cadence time mismatch")
    cadence_check(pp,coords)
    return files,pp
def select(lane):
    started=time.perf_counter();context=load("exclusions.json");parent=load("parent_manifest.json")
    earlier=load("same_receiver_manifest.json")["selected"] if lane=="uwl" else []
    names,days,ns,urls,coords=exclude(context,earlier)
    idx=set(scout.walk(json.loads((O/"target-index.json").read_bytes())))
    aliases={key(x[:-2]) for x in idx if re.fullmatch(r"[A-Za-z0-9_-]{2,48}_S",x,re.I)}
    fixed=lambda t:hashlib.sha256((C["seed"]+"|"+lane+"|"+t).encode()).hexdigest()
    if lane=="same_receiver":
        old={x["target"]:x for x in parent["attempts"]}
        order=sorted((set(old)|aliases)-names,key=lambda t:(not bool(old.get(t,{}).get("cadences",0)),fixed(t)))[:C["max_targets"]]
    else:
        pub=set(S["published_uwl"]["targets"]);tic={key(x) for x in idx if re.fullmatch(r"TIC\d+(?:_[SR])?",x,re.I)}
        tic={re.sub(r"_[SR]$","",x) for x in tic}
        order=sorted((pub|tic)-names,key=lambda t:(t not in pub,fixed(t)))[:100]
    order=[t for t in order if "ATLAS" not in t and t not in names]
    manifest={"lane":lane,"protocol_sha256":c.sha(P.read_bytes()),"parent_manifest_sha256":S["source"]["manifest_sha256"],"candidate_order":order,"selected":[],"attempts":[]}
    save(lane+"_manifest.json",manifest)
    # Bounded parallel metadata query batches; header qualification deterministic and cached per file.
    cache={}
    with concurrent.futures.ThreadPoolExecutor(max_workers=C["workers"]) as pool:
        for k in range(0,len(order),12):
            batch=order[k:k+12];results=list(pool.map(query,batch))
            for target,(rr,errors) in zip(batch,results):
                if len(manifest["selected"])>=4:break
                if target in names:continue
                rec={"target":target,"queries":len(rr),"index_sha256":[x["index_sha256"] for x in rr],"query_errors":errors,"errors":[]}
                if lane=="same_receiver":
                    fs=qual.entries(rr,"Parkes",target,days,ns,urls);cand=qual.cadence(fs);rec.update(files=len(fs),cadences=len(cand))
                else:
                    fs=uwl_scan_candidates(rr,target,days,ns,urls);cand=uwl_cadences(fs);rec.update(scans=len(fs),cadences=len(cand))
                for group in cand[:C["header_attempts"]]:
                    try:
                        day=group[0]["day"];require(day not in days,"previous observation day")
                        if lane=="same_receiver":
                            hs,pp=headers_same(target,group,coords,context["reference"],cache)
                            files=[{**v,"header":hd,"role":role} for v,hd,role in zip(group,hs,ROLES)]
                        else:files,pp=qualify_uwl(target,group,coords,cache)
                        g={"telescope":"Parkes","receiver_lane":lane,"target":target,"day":day,"split":"holdout","pointing":pp[0],"files":files}
                        manifest["selected"].append(g);names.add(target);days.add(day);coords.append({"target":target,"pointing":pp[0]})
                        for ff in files:
                            for f in ff.get("subbands",[ff]):ns.add(f["name"]);urls.add(f["url"])
                        emit("R40N10B_SELECTED",{"lane":lane,"target":target,"day":day,"count":len(manifest["selected"])})
                        break
                    except Exception as e:rec["errors"].append({"day":group[0]["day"],"error":msg(e)})
                manifest["attempts"].append(rec)
                save(lane+"_manifest.json",manifest)
            if (k//12)%3==0:emit("R40N10B_SCOUT",{"lane":lane,"processed":len(manifest["attempts"]),"total":len(order),"qualified":len(manifest["selected"])})
            if len(manifest["selected"])>=4:break
    manifest["eligible"]=len(manifest["selected"])>=2
    save(lane+"_manifest.json",manifest)
    cp("01-"+lane,{"stage":"SELECTION","lane":lane,"qualified":len(manifest["selected"]),"eligible":manifest["eligible"],"attempts":len(manifest["attempts"]),"cadences":sum(x["cadences"] for x in manifest["attempts"]),"rejections":dict(collections.Counter(y["error"] for x in manifest["attempts"] for y in x["errors"])),"manifest_sha256":c.sha((O/(lane+"_manifest.json")).read_bytes()),"seconds":time.perf_counter()-started})
def freeze():
    context=load("exclusions.json");lanes={k:load(k+"_manifest.json") for k in ("same_receiver","uwl")}
    allg=sum([m["selected"] for m in lanes.values()],[])
    require(len({g["day"] for g in allg})==len(allg) and len({g["target"] for g in allg})==len(allg),"cross lane overlap")
    fresh=sum([m["selected"] if m["eligible"] else [] for m in lanes.values()],[])
    m={"status":"FROZEN_METADATA_ONLY","source_n10_lock_sha256":S["source"]["lock_sha256"],"lanes":{k:{"manifest_sha256":c.sha((O/(k+"_manifest.json")).read_bytes()),"qualified":len(v["selected"]),"eligible":v["eligible"],"targets":[{"target":g["target"],"day":g["day"],"pointing":g["pointing"]} for g in v["selected"]]} for k,v in lanes.items()},"selected":fresh}
    save("manifest.json",m)
    plan={"protocol_sha256":c.sha(P.read_bytes()),"manifest_sha256":c.sha((O/"manifest.json").read_bytes()),"development":[],"fresh":fresh,"excluded_raw_hashes":context["excluded_raw_hashes"],"fresh_allowed":bool(fresh)}
    save("plan.json",plan)
    require(c.sha((O/"lock.json").read_bytes())==S["source"]["lock_sha256"],"model changed")
    save("freeze_receipt.json",{"lock_sha256":S["source"]["lock_sha256"],"plan_sha256":c.sha((O/"plan.json").read_bytes())})
    cp("02-freeze",{"stage":"FROZEN","lanes":m["lanes"],"total_groups":len(fresh),"lock_sha256":S["source"]["lock_sha256"],"manifest_sha256":plan["manifest_sha256"],"no_training":True})
def run_group(lane,slot):
    plan=a.verified_plan();gs=[(i,g) for i,g in enumerate(plan["fresh"]) if g["receiver_lane"]==lane]
    label=lane+"-"+str(slot+1).zfill(2)
    if slot>=len(gs):
        cp(label,{"stage":"SKIP","lane":lane,"slot":slot+1,"reason":"not qualified before freeze"});return
    i,g=gs[slot];original=n.ref.tiles;full=[]
    def uwl_tiles(f):
        if "subbands" not in f:return original(f)
        spec=eng.SPEC;out=[]
        try:
            for x in f["subbands"]:
                eng.SPEC={**spec,"window_centers_mhz":[x["window_mhz"]],"limits":{**spec["limits"],"spectrum_budget_per_file_mib":80}}
                rr=original(x);require(len(rr)==1,"one UWL native window")
                md,X=rr[0];full.append(md["slice_sha256"]);require(X.shape[0]>=16,"UWL native time samples")
                Y=np.ascontiguousarray(X[:16]);out.append(({**md,"original_shape":list(X.shape),"original_slice_sha256":md["slice_sha256"],"shape":list(Y.shape),"slice_sha256":c.sha(Y.tobytes()),"native_row_selection":[0,16]},Y))
        finally:eng.SPEC=spec
        return out
    if lane=="uwl":n.ref.tiles=uwl_tiles
    try:
        a.group("fresh",i)
        path=PROG/("fresh-"+str(i+1).zfill(2)+".json")
        rec=json.loads(path.read_text());rec.update(lane=lane,lane_slot=slot+1,stage="MEASURED",full_native_tile_hashes=full)
        if full:
            require(len(full)==len(set(full))==16,"UWL full native hash count")
            seen=set(plan["excluded_raw_hashes"])
            for fn in PROG.glob("*-*.json"):
                if fn.name==path.name:continue
                d=json.loads(fn.read_text())
                if d.get("stage")=="MEASURED":seen.update(d.get("raw_tile_hashes",[]));seen.update(d.get("full_native_tile_hashes",[]))
            require(not set(full)&seen and not set(rec["raw_tile_hashes"])&seen,"UWL full or cropped prior hash overlap")
        cp(label,rec)
    except Exception as e:
        traceback.print_exc()
        cp(label,{"stage":"FAILED","lane":lane,"slot":slot+1,"target":g["target"],"day":g["day"],"error":msg(e),"paired_gate":False})
        # No power-based replacements. Quarantine any partial outputs; failed lane never passes.
    finally:n.ref.tiles=original
    require(c.sha((O/"lock.json").read_bytes())==S["source"]["lock_sha256"],"frozen weights changed")
def summarize():
    plan=a.verified_plan();m=load("manifest.json");results={}
    for lane in ("same_receiver","uwl"):
        recs=[json.loads((PROG/(lane+"-"+str(i+1).zfill(2)+".json")).read_text()) for i in range(4)]
        good=[d for d in recs if d["stage"]=="MEASURED"];bad=[d for d in recs if d["stage"]=="FAILED"]
        expected=sum(g["receiver_lane"]==lane for g in plan["fresh"])
        out={"status":"HOLD_INSUFFICIENT_QUALIFIED_GROUPS" if expected<2 else "HOLD_ACQUISITION_OR_PARITY_FAILURE" if bad or len(good)!=expected else "FROZEN_TRANSFER_MEASURED","qualified":m["lanes"][lane]["qualified"],"planned":expected,"measured":len(good),"failures":bad,"controls":128*len(good),"groups":good}
        if good:
            def total(kind):
                rows=[r["metrics"] for g in good for r in g["rows"] if r["kind"]==kind]
                return {a0:{k:sum(x[a0][k] for x in rows) for k in rows[0][a0]} for a0 in rows[0]}
            ev=total("eval");sh=total("shared");nc=128*len(good);base=ev["learned05"]["control_alerts"]
            gates={k:out["status"]=="FROZEN_TRANSFER_MEASURED" and ev[k]["lost_parent"]==ev[k]["lost_joint"]==sh[k]["lost_parent"]==sh[k]["lost_joint"]==0 and ev[k]["control_alerts"]*640<=35*nc and ev[k]["control_alerts"]*8<=7*base and ev[k]["repeat_stress_lost"]==sh[k]["repeat_stress_lost"]==0 and ev[k]["repeat_stress_cases"]+sh[k]["repeat_stress_cases"]>0 for k in ev}
            fast=sum(g["turbo"]["paired_turbo_seconds"] for g in good);ref=sum(g["turbo"]["paired_reference_seconds"] for g in good)
            out.update(eval=ev,shared=sh,gates=gates,paired_speedup=ref/fast if fast else None)
        results[lane]=out
    require(c.sha((O/"lock.json").read_bytes())==S["source"]["lock_sha256"],"frozen weights changed")
    result={"status":"FROZEN_TWO_LANE_REPORT","source_n10_lock_sha256":S["source"]["lock_sha256"],"manifest_sha256":plan["manifest_sha256"],"no_training":True,"production_adoption":False,"discovery_relaunch_allowed":False,"lanes":results}
    save("result.json",result);emit("R40N10B_SUMMARY",{**result,"lanes":{k:{x:y for x,y in v.items() if x!="groups"} for k,v in results.items()}})
    cp("99-summary",{"stage":"COMPLETE","lanes":{k:{"status":v["status"],"measured":v["measured"],"passed_arms":[x for x,y in v.get("gates",{}).items() if y]} for k,v in results.items()},"lock_sha256":S["source"]["lock_sha256"]})
if __name__=="__main__":
    if ARGS[0]=="prepare":prepare()
    elif ARGS[0]=="select":select(ARGS[1])
    elif ARGS[0]=="freeze":freeze()
    elif ARGS[0]=="group":run_group(ARGS[1],int(ARGS[2]))
    elif ARGS[0]=="summary":summarize()
    else:raise ValueError(ARGS)
