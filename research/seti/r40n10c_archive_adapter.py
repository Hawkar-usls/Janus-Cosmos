#!/usr/bin/env python3
"""R40N10C: reproducible archive adapter audit, no model fitting."""
import sys
ARGS=sys.argv[1:]
import collections,hashlib,json,math,pathlib,re,time,traceback,urllib.parse
import requests
from bs4 import BeautifulSoup
import r40n10b_expansion as prev
a=prev.a;n=a.n;c=a.c;qual=prev.qual;eng=prev.eng;scout=prev.scout
P=pathlib.Path("research/seti/r40n10c_archive_adapter_protocol.json");S=json.loads(P.read_text());C=S["selection"]
O=pathlib.Path("runtime/r40n10c");O.mkdir(parents=True,exist_ok=True);CACHE=O/"cache";CACHE.mkdir(exist_ok=True);PROG=O/"progress";PROG.mkdir(exist_ok=True)
for m in (prev,a):
    m.P=P;m.S=S;m.C=C;m.O=O;m.CACHE=CACHE;m.PROG=PROG
n.O=O;n.k.O=O;n.l.O=O;eng.OUT=O
def emit(k,d):print(k.replace("R40N10B_","R40N10C_").replace("R40N10_","R40N10C_"),json.dumps(d,allow_nan=False),flush=True)
prev.emit=emit;a.emit=emit
save=prev.save;load=prev.load;cp=prev.cp;require=prev.require
LANES=("same_receiver","uwl","gbt_s")
def error(e):
    tb=traceback.extract_tb(e.__traceback__)
    return {"type":type(e).__name__,"message":str(e)[:300],"frames":[{"file":pathlib.Path(x.filename).name,"function":x.name,"line":x.lineno,"code":(x.line or "")[:240]} for x in tb[-5:]]}
def header(f,cache):
    if f["url"] not in cache:cache[f["url"]]=eng.hdr(f)
    return cache[f["url"]]
def raw_parent():
    z,nn=c.artifact(S["parent"],128)
    for k,v in S["parent"]["manifests"].items():
        require(c.sha(z.read(nn[k+"_manifest.json"]))==v,"N10B manifest digest "+k)
    return z,nn
def prepare():
    prev.prepare()
    z,nn=raw_parent()
    for k in ("same_receiver","uwl"):
        (O/("old_"+k+".json")).write_bytes(z.read(nn[k+"_manifest.json"]))
    old=load("old_same_receiver.json")
    require(sum(e["error"]=="AssertionError: " for q in old["attempts"] for e in q["errors"])==80,"expected 80 old assertion records")
    z0,nn0,p=prev.old_source()
    g=p["development"][0]["group"];cache={}
    hd,pp=prev.headers_same(g["target"],g["files"],[],load("exclusions.json")["reference"],cache)
    require(hd==[f["header"] for f in g["files"]],"known N10 header parity")
    cp("01-positive-parkes",{"stage":"KNOWN_PARKES_CONTROL","target":g["target"],"day":g["day"],"headers":hd,"pointing":pp,"pass":True,"eligible_for_fresh":False})
def summarize_error(e):
    f=e["frames"][-1] if e["frames"] else {}
    return "|".join([e["type"],e["message"] or f.get("code",""),f.get("file",""),str(f.get("line",""))])
def parkes():
    t0=time.perf_counter();old=load("old_same_receiver.json");ctx=load("exclusions.json");names,days,ns,urls,coords=prev.exclude(ctx,[])
    manifest={"lane":"same_receiver","selected":[],"attempts":[],"old_assertions":80,"old_manifest_sha256":S["parent"]["manifests"]["same_receiver"]}
    cache={};resolved=0
    for rec in old["attempts"]:
        if not rec.get("cadences"):continue
        target=rec["target"];rr,qe=prev.query(target)
        # Reconstruct original N10B candidates against the unchanged historical exclusions,
        # before adding any newly chosen target/day. Match old records by day AND ordinal.
        fs=qual.entries(rr,"Parkes",target,set(ctx["days"]),set(ctx["files"]),set(ctx["urls"]))
        cand=qual.cadence(fs)[:C["header_attempts"]]
        item={"target":target,"old_cadences":rec["cadences"],"current_cadences":len(cand),"query_errors":qe,"records":[]}
        for j,old_err in enumerate(rec["errors"]):
            e={"ordinal":j,"old":old_err,"matched":False}
            if len(cand)!=rec["cadences"] or j>=len(cand) or cand[j][0]["day"]!=old_err["day"]:
                e["status"]="UNRESOLVED_METADATA_CHANGED"
                item["records"].append(e);continue
            g=cand[j];e.update(matched=True,files=[f["name"] for f in g],day=g[0]["day"])
            try:
                hs,pp=prev.headers_same(target,g,ctx["coords"],ctx["reference"],cache)
                e["status"]="NOW_QUALIFIED_METADATA"
                if target not in names and g[0]["day"] not in days and len(manifest["selected"])<4 and all(qual.history_distance(pp[0],v)>=qual.A["pointing_separation_degrees"] for v in coords):
                    files=[{**f,"header":hd,"role":role} for f,hd,role in zip(g,hs,prev.ROLES)]
                    rec0={"telescope":"Parkes","receiver_lane":"same_receiver","target":target,"day":g[0]["day"],"pointing":pp[0],"files":files,"split":"holdout"}
                    manifest["selected"].append(rec0);names.add(target);days.add(g[0]["day"]);coords.append({"target":target,"pointing":pp[0]})
            except Exception as ex:
                e.update(status="REJECTED_WITH_TRACE",error=error(ex))
            # Diagnose specific measured header/metadata values without reading power.
            e["headers"]=[{"name":f["name"],"header":cache.get(f["url"]),"metadata":{k:f["metadata"].get(k) for k in ("target","mjd","ra","dec","decl","declination","center_freq")}} for f in g]
            if old_err["error"]=="AssertionError: " and e["matched"]:resolved+=1
            item["records"].append(e)
        manifest["attempts"].append(item);save("same_receiver_manifest.json",manifest)
        emit("R40N10C_PARKES_AUDIT",{"target":target,"resolved_old_assertions":resolved,"total_old_assertions":80,"qualified":len(manifest["selected"])})
    manifest["eligible"]=len(manifest["selected"])>=2
    rr=[e for q in manifest["attempts"] for e in q["records"] if e["old"]["error"]=="AssertionError: "]
    require(len(rr)==80,"all original assertion records retained")
    manifest["diagnosis"]={"matched":sum(e["matched"] for e in rr),"unresolved":sum(not e["matched"] for e in rr),"reasons":dict(collections.Counter(summarize_error(e["error"]) if "error" in e else e["status"] for e in rr))}
    save("same_receiver_manifest.json",manifest)
    cp("02-parkes",{"stage":"PARKES_AUDIT","diagnosis":manifest["diagnosis"],"qualified":len(manifest["selected"]),"seconds":time.perf_counter()-t0})
def uwl_control():
    # A known public dataset, never included in new science. Bounded directory traversal.
    roots=C["uwl_control_dirs"];queue=list(roots);visited=set();files=[];errors=[];seen=set()
    while queue and len(visited)<24 and len(files)<48:
        u=queue.pop(0)
        if u in visited:continue
        visited.add(u)
        try:
            with requests.get(u,stream=True,timeout=(15,45)) as p:
                p.raise_for_status()
                require(urllib.parse.urlparse(p.url).hostname=="bldata.berkeley.edu" and any(p.url.startswith(x) for x in roots),"control redirect outside known public root")
                buf=bytearray()
                for x in p.iter_content(65536):
                    buf.extend(x);require(len(buf)<=4*1024**2,"directory listing budget")
            soup=BeautifulSoup(bytes(buf),"html.parser")
            for el in soup.find_all("a",href=True):
                url=urllib.parse.urljoin(u,el["href"]);v=urllib.parse.urlparse(url)
                if v.hostname!="bldata.berkeley.edu" or not any(url.startswith(x) for x in roots) or v.query or v.fragment:continue
                if url.endswith("/") and url not in visited:queue.append(url)
                elif url.lower().endswith(".h5") and url not in seen:files.append(url);seen.add(url)
        except Exception as e:errors.append({"url":u,"error":error(e)})
    rows=[]
    for u in sorted(files)[:16]:
        try:
            h=eng.hdr({"url":u});rows.append({"url":u,"header":h,"high_resolution":0<abs(h["df"])*1e6<=10,"old_name_filter":bool(re.search(r"(?:fine|rawspec\.0000|seti)",u,re.I))})
        except Exception as e:rows.append({"url":u,"error":error(e)})
    d={"stage":"UWL_KNOWN_ARCHIVE_CONTROL","directory_count":len(visited),"listed_hdf5":len(files),"headers":rows,"errors":errors,"readable_highres_headers":sum(x.get("high_resolution",False) for x in rows),"fresh_eligible":False}
    save("uwl_control.json",d);cp("03-uwl-control",d)
def role(x,target):
    # Separators vary in SIGPROC source_name and catalogue metadata.
    # Exact base ID is still mandatory; no fuzzy name or position matching.
    m=re.fullmatch(r"(.+?)[ _-]+([SR])",str(x).strip(),re.I)
    if not m:return None
    canon=lambda z:re.sub(r"[ _-]+","",str(z)).upper()
    return m[2].upper() if canon(m[1])==canon(target) else None
def uwl():
    t0=time.perf_counter();ctx=load("exclusions.json");prev_groups=load("same_receiver_manifest.json")["selected"]
    names,days,ns,urls,coords=prev.exclude(ctx,prev_groups)
    order=list(dict.fromkeys(S["published_uwl"]["targets"]+C["uwl_aliases"]))
    require(len(order)<=38,"bounded candidate list")
    manifest={"lane":"uwl","selected":[],"attempts":[]};cache={}
    for target in order:
        if len(manifest["selected"])>=4:break
        canon=re.sub(r"[ _-]+","",target).upper()
        if canon in names:continue
        rr,qe=prev.query(target)
        allf={f["url"]:f for q in rr for f in q["files"]}
        count=collections.Counter();scans={};sample=[];errs=[]
        # Record API schema and prefilter effects; do not assume a name or center-frequency field.
        keys=collections.Counter(k for f in allf.values() for k in f["metadata"])
        candidates=[]
        for f in sorted(allf.values(),key=lambda f:(f["day"],f["mjd"],f["name"],f["url"])):
            count["hdf5"]+=1
            if f["day"] in days or f["name"] in ns or f["url"] in urls:count["excluded_prior"]+=1;continue
            mt=f["metadata"];pos=role(mt.get("target",""),target)
            old=bool(re.search(r"(?:fine|rawspec\.0000|seti)",f["name"],re.I))
            count["old_name_filter"]+=int(old)
            count["old_center_freq_present"]+=int(mt.get("center_freq") is not None)
            count["exact_role"]+=int(pos is not None)
            if len(sample)<6:sample.append({"name":f["name"],"metadata":mt,"old_name_filter":old,"exact_role":pos})
            if pos is None:continue
            # Metadata center frequency is only a deterministic priority hint, never a rejection.
            try:
                fc=float(mt.get("center_freq"));fc=fc/1e6 if fc>1e6 else fc
                dist=min(abs(fc-w) for w in C["windows_mhz"])
                if not math.isfinite(dist):dist=1e9
            except (ValueError,TypeError):dist=1e9
            candidates.append((dist,f["mjd"],f["name"],f,pos))
        for _,_,_,f,pos in sorted(candidates,key=lambda z:z[:3])[:C["max_header_per_uwl_target"]]:
            count["header_attempted"]+=1
            try:
                hd=header(f,cache);count["header_read"]+=1
                require(role(hd["source"],target)==pos,"header and catalogue exact role disagree")
                require(len(hd["shape"])==3 and 16<=hd["shape"][0]<=256 and hd["shape"][1]==1,"UWL native single-intensity rows")
                require(0<abs(hd["df"])*1e6<=10 and 0<hd["ts"]<=60,"UWL native time/frequency resolution")
                require(hd["bounds"][1]-hd["bounds"][0]<=150,"not a UWL subband")
                require(abs(hd["mjd"]-f["mjd"])*86400<=5,"header and index time mismatch")
                require(int(hd["mjd"])==f["day"],"header and index day mismatch")
                pt=qual.pointing(hd,f["metadata"])
                found=[w for w in C["windows_mhz"] if prev.contains(hd,w)]
                require(len(found)==1,"expected one native window per subband")
                k=(f["day"],hd["mjd"],pos);scans.setdefault(k,{}).setdefault(found[0],[]).append({**f,"header":hd,"position":pos,"pointing":pt})
                count["usable_subband"]+=1
            except Exception as e:errs.append({"name":f["name"],"error":error(e)})
        series=[]
        for (day,mjd,pos),bywin in sorted(scans.items()):
            if all(w in bywin for w in C["windows_mhz"]):
                series.append({"day":day,"mjd":mjd,"position":pos,"kind":"ON" if pos=="S" else "OFF","candidates":[sorted(bywin[w],key=lambda f:(f["name"],f["url"]))[:3] for w in C["windows_mhz"]]})
        # No heuristic source-name changes: give the previous qualifier equivalent canonical
        # source spelling only AFTER exact role agreement, retain raw header in the ledger.
        oldhdr=eng.hdr
        try:
            def canonical_hdr(f):
                hd=header(f,cache);pos=role(hd["source"],target)
                return {**hd,"source":prev.key(target)+"_"+pos} if pos is not None else hd
            eng.hdr=canonical_hdr
            cand=prev.uwl_cadences(series);count["complete_scans"]=len(series);count["cadences"]=len(cand)
            for g in cand[:C["header_attempts"]]:
                try:
                    fs,pp=prev.qualify_uwl(prev.key(target),g,coords,{})
                    require(g[0]["day"] not in days,"prior day")
                    # Preserve original source_name in frozen manifest: no invented header metadata.
                    for f in fs:
                        f["header"]=cache[f["url"]]
                        for sub in f["subbands"]:sub["header"]=cache[sub["url"]]
                    group={"telescope":"Parkes","receiver_lane":"uwl","target":canon,"day":g[0]["day"],"split":"holdout","pointing":pp[0],"files":fs}
                    manifest["selected"].append(group);names.add(canon);days.add(g[0]["day"]);coords.append({"target":canon,"pointing":pp[0]})
                    for f in fs:
                        for sub in f["subbands"]:ns.add(sub["name"]);urls.add(sub["url"])
                    break
                except Exception as e:errs.append({"day":g[0]["day"],"error":error(e)})
        finally:eng.hdr=oldhdr
        manifest["attempts"].append({"target":target,"counts":dict(count),"metadata_keys":dict(keys),"metadata_examples":sample,"query_errors":qe,"rejections":errs,"index_sha256":[r["index_sha256"] for r in rr]})
        save("uwl_manifest.json",manifest)
        emit("R40N10C_UWL",{"target":target,"counts":dict(count),"qualified":len(manifest["selected"])})
    manifest["eligible"]=len(manifest["selected"])>=2
    totals=collections.Counter()
    for x in manifest["attempts"]:totals.update(x["counts"])
    manifest["totals"]=dict(totals);save("uwl_manifest.json",manifest)
    cp("04-uwl",{"stage":"UWL_ADAPTER_AUDIT","targets":len(manifest["attempts"]),"totals":dict(totals),"qualified":len(manifest["selected"]),"seconds":time.perf_counter()-t0})
def gbt():
    t0=time.perf_counter();ctx=load("exclusions.json");prev_groups=load("same_receiver_manifest.json")["selected"]+load("uwl_manifest.json")["selected"]
    names,days,ns,urls,coords=prev.exclude(ctx,prev_groups)
    strings=set(scout.walk(json.loads((O/"target-index.json").read_bytes())))
    pool={x[:-4].upper() for x in strings if re.fullmatch(r"[A-Za-z0-9_-]{2,48}_OFF",x,re.I)}
    order=sorted(pool-names,key=lambda x:c.sha((C["seed"]+"|GBT-S|"+x).encode()))[:C["max_gbt_targets"]]
    qual.S["gbt"]={**qual.S["gbt"],"windows_mhz":C["gbt_windows_mhz"],"half_hz":C["half_hz"]}
    manifest={"lane":"gbt_s","selected":[],"candidate_order":order,"attempts":[]}
    for target in order:
        if len(manifest["selected"])>=4:break
        if target in names or "ATLAS" in target:continue
        rec={"target":target,"errors":[]}
        try:
            rr=[a.index_query("GBT",target,q) for q in (target,target+"_OFF")]
            fs=qual.entries(rr,"GBT",target,days,ns,urls);cand=qual.cadence(fs)
            rec.update(files=len(fs),cadences=len(cand),index_sha256=[x["index_sha256"] for x in rr])
            for g in cand[:C["header_attempts"]]:
                try:
                    hs,pp=qual.qualify("GBT",target,g,coords)
                    require(all(x["shape"][1]==1 for x in hs),"GBT single intensity plane")
                    group={"telescope":"GBT","receiver_lane":"gbt_s","target":target,"day":g[0]["day"],"pointing":pp[0],"split":"holdout","files":[{**f,"header":hd,"role":role0} for f,hd,role0 in zip(g,hs,prev.ROLES)]}
                    manifest["selected"].append(group);names.add(target);days.add(group["day"]);coords.append({"target":target,"pointing":pp[0]});ns.update(f["name"] for f in g);urls.update(f["url"] for f in g);break
                except Exception as e:rec["errors"].append({"day":g[0]["day"],"error":error(e)})
        except Exception as e:rec["errors"].append({"error":error(e)})
        manifest["attempts"].append(rec);save("gbt_s_manifest.json",manifest)
        if len(manifest["attempts"])%10==0:emit("R40N10C_GBT",{"targets":len(manifest["attempts"]),"qualified":len(manifest["selected"])})
    manifest["eligible"]=len(manifest["selected"])>=2;save("gbt_s_manifest.json",manifest)
    cp("05-gbt",{"stage":"GBT_S_METADATA_TRANSFER","targets":len(manifest["attempts"]),"cadences":sum(x.get("cadences",0) for x in manifest["attempts"]),"qualified":len(manifest["selected"]),"reasons":dict(collections.Counter(summarize_error(e["error"]) for x in manifest["attempts"] for e in x["errors"])),"seconds":time.perf_counter()-t0})
def freeze():
    ctx=load("exclusions.json");lanes={k:load(k+"_manifest.json") for k in LANES};allg=sum([m["selected"] for m in lanes.values()],[])
    require(len({g["day"] for g in allg})==len(allg) and len({g["target"] for g in allg})==len(allg),"overlap across lanes")
    fresh=sum([m["selected"] if m["eligible"] else [] for m in lanes.values()],[])
    m={"status":"FROZEN_METADATA_ONLY","source_n10_lock_sha256":S["source"]["lock_sha256"],"lanes":{k:{"manifest_sha256":c.sha((O/(k+"_manifest.json")).read_bytes()),"qualified":len(v["selected"]),"eligible":v["eligible"],"targets":[{"target":g["target"],"day":g["day"]} for g in v["selected"]]} for k,v in lanes.items()},"selected":fresh}
    save("manifest.json",m)
    plan={"protocol_sha256":c.sha(P.read_bytes()),"manifest_sha256":c.sha((O/"manifest.json").read_bytes()),"development":[],"fresh":fresh,"excluded_raw_hashes":ctx["excluded_raw_hashes"],"fresh_allowed":bool(fresh)}
    save("plan.json",plan)
    require(c.sha((O/"lock.json").read_bytes())==S["source"]["lock_sha256"],"frozen model")
    save("freeze_receipt.json",{"lock_sha256":S["source"]["lock_sha256"],"plan_sha256":c.sha((O/"plan.json").read_bytes())})
    cp("06-freeze",{"stage":"FROZEN","lanes":m["lanes"],"total_groups":len(fresh),"lock_sha256":S["source"]["lock_sha256"],"manifest_sha256":plan["manifest_sha256"],"no_training":True})
def run_group(lane,slot):
    if lane=="gbt_s":a.C={**C,"windows_mhz":C["gbt_windows_mhz"]}
    prev.run_group(lane,slot)
def summary():
    plan=a.verified_plan();m=load("manifest.json");results={}
    for lane in LANES:
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
            out.update(eval=ev,shared=sh,gates=gates)
        results[lane]=out
    require(c.sha((O/"lock.json").read_bytes())==S["source"]["lock_sha256"],"model changed")
    result={"status":"ARCHIVE_ADAPTER_AUDIT_FROZEN_TRANSFER","source_n10_lock_sha256":S["source"]["lock_sha256"],"manifest_sha256":plan["manifest_sha256"],"no_training":True,"production_adoption":False,"discovery_relaunch_allowed":False,"lanes":results,"parkes_diagnosis":load("same_receiver_manifest.json")["diagnosis"],"uwl_totals":load("uwl_manifest.json")["totals"],"uwl_positive_control_readable":load("uwl_control.json")["readable_highres_headers"]}
    save("result.json",result);emit("R40N10C_SUMMARY",{**result,"lanes":{k:{x:y for x,y in v.items() if x!="groups"} for k,v in results.items()}})
    cp("99-summary",{"stage":"COMPLETE","parkes_diagnosis":result["parkes_diagnosis"],"uwl_totals":result["uwl_totals"],"uwl_positive_control_readable":result["uwl_positive_control_readable"],"lanes":{k:{"status":v["status"],"measured":v["measured"],"passed_arms":[x for x,y in v.get("gates",{}).items() if y]} for k,v in results.items()},"lock_sha256":S["source"]["lock_sha256"]})
if __name__=="__main__":
    actions={"prepare":prepare,"parkes":parkes,"uwl-control":uwl_control,"uwl":uwl,"gbt":gbt,"freeze":freeze,"summary":summary}
    if ARGS[0]=="group":run_group(ARGS[1],int(ARGS[2]))
    else:actions[ARGS[0]]()
