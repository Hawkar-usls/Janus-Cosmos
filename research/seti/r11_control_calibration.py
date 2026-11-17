#!/usr/bin/env python3
"""R11 prospective control occupancy, reverse cadence and full-window injections."""
import collections,hashlib,json,math,pathlib,re,sys
import numpy as np
import r9_expanded_search as engine
import r4_chunk_replay as base

PATH=pathlib.Path("research/seti/r11_preregistered_search.json")
S=json.loads(PATH.read_text())
OLD=json.loads(pathlib.Path("research/seti/r9_preregistered_search.json").read_text())
for k in ("drift_range_hz_s","anchor_min_score","anchor_min_row_excess_fraction","top_candidates_per_tile","off_and_repeat_threshold_score","frequency_guard_hz"):
    assert S[k]==OLD[k],k
OUT=pathlib.Path("runtime/r11");OUT.mkdir(parents=True,exist_ok=True)
engine.SPEC=S
engine.OUT=OUT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    (OUT/name).write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n")
def emit(k,x):print(k,json.dumps(x,ensure_ascii=False),flush=True)

def names(obj):
    if isinstance(obj,str):
        if re.fullmatch(r"HIP\d+_S",obj,re.I):yield obj.upper()[:-2]
    elif isinstance(obj,list):
        for v in obj:yield from names(v)
    elif isinstance(obj,dict):
        for v in obj.values():yield from names(v)

def freeze():
    manifest={"schema":"janus.seti.r11.manifest.v1","spec_sha256":sha(PATH),"selected":[],"attempts":[],"metadata_only":True}
    try:
        raw,url=base.fetch(S["selection"]["index"],limit=32*base.MIB)
        (OUT/"target-index.json").write_bytes(raw)
        available=set(names(json.loads(raw)))-set(S["exclude_prior_targets"])
        ordered=sorted(available,key=lambda x:hashlib.sha256(("R11-20261117|"+x).encode()).hexdigest())
        manifest["target_index"]={"url":url,"sha256":hashlib.sha256(raw).hexdigest(),"eligible_name_count":len(ordered),"candidate_order":ordered[:32]}
        for target in ordered[:S["limits"]["max_index_queries"]]:
            attempt={"target":target,"cadence_attempts":[]};manifest["attempts"].append(attempt)
            try:
                fs,q=engine.query(target);attempt.update(q)
                fs.sort(key=lambda f:(f["day"],f["seconds"],f["url"]))
                cad=[]
                for i in range(len(fs)-3):
                    a=fs[i:i+4]
                    if [x["position"] for x in a]!=["S","R","S","R"] or len({x["day"] for x in a})!=1:continue
                    if not all(0<y["seconds"]-x["seconds"]<=600 for x,y in zip(a,a[1:])):continue
                    cad.append(a)
                cad.sort(key=lambda a:(a[0]["day"]<58575,a[0]["day"],a[0]["seconds"],a[0]["url"]))
                for a in cad[:5]:
                    try:
                        hs=[engine.hdr(f) for f in a]
                        if not all(engine.valid(h,f,target) for h,f in zip(hs,a)):raise ValueError("header ineligible")
                        if len({(h["shape"][0],round(h["ts"],7),round(h["df"],12)) for h in hs})!=1:raise ValueError("resolution mismatch")
                        coords=[(h["src_raj"],h["src_dej"]) for h in hs]
                        if coords[0]!=coords[2] or coords[1]!=coords[3] or coords[0]==coords[1]:raise ValueError("pointing mismatch")
                        manifest["selected"].append({"target":target,"cohort":"post_2019_march" if a[0]["day"]>=58575 else "archival_reanalysis","files":[{**f,"header":h,"role":role} for f,h,role in zip(a,hs,("ON1","OFF1","ON2","OFF2"))]})
                        attempt["cadence_attempts"].append({"day":a[0]["day"],"seconds":a[0]["seconds"],"status":"SELECTED"});break
                    except Exception as e:attempt["cadence_attempts"].append({"day":a[0]["day"],"seconds":a[0]["seconds"],"error":str(e)})
            except Exception as e:attempt["error"]=str(e)
            if len(manifest["selected"])==S["limits"]["targets"]:break
    except Exception as e:manifest["index_error"]=str(e)
    save("manifest.json",manifest)
    emit("R11_MANIFEST",{"target_index":manifest.get("target_index"),"selected":[{"target":x["target"],"cohort":x["cohort"],"files":[{"name":f["name"],"role":f["role"],"header":f["header"]} for f in x["files"]]} for x in manifest["selected"]],"attempts":manifest["attempts"],"index_error":manifest.get("index_error")})

def top(y,t):
    return engine.top(y,t["df_hz"],t["ts"],S["top_candidates_per_tile"],S["anchor_min_score"])

def track_match(y,t,f,dr):
    return engine.matched(y,t,f,dr)

def candidate(h,observations,ts,ys,roles,target,fc):
    a,b,c=ts; y=ys[0]
    f=a["channel0_mhz"]+h["bin"]*a["df_hz"]/1e6;dr=h["drift_hz_s"]
    item={"target":target,"center_mhz":fc,"frequency_mhz":f,**h,"cadence":roles,"comparisons":[],"status":"HOLD_NO_DIRECTIONAL_REPEAT"}
    ds=np.rint(np.arange(len(y))*h["shift"]/(len(y)-1)).astype(int)
    amp=float(np.median(y[np.arange(len(y)),h["bin"]+ds]))
    for o,t,z in zip(observations[1:],ts[1:],ys[1:]):
        dt=(o["header"]["mjd"]-observations[0]["header"]["mjd"])*86400
        pred=f+dr*dt/1e6
        cv=track_match(z,t,pred,dr);fv=track_match(z,t,f,dr)
        injected=[]
        for offset in (-70,70):
            ft=pred+offset*abs(t["df_hz"])/1e6;k=int(round((ft-t["channel0_mhz"])*1e6/t["df_hz"]))
            s=int(round(dr*(len(z)-1)*t["ts"]/t["df_hz"]))
            dsi=np.rint(np.arange(len(z))*s/(len(z)-1)).astype(int)
            if k+dsi.min()<0 or k+dsi.max()>=z.shape[1]:continue
            zz=z.copy();zz[np.arange(len(z)),k+dsi]+=amp
            injected.append(engine.present(track_match(zz,t,ft,dr)))
        item["comparisons"].append({"role":o["role"],"dt_seconds":dt,"continuous":cv,"fixed":fv,"injections_passed":sum(injected),"injections":len(injected)})
    off,rep=item["comparisons"]
    if engine.present(off["continuous"]) or engine.present(off["fixed"]):item["status"]="CONTROL_TRACK_PRESENT_RFI_CONCERN"
    elif off["continuous"] is not None and off["fixed"] is not None and off["injections"]==2 and off["injections_passed"]==2 and engine.present(rep["continuous"]):
        item["status"]="PRELIMINARY_DIRECTIONAL_HOLD"
    return item

def calibration(x,t,baseline,target,role,fc):
    n,m=x.shape;df=t["df_hz"];ts=t["ts"]
    seed=int.from_bytes(hashlib.sha256(("R11-20261117|"+target+"|"+role+"|"+str(fc)).encode()).digest()[:8],"big")
    rng=np.random.default_rng(seed);lim=int(math.floor(4*(n-1)*ts/abs(df)))
    shift=int(rng.integers(-lim,lim+1));ds=np.rint(np.arange(n)*shift/(n-1)).astype(int)
    lo=max(32,32-int(ds.min()));hi=min(m-32,m-32-int(ds.max()))
    if hi<=lo:raise ValueError("no legal injection")
    c=int(rng.integers(lo,hi))
    def hit(hs):
        return next((h for h in hs if abs(h["bin"]-c)<=2 and abs(h["bin"]+h["shift"]-c-shift)<=2),None)
    baseline_coincidence=bool(hit(baseline))
    scale=1.4826*np.median(np.abs(x-np.median(x,axis=1)[:,None]),axis=1)
    result={"target":target,"role":role,"center_mhz":fc,"shape":list(x.shape),"seed":seed,"start_bin":c,"shift":shift,"drift_hz_s":shift*df/((n-1)*ts),"baseline_candidates":len(baseline),"baseline_coincidence":baseline_coincidence,"levels":[]}
    for amp in S["calibration"]["amplitudes"]:
        z=x.copy();z[np.arange(n),c+ds]+=amp*scale
        hs,_=top(engine.normalized(z),t);h=hit(hs)
        result["levels"].append({"amplitude":amp,"recovered":bool(h),"new_recovered":bool(h) and not baseline_coincidence,"score":h["score"] if h else None,"candidate_count":len(hs)})
    emit("R11_CALIBRATION",result)
    return result

def summarize(manifest,results,tiles,cal,forward,reverse):
    off=[r for t in tiles for r in t["scans"] if r["role"].startswith("OFF")]
    n=len(off);cnt=sum(r["candidate_count"] for r in off)
    on=[r for t in tiles for r in t["scans"] if r["role"]=="ON1"]
    levels=[]
    for amp in S["calibration"]["amplitudes"]:
        a=[l for t in cal for l in t["levels"] if l["amplitude"]==amp]
        k=sum(l["recovered"] for l in a);new=sum(l["new_recovered"] for l in a)
        levels.append({"amplitude":amp,"recovered":k,"new_recovered":new,"trials":len(a),"fraction":k/len(a) if a else None})
    return {"schema":"janus.seti.r11.result.v1","manifest_sha256":sha(OUT/"manifest.json"),"selected":len(manifest["selected"]),"completed":sum(r["status"]=="COMPLETE" for r in results),"completed_tiles":len(tiles),"on_anchor_candidate_count":len(forward),"on_anchor_tiles":len(on),"off_tiles":n,"off_tiles_with_alert":sum(r["candidate_count"]>0 for r in off),"off_returned_alerts":cnt,"off_tiles_saturated_at_top5":sum(r["candidate_count"]==5 for r in off),"off_template_upper_bound":sum(r["counts"]["upper_bound_templates"] for r in off),"off_alerts_per_mhz_scan":cnt/(n*.1) if n else None,"forward_status_counts":dict(collections.Counter(h["status"] for h in forward)),"reverse_status_counts":dict(collections.Counter(h["status"] for h in reverse)),"reverse_anchor_tiles":len(tiles),"calibration":{"paired_tracks":len(cal),"baseline_coincidences":sum(x["baseline_coincidence"] for x in cal),"levels":levels},"confirmed_technosignatures":0,"independently_confirmed":0,"limitations":["OFF alert incidence and reverse cadence patterns are empirical interference proxies, not a known-null false-alarm probability.","Top-five cap and per-drift peak cap censor counts. Repeated windows/scans and paired injections are correlated.","Full 100 kHz windows, not a full receiver band. One discrete channel constant-power track per sampled OFF tile; no absolute Jy/EIRP or beam calibration.","Archival controls, not new independent observations. Any directional screen pass remains HOLD."]}

def replay():
    engine.selftest()
    manifest=json.loads((OUT/"manifest.json").read_text())
    assert manifest["spec_sha256"]==sha(PATH)
    results=[];tiles=[];cal=[];forward=[];reverse=[]
    for q in manifest["selected"]:
        r={"target":q["target"],"status":"HOLD","observations":[]};results.append(r)
        try:
            for f in q["files"]:r["observations"].append({**f,"tiles":engine.get_tiles(f)})
            obs=r["observations"]
            for i,fc in enumerate(S["window_centers_mhz"]):
                ts=[o["tiles"][i] for o in obs];xs=[]
                for t in ts:
                    x=np.load(OUT/t["artifact"],allow_pickle=False)
                    assert list(x.shape)==t["shape"] and hashlib.sha256(x.tobytes()).hexdigest()==t["slice_sha256"]
                    assert x.shape[1]*abs(t["df_hz"])>=100000
                    xs.append(x)
                ys=[engine.normalized(x) for x in xs]
                hits=[];records=[]
                for o,t,y in zip(obs,ts,ys):
                    hs,counts=top(y,t);hits.append(hs);records.append({"role":o["role"],"candidate_count":len(hs),"counts":counts,"candidates":hs})
                fwd=[candidate(h,obs[:3],ts[:3],ys[:3],"ON_OFF_ON",q["target"],fc) for h in hits[0]]
                rev=[candidate(h,obs[1:],ts[1:],ys[1:],"OFF_ON_OFF",q["target"],fc) for h in hits[1]]
                cc=[calibration(xs[j],ts[j],hits[j],q["target"],obs[j]["role"],fc) for j in (1,3)]
                tile={"target":q["target"],"center_mhz":fc,"scans":records}
                tiles.append(tile);forward.extend(fwd);reverse.extend(rev);cal.extend(cc)
                for h in fwd:emit("R11_FORWARD",h)
                for h in rev:emit("R11_REVERSE",h)
                emit("R11_TILE",{"target":q["target"],"center_mhz":fc,"counts":{a["role"]:a["candidate_count"] for a in records},"off_alerts":records[1]["candidate_count"]+records[3]["candidate_count"]})
                save("partial.json",{"results":results,"tiles":tiles,"calibration_trials":cal,"forward":forward,"reverse":reverse})
            r["status"]="COMPLETE"
        except Exception as e:r["error"]=str(e);emit("R11_HOLD",{"target":q["target"],"error":str(e)})
    s=summarize(manifest,results,tiles,cal,forward,reverse)
    save("result.json",{**s,"results":results,"tiles":tiles,"calibration_trials":cal,"forward":forward,"reverse":reverse})
    emit("R11_SUMMARY",s)

if __name__=="__main__":
    {"freeze":freeze,"scan":replay}[sys.argv[1]]()
