#!/usr/bin/env python3
"""Qualify source candidates metadata-only; independent aliases remain review-required."""
import hashlib,io,itertools,json,math,os,pathlib,re,zipfile
import requests
import r9_expanded_search as eng
P=pathlib.Path("research/seti/r37b_source_header_protocol.json");S=json.loads(P.read_text())
H=json.loads(pathlib.Path(S["history"]).read_text());A=S["selection"]
O=pathlib.Path("runtime/r37b");O.mkdir(parents=True,exist_ok=True)
def sha(x):return hashlib.sha256(x).hexdigest()
def save(n,v):(O/n).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def artifact(i,d):
    r=requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{i}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=90);r.raise_for_status()
    b=r.content;assert len(b)<=128*1024**2 and sha(b)==d
    z=zipfile.ZipFile(io.BytesIO(b));assert sum(i.file_size for i in z.infolist())<=512*1024**2
    names={pathlib.PurePosixPath(i.filename).name:i.filename for i in z.infolist() if not i.is_dir()}
    assert len(names)==len([i for i in z.infolist() if not i.is_dir()])
    return z,names
def angle(v,ra):
    s=str(v).strip();sign=-1 if s.startswith("-") else 1;s=s.lstrip("+-")
    if ":" in s:
        a,b,c=map(float,s.split(":"))
    else:
        q=float(s);a=int(q//10000);b=int(q//100)%100;c=q%100
    assert 0<=b<60 and 0<=c<60
    u=(a+b/60+c/3600)*(15 if ra else sign)
    assert (0<=u<360 if ra else abs(u)<=90)
    return u
def pointing(h):
    r=angle(h["src_raj"],True);d=angle(h["src_dej"],False);assert r!=0 or d!=0
    return r,d
def sep(a,b):
    r,d=map(math.radians,a);v,w=map(math.radians,b)
    x=math.sin((d-w)/2)**2+math.cos(d)*math.cos(w)*math.sin((r-v)/2)**2
    return math.degrees(2*math.asin(math.sqrt(min(1,max(0,x)))))
def history():
    targets=set(H["exclude_extra"])|{"HIP93069"};days=set(H["shield"]["exclude_days"])|{58000};ns=set();urls=set();coords=[]
    for a in H["r36f2"]["prior_manifests"]:
        z,n=artifact(a["id"],a["digest"]);d=json.loads(z.read(n["manifest.json"]));z.close()
        for q in d["selected"]:
            targets.add(q["target"]);days.add(int(q.get("day",q["files"][0]["day"])))
            hh=q["files"][0]["header"]; mm=q["files"][0].get("metadata",{})
            print("R37B_HISTORY_COORD",json.dumps({"target":q["target"],"stage":a["stage"],"src_raj":hh.get("src_raj"),"src_dej":hh.get("src_dej"),"archive_coordinates":{k:v for k,v in mm.items() if k.lower() in ("ra","dec","raj","dej","src_raj","src_dej","ra_deg","dec_deg")}}),flush=True)
            coords.append({"target":q["target"],"pointing":pointing(hh)})
            for v in q["files"]:ns.add(v["name"]);urls.add(v["url"])
    return targets,days,ns,urls,coords
def entries(rows,tel,target,days,ns,urls):
    fs={}
    for r in rows:
        if r["telescope"]!=tel or r["target"]!=target:continue
        for v in r.get("files",[]):
            u=v["url"];n=v["name"];day=v["day"]
            if day in days or n in ns or u in urls or not eng.allowed(u):continue
            if tel=="Parkes":
                m=re.fullmatch(r"Parkes_(\d{5})_(\d{5})_"+re.escape(target)+r"_([SR])_fine\.h5",n,re.I)
                if not m or int(m[1])!=day:continue
                sec=int(m[2]);role="ON" if m[3].upper()=="S" else "OFF";backend="parkes"
            else:
                m=re.search(r"(?P<backend>.*?)guppi_(?P<day>\d{5})_(?P<sec>\d{5})_",n,re.I)
                if not m or int(m["day"])!=day:continue
                sec=int(m["sec"]);backend=m["backend"].lower()
                # Filename must carry the exact target token and role.
                k=re.search(r"(?<![A-Za-z0-9])"+re.escape(target)+r"(_OFF)?(?=[_.-]|$)",n,re.I)
                if not k:continue
                role="OFF" if k[1] else "ON"
                if r["query"].upper()!=(target+"_OFF" if role=="OFF" else target).upper():continue
            if not 0<=sec<86400:continue
            fs[u]={"url":u,"name":n,"day":day,"seconds":sec,"position":("S" if role=="ON" else "R"),"kind":role,"backend":backend,"metadata":v["metadata"]}
    return list(fs.values())
def cadence(fs):
    by={}
    for v in fs:by.setdefault((v["day"],v["backend"]),[]).append(v)
    out=[]
    for key,vs in by.items():
        vs.sort(key=lambda v:(v["seconds"],v["url"]))
        for i in range(len(vs)-3):
            q=vs[i:i+4]
            if [v["kind"] for v in q]!=["ON","OFF","ON","OFF"]:continue
            if not all(0<b["seconds"]-a["seconds"]<=A["max_gap_seconds"] for a,b in zip(q,q[1:])):continue
            if len({v["url"] for v in q})==4 and len({v["name"] for v in q})==4:out.append(q)
    return sorted(out,key=lambda x:(x[0]["day"],x[0]["seconds"],x[0]["backend"],x[0]["url"]))
def qualify(tel,target,group,coords):
    hs=[eng.hdr(f) for f in group]
    spec=S["parkes" if tel=="Parkes" else "gbt"]
    for h,v in zip(hs,group):
        assert len(h["shape"])==3 and h["shape"][0]==16 and h["shape"][1]>=1
        assert 0<abs(h["df"])*1e6<=10 and 0<h["ts"]<=60
        assert abs((h["mjd"]-v["day"])*86400-v["seconds"])<=3
        assert all(h["bounds"][0]<x-(spec["half_hz"]/1e6+.001) and h["bounds"][1]>x+(spec["half_hz"]/1e6+.001) for x in spec["windows_mhz"])
        expected=target+"_"+v["position"] if tel=="Parkes" else target+("_OFF" if v["kind"]=="OFF" else "")
        assert h["source"].upper()==expected.upper()
    assert len({(h["shape"][0],round(h["ts"],7),round(h["df"],12)) for h in hs})==1
    df=hs[0]["df"];assert all(abs((h["f0"]-hs[0]["f0"])/df-round((h["f0"]-hs[0]["f0"])/df))<.001 for h in hs)
    pp=[pointing(h) for h in hs];assert sep(pp[0],pp[2])<.001
    assert sep(pp[0],pp[1])>=.01 and sep(pp[0],pp[3])>=.01
    if tel=="Parkes":assert pp[1]==pp[3]
    assert all(sep(pp[0],c["pointing"])>=A["pointing_separation_degrees"] for c in coords),"Prior or duplicate sky pointing"
    return hs,pp
def main():
    q=S["source"];z,n=artifact(q["id"],q["digest"]);rr=json.loads(z.read(n["queries.json"]));z.close()
    targets,days,ns,urls,coords=history();selected={"Parkes":[],"GBT":[]};attempts=[]
    for tel in ("Parkes","GBT"):
        order=sorted({r["target"] for r in rr if r["telescope"]==tel and r["target"] not in targets},key=lambda x:sha(("R37B|"+tel+"|"+x).encode()))
        for target in order:
            fs=entries(rr,tel,target,days,ns,urls);cand=cadence(fs)
            a={"telescope":tel,"target":target,"eligible_files":len(fs),"cadence_candidates":len(cand),"header_attempts":[]}
            for group in cand[:A["header_attempts_per_target"]]:
                try:
                    hs,pp=qualify(tel,target,group,coords);i=len(selected[tel]);day=group[0]["day"]
                    split="component_fit" if i<2 else "null_calibration" if i<4 else "evaluation"
                    selected[tel].append({"telescope":tel,"target":target,"day":day,"split":split,"pointing":pp[0],"files":[{**f,"header":h,"role":role} for f,h,role in zip(group,hs,A["required_roles"])]})
                    days.add(day);targets.add(target);coords.append({"target":target,"pointing":pp[0]});ns.update(f["name"] for f in group);urls.update(f["url"] for f in group)
                    a["header_attempts"].append({"day":day,"status":"HEADER_QUALIFIED"});break
                except Exception as e:a["header_attempts"].append({"day":group[0]["day"],"names":[v["name"] for v in group],"error":type(e).__name__+": "+str(e)[:200]})
            attempts.append(a);save("partial.json",{"selected":selected,"attempts":attempts})
            if len(selected[tel])==A["per_telescope"]:break
    flat=selected["Parkes"]+selected["GBT"];assert len({r["day"] for r in flat})==len(flat) and len({r["target"] for r in flat})==len(flat)
    d={"status":"METADATA_QUALIFICATION_COMPLETE","protocol_sha256":sha(P.read_bytes()),"source":S["source"],"selected":selected,"attempts":attempts,"parkes_count":len(selected["Parkes"]),"gbt_count":len(selected["GBT"]),"required_per_telescope":8,"parkes_gate":"READY_HEADER_ONLY" if len(selected["Parkes"])==8 else "HOLD","gbt_gate":"READY_HEADER_ONLY" if len(selected["GBT"])==8 else "HOLD","new_obs_power_read":False,"source_origin_labels":False,"exam_started":False,"production_adoption":False}
    save("manifest.json",d);summary={k:v for k,v in d.items() if k not in ("selected","attempts")};summary["groups"]={k:[{i:v[i] for i in ("target","day","split","pointing")} for v in vv] for k,vv in selected.items()};summary["attempts"]=len(attempts);summary["header_attempts"]=sum(len(a["header_attempts"]) for a in attempts)
    print("R37B_SUMMARY",json.dumps(summary),flush=True)
if __name__=="__main__":main()
