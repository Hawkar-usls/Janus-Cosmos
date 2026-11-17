#!/usr/bin/env python3
"""Exact candidate relation filtering, evidence planning and provenance catalog."""
import hashlib,io,itertools,json,math,os,pathlib,statistics,time,zipfile
import numpy as np
import requests
import donor_eye_solver as eye
P=pathlib.Path("research/seti/r19bcd_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r19bcd");OUT.mkdir(parents=True,exist_ok=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
def emit(k,d):print(k,json.dumps(d),flush=True)
def archive(key):
    t=time.perf_counter();a=S["artifacts"][key];h=hashlib.sha256();p=OUT/(key+".zip");n=0
    with requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['id']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},stream=True,timeout=(15,120)) as r:
        r.raise_for_status()
        with p.open("wb") as f:
            for b in r.iter_content(1024*1024):
                n+=len(b);assert n<=500*1024**2;h.update(b);f.write(b)
    assert h.hexdigest()==a["digest"]
    with zipfile.ZipFile(p) as z:
        names=[v for v in z.namelist() if pathlib.PurePosixPath(v).name==("manifest.json" if key=="manifest" else "result.json")]
        assert len(names)==1 and z.getinfo(names[0]).file_size<32*1024**2
        d=json.loads(z.read(names[0]))
    p.unlink();return d,{"bytes":n,"seconds":time.perf_counter()-t}
def relation(a,b,ft,dt):
    # columns are frequency projected to common time, drift in Hz/s.
    return (np.abs(a[:,None,0]-b[None,:,0])<=ft)&(np.abs(a[:,None,1]-b[None,:,1])<=dt)
def mats(A,ft,dt):
    return relation(A[0],A[1],ft,dt),relation(A[0],A[2],ft,dt),relation(A[1],A[2],ft,dt)
def dense(A,ft,dt):
    x,y,z=mats(A,ft,dt)
    return np.argwhere(x[:,:,None]&y[:,None,:]&z[None,:,:]),int(np.prod([len(a) for a in A]))
def semi(A,ft,dt):
    x,y,z=mats(A,ft,dt)
    # Intersect common anchor supports. Filter pair neighborhoods by the other support before exact join.
    anchor=np.flatnonzero(x.any(1)&y.any(1));out=[];expanded=0
    for i in anchor:
        b=np.flatnonzero(x[i]&z.any(1));c=np.flatnonzero(y[i]&z.any(0))
        expanded+=len(b)*len(c)
        for k,l in np.argwhere(z[np.ix_(b,c)]):out.append((int(i),int(b[k]),int(c[l])))
    return np.asarray(out,dtype=np.int64).reshape(-1,3),expanded
def bench(A,ft,dt):
    x,xc=dense(A,ft,dt);y,yc=semi(A,ft,dt)
    assert np.array_equal(x,y)
    timing={}
    for method,fn in (("dense",dense),("semijoin",semi)):
        ts=[]
        for _ in range(3):
            t=time.perf_counter();ans,count=fn(A,ft,dt);ts.append(time.perf_counter()-t)
            assert np.array_equal(x,ans)
        timing[method]=statistics.median(ts)
    return {"dense_combinations":xc,"reduced_combinations":yc,"matches":len(x),"seconds":timing}
def stress():
    rows=[];rng=np.random.default_rng(S["stress"]["seed"])
    for n in S["stress"]["sizes"]:
        for density in S["stress"]["density"]:
            for rep in range(S["stress"]["replicates"]):
                span=10000 if density=="sparse" else 10 if density=="dense" else 1000
                a=[np.column_stack((rng.uniform(0,span,n),rng.normal(0,1,n))) for _ in range(3)]
                if density=="hostile":a[2][:,0]+=5000
                # known match except in hostile control.
                if density!="hostile":
                    for v in a:v[0]=[5,.1]
                r=bench(a,5,.1)
                if density!="hostile":assert r["matches"]>=1
                else:assert r["matches"]==0
                rows.append({"n":n,"density":density,**r})
    return rows
def eye_run(real_count):
    raw=pathlib.Path("research/seti/donor_eye_solver.py").read_bytes()
    assert hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()==S["sources"]["eye"]["blob"]
    causes=["sky","rfi","instrument"]
    allpairs=[list(x) for x in itertools.combinations(causes,2)]
    def solve(tests,k):
        c,t,r,a=eye.parse({"cause_classes":causes,"tests":tests,"requirements":{"required_independent_separators_per_pair":k}})
        z=eye.search(c,t,r,a)
        return {"status":z["status"],"exact":z["exact"],"solutions":[sorted(x) for x in z["solutions"]],"cost":min((eye.cost(x,t,a) for x in z["solutions"]),default=None),"examined":z["examined"],"unresolved":[list(x) for x in sorted(eye.evaluate(z["best_partial"],z["pairs"],t,r)["unresolved"])]}
    synthetic=[{"id":"a","cost":2,"distinguishes":allpairs,"failure_domains":["receiver_A"]},{"id":"copy_a","cost":.1,"distinguishes":allpairs,"failure_domains":["receiver_A"]},{"id":"b","cost":3,"distinguishes":allpairs,"failure_domains":["receiver_B"]},{"id":"unknown","cost":0,"distinguishes":allpairs,"failure_domains":[]}]
    a=solve(synthetic,2);assert a["solutions"]==[["b","copy_a"]] and a["cost"]==3.1
    b=solve(synthetic[:2]+synthetic[3:],2);assert b["status"]=="NON_IDENTIFIABLE_UNDER_CURRENT_MEASUREMENT_MODEL"
    # Same-day archive cannot prove discriminatory capability for cause pairs.
    real=[{"id":"archive_onoff","cost":0,"distinguishes":[],"failure_domains":["GBT_61027"]},{"id":"archive_reverse","cost":0,"distinguishes":[],"failure_domains":["GBT_61027"]},{"id":"archive_time_sham","cost":0,"distinguishes":[],"failure_domains":["GBT_61027"]},{"id":"independent_receiver","cost":1,"distinguishes":[],"failure_domains":["UNVERIFIED_FUTURE"],"available":False}]
    c=solve(real,2);assert c["status"]=="NON_IDENTIFIABLE_UNDER_CURRENT_MEASUREMENT_MODEL"
    return {"synthetic_known_minimum":a,"shared_domain_control":b,"real_archive":c,"real_catalog_count":real_count,"scope":"No empirical separation edges asserted; real outcome is missing validated evidence, not rejection of sky."}
def fingerprint(record):
    return hashlib.sha256(json.dumps(record,sort_keys=True,separators=(",",":")).encode()).hexdigest()
def dedup(records):
    # TOPA-inspired lineage and exact duplicate handling; never merge across provenance.
    seen={};duplicates=[];conflicts=[];unique=[]
    for r in records:
        key=(r["source_family"],r["record_id"])
        h=fingerprint(r["payload"])
        if key in seen:
            if seen[key]==h:duplicates.append(key);continue
            conflicts.append(key)
        else:seen[key]=h
        unique.append(r)
    return unique,duplicates,conflicts
def catalog(windows,manifest):
    by={x["center_mhz"]:x for x in manifest["selected"]};records=[];joins=[]
    for w in windows:
        if w["status"]!="COMPLETE":continue
        fc=w["center_mhz"];files=by[fc]["files"];rows=w["rows"]
        tref=files[0]["header"]["mjd"];df=abs(rows[0]["tile"]["df_hz"]);dur=(rows[0]["tile"]["shape"][0]-1)*rows[0]["tile"]["ts"]
        arr=[]
        for i,r in enumerate(rows):
            t=r["tile"];delta=(files[i]["header"]["mjd"]-tref)*86400
            hits=[]
            for k,h in enumerate(r["hits"]):
                f=t["channel0_mhz"]*1e6+h["bin"]*t["df_hz"];d=h["drift_hz_s"]
                projected=f-d*delta;hits.append([projected,d])
                records.append({"source_family":"GBT_MJD61027","record_id":f'{fc}:{r["name"]}:{k}',"payload":{"window_mhz":fc,"role":r["role"],"name":r["name"],"tile_sha256":t["slice_sha256"],"frequency_hz":f,"projected_hz":projected,"drift":d,"score":h["score"]}})
            arr.append(np.asarray(hits,float).reshape(-1,2))
        for label,idx in (("ON",[0,2,4]),("OFF",[1,3,5])):
            joins.append({"window":fc,"cadence":label,**bench([arr[i] for i in idx],S["match"]["freq_tolerance_bins"]*df,S["match"]["drift_tolerance_bins_per_scan"]*df/dur)})
    unique,dup,conf=dedup(records)
    assert len(unique)==len(records) and not dup and not conf
    # Stress exact repeat; retain distinct archive IDs and contradictory payloads.
    r=records[0];copy=dict(r);other={**r,"source_family":"SYNTHETIC_OTHER_DOMAIN"};bad={**r,"payload":{**r["payload"],"score":-123}}
    u,d,c=dedup([r,copy,other,bad]);assert len(u)==3 and len(d)==1 and len(c)==1
    # Advisory compatibility edges, no transitive or cause equivalence claim.
    edges=[]
    for i,a in enumerate(records):
        for j,b in enumerate(records[:i]):
            x=a["payload"];y=b["payload"]
            if x["window_mhz"]==y["window_mhz"] and abs(x["projected_hz"]-y["projected_hz"])<=25 and abs(x["drift"]-y["drift"])<=.03:
                edges.append({"a":a["record_id"],"b":b["record_id"],"same_archive_domain":True})
    save("rfi_catalog.json",{"records":records,"similarity_edges":edges,"scope":"Advisory only; all raw hits retained; same receiver/day is one failure domain, not independently confirmed RFI families."})
    return joins,{"raw_records":len(records),"unique_records":len(unique),"exact_duplicates":len(dup),"conflicts":len(conf),"similarity_edges":len(edges),"source_domains":1,"synthetic_duplicate_control":"PASS","scope":"TOPA-inspired numeric adaptation; no classifier or detection authority"}
def main():
    t=time.perf_counter();d,io1=archive("r16c");m,io2=archive("manifest")
    joins,topa=catalog(d["windows"],m);st=stress();e=eye_run(topa["raw_records"])
    def aggregate(rows):
        return {"cases":len(rows),"dense_combinations":sum(x["dense_combinations"] for x in rows),"reduced_combinations":sum(x["reduced_combinations"] for x in rows),"matches":sum(x["matches"] for x in rows),"seconds":{k:sum(x["seconds"][k] for x in rows) for k in ("dense","semijoin")},"mismatches":0}
    out={"eye":e,"semijoin":{"real":aggregate(joins),"stress":{k:aggregate([x for x in st if x["density"]==k]) for k in S["stress"]["density"]},"per_window":joins,"stress_rows":st},"topa":topa,"archive_io":[io1,io2],"run_seconds":time.perf_counter()-t,"confirmed_technosignatures":0,"production_adoption":False}
    save("result.json",out);emit("R19BCD_SUMMARY",{**out,"semijoin":{k:v for k,v in out["semijoin"].items() if k not in ("per_window","stress_rows")}})
if __name__=="__main__":main()
