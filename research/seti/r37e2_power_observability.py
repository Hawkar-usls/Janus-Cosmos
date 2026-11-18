#!/usr/bin/env python3
"""R37E2 aggregate power observation; no device identity models or per-device descriptors."""
import collections,hashlib,io,json,math,os,pathlib,urllib.parse,zipfile
import numpy as np
import requests
P=pathlib.Path("research/seti/r37e2_power_observability_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r37e2");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,v):(O/n).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def download(cfg):
    name=cfg["name"];p=O/name;h=hashlib.md5();s=hashlib.sha256();n=0
    with requests.get(f"https://zenodo.org/records/{S['source']['record']}/files/{name}?download=1",stream=True,timeout=(20,120)) as r:
        r.raise_for_status();assert urllib.parse.urlparse(r.url).hostname=="zenodo.org"
        with p.open("wb") as f:
            for b in r.iter_content(1048576):
                n+=len(b);assert n<=cfg["bytes"];h.update(b);s.update(b);f.write(b)
    assert n==cfg["bytes"] and h.hexdigest()==cfg["md5"] and s.hexdigest()==cfg["sha256"]
    return p
def inventory():
    cfg=S["parent"];u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{cfg['artifact']}/zip"
    r=requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=(15,90))
    r.raise_for_status();raw=r.content;assert len(raw)<32*1024**2 and sha(raw)==cfg["sha256"]
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        by={}
        for c in S["source"]["archives"]:
            env=c["environment"];n=[i for i in z.infolist() if pathlib.PurePosixPath(i.filename).name==env+"_inventory.json"]
            assert len(n)==1 and n[0].file_size<=8*1024**2
            rows=json.loads(z.read(n[0]));by[env]=rows
        return by
def iq(raw,typ):
    if typ in ("cf32_le","cf32_be","cf64_le","cf64_be"):
        dt={"cf32_le":"<c8","cf32_be":">c8","cf64_le":"<c16","cf64_be":">c16"}[typ]
        return np.frombuffer(raw,dt).astype(np.complex128)
    dt={"ci16_le":"<i2","ci16_be":">i2","ci8":"i1","cu8":"u1"}[typ]
    q=np.frombuffer(raw,dt).astype(float).reshape(-1,2)
    if typ=="cu8":q-=127.5
    return q[:,0]+1j*q[:,1]
def loss(p,k):
    assert len(p)%k==0
    c=p.reshape(-1,k).sum(1);r=np.repeat(c/k,k)
    assert abs(c.sum()-p.sum())<=1e-12
    return float(.5*np.sum(abs(p-r)))
def spectrum(x,n):
    blocks=len(x)//n
    if blocks<1:return None
    energy=np.zeros(n);sum_t=0.
    for i in range(0,blocks,64):
        z=x[i*n:min(blocks,i+64)*n].reshape(-1,n)
        energy+=np.sum(abs(np.fft.fft(z,axis=1))**2,axis=0)/n
        sum_t+=float(np.sum(abs(z)**2))
    rel=abs(energy.sum()-sum_t)/max(1.,sum_t);assert rel<=1e-10
    if energy.sum()<=0:return None
    p=energy/energy.sum()
    return {"blocks":blocks,"parseval_relative_error":float(rel),"coarse_total_variation":{str(k):loss(p,k) for k in S["method"]["coarsening_factors"]}}
def stats(v):
    if not v:return None
    return {"count":len(v),"minimum":float(min(v)),"p10":float(np.quantile(v,.1)),"median":float(np.median(v)),"p90":float(np.quantile(v,.9)),"maximum":float(max(v))}
def main():
    # Deterministic canaries: conservation, exact flat spectrum, loss for a delta.
    assert loss(np.ones(1024)/1024,64)==0
    v=np.zeros(1024);v[0]=1;assert abs(loss(v,64)-63/64)<1e-12
    rng=np.random.default_rng(3702);x=rng.normal(size=2048)+1j*rng.normal(size=2048);assert spectrum(x,1024)["parseval_relative_error"]<1e-10
    inv=inventory();seen={};cohorts=collections.defaultdict(list);sessions=set();duplicate=0
    rates=S["method"]["sample_rate_hz"];sizes=[S["method"]["fft_size"]]+S["method"]["additional_fft_sizes"]
    receipt=[]
    for cfg in S["source"]["archives"]:
        env=cfg["environment"];path=download(cfg)
        try:
            with zipfile.ZipFile(path) as z:
                members={i.filename:i for i in z.infolist() if not i.is_dir()}
                rows=inv[env];assert len(rows)==sum(n.endswith(".sigmf-data") for n in members)
                for r in rows:
                    day=r["day"];sessions.add((env,r["session"]))
                    m=members[r["meta"]];d=members[r["data"]]
                    assert m.file_size<=65536 and d.file_size==r["data_bytes"] and d.file_size<=128*1024**2
                    mr=z.read(m);assert sha(mr)==r["meta_sha256"]
                    md=json.loads(mr);assert md["global"]["core:sample_rate"]==rates and all(c["core:frequency"]==S["method"]["receiver_center_hz"] for c in md["captures"])
                    raw=z.read(d);h=sha(raw)
                    if h in seen:duplicate+=1;continue
                    seen[h]=day
                    x=iq(raw,r["datatype"]);assert len(x)>0 and np.isfinite(x).all()
                    # No label, frequency-offset extraction, demodulation or individual descriptor output.
                    record={"duration_s":len(x)/rates,"sizes":{str(n):spectrum(x,n) for n in sizes}}
                    cohorts[day].append(record)
            receipt.append({"environment":env,"archive_sha256":cfg["sha256"],"n":len(rows)})
        finally:path.unlink(missing_ok=True)
    assert len(seen)+duplicate==S["source"]["expected_recordings"] and len(sessions)==S["source"]["expected_sessions"]
    assert duplicate==0,"duplicate sample hashes"
    assert {k:len(v) for k,v in cohorts.items()}==S["source"]["expected_days"]
    days={}
    for day,rows in sorted(cohorts.items()):
        dur=[r["duration_s"] for r in rows]
        a={"recordings":len(rows),"duration_seconds":stats(dur),"best_possible_full_record_bin_width_hz":stats([1/v for v in dur]),"long_enough_for_telescope_interval":sum(v>=min(S["telescope_reference"]["integration_s"]) for v in dur),"sizes":{}}
        for n in sizes:
            q=[r["sizes"][str(n)] for r in rows if r["sizes"][str(n)] is not None]
            a["sizes"][str(n)]={"supported":len(q),"bin_width_hz":rates/n,"blocks":stats([v["blocks"] for v in q]),"max_parseval_relative_error":max((v["parseval_relative_error"] for v in q),default=None),"coarse_total_variation":{str(k):stats([v["coarse_total_variation"][str(k)] for v in q]) for k in S["method"]["coarsening_factors"]}}
        days[day]=a
    result={"status":"PASS_AGGREGATE_OBSERVABILITY_ONLY","recordings":len(seen),"sessions":len(sessions),"duplicate_sample_hashes":duplicate,"archive_receipts":receipt,"days":days,"telescope_gate":"HOLD_INSUFFICIENT_OBSERVATION" if not any(v["long_enough_for_telescope_interval"] for v in days.values()) else "HOLD_REQUIRES_PHYSICAL_RESOLUTION_QUALIFICATION","identification_gate":"HOLD","scope":"Pooled scalar power representation loss; no hardware grouping or classification. Different dates and environments confounded. Frequency coarsening only at native available discrete resolution; no fake fine telescope channels or time coverage.","identity_classifier_trained":False,"exam_started":False,"production_adoption":False}
    save("result.json",result);print("R37E2_POWER_SUMMARY",json.dumps(result),flush=True)
if __name__=="__main__":main()
