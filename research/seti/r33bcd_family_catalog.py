#!/usr/bin/env python3
"""R33B/C/D: source catalog, raw-replay sidecar, held-scan diagnostic and protection."""
import os,json,hashlib,pathlib,time,datetime,math
import numpy as np
import h5py,hdf5plugin
import r30c_hannah_memory as mem
import r26c_reference_alignment as ref
import r33_verified_families as a
import r9_expanded_search as eng
import r26f_voyager as voy
import r32_hannah_adaptive as adapt
P=pathlib.Path("research/seti/r33bcd_family_catalog_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r33bcd");OUT.mkdir(parents=True,exist_ok=True);mem.OUT=OUT;ref.BLC_OUT=OUT
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,v):mem.save(n,v)
def emit(k,v):mem.emit(k,v)
def archive(key,cfg):
    ref.BLC_P["artifacts"][key]={"id":cfg["id"],"digest":cfg["digest"]}
    return ref.archive(key)
def revalidate():
    c=S["r33a"];z,n,p=archive("r33a",c)
    raw=z.read(n["external_label_candidates.json"]);assert sha(raw)==c["candidate_sha256"]
    assert sha(z.read(n["manifest.json"]))==c["manifest_sha256"]
    source=ref.get(z,n,"raw_receipts.json");rows={r["tile"]["slice_sha256"]:r for r in source}
    data=json.loads(raw);out=[]
    for e in data["entries"]:
        assert e["evidence_sha256"]==mem.digest({k:v for k,v in e.items() if k!="evidence_sha256"})
        pr=e["passport"];mem.validate(pr);r=rows[pr["provenance"]["source_tile_sha256"]]
        x=ref.tile(z,n,r["tile"],r["artifact"]);path,peak=a.peak(x,r["tile"],e["published_frequency_mhz"],4)
        assert path is not None and peak==e["measured_peak"]
        again=mem.make_record(x,path,pr["provenance"]);assert again==pr
        matches=[rows[h] for h in e["paired_tile_sha256"]]
        assert len(matches)==3 and {q["role"] for q in matches}=={"ON","OFF"}
        on=next(q for q in matches if q["role"]=="ON")
        assert all((q["header"]["src_raj"],q["header"]["src_dej"])!=(on["header"]["src_raj"],on["header"]["src_dej"]) for q in matches if q["role"]=="OFF")
        out.append({"scope":"PUBLISHED_CENTRAL_INTERFERER_OCCURRENCE","family_group":e["family"],"family_member_kind":"STATIONARY_CENTRAL_ANCHOR_NOT_LOOKALIKE","origin_label":"VERIFIED_RFI_EXTERNAL_SOURCE_REPLAY","publication_doi":S["doi"],"publication_scope":"Section 4.4.2; central frequency only; no specific hardware attribution","raw_evidence":e,"scan_role":r["role"],"source_artifact":c["id"],"source_zip_sha256":c["digest"]})
    z.close();p.unlink()
    assert len(out)==6 and len({x["raw_evidence"]["passport"]["id"] for x in out})==6
    return out
def timestamp(utc):
    return datetime.datetime.fromisoformat("2019-04-29T"+utc).replace(tzinfo=datetime.timezone.utc).timestamp()/86400+40587
def track_table(x,t,h,row):
    y=mem.norm(x.astype(float));n,w=y.shape;ts=t["ts"];df=t["df_hz"]
    start=int(round((timestamp(row["utc"])-h["tstart"])*86400/ts))
    freq=t["channel0_mhz"]*1e6+np.arange(w)*df
    center=row["frequency_mhz"]*1e6
    jj=np.flatnonzero(abs(freq-center)<=8.000001)
    if not len(jj):return None
    lim=max(.008,row["drift_hz_s"]*row["drift_error_proxy_percent"]/100)
    rates=np.linspace(max(0.,row["drift_hz_s"]-lim),row["drift_hz_s"]+lim,25)
    best=None
    for st in range(max(0,start-1),min(n-15,start+2)):
        for length in sorted(set(min(n-st,k) for k in (16,32,n-st))):
            tt=np.arange(length)*ts
            for d in rates:
                off=np.rint(d*tt/df).astype(int);lo=max(4,4-int(off.min()));hi=min(w-4,w-4-int(off.max()))
                if hi<=lo:continue
                score=sum(y[st+i,lo+k:hi+k] for i,k in enumerate(off))
                med=np.median(score);sd=1.4826*np.median(abs(score-med))
                if sd<=0:continue
                for j in jj:
                    if not lo<=j<hi:continue
                    q=int(j)+off;occ=float(np.mean(y[st+np.arange(length),q]>3));sn=float((score[j-lo]-med)/sd)
                    if sn<10 or occ<.5:continue
                    o={"start":st,"rows":length,"path":q.tolist(),"frequency_hz":float(freq[j]),"searched_drift_hz_s":float(d),"robust_integrated_score":sn,"occupancy":occ}
                    if best is None or sn>best["robust_integrated_score"]:best=o
    return best
def table_join():
    z,n,p=archive("r33b_blc",S["blc"]);m,mn,mp=archive("r33b_blc_manifest",S["blc_manifest"])
    mf=ref.get(m,mn,"manifest.json");rr=ref.get(z,n,"result.json")["rows"];by={q["name"]:q for q in rr if "tile" in q};hs={q["name"]:q for q in mf["selected"]}
    results=[]
    for row in S["table1"]:
        t0=timestamp(row["utc"])
        fs=[q for q in hs.values() if q["name"].startswith("blc03_ProxCen_S_") and q["header"]["tstart"]-q["header"]["tsamp"]/86400<=t0<q["header"]["tstart"]+q["header"]["rows"]*q["header"]["tsamp"]/86400]
        rec={"published":row,"matches":[q["name"] for q in fs],"status":"HOLD_NO_UNIQUE_RAW_JOIN"}
        if len(fs)==1 and fs[0]["name"] in by:
            f=fs[0];r=by[f["name"]];t=r["tile"];h=f["header"];x=ref.tile(z,n,t,r["artifact"]);v=track_table(x,t,h,row)
            rec.update(scan=f["name"],tile=t,header=h,archive_url=f["url"])
            if v is not None:
                st=v["start"];nr=v["rows"];q=np.array(v["path"],int)
                meta=mem.meta(t,f["name"],S["blc"],st);meta.update({"source_artifact_role":"VERIFIED_RAW_TILE_ARCHIVE","archive_url":f["url"],"day":58602,"role":"ON","source_telescope":"Parkes/Murriyang"})
                passport=mem.make_record(x[st:st+nr],q,meta,verify_full=False)
                # Full raw array was checked above; separately retain exact row-slice digest.
                rec.update(status="EXISTING_ANCHOR_REPLAY" if row["id"]=="BLC1-T1-2" else "SOURCE_LINKED_CANDIDATE_REVIEW_REQUIRED",track=v,passport=passport,row_slice_sha256=sha(x[st:st+nr].tobytes()))
            else:rec["status"]="HOLD_TRACK_NOT_RECOVERED"
        results.append(rec);emit("R33B_TABLE1",{k:v for k,v in rec.items() if k not in ("passport","header","tile","track")})
    z.close();m.close();p.unlink();mp.unlink()
    return results
def recognition(side):
    # Same-campaign held scans; cannot establish transfer to unseen RFI families.
    test=[];refs=[q for q in side if q["scan_role"]=="ON"]
    for q in side:
        if q["scan_role"]!="OFF":continue
        p=q["raw_evidence"]["passport"];rank=[]
        for e in refs:
            r=e["raw_evidence"]["passport"]
            assert p["provenance"]["source_scan"]!=r["provenance"]["source_scan"]
            assert p["provenance"]["source_tile_sha256"]!=r["provenance"]["source_tile_sha256"]
            rank.append({"group":e["family_group"],"full":mem.compare(p,r)["shape_similarity"],"adaptive":adapt.family(p["features"],r["features"])})
        byfull=sorted(rank,key=lambda x:(-x["full"],x["group"]));byadapt=sorted(rank,key=lambda x:(-x["adaptive"],x["group"]))
        test.append({"scan":p["provenance"]["source_scan"],"group":q["family_group"],"scores":rank,"full_top1":byfull[0]["group"]==q["family_group"] and byfull[0]["full"]>byfull[1]["full"],"adaptive_top1":byadapt[0]["group"]==q["family_group"] and byadapt[0]["adaptive"]>byadapt[1]["adaptive"],"correct_group_match_0_9":any(r["group"]==q["family_group"] and r["adaptive"]>=.9 for r in rank),"wrong_group_match_0_9":any(r["group"]!=q["family_group"] and r["adaptive"]>=.9 for r in rank)})
    return {"status":"SAME_CAMPAIGN_DIAGNOSTIC_ONLY","held_scans":len(test),"reference_scans":len(refs),"full_top1":sum(x["full_top1"] for x in test),"adaptive_top1":sum(x["adaptive_top1"] for x in test),"correct_group_matches":sum(x["correct_group_match_0_9"] for x in test),"wrong_group_matches":sum(x["wrong_group_match_0_9"] for x in test),"rows":test,"independent_family_generalization":"NOT_TESTED","scope":"Stationary central anchors, not drifting TF/SF lookalike family members."}
def protection(side):
    # Official public spacecraft signal, actual data, exact SHA already pinned in R26F.
    path,secs=voy.download();expected=voy.P["expected"]
    with h5py.File(path,"r") as h:
        d=h["data"];a=d.attrs;f0=float(a["fch1"]);df=float(a["foff"])*1e6;ts=float(a["tsamp"]);x=np.asarray(d[:,0,:],np.float32)
        name=a["source_name"];name=name.decode() if isinstance(name,bytes) else str(name)
        assert name.lower()=="voyager1" and x.shape==(16,1048576)
    j0=round((expected["frequency_mhz"]-f0)*1e6/df);lo=j0-256;hi=j0+257;assert lo>=0 and hi<=x.shape[1]
    raw=x[:,lo:hi].copy();del x;y=mem.norm(raw.astype(float));n,w=y.shape;freq=f0*1e6+(lo+np.arange(w))*df
    best=None
    for end in range(-150,151):
        dr=end*df/((n-1)*ts)
        if abs(dr-expected["drift_hz_s"])>expected["drift_tolerance"]:continue
        off=np.rint(np.arange(n)*end/(n-1)).astype(int)
        for j in np.flatnonzero(abs(freq-expected["frequency_mhz"]*1e6)<=expected["frequency_tolerance_hz"]):
            p=int(j)+off
            if p.min()<4 or p.max()>=w-4:continue
            score=float(np.sum(y[np.arange(n),p])/np.sqrt(n))
            if best is None or score>best["score"]:best={"score":score,"frequency_mhz":float(freq[j]/1e6),"drift_hz_s":float(dr),"track":p.tolist()}
    assert best is not None
    p=np.array(best["track"],int);f=mem.features(raw,p,df,ts)
    scores=[{"family":r["family_group"],"grid_compatible":math.isclose(f["channel_width_hz"],r["raw_evidence"]["passport"]["features"]["channel_width_hz"],rel_tol=.02) and math.isclose(f["sample_interval_s"],r["raw_evidence"]["passport"]["features"]["sample_interval_s"],rel_tol=.02),"adaptive":adapt.family(f,r["raw_evidence"]["passport"]["features"])} for r in side]
    rec={"status":"KNOWN_PROTECTED_CARRIER_RECOVERED" if best["score"]>=10 else "HOLD_CARRIER_NOT_RECOVERED","source_sha256":voy.P["sha256"],"source":"Voyager1","anchor":best,"native_grid_compatible":sum(r["grid_compatible"] for r in scores),"matched_at_0_9":sum(r["adaptive"]>=.9 for r in scores),"scores":scores,"scope":"Protection challenge only. Grid mismatch/abstention is not learned discrimination or evidence of sensitivity improvement.","production_adoption":False}
    path.unlink();return rec
def main():
    raw=pathlib.Path(S["parent"]["path"]).read_bytes();assert sha(raw)==S["parent"]["sha256"];parent,_=mem.append(json.loads(raw),[]);assert len(parent["entries"])==169
    side=revalidate();assert not {e["raw_evidence"]["passport"]["id"] for e in side}&{p["id"] for p in parent["entries"]}
    catalog={"schema":"janus.hannah.published_catalog.v1","doi":S["doi"],"table1":S["table1"],"direct":[{"figure":7,"index":i+1,"rounded_mhz":v,"direction":"DIRECT","raw_join":"OPEN"} for i,v in enumerate(S["figure7_direct_mhz"])],"mirrored":[{"figure":9,"index":i+1,"rounded_mhz":v,"direction":"MIRRORED","raw_join":"OPEN"} for i,v in enumerate(S["figure9_mirrored_mhz"])],"harmonics":S["table2_harmonic_mhz"],"family_metadata":S["family_catalog"],"scope":S["catalog_scope"]}
    save("published_catalog.json",catalog);save("verified_central_sidecar.json",{"schema":"janus.hannah.verified_occurrence_sidecar.v1","entries":side,"parent_sha256":S["parent"]["sha256"],"authority":"RETRIEVAL_REFERENCE_ONLY_NO_VETO"})
    table=table_join();save("table1_replays.json",table)
    diag=recognition(side);save("recognition.json",diag);emit("R33C_SUMMARY",{k:v for k,v in diag.items() if k!="rows"})
    try:protect=protection(side)
    except Exception as e:protect={"status":"HOLD_PROTECTED_CONTROL_BLOCKED","reason":str(e),"production_adoption":False}
    save("protection.json",protect);emit("R33D_SUMMARY",{k:v for k,v in protect.items() if k!="scores"})
    res={"stage":"R33B","status":"PASS_DOCUMENTED_CATALOG_AND_CENTRAL_REPLAY","parent_entries":169,"verified_main_ledger":1,"verified_sidecar_occurrences":len(side),"central_groups":len({e["family_group"] for e in side}),"direct_catalog":len(catalog["direct"]),"mirrored_catalog":len(catalog["mirrored"]),"table1_linked_candidates":sum(q["status"]=="SOURCE_LINKED_CANDIDATE_REVIEW_REQUIRED" for q in table),"table1_anchor_replays":sum(q["status"]=="EXISTING_ANCHOR_REPLAY" for q in table),"table1_holds":sum(q["status"].startswith("HOLD") for q in table),"catalog_sha256":sha((OUT/"published_catalog.json").read_bytes()),"sidecar_sha256":sha((OUT/"verified_central_sidecar.json").read_bytes()),"recognition":diag["status"],"protection":protect["status"],"r34":"HOLD_FAMILY_DISJOINT_VALIDATION_REQUIRED","production_adoption":False}
    save("result.json",res);emit("R33B_SUMMARY",res)
    # bounded exact sidecar export for durable publication after review
    import base64
    for name in ("verified_central_sidecar.json","published_catalog.json"):
        payload=(OUT/name).read_bytes();b64=base64.b64encode(payload).decode()
        for i in range(0,len(b64),3000):print("R33B_EXPORT",name,i//3000,b64[i:i+3000],flush=True)
if __name__=="__main__":main()
