"""R42A4A: frozen NOISE raw score coverage beyond Minimal candidate gate; historical pilot only."""
import json,hashlib,pathlib,io
import numpy as np
import r41s1f_unclipped_joint as f
import r42a3b_full_score_contract as full
base=f.base;a=f.a;n=f.n
P=pathlib.Path("research/seti/r42a4a_noise_coverage_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r42a4a");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,v):(O/n).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def frozen_score(features,model):
    if features is None:return None
    keys=model["features"];v=np.array([np.nan if features.get(k) is None else features[k] for k in keys],float)
    drop=np.asarray(model["dropped_all_missing"],bool);im=np.asarray(model["imputation"],float)
    assert len(v)==len(drop)==len(im) and np.isfinite(im).all()
    x=np.where(np.isfinite(v),v,im)[~drop]
    mu=np.asarray(model["scaler_mean"],float);scale=np.asarray(model["scaler_scale"],float)
    coef=np.asarray(model["coefficients"],float);intercept=float(model["intercept"])
    assert len(x)==len(mu)==len(scale)==len(coef) and (scale>0).all() and np.isfinite(x).all()
    z=float(((x-mu)/scale)@coef+intercept)
    return z
def prepare():
    base.prepare()
    original,rd=f.download(f.L["source"])
    assert sha(rd("result.json"))==S["reference"]["noise_result_sha256"]
    assert original["baseline"]=={"controls":53,"localized":318,"weak":69,"ordinary":205,"shared":113}
    assert original["models"]["raw_noise"]["hypothetical"]["controls"]==51
    models=original["models"]["raw_noise"]
    for lane in ("parent","rescue"):
        m=models[lane];assert m["status"]=="HISTORICAL_ONLY"
        assert len(m["features"])==len(m["dropped_all_missing"])==len(m["imputation"])
    common,_=full.unpack()
    assert len(common["records"])==22 and sum(r["rows"] for r in common["records"])==3168
    save("preparation.json",{"status":"AUTHENTICATED_FROZEN_ARCHIVES","protocol_sha256":sha(P.read_bytes()),"original_fit_only":True,"baseline_reference":original["baseline"],"noise_hypothetical":models["hypothetical"]})
    print("R42A4A_PREPARED",flush=True)
def group(i):
    assert (i+1) in S["pilot_groups"]
    original,rd=f.download(f.L["source"])
    assert sha(rd("result.json"))==S["reference"]["noise_result_sha256"]
    models=original["models"]["raw_noise"];common,rc=full.unpack()
    plan=base.plan();g=plan["groups"][i];num=f"{i+1:02d}"
    seed={"n6":a.v.C["seed"],"n7":base.w.C["seed"],"n9":base.old9.C["seed"]}[g["source"]]
    n.set_config(seed,a.C["windows_mhz"],a.C["half_hz"],a.C["max_mib_per_scan"]);base.b.taylor.init_tree()
    m,LL,old=base.locks("minimal",plan)
    rr,R,ch,hs,meta,cm,pre=n.gather(g["group"],"r42a4a_"+num)
    assert sorted(hs)==g["raw_tile_hashes"]
    old_group=json.loads(rd("group-"+num+".json"))
    assert old_group["group"]==i+1 and old_group["fit_viewed"]==g["fit"]
    old_rows={(r["kind"],r["row"]):r for r in old_group["records"]}
    archived_pred={}
    if not g["fit"]:
        for lane in ("parent","rescue"):
            for r in models[lane]["predictions"]:
                if r["group"]==i+1:
                    key=(r["kind"],r["row"]);assert key not in archived_pred
                    archived_pred[key]=r
    outputs=[];allrows=[]
    for kind,rows,refs,mm in (("eval",rr,R,meta),("shared",[x[0] for x in ch],np.stack([x[1] for x in ch]),cm)):
        d=base.infer("minimal",rows,refs,mm,pre,LL,old,m,None,kind)
        oldarr=a.read_npz(base.C/("dev_"+num+"_"+kind+".npz"));a.verify_rows(oldarr,rows);base.check(d,oldarr,m)
        name=num+"-"+kind+".npz";rec=next(x for x in common["records"] if x["name"]==name)
        raw=rc(name);assert sha(raw)==rec["sha256"]
        with np.load(io.BytesIO(raw),allow_pickle=False) as q:
            X=n.ref.norm(np.stack([v[0] for v in rows]));Y,_,_,_=n.ref.register(X,n.ref.norm(refs))
            assert sha(Y.astype(np.float32).tobytes())==rec["array_sha256"]
            assert np.array_equal(q["minimal_alert"],d["alert"]) and np.array_equal(q["minimal_hit"],d["hit"])
        scores=np.full(len(rows),np.nan,float);valid=np.zeros(len(rows),bool);eligible=np.zeros(len(rows),bool);causes=[]
        for j in range(len(rows)):
            key=(kind,j);row=old_rows.get(key)
            base_accept=bool(d["alert"][j]);lane="parent" if d["parent_alert"][j] else "rescue"
            assert (row is not None)==base_accept
            if row is not None:assert row["lane"]==lane and bool(row["localized_evaluation_only"])==bool(d["hit"][j])
            path=d["parent_path"][j] if lane=="parent" else d["p"][j]
            prov=f.s.provenance(g["group"]["files"],pre,mm[j],str(d["source_tile_sha"][j]))
            if prov is None:reason="NO_AUTHENTIC_PHYSICAL_MAPPING";features=None
            else:
                blocks,receipts=prov;slot=mm[j][0];oth=2-slot
                au=blocks[slot].copy();af=blocks[slot+1].copy()
                au[:,f.HALO:f.HALO+128]=rows[j][0];af[:,f.HALO:f.HALO+128]=refs[j]
                diag=f.raw_diagnose(np.stack([au,af,blocks[oth],blocks[oth+1]]),np.asarray(path,int)+f.HALO)
                features=None if diag is None else {k:v for k,v in diag["features"].items() if k.startswith("noise_")}
                reason="RAW_SCALE_UNSUPPORTED" if features is None else "SCORED"
            if features is not None:
                scores[j]=frozen_score(features,models[lane]);assert np.isfinite(scores[j])
                valid[j]=True
            if row is not None:
                # Parity with original advanced raw feature extraction and historical models.
                original_features=f.arm(row,"raw_noise")
                if original_features is not None:
                    assert features is not None,"Original qualified row failed physical score reconstruction"
                    assert sorted(original_features)==sorted(features)
                    for k,v in original_features.items():
                        u=features[k]
                        if v is None:assert u is None
                        else:assert u is not None and np.isclose(float(u),float(v),rtol=1e-10,atol=1e-10)
                    eligible[j]=True
                pred=archived_pred.get(key)
                if pred is not None:
                    assert eligible[j] and np.isclose(scores[j],pred["score"],rtol=1e-9,atol=1e-9)
                    cut=float(models[lane]["train_positive_threshold"])
                    assert bool(scores[j]>=cut)==bool(pred["would_keep"])
                elif not g["fit"]:assert not eligible[j]
            causes.append(reason)
            allrows.append({"group":i+1,"kind":kind,"row":j,"baseline_alert":base_accept,"lane":lane,"supported":bool(valid[j]),"frozen_eligible":bool(eligible[j]),"reason":reason,"score":None if not valid[j] else float(scores[j]),"archived_prediction":key in archived_pred})
        np.savez_compressed(O/(num+"-"+kind+"-noise-all.npz"),score=scores,supported=valid,frozen_eligible=eligible,original_minimal_alert=d["alert"],original_minimal_hit=d["hit"],parent_route=d["parent_alert"],parent_path=d["parent_path"],rescue_path=d["p"],tile_sha=d["source_tile_sha"])
        x={"kind":kind,"rows":len(rows),"score_supported":int(valid.sum()),"original_alerts":int(d["alert"].sum()),"original_qualified":int(eligible.sum()),"rejected_scored":int((valid&~d["alert"]).sum()),"rejected_total":int((~d["alert"]).sum()),"unsupported_reasons":{r:causes.count(r) for r in sorted(set(causes))},"file_sha256":sha((O/(num+"-"+kind+"-noise-all.npz")).read_bytes())}
        outputs.append(x);print("R42A4A_KIND",json.dumps({"group":i+1,**x}),flush=True)
    result={"stage":"R42A4A","group":i+1,"fit":g["fit"],"rows":sum(r["rows"] for r in outputs),"original_alerts":sum(r["original_alerts"] for r in outputs),"original_qualified":sum(r["original_qualified"] for r in outputs),"score_supported":sum(r["score_supported"] for r in outputs),"rejected_scored":sum(r["rejected_scored"] for r in outputs),"records":outputs,"production_adoption":False,"external_validation":False}
    assert result["rows"]==288
    if not g["fit"]:assert sum(r["archived_prediction"] for r in allrows)==sum(r["original_qualified"] for r in outputs)
    save("group-"+num+".json",result);save("rows-"+num+".json",allrows)
    print("R42A4A_GROUP",json.dumps(result),flush=True)
def summary():
    gr=[json.loads((O/f"group-{i:02d}.json").read_text()) for i in S["pilot_groups"]]
    assert {x["group"] for x in gr}==set(S["pilot_groups"]) and all(x["rows"]==288 for x in gr)
    d={"status":"PASS_ALL_CANDIDATE_NOISE_COVERAGE_PILOT__HOLD_GLOBAL_SCORE_AND_PROSPECTIVE_DATA","pilot_groups":S["pilot_groups"],"rows":sum(x["rows"] for x in gr),"score_supported":sum(x["score_supported"] for x in gr),"original_qualified":sum(x["original_qualified"] for x in gr),"rejected_scored":sum(x["rejected_scored"] for x in gr),"per_group":gr,"production_adoption":False,"external_validation":False,"full_pipeline_roc":False}
    save("result.json",d);print("R42A4A_SUMMARY",json.dumps({k:v for k,v in d.items() if k!="per_group"}),flush=True)
