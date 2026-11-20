"""R42A CPU native engines on the same audited historical R41 registered arrays."""
import sys,os,json,hashlib,pathlib,time,resource,io,zipfile
import numpy as np,requests
import r35g_native_compare as native
MODE=sys.argv[1];O=pathlib.Path("runtime/r42a-native")/MODE;O.mkdir(parents=True,exist_ok=True)
native.OUT=O
def save(k,d):(O/k).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
    aid=os.environ["INPUT_ID"];digest=os.environ["INPUT_DIGEST"].removeprefix("sha256:");assert aid.isdigit() and len(digest)==64
    url=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{aid}/zip"
    h={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"}
    r=requests.get(url,headers=h,timeout=120);r.raise_for_status();raw=r.content;assert len(raw)<128*1024**2 and sha(raw)==digest
    z=zipfile.ZipFile(io.BytesIO(raw));files=[x for x in z.infolist() if not x.is_dir()]
    assert sum(x.file_size for x in files)<512*1024**2
    names={pathlib.PurePosixPath(x.filename).name:x.filename for x in files};assert len(names)==len(files)
    rec=json.loads(z.read(names["receipt.json"]));assert rec["status"]=="PASS_HISTORICAL_COMMON_INPUT" and len(rec["records"])==22
    fn,init=native.init()
    def load(row):
        b=z.read(names[row["name"]]);assert sha(b)==row["sha256"]
        q=np.load(io.BytesIO(b),allow_pickle=False);assert sha(q["V"].tobytes())==row["array_sha256"];return q
    cal=[];calbase=0;cal_rows=[];cal_pos=0
    for row in rec["records"]:
        if not row["fit"]:continue
        with load(row) as q:
            s,p,t=fn(q["V"],row["name"][:-4]);assert np.isfinite(s).all()
            null=q["y"]==0;cal.append(s[null]);calbase+=int((q["minimal_alert"]&null).sum())
            cal_pos+=int((~null).sum());cal_rows.append({"name":row["name"],"timing":t})
    c=np.concatenate(cal);assert len(c)==rec["fit"]["negatives"] and calbase==rec["fit"]["controls"]
    assert len(c)>0
    def thresh(k):
        return float(np.nextafter(np.sort(c)[::-1][k],np.inf)) if k<len(c) else float(np.min(c))
    budgets={"five_percent":int(len(c)*.05),"matched_minimal_fit":calbase}
    ts={k:thresh(v) for k,v in budgets.items()}
    assert all(int((c>=v).sum())<=budgets[k] for k,v in ts.items())
    save("lock.json",{"thresholds":ts,"budgets":budgets,"fit_controls":len(c),"achieved":{k:int((c>=v).sum()) for k,v in ts.items()},"fit_positives":cal_pos,"input_id":aid,"input_digest":digest,"calibration_timing":cal_rows})
    results=[]
    for row in rec["records"]:
        if row["fit"]:continue
        with load(row) as q:
            s,p,t=fn(q["V"],row["name"][:-4]);assert np.isfinite(s).all()
            null=q["y"]==0;weak=q["weak"];loc=((abs(p-q["truth"])<=2)&q["active"]).sum(1)>=.75*q["active"].sum(1)
            d={"group":row["group"],"kind":row["kind"],"target":row["target"],"positives":int((~null).sum()),"negatives":int(null.sum()),"timing":t,"arms":{}}
            arr={"score":s,"path":p}
            for k,th in ts.items():
                alert=s>=th;hit=alert&loc&(~null)
                v={"localized":int(hit.sum()),"weak":int((hit&weak).sum()),"controls":int((alert&null).sum())}
                for b in ("minimal","noise"):
                    ref=q[b+"_hit"];v["added_vs_"+b]=int((hit&~ref).sum());v["lost_vs_"+b]=int((ref&~hit).sum());v["lost_weak_vs_"+b]=int((ref&~hit&weak).sum())
                d["arms"][k]=v;arr[k+"_alert"]=alert;arr[k+"_hit"]=hit
            np.savez_compressed(O/(row["name"][:-4]+"-scores.npz"),**arr);results.append(d);save("partial.json",results)
    assert len(results)==10
    allarms={}
    for k in ts:
        ds=[d["arms"][k] for d in results];v={n:sum(d[n] for d in ds) for n in ds[0]}
        v["ordinary"]=sum(d["arms"][k]["localized"] for d in results if d["kind"]=="eval");v["shared"]=sum(d["arms"][k]["localized"] for d in results if d["kind"]=="shared")
        v["positives"]=sum(d["positives"] for d in results);v["negatives"]=sum(d["negatives"] for d in results)
        v["control_fraction"]=v["controls"]/v["negatives"];allarms[k]=v
    result={"status":"PASS_HISTORICAL_NATIVE_COMPARISON","method":MODE,"arms":allarms,"janus_reference":rec["audit"],"per_group":results,"compile_setup_seconds":init,"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"input_digest":digest,"input_id":aid,"historical_only":True,"production_adoption":False,"scope":"registered single-view CPU comparator vs R41 richer paired evidence; not full origin pipeline or independent exam"}
    save("result.json",result);print("R42A_NATIVE",json.dumps({k:v for k,v in result.items() if k!="per_group"}),flush=True)
if __name__=="__main__":main()
