#!/usr/bin/env python3
"""Revalidate and export exact R33B2 candidate bytes; no origin promotion."""
import base64,hashlib,json,pathlib
import numpy as np
import r26c_reference_alignment as ref
import r33b2_drift_replay as r
import r30c_hannah_memory as mem
OUT=pathlib.Path("runtime/r33b3");OUT.mkdir(parents=True,exist_ok=True);ref.BLC_OUT=OUT
ART={"id":14414154131,"digest":"fe97f959ec85a56457002cdf17f07a19e78709ee363d1fa75b7ff92b44200074"}
CAND="9593bce5c10af104dc9f3556aa3869909d849a06418bdc3c5f139a6826abc648"
MAN="c8f7e84042ddb756c5e9b2d385708b8792a016fb77b4a6d68ad956c9126bf786"
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
    ref.BLC_P["artifacts"]["r33b2"]=ART
    z,n,p=ref.archive("r33b2")
    raw=z.read(n["drift_candidates.json"]);assert sha(raw)==CAND
    assert sha(z.read(n["manifest.json"]))==MAN
    parent=pathlib.Path(r.S["parent"]["path"]).read_bytes();assert sha(parent)==r.S["parent"]["sha256"]
    rows=ref.get(z,n,"raw_receipts.json");assert len(rows)==36
    cand=json.loads(raw);assert len(cand["entries"])==2
    assert {e["published"]["id"] for e in cand["entries"]}=={"F7-9","F9-16"}
    results=[];ids=set();hashes=set()
    for e in cand["entries"]:
        assert e["evidence_sha256"]==mem.digest({k:v for k,v in e.items() if k!="evidence_sha256"})
        row=next(v for v in rows if v["published"]["id"]==e["published"]["id"])
        assert row["status"]=="SOURCE_LINKED_DRIFT_CANDIDATE_REVIEW_REQUIRED"
        pair=[(e["passport"],e["on_track"])] + [(o["passport"],o["measured_track"]) for o in e["off_matches"]]
        assert len(pair)==2
        for pr,hit in pair:
            mem.validate(pr);assert pr["label"]=="UNKNOWN"
            meta=pr["provenance"];h=meta["source_tile_sha256"]
            rr=[v for v in row["measurements"] if v.get("tile",{}).get("slice_sha256")==h and v["scan"]==meta["source_scan"]]
            assert len(rr)==1;v=rr[0];x=ref.tile(z,n,v["tile"],v["raw_artifact"])
            st=hit["start"];length=hit["rows"];path=np.array(hit["path"],int)
            assert len(path)==length and sha(x[st:st+length].tobytes())==meta["row_slice_sha256"]
            assert mem.make_record(x[st:st+length],path,meta,verify_full=False)==pr
            assert v["search"]["unique"] and v["search"]["best"]==hit
            assert pr["id"] not in ids and h not in hashes;ids.add(pr["id"]);hashes.add(h)
            assert pr["features"]["drift_hz_s"]*e["published"]["sign"]>0
            if v["role"]=="ON":
                rate=row["expected_drift_hz_s"];tol=abs(rate)*r.S["relative_rate_tolerance"]
                center=e["published"]["rounded_mhz"]*1e6;rad=r.A["start_frequency_within_hz"]
            else:
                rate=e["on_track"]["drift_hz_s"];tol=abs(rate)*r.G["off_drift_tolerance_fraction"]
                center=v["predicted_frequency_hz"];rad=r.G["off_prediction_tolerance_hz"]
            check=r.search(x,v["tile"],rate,tol,center,rad)
            assert check==v["search"]
            # Quantization permits endpoint estimate to differ from a grid trial.
            endpoint_resolution=abs(meta["df_hz"])/((length-1)*meta["ts"])
            assert abs(hit["drift_hz_s"]-rate)<=tol+endpoint_resolution+1e-12
            results.append({"id":pr["id"],"figure":e["published"]["id"],"role":v["role"],"scan":v["scan"],"raw_hash":h,"rows":length,"frequency_hz":pr["frequency_hz"],"drift_hz_s":pr["features"]["drift_hz_s"],"endpoint_resolution_hz_s":endpoint_resolution,"expected_rate_hz_s":rate,"score":hit["score"],"occupancy":hit["occupancy"]})
        on,off=results[-2:];assert on["role"]=="ON" and off["role"]=="OFF"
        assert abs(e["off_matches"][0]["frequency_residual_hz"])<=150.00001
    z.close();p.unlink()
    assert len(ids)==len(hashes)==4
    (OUT/"drift_candidates.json").write_bytes(raw)
    res={"status":"PASS_EXACT_MEASURED_DRIFT_REVALIDATION","candidate_count":2,"passports":4,"distinct_raw_tiles":4,"sha256":CAND,"bytes":len(raw),"results":results,"label":"REVIEW_REQUIRED","verified_promoted":0,"production_adoption":False}
    (OUT/"verification.json").write_text(json.dumps(res,indent=2)+"\n")
    enc=base64.b64encode(raw).decode()
    for i in range(0,len(enc),3000):print("R33B3_EXPORT",i//3000,enc[i:i+3000],flush=True)
    print("R33B3_SUMMARY",json.dumps(res),flush=True)
if __name__=="__main__":main()
