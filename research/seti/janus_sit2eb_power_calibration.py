"""SIT2E-B July development / September calibration; October power never opened.

Research-only. Genuine archive power is decoded and SHA-checked before any score.
Injection truth is separate from scorer inputs. Sept p95 is in-sample calibration,
NOT an externally validated 5% astronomical false-alarm guarantee.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import pathlib
import sys
import time
import numpy as np
import h5py
import hdf5plugin  # Register HDF5 filters before opening archives.

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "research/seti"))
import janus_sit2eb_full_integrity as integrity
import janus_sit2e_parkes_range_headers as headers
import janus_sit2b_temporal_recovery as sit2b

OUT = ROOT / "runtime/sit2eb_power"
DATES = ("2025-07-31", "2025-09-12")
# Recorded in successful run 45061934842; the first complete file hashes.
SOURCE_SHA = {
 "guppi_60887_50760_6690867462_3IATLAS_S_0001.0002.fil": "7aaa17698201c4eb84694432b0ea193b02ad9604fa2e91ede320f60ea1b399eb",
 "guppi_60887_51102_6690868766_3IATLAS_R_0001.0002.fil": "3c925ea095ed0bfaeae418136550d7990a8776efc465bf0cfbd3db01675aa5ff",
 "guppi_60930_18846_6704918083_3IATLAS_S_0001.0002.h5": "eb0418945b33ddf6ade273534290bb5080d524766ea3838688cb73d128f7c69b",
 "guppi_60930_19193_6704919406_3IATLAS_R_0001.0002.h5": "323fda755f581204b881de151ca039c21f6dd3c8953b1fbccafbabc9c7324962",
}
FRACTIONS = (0.20, 0.40, 0.60, 0.80)
SEED = 22061122
METHODS = ("JANUS_NATIVE_RIDGE", "DIRECT_672_MATCHED", "DIRECT_672_SIDEBAND")
CONTROL_CLASSES = ("NO_INJECTION", "SYNTHETIC_RFI")
POS_CLASSES = ("WEAK_LINE", "STRONG_LINE")
ALL_CLASSES = CONTROL_CLASSES + POS_CLASSES

def sha_file(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(4*1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()

def assert_file(path, spec):
    if not path.is_file() or path.stat().st_size != spec["bytes"]:
        raise ValueError("SOURCE_MISSING_OR_SIZE_MISMATCH:" + spec["filename"])
    actual = sha_file(path)
    if SOURCE_SHA.get(spec["filename"]) != actual:
        raise ValueError("FULL_SHA256_MISMATCH:" + spec["filename"])
    return actual

def read_panels(path, date, spec):
    # Mean aggregation of native 300x128 data into fixed 16x128 panels.
    # Frequency is NEVER rebinned or interpolated. IF=0 is explicitly fixed.
    def reduce(get_block, ntime, nfreq, f0, foff):
        if ntime != 300 or nfreq != 64000 or abs(foff-0.002) > 1e-9 or abs(f0-704) > 1e-6:
            raise ValueError("UNEXPECTED_MEASUREMENT_GRID")
        out=[]
        for fraction in FRACTIONS:
            start = int((nfreq-128)*fraction)
            values = np.asarray(get_block(start, start+128),dtype=np.float64)
            if values.shape != (300,128) or not np.isfinite(values).all():
                raise ValueError("INVALID_POWER_VALUES")
            panel = np.stack([x.mean(axis=0) for x in np.array_split(values,16,axis=0)]).astype(np.float32)
            if not np.isfinite(panel).all() or np.any(panel < 0):
                raise ValueError("INVALID_AGGREGATED_POWER")
            out.append({"panel":panel, "start_bin":start,
                        "frequency_start_mhz":f0 + start*foff,
                        "frequency_end_mhz":f0 + (start+127)*foff,
                        "panel_sha256":hashlib.sha256(panel.tobytes()).hexdigest()})
        return out
    if date == "2025-07-31":
        with path.open("rb") as fh:
            raw = fh.read(8192)
        h = headers.parse_sigproc_header(raw)
        if h["nifs"] != 4 or h["nbits"] != 32:
            raise ValueError("UNEXPECTED_SIGPROC_IF_OR_BITS")
        size = 300*4*64000
        if path.stat().st_size != h["header_bytes"]+size*4:
            raise ValueError("SIGPROC_HEADER_OFFSET_OR_LENGTH_MISMATCH")
        mm = np.memmap(path, mode="r", dtype="<f4", offset=h["header_bytes"],shape=(300,4,64000))
        try:
            return reduce(lambda a,b:mm[:,0,a:b],300,64000,h["fch1"],h["foff"])
        finally:
            del mm
    with h5py.File(path,"r") as h:
        d=h["data"]
        if tuple(d.shape)!=(300,4,64000):
            raise ValueError("HDF5_SHAPE_MISMATCH")
        attrs=dict(h.attrs)
        attrs.update(dict(d.attrs))
        f0=float(attrs.get("fch1",attrs.get("FCH1")))
        foff=float(attrs.get("foff",attrs.get("FOFF")))
        return reduce(lambda a,b:d[:,0,a:b],300,64000,f0,foff)

def score(panel, model):
    native_path, shadow, ridge = sit2b.SIT1.score_one(panel,model)
    norm, _ = sit2b.robust_normalize(panel)
    full_path, full_score, _ = sit2b.direct_matched(panel,sideband=False)
    side_path, side_score, _ = sit2b.direct_matched(panel,sideband=True)
    if not np.array_equal(native_path,full_path) or not np.isclose(ridge,full_score,atol=1e-9,rtol=0):
        raise AssertionError("NATIVE_AND_MATCHED_SAME_OBJECTIVE_DIVERGE")
    return {
        METHODS[0]: (native_path,ridge),
        METHODS[1]: (full_path,full_score),
        METHODS[2]: (side_path,side_score)
    }

def create_overlay(panel, category, rng):
    # Existing SIT2B generator; no ground truth sent into score().
    return sit2b.generate_overlay(panel, sit2b.CATEGORIES[ALL_CLASSES.index(category)], rng)

def run():
    manifest_bytes=integrity.MANIFEST.read_bytes()
    manifest=json.loads(manifest_bytes)
    selected=integrity.selected_files(manifest)
    if len(selected)!=4 or any(d not in DATES for d,_,_ in selected):
        raise ValueError("INVALID_FROZEN_SPLIT")
    OUT.mkdir(parents=True,exist_ok=True)
    source_dir=OUT/"sources"
    source_dir.mkdir(exist_ok=True)
    report={"schema":"janus.sit2eb.parkes_power_calibration.v1",
            "status":"POWER_CALIBRATION_HOLD", "october_power_read":False,
            "october_access_forbidden":True,"production_changed":False,
            "source_manifest_sha256":hashlib.sha256(manifest_bytes).hexdigest(),
            "files":[],"panels":[],"errors":[],"results":{},
            "numeric_thresholds_frozen":False,
            "scientific_gain_confirmed":False,
            "no_independent_telescope_pointing":True,
            "threshold_origin":"September synthetic/control scores only"}
    file_panels={}
    try:
        with __import__("requests").Session() as sess:
            for date,role,spec in selected:
                # Re-download each file exactly once and verify complete byte SHA before decoding.
                response = sess.get(spec["url"],stream=True,timeout=(20,100),
                                    headers={"Accept-Encoding":"identity"},allow_redirects=False)
                try:
                    response.raise_for_status()
                    if response.status_code != 200:
                        raise ValueError("RANGE_OR_REDIRECT_NOT_ALLOWED")
                    dest=source_dir/spec["filename"]
                    total=0
                    with dest.open("wb") as fh:
                        for block in response.iter_content(chunk_size=4*1024*1024):
                            if not block:continue
                            total+=len(block)
                            if total>spec["bytes"]:raise ValueError("OVERLONG_DOWNLOAD")
                            fh.write(block)
                    digest=assert_file(dest,spec)
                    panels=read_panels(dest,date,spec)
                    file_panels[(date,role)]=panels
                    report["files"].append({"date":date,"role":role,"sha256":digest,
                                            "bytes":total,"filename":spec["filename"]})
                    report["panels"].extend([{"date":date,"role":role,"start_bin":p["start_bin"],
                        "frequency_start_mhz":p["frequency_start_mhz"],
                        "frequency_end_mhz":p["frequency_end_mhz"],
                        "panel_sha256":p["panel_sha256"]} for p in panels])
                finally:
                    response.close()
                    (source_dir/spec["filename"]).unlink(missing_ok=True)
        model=sit2b.SIT1.cosmos.frozen_model()
        rng=np.random.default_rng(SEED)
        records=[]
        # July is development. September controls calibrate; no October path exists.
        for date in DATES:
            src=file_panels[(date,"source")]
            ref=file_panels[(date,"reference")]
            for idx in range(4):
                assert src[idx]["start_bin"]==ref[idx]["start_bin"]
                for role,panel_entry in (("source",src[idx]),("reference",ref[idx])):
                    base=panel_entry["panel"]
                    for rep in range(12):
                        for category in ALL_CLASSES:
                            values,truth=create_overlay(base,category,rng)
                            scored=score(values,model)
                            records.append({"date":date,"role":role,"band":idx,
                                            "rep":rep,"kind":category,
                                            "scores":{k:float(v[1]) for k,v in scored.items()},
                                            "located":{k:sit2b.localizes(v[0],truth)
                                                       for k,v in scored.items()}})
        july=[r for r in records if r["date"]==DATES[0]]
        september=[r for r in records if r["date"]==DATES[1]]
        assert len(july)==len(september)==384
        results={}
        for method in METHODS:
            negative=[r["scores"][method] for r in september if r["kind"] in CONTROL_CLASSES]
            assert len(negative)==192
            threshold=float(np.quantile(negative,.95,method="higher"))
            def group_eval(rows):
                controls=[r for r in rows if r["kind"] in CONTROL_CLASSES]
                strong=[r for r in rows if r["kind"]=="STRONG_LINE"]
                weak=[r for r in rows if r["kind"]=="WEAK_LINE"]
                return {"control_exceedances":sum(r["scores"][method]>threshold for r in controls),
                        "control_count":len(controls),
                        "strong_detected":sum(r["scores"][method]>threshold and r["located"][method] for r in strong),
                        "strong_count":len(strong),
                        "weak_detected":sum(r["scores"][method]>threshold and r["located"][method] for r in weak),
                        "weak_count":len(weak)}
            results[method]={"threshold_september_p95":threshold,
                             "july_retrospective_only":group_eval(july),
                             "september_in_sample_calibration":group_eval(september)}
        assert np.isclose(results[METHODS[0]]["threshold_september_p95"],
                          results[METHODS[1]]["threshold_september_p95"],atol=1e-9)
        report["results"]=results
        report["records_total"]=len(records)
        report["numeric_thresholds_frozen"]=True
        report["status"]="SEPTEMBER_NUMERIC_P95_CALIBRATED__IN_SAMPLE_ONLY__OCTOBER_HOLD"
        report["calibration_lock_sha256"]=hashlib.sha256(json.dumps({
            "source_manifest_sha256":report["source_manifest_sha256"],
            "files":report["files"],"thresholds":{k:v["threshold_september_p95"] for k,v in results.items()},
            "seed":SEED,"methods":METHODS},sort_keys=True).encode()).hexdigest()
    except Exception as exc:
        report["errors"].append({"type":type(exc).__name__,"message":str(exc)[:400]})
    finally:
        (OUT/"sit2eb_calibration.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("SIT2EB_POWER_RESULT",json.dumps({k:report.get(k) for k in
          ("status","files","results","errors","calibration_lock_sha256",
           "october_power_read")},sort_keys=True),flush=True)
    if report["errors"]:
        raise RuntimeError(report["errors"][0]["message"])
    return report

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--calibrate",action="store_true",required=True)
    ap.parse_args()
    run()
