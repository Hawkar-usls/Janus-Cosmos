#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, math
from collections import defaultdict
from pathlib import Path
import numpy as np
from janus_cosmos.luci import read_luci_fits_image
from janus_cosmos.luci_psf import _inject_gaussian, robust_background
from janus_cosmos.luci_psf_r1 import measure_psf_at, psf_relative_injection_recovery_gate
from janus_cosmos.pipeline import EventWriter, download_source, sha256_file
from experiments.luci.run_palomar_2f_d import counterpart_with_matched_controls
from experiments.luci.run_palomar_2f_e import local_coordinate_injection

PARENT_ARTIFACT_ID=9248259268
PARENT_ZIP_SHA256="87ce38fe9c6d0977a009efbd2df67c668656a4a82889128bbcf962d666aecff7"
PARENT_EXACT_SHA256="53cb2d94566eabb35b431f027d5958362f0e3fcd0640efb61593eb06c0068c22"
PARENT_REP_SHA256="b5d81da2ddb4a0bf21370184641b1616693af300458b00117494b5f84f4df52b"
EXPECTED_PARENT_NEGATIVE=92
EXPECTED_PROBLEM_SOURCES=(
"tile_RA122.762_DECp56.240:1372","tile_RA160.223_DECp32.041:1277",
"tile_RA188.616_DECp45.284:633","tile_RA217.768_DECp44.372:762",
"tile_RA264.754_DECp26.740:2753","tile_RA307.409_DECp40.547:701",
"tile_RA41.145_DECp27.430:1698","tile_RA65.091_DECp15.163:1317")
EXPECTED_UNTESTED_PAIRS=100
EXPECTED_UNTESTED_FILES=100
MIN_RETAINED_PSF_FLUX=0.995
SNR_GRID=(8.0,12.0)

def _rows(path):
    return list(csv.DictReader(Path(path).open("r",encoding="utf-8",newline="")))

def derive_problem_sources(parent_receipt):
    by=defaultdict(list)
    for r in parent_receipt["pixel_replay"]["results"]: by[r["src_id"]].append(r)
    problem=[]; neg=0
    for sid,xs in by.items():
        q=sum(x.get("counterpart_test",{}).get("counterpart_present") is False for x in xs); neg+=q
        if q < len(xs): problem.append(sid)
    if neg != EXPECTED_PARENT_NEGATIVE: raise RuntimeError(f"parent qualified-negative count changed: {neg}")
    if tuple(sorted(problem)) != tuple(sorted(EXPECTED_PROBLEM_SOURCES)): raise RuntimeError(f"problem-source set changed: {sorted(problem)}")
    return tuple(sorted(problem))

def build_untested_rows(exact,reps,problem_sources):
    used={(r["src_id"],r["file_name"]) for r in reps}
    out=[r for r in exact if r["src_id"] in problem_sources and (r["src_id"],r["file_name"]) not in used]
    out.sort(key=lambda r:(r["src_id"],r.get("date_obs",""),r["file_name"]))
    if len(out)!=EXPECTED_UNTESTED_PAIRS or len({r["file_name"] for r in out})!=EXPECTED_UNTESTED_FILES: raise RuntimeError("untested recovery cardinality changed")
    return out

def gaussian_retained_fraction(shape,y,x,fwhm_px):
    h,w=shape; sig=float(fwhm_px)/2.354820045
    def cdf(z): return 0.5*(1.0+math.erf(z/math.sqrt(2.0)))
    return float((cdf((w-0.5-x)/sig)-cdf((-0.5-x)/sig))*(cdf((h-0.5-y)/sig)-cdf((-0.5-y)/sig)))

def edge_aware_local_injection(image,x,y,fwhm_px):
    a=np.asarray(image,dtype=float); med,sigma=robust_background(a); fwhm=max(2.0,min(8.0,float(fwhm_px)))
    retained=gaussian_retained_fraction(a.shape,y,x,fwhm)
    if retained < MIN_RETAINED_PSF_FLUX: return {"passed":False,"reason":"EDGE_PSF_FLUX_RETENTION_BELOW_THRESHOLD","retained_psf_flux_fraction":retained}
    trials=[]
    for snr in SNR_GRID:
        z=np.array(a,copy=True); _inject_gaussian(z,y,x,fwhm,snr*sigma); trials.append({"snr":snr,"recovered":measure_psf_at(z,y,x) is not None})
    return {"passed":all(t["recovered"] for t in trials),"background_median":med,"background_sigma":sigma,"injection_fwhm_px":fwhm,"retained_psf_flux_fraction":retained,"minimum_retained_psf_flux_fraction":MIN_RETAINED_PSF_FLUX,"trials":trials}

def _hdu_binding_ok(row,meta):
    return int(row["exact_hdu"])==int(meta["selected_hdu"]) and [int(row["exact_naxis2"]),int(row["exact_naxis1"])]==list(meta["native_shape"])

def replay_row(row,image,meta,seed):
    if not _hdu_binding_ok(row,meta): return {"status":"BLOCKED_HDU_BINDING"}
    gate=psf_relative_injection_recovery_gate(image,seed=seed)
    if not gate.get("passed"): return {"status":"BLOCKED_BY_OVERLAP_FRAME_R1_GATE","overlap_frame_r1_gate":gate}
    x=float(row["exact_x"]); y=float(row["exact_y"]); target=measure_psf_at(image,y,x)
    if target is not None: return {"status":"COUNTERPART_CANDIDATE","overlap_frame_r1_gate":gate,"counterpart_test":counterpart_with_matched_controls(image,x,y)}
    strict=local_coordinate_injection(image,x,y,gate["injection_base_fwhm_px"])
    if strict.get("passed"): return {"status":"QUALIFIED_NO_COUNTERPART_INHERITED_R1","overlap_frame_r1_gate":gate,"strict_local":strict}
    edge=None
    if strict.get("reason")=="COORDINATE_TOO_CLOSE_TO_EDGE":
        edge=edge_aware_local_injection(image,x,y,gate["injection_base_fwhm_px"])
        if edge.get("passed"): return {"status":"QUALIFIED_NO_COUNTERPART_EDGE_AWARE_R1","overlap_frame_r1_gate":gate,"strict_local":strict,"edge_local":edge}
    return {"status":"BLOCKED_LOCAL_SENSITIVITY","overlap_frame_r1_gate":gate,"strict_local":strict,"edge_local":edge}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--parent-artifact-dir",required=True); ap.add_argument("--output-dir",default="results/luci_palomar_2f_f"); ap.add_argument("--cache-dir",default=".cache/luci_palomar_2f_f"); a=ap.parse_args()
    p=Path(a.parent_artifact_dir); out=Path(a.output_dir); out.mkdir(parents=True,exist_ok=True)
    exactp=p/"frozen_exact_fits_wcs_pairs.csv"; repp=p/"frozen_temporal_representatives.csv"; recp=p/"receipt.json"
    if sha256_file(exactp)!=PARENT_EXACT_SHA256 or sha256_file(repp)!=PARENT_REP_SHA256: raise RuntimeError("parent frozen SHA mismatch")
    parent=json.loads(recp.read_text(encoding="utf-8")); problem=derive_problem_sources(parent); exact=_rows(exactp); reps=_rows(repp); untested=build_untested_rows(exact,reps,problem)
    fields=list(untested[0].keys()); targetp=out/"frozen_untested_recovery_pairs.csv"
    with targetp.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(untested)
    target_sha=sha256_file(targetp); events=EventWriter(out/"events.jsonl"); results=[]; cache=Path(a.cache_dir)
    for i,row in enumerate(untested):
        path,_=download_source(row["file_url"],cache,events,target="PALOMAR_2F_F_UNTESTED",filter_name=row.get("filters","UNKNOWN")); image,meta=read_luci_fits_image(path,require_imaging=True,expected_instrument=row.get("instrument")); results.append({**row,"image_meta":meta,"recovery":replay_row(row,image,meta,20262200+i)})
    parent_edge=[]
    for r in parent["pixel_replay"]["results"]:
        if r["src_id"] not in problem or r.get("counterpart_test",{}).get("status")!="NO_SOURCE_RESULT_BLOCKED_BY_LOCAL_SENSITIVITY" or r.get("local_coordinate_injection",{}).get("reason")!="COORDINATE_TOO_CLOSE_TO_EDGE": continue
        path,_=download_source(r["file_url"],cache,events,target="PALOMAR_2F_F_PARENT_EDGE_ONLY",filter_name=r.get("filters","UNKNOWN")); image,meta=read_luci_fits_image(path,require_imaging=True,expected_instrument=r.get("instrument")); edge=edge_aware_local_injection(image,float(r["exact_x"]),float(r["exact_y"]),r["overlap_frame_r1_gate"]["injection_base_fwhm_px"]); parent_edge.append({"src_id":r["src_id"],"file_name":r["file_name"],"edge_local":edge,"hdu_binding_ok":_hdu_binding_ok(r,meta)})
    source_evidence={sid:{"parent_qualified":0,"new_qualified":0,"edge_recovered":0,"counterpart_candidates":0} for sid in problem}
    for r in parent["pixel_replay"]["results"]:
        if r["src_id"] in source_evidence and r.get("counterpart_test",{}).get("counterpart_present") is False: source_evidence[r["src_id"]]["parent_qualified"]+=1
    for r in results:
        s=r["recovery"]["status"]; d=source_evidence[r["src_id"]]
        if s.startswith("QUALIFIED_NO_COUNTERPART"): d["new_qualified"]+=1
        if s=="COUNTERPART_CANDIDATE": d["counterpart_candidates"]+=1
    for r in parent_edge:
        if r["hdu_binding_ok"] and r["edge_local"].get("passed"): source_evidence[r["src_id"]]["edge_recovered"]+=1
    unresolved=[sid for sid,d in source_evidence.items() if d["parent_qualified"]+d["new_qualified"]+d["edge_recovered"]==0 and d["counterpart_candidates"]==0]; candidates=[sid for sid,d in source_evidence.items() if d["counterpart_candidates"]>0]
    status="PASS" if not unresolved and not candidates else "BLOCKED"; scientific_status="ALL_42_SOURCES_HAVE_AT_LEAST_ONE_SENSITIVITY_QUALIFIED_NO_COUNTERPART_EPOCH" if status=="PASS" else "TARGETED_RECOVERY_EXECUTED__SOME_SOURCES_REMAIN_UNRESOLVED_OR_HAVE_CANDIDATES"
    receipt={"schema":"janus.cosmos.luci_palomar.jpfm_2f_f.receipt.v1","experiment_id":"LUCI-PALOMAR-JPFM-2F-F-UNRESOLVED-SENSITIVITY-RECOVERY","status":status,"scientific_status":scientific_status,"parent":{"artifact_id":PARENT_ARTIFACT_ID,"artifact_zip_sha256":PARENT_ZIP_SHA256,"qualified_negative_pairs_preserved_untouched":92},"problem_sources":list(problem),"frozen_untested_recovery":{"pair_count":len(untested),"unique_files":len({r["file_name"] for r in untested}),"sha256":target_sha},"edge_recovery_rule":{"minimum_retained_psf_flux_fraction":MIN_RETAINED_PSF_FLUX,"snr_grid":list(SNR_GRID),"detector":"unchanged measure_psf_at"},"untested_results":results,"parent_edge_recovery":parent_edge,"source_evidence":source_evidence,"unresolved_sources":unresolved,"counterpart_candidate_sources":candidates,"claim_ceiling":"TARGETED_NEAR_IR_SENSITIVITY_RECOVERY_ONLY__NO_ANOMALY_OR_UAP_ORIGIN_CLAIM__NO_CAUSALITY"}
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8"); print(json.dumps({"status":status,"scientific_status":scientific_status,"unresolved_sources":unresolved,"counterpart_candidate_sources":candidates,"source_evidence":source_evidence},indent=2)); return 0 if status=="PASS" else 3
if __name__=="__main__": raise SystemExit(main())
