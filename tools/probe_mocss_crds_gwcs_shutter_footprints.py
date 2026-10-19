#!/usr/bin/env python3
"""Original-MSA-selected small WCS detector illumination footprint proof.

Stage 1: SHA-verify original RATE & MSA, select 3 physically OPEN consecutive
shutters along dispersion from ORIGINAL MSA by (MSAMETID,PATT_NUM).
Stage 2: construct a TEST copy containing nine independently labeled
BACKGROUND shutter entries (3x3 around source), preserving exact original
shutter coordinates; run PINNED official JWST AssignWcsStep CRDS.
Stage 3: for each selected slit, ask official NIRSpec GWCS for a
slit_frame -> detector transform at 3 frozen physical wavelengths.
Record 9-point rectangle corners/center and bounding envelope; never
substitute an optical envelope for an illumination-weight map.

This is a CALIBRATION FEASIBILITY PROBE, not a flux/extraction pipeline.
Original RATE/MSA files are not modified. Source original MPT identity
is not asserted. Slitlet naming by copied MSA is a diagnostic device.
"""
import argparse, hashlib, json, math, os, traceback
from pathlib import Path
from astropy.io import fits
from fetch_jwst_mocss_mast_subset import read_and_store, BASE, ROOT_URI
from prepare_mocss_wcs_probe import create_probe
from urllib.parse import quote
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"domains/space/MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json"
MSAS={
"A":("jw03290001001_05_msa.fits","fb6239b287d25c478d2e22737ff956850a9bc1e0f365500cfc5661abbfe17a44"),
"B":("jw03290001001_05_msa.fits","fb6239b287d25c478d2e22737ff956850a9bc1e0f365500cfc5661abbfe17a44"),
"C":("jw03290002001_01_msa.fits","25a061a226521a8bbf2f30d3f6922a22fb13437800047f12c8511e354820c91e"),
"D":("jw03290002001_01_msa.fits","25a061a226521a8bbf2f30d3f6922a22fb13437800047f12c8511e354820c91e")
}
WAVES={"Pa_gamma":13656.7483309,"OIII_4960":15850.7932988,"primary":16004.}
FRAC=[-.5,0,.5]

def hash_file(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for c in iter(lambda:f.read(1048576),b""):h.update(c)
    return h.hexdigest()

def fetch(originaldir,name,expected):
    p=originaldir/name
    if not p.exists():
        read_and_store(BASE+quote(ROOT_URI+name,safe=""),p)
    assert hash_file(p)==expected,(name,"ARCHIVE_SHA_FAIL")
    return p

def finite_transform(f,ww):
    # Official WCS slit_frame to detector uses wavelengths IN METRES.
    XX,YY=np.meshgrid(FRAC,FRAC)
    xx,yy=f(XX,YY,np.ones_like(XX)*ww*1e-10)
    xx=np.asarray(xx,dtype=float);yy=np.asarray(yy,dtype=float)
    mask=np.isfinite(xx)&np.isfinite(yy)
    if not np.all(mask):return {"status":"NONFINITE_TRANSFORM_POINTS",
                                "valid":int(mask.sum()),"points":9}
    bbox={"x_min":float(xx.min()),"x_max":float(xx.max()),
          "y_min":float(yy.min()),"y_max":float(yy.max())}
    return {"status":"FINITE_GWCS_SAMPLE","detector_box":bbox,
            "center_x_px":float(xx[1,1]),"center_y_px":float(yy[1,1]),
            "grid_detector_x_px":xx.tolist(),"grid_detector_y_px":yy.tolist(),
            "on_detector_center":bool(0<=xx[1,1]<2048 and 0<=yy[1,1]<2048)}

def process(entry,dir):
    from jwst.assign_wcs import AssignWcsStep,nrs_wcs_set_input
    group=entry["group"];det=entry["expected_science_detector_from_WCS"]
    orig=entry["detectors"][det]
    rate=fetch(dir/"original",orig["rate_file"],orig["sha256"])
    msa_name,msa_sha=MSAS[group]
    msa=fetch(dir/"original",msa_name,msa_sha)
    proposal=orig["selected_virtual_shutter"]
    testdir=dir/("WCS_"+group+"_"+str(entry["exposure_index"]))
    probe=create_probe(rate,msa,testdir,grid=1,
        target={key:proposal[key] for key in ("quadrant","row","column")},half_width=1)
    with fits.open(probe["probe_msa"],memmap=False) as h:
        data=h["SHUTTER_INFO"].data
        desired={(proposal["quadrant"],proposal["row"]+d,proposal["column"]) for d in (-1,0,1)}
        found={}
        for i,z in enumerate(data):
            pos=(int(z["shutter_quadrant"]),int(z["shutter_row"]),int(z["shutter_column"]))
            if pos in desired:
                assert str(z["shutter_state"]).strip()=="OPEN"
                found[pos]=str(int(z["slitlet_id"]))
        assert set(found)==desired,(found,desired)
    result={
        "group":group,"exposure_index":entry["exposure_index"],
        "original_rate_file":rate.name,"original_rate_sha256":hash_file(rate),
        "original_msa_file":msa.name,"original_msa_sha256":hash_file(msa),
        "original_signed_MSA_id":int(entry["source_msametid_signed16"]),
        "dither_PATT_NUM":int(probe["dither"]),
        "test_virtual_MSA_sha256":hash_file(Path(probe["probe_msa"])),
        "original_source_nearest_shutter":proposal,
        "three_rows_open_in_original_MSA":True,
        "test_source_id":0,"original_source_MPT_identity_claim":False,
        "stage2_absolute_flux":None,"spectral_background_calibrated":False,
        "samples_by_shutter":[]}
    try:
        # assign_wcs alone, NOT an all-rate calwebb_spec2 complete reduction
        cal=AssignWcsStep.call(probe["probe_rate"],save_results=False)
        gwcs=cal.meta.wcs
        assert gwcs is not None
        opens=gwcs.get_transform("gwa","slit_frame").slits
        result["official_assign_wcs_cal_step"]=str(cal.meta.cal_step.assign_wcs)
        result["official_gwcs_open_slits"]=len(opens)
        slitnames={str(s.name) for s in opens}
        for coord,name in sorted(found.items(),key=lambda x:x[0][1]):
            out={"quadrant":coord[0],"row":coord[1],"column":coord[2],
                 "test_virtual_slit_name":name,"wcs_slit_exists":name in slitnames,
                 "lines":{}}
            if name in slitnames:
                sw=nrs_wcs_set_input(cal,name)
                mapper=sw.get_transform("slit_frame","detector")
                for label,w in WAVES.items():
                    try:out["lines"][label]=finite_transform(mapper,w)
                    except Exception as exc:
                        out["lines"][label]={"status":"EXCEPTION","exception":type(exc).__name__,"details":str(exc)[:600]}
            result["samples_by_shutter"].append(out)
        result["status"]="OFFICIAL_CRDS_ASSIGN_WCS_CALL_RETURNED"
        cal.close()
    except Exception as exc:
        result["status"]="WCS_EXECUTION_FAILURE"
        result["error_type"]=type(exc).__name__
        result["error"]=str(exc)[:1000]
        result["traceback_tail"]=traceback.format_exc()[-3000:]
    print("WCS",group,result["status"],
          [(q["row"],q["wcs_slit_exists"]) for q in result["samples_by_shutter"]],flush=True)
    return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--data-dir",type=Path,default=Path("/tmp/janus_mocss_crds_open_shutter_footprints"))
    ap.add_argument("--output",type=Path,default=Path("/tmp/janus_mocss_crds_open_shutter_footprints.json"))
    ap.add_argument("--indices",default="1,8,10,13")
    a=ap.parse_args()
    import jwst
    assert jwst.__version__=="1.19.1",jwst.__version__
    assert os.environ.get("CRDS_CONTEXT")=="jwst_1413.pmap"
    a.data_dir.mkdir(parents=True,exist_ok=True);a.output.parent.mkdir(parents=True,exist_ok=True)
    all_e=json.loads(RECEIPT.read_text())["exposures"]
    selected=[int(x) for x in a.indices.split(",")]
    assert len(selected)==len(set(selected)) and all(1<=i<=16 for i in selected)
    obj={"schema":"janus-mocss-original-MSA-neighbor-CRDS-GWCS-pixel-footprint-pilot-v1",
         "jwst_version":jwst.__version__,"crds_context":"jwst_1413.pmap",
         "wavelengths_A":WAVES,
         "test_MSA_inputs_are_synthetic_background_source_zero":True,
         "real_MSA_shutter_coordinates_and_open_states":True,
         "test_subsampled_9_slit_WCS_not_all_original_open_shutters":True,
         "exact_full_detector_illumination_for_all_open_shutters":False,
         "STPSF_weighted_detector_pixel_throughput":False,
         "calibrated_physical_emission_line_flux":None,
         "source_identity":None,"redshift":None,"gas_trail":None,
         "records":[]}
    a.output.write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")
    for i in selected:
        row=process(all_e[i-1],a.data_dir)
        obj["records"].append(row)
        a.output.write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")
    if any(q["status"]=="WCS_EXECUTION_FAILURE" for q in obj["records"]):
        raise SystemExit("Official JWST assign_wcs failed; consult artifact and fix actual blocker")
    print("GWCS_APERTURE_WAVELENGTH_FOOTPRINT_ONLY; FULL_DETECTOR_ILLUMINATION_UNCALIBRATED",flush=True)
if __name__=="__main__":main()
