#!/usr/bin/env python3
"""Historical natural galaxy mask provenance and independent Tsys arithmetic."""
import hashlib,importlib.metadata,inspect,json,pathlib,subprocess,sys
import numpy as np
from astropy.io import fits
from astropy.utils import iers
from astropy.coordinates import solar_system_ephemeris
from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.spectra import core,scan
P=pathlib.Path("research/seti/r37f3d1bb2a_vegas_protocol.json")
S=json.loads(P.read_text());M=S["method"];O=pathlib.Path("runtime/r37f3d1bb2a");O.mkdir(parents=True,exist_ok=True)
for n,h in S["pinned_blobs"].items():
    b=(P.parent/n).read_bytes()
    assert hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==h,n
import r37f3d1bb1_mask as b1
prev=b1.prev
def assertclose(x,y,atol,rtol,label):
    a=np.asarray(x,float);b=np.asarray(y,float)
    assert a.shape==b.shape and a.size and np.isfinite(a).all() and np.isfinite(b).all(),label
    err=float(np.max(np.abs(a-b)))
    assert np.allclose(a,b,atol=atol,rtol=rtol),f"{label} max abs error {err}"
    return err
def measure(sdf,p,q,row):
    block=sdf.getps(scan=152,ifnum=M["ifnum"],plnum=q,fdnum=0)
    assert len(block)==1
    ps=block[0];assert ps.nint==M["nint"] and ps.nchan==M["nchan"]
    assert not ps._nocal and ps._smoothref==1 and ps.tscale.lower()=="ta"
    defs=[("sigcalon","sigonrows",152,"T"),("sigcaloff","sigoffrows",152,"F"),("refcalon","refonrows",153,"T"),("refcaloff","refoffrows",153,"F")]
    arrays=[];flags=[];masks=[];predicted=[];exp=[];df=[];dates=[];tcal=[];raw_meta=[]
    with fits.open(p,memmap=True) as hd:
        d=next(h.data for h in hd if isinstance(h,fits.BinTableHDU) and h.header.get("EXTNAME")=="SINGLE DISH")
        row["fits_flags_column"]="FLAGS" in d.names
        for a,r,s,phase in defs:
            ix=np.array(getattr(ps,"_"+r),dtype=int)
            assert len(ix)==M["nint"] and len(set(ix.tolist()))==len(ix)
            assert np.all(d["SCAN"][ix]==s) and np.all(d["IFNUM"][ix]==M["ifnum"]) and np.all(d["PLNUM"][ix]==q) and np.all(d["FDNUM"][ix]==0)
            assert set(np.char.strip(d["CAL"][ix].astype(str)))=={phase}
            dates.append(np.char.strip(d["DATE-OBS"][ix].astype(str)))
            raw=np.stack([np.asarray(d["DATA"][i],float).ravel() for i in ix])
            assert raw.shape==(M["nint"],M["nchan"]) and np.isfinite(raw).all()
            fl=np.stack([np.asarray(d["FLAGS"][i],bool).ravel() for i in ix]) if row["fits_flags_column"] else np.zeros(raw.shape,bool)
            assert fl.shape==raw.shape
            lib=np.ma.asarray(getattr(ps,"_"+a))
            assert lib.shape==raw.shape
            mask=np.ma.getmaskarray(lib)
            assertclose(np.asarray(lib.data)[~mask],raw[~mask],0,0,"raw archive read")
            arrays.append(raw);masks.append(mask);flags.append(fl)
            assert set(np.char.strip(d["BACKEND"][ix].astype(str)))=={"VEGAS"}
            v,pixel,step=[np.asarray(d[k][ix],float) for k in ("VSPRVAL","VSPRPIX","VSPDELT")]
            assert np.isfinite(v).all() and np.isfinite(pixel).all() and np.isfinite(step).all()
            pred=np.zeros_like(mask)
            for k in range(len(ix)):
                for j in range(32):
                    if j==16:continue
                    ch=int(j*step[k]-v[k]*step[k]+pixel[k]-1)
                    if 0<=ch<M["nchan"]:pred[k,ch]=True
            predicted.append(pred)
            exp.append(np.asarray(d["EXPOSURE"][ix],float));df.append(np.asarray(d["CDELT1"][ix],float))
            tcal.append(np.asarray(d["TCAL"][ix],float))
            raw_meta.append({k:np.asarray(d[k][ix],float) for k in ("CRVAL1","CRPIX1","CDELT1")})
    assert np.array_equal(dates[0],dates[1]) and np.array_equal(dates[2],dates[3])
    assert len(set(dates[0]))==len(set(dates[2]))==M["nint"]
    row["input_group_flags"]=[{"group":a,"fits_flagged":int(f.sum()),"dysh_masked":int(m.sum()),"exact_match":bool(np.array_equal(f,m))} for (a,*_),f,m in zip(defs,flags,masks)]
    expected_input_mask=np.logical_or.reduce(masks)
    archive_union=np.logical_or.reduce(flags)
    sig=0.5*(arrays[0]+arrays[1]);ref=0.5*(arrays[2]+arrays[3])
    denom_bad=(ref==0)|~np.isfinite(ref)
    y=[];g=[];axes=[]
    for i in range(M["nint"]):
        sp=ps.getspec(i,use_wcs=True);v,f,good,_=b1.unpack(sp)
        y.append(v);axes.append(f);g.append(good)
        assert prev.text(sp.meta["DATE-OBS"])==dates[1][i]
    y=np.stack(y);g=np.stack(g);axes=np.stack(axes)
    invalid=~g;expected=expected_input_mask|denom_bad
    row.update(calibrated_invalid_count=int(invalid.sum()),input_mask_union_count=int(expected_input_mask.sum()),archive_flag_union_count=int(archive_union.sum()),zero_or_nonfinite_reference_count=int(denom_bad.sum()),calibrated_mask_matches_inputs_and_denominator=bool(np.array_equal(invalid,expected)),fits_flags_match_all_inputs=all(x["exact_match"] for x in row["input_group_flags"]))
    vegas_union=np.logical_or.reduce(predicted)
    row["vegas_group_exact_matches"]=[bool(np.array_equal(p,m)) for p,m in zip(predicted,masks)]
    row["vegas_predicted_union_count"]=int(vegas_union.sum())
    if not archive_union.any() and not denom_bad.any() and invalid.any() and all(row["vegas_group_exact_matches"]) and np.array_equal(invalid,vegas_union):
        origin="VEGAS_METADATA_RULE_PROPAGATED"
    elif np.array_equal(invalid,archive_union) and row["fits_flags_match_all_inputs"] and not denom_bad.any() and invalid.any():
        origin="ARCHIVE_FLAGS_PROPAGATED_THROUGH_CALIBRATION"
    elif np.array_equal(invalid,denom_bad) and not expected_input_mask.any() and denom_bad.any():
        origin="INVALID_REFERENCE_DENOMINATOR"
    else:origin="UNRESOLVED_OR_MIXED"
    row["mask_origin"]=origin
    row["physical_reason_for_archive_flag_known"]=False
    row["software_mask_origin_resolved"]=origin!="UNRESOLVED_OR_MIXED"
    assert origin!="UNRESOLVED_OR_MIXED","mask origin not uniquely established"
    # Independent noise diode system temperature using the documented central range.
    tc=np.unique(np.concatenate(tcal[2:]));assert len(tc)==1 and tc[0]>0
    edge=int(M["nchan"]*0.1);sl=slice(edge,-(edge-1))
    ts=[]
    for i in range(M["nint"]):
        off=arrays[3][i,sl];on=arrays[2][i,sl]
        # Distinct mean supports follow mean_tsys mode=0.
        mo=~masks[3][i,sl];md=mo&~masks[2][i,sl]
        assert mo.any() and md.any()
        diff=np.mean((on-off)[md]);assert diff!=0
        ts.append(float(tc[0]*(np.mean(off[mo])/diff+0.5)))
    ts=np.asarray(ts)
    assert np.isfinite(ts).all() and (ts>0).all()
    row["tsys_max_abs_error_K"]=assertclose(ps.tsys,ts,1e-9,1e-9,"independent Tsys")
    row["tsys_min_K"]=float(ts.min());row["tsys_max_K"]=float(ts.max())
    sig_exp=exp[0]+exp[1];ref_exp=exp[2]+exp[3]
    effective=sig_exp*ref_exp/(sig_exp+ref_exp)
    width=np.abs(sum(df)/4.0)
    assert (effective>0).all() and (width>0).all()
    row["exposure_max_abs_error_s"]=assertclose(ps.exposure,effective,1e-10,1e-10,"effective exposure")
    row["channel_width_max_abs_error_hz"]=assertclose(np.abs(ps.delta_freq),width,1e-10,1e-10,"channel width")
    w=effective*width/(ts*ts);assert (w>0).all() and np.isfinite(w).all()
    row["weight_max_abs_error"]=assertclose(ps.tsys_weight,w,1e-10,M["weight_rtol"],"Tsys weights")
    row["weight_min"]=float(w.min());row["weight_max"]=float(w.max())
    # Only legitimate, unmasked measured values enter arithmetic. No invalid values filled.
    rr=np.broadcast_to(ts[:,None],y.shape)
    manual=rr[g]*(sig[g]-ref[g])/ref[g]
    row["independent_calibration_max_abs_K"]=assertclose(y[g],manual,M["calibration_atol_K"],M["calibration_rtol"],"per-integration calibration")
    row["calibration_valid_samples"]=int(g.sum())
    fully=g.all(0);anyvalid=g.any(0);assert fully.sum()==32767 and np.array_equal(fully,anyvalid)
    before=np.ma.array(ps._calibrated,copy=True)
    avg,axis,agood,_=b1.unpack(ps.timeaverage(weights="tsys",use_wcs=True))
    avg2,axis2,agood2,_=b1.unpack(block.timeaverage(weights="tsys",use_wcs=True))
    assert np.array_equal(agood,fully) and np.array_equal(agood2,fully)
    truth=np.sum(y[:,fully]*w[:,None],axis=0)/np.sum(w)
    row["weighted_scan_max_abs_K"]=assertclose(avg[fully],truth,M["average_atol_K"],M["average_rtol"],"weighted scan")
    row["weighted_block_max_abs_K"]=assertclose(avg2[fully],truth,M["average_atol_K"],M["average_rtol"],"weighted block")
    row["weighted_vs_uniform_rms_K"]=b1.rms(truth-y[:,fully].mean(0))
    row["average_first_axis_max_error_hz"]=assertclose(axis,axes[0],.01,0,"first WCS")
    assertclose(axis2,axes[0],.01,0,"block WCS")
    m=raw_meta[1];grid=m["CRVAL1"][:,None]+(np.arange(M["nchan"])[None,:]+1-m["CRPIX1"][:,None])*m["CDELT1"][:,None]
    row["per_integration_wcs_max_error_hz"]=assertclose(axes,grid,.01,0,"native per-integration WCS")
    after=np.ma.asarray(ps._calibrated)
    row["calibrated_integrations_unchanged"]=bool(np.array_equal(np.ma.getmaskarray(before),np.ma.getmaskarray(after)) and np.array_equal(before.data[g],after.data[g]))
    assert row["calibrated_integrations_unchanged"]
    row.update(valid_average_channels=int(fully.sum()),masked_average_channels=int((~fully).sum()),imputed_samples=0,interpolated_samples=0,official_tsys_arithmetic_pass=True)
def main():
    iers.conf.auto_download=False;iers.conf.iers_degraded_accuracy="error";solar_system_ephemeris.set("builtin")
    p=O/"original.fits";d={"status":"HOLD_TSYS_AUDIT","source":"NGC2415","production_adoption":False,"detector_tested":False,"lossless_average_proven":False,"absolute_astrophysical_accuracy_proven":False,"full_five_window_calibration":False,"rows":[]}
    try:
        assert importlib.metadata.version("dysh")=="1.0.1"
        d["code_receipts"]={n:hashlib.sha256(inspect.getsource(inspect.unwrap(fn)).encode()).hexdigest() for n,fn in (("calibrate",scan.PSScan.calibrate),("average",scan.PSScan.timeaverage),("tsys",core.mean_tsys),("weights",core.tsys_weight))}
        d["integrity"]=prev.meta.prev.download(p)
        assert d["integrity"]["sha256"]==prev.S["archive_sha256"] and d["integrity"]["bytes"]==prev.S["archive_bytes"]
        p.chmod(0o444);m={};prev.meta.inspect(p,m);assert len(m["gates"])==14 and all(m["gates"].values())
        d["spatial_off_verified"]=True
        sdf=GBTFITSLoad(p)
        for q in M["pols"]:
            row={"pol":q,"status":"HOLD"}
            try:measure(sdf,p,q,row);row["status"]="PASS_TSYS_ARITHMETIC_SUPPORT"
            except Exception as e:
                import traceback
                traceback.print_exc();row["error"]=type(e).__name__+": "+str(e)[:500]
            d["rows"].append(row);print("R37F3D1BB2A_POL",json.dumps(row),flush=True)
        d["original_unchanged"]=prev.hashfile(p)==prev.S["archive_sha256"]
        d["status"]="PASS_ISOLATED_TSYS_CALIBRATION_AUDIT" if d["original_unchanged"] and len(d["rows"])==2 and all(r["status"]=="PASS_TSYS_ARITHMETIC_SUPPORT" for r in d["rows"]) else "HOLD_TSYS_AUDIT"
    except Exception as e:
        import traceback;traceback.print_exc();d["error"]=type(e).__name__+": "+str(e)[:500]
    finally:p.unlink(missing_ok=True)
    (O/"result.json").write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    (O/"pip-freeze.txt").write_text(subprocess.check_output([sys.executable,"-m","pip","freeze"],text=True))
    print("R37F3D1BB2A_SUMMARY",json.dumps(d,allow_nan=False),flush=True)
if __name__=="__main__":main()
