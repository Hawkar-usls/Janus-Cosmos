#!/usr/bin/env python3
"""R42A7B2: authentic *archive* ON/OFF bytes/headers; never pretend target catalog
coordinates are an independently corroborated telescope pointing. No model run.
"""
import hashlib,json,math,pathlib,requests,os
import numpy as np
import h5py
try:
 import hdf5plugin  # noqa: F401 -- registers compression filters for BL archives
except ImportError:
 pass
P=json.loads(pathlib.Path("research/seti/r42a7b2_frozen_raw_pair.json").read_text())
OUT=pathlib.Path("runtime/r42a7b2");OUT.mkdir(parents=True,exist_ok=True)
TMP=OUT/"private_originals";TMP.mkdir(parents=True,exist_ok=True)
def save(name,data):(OUT/name).write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def sha(b):return hashlib.sha256(b).hexdigest()
def value(v):
 if isinstance(v,bytes):return v.decode("utf8","replace")
 if isinstance(v,np.ndarray):
  if v.size>32:return {"array_shape":list(v.shape),"dtype":str(v.dtype)}
  return [value(q) for q in v.tolist()]
 if isinstance(v,np.generic):return value(v.item())
 if isinstance(v,(float,int,str,bool,type(None))):
  return v if not isinstance(v,float) or math.isfinite(v) else str(v)
 return str(v)
def download(spec):
 dst=TMP/(spec["role"]+".h5")
 limit=P["limits"]["max_file_bytes"]
 h=hashlib.md5();sha256=hashlib.sha256();count=0
 rec={"role":spec["role"],"target":spec["target"],"archive_url":spec["url"],
       "archive_mjd":spec["mjd"],"archive_ra_deg":spec["archive_ra_deg"],
       "archive_dec_deg":spec["archive_dec_deg"],"expected_md5":spec["md5"],
       "expected_bytes":spec["expected_bytes"],"downloaded":False}
 try:
  with requests.get(spec["url"],stream=True,timeout=(10,45),
        headers={"User-Agent":"JANUS-R42A7B2 raw-readonly experiment/1.0"}) as r:
   r.raise_for_status()
   if r.headers.get("Content-Length") and int(r.headers["Content-Length"])>limit:
    raise ValueError("DECLARED_CONTENT_LENGTH_EXCEEDS_CAP")
   with dst.open("wb") as f:
    for chunk in r.iter_content(chunk_size=2*1024*1024):
     if not chunk:continue
     count+=len(chunk)
     if count>limit:raise ValueError("DOWNLOAD_LIMIT_EXCEEDED")
     h.update(chunk);sha256.update(chunk);f.write(chunk)
  rec.update({"bytes":count,"md5":h.hexdigest(),"sha256":sha256.hexdigest(),
              "size_match":count==spec["expected_bytes"],"md5_match":h.hexdigest()==spec["md5"]})
  if not rec["size_match"] or not rec["md5_match"]:raise ValueError("OFFICIAL_ARCHIVE_BYTE_INTEGRITY_MISMATCH")
  rec["downloaded"]=True
  print("R42A7B_RAW_FILE",json.dumps({k:v for k,v in rec.items() if k!="archive_url"}),flush=True)
 except Exception as e:
  rec["error"]=str(e)[:400]
  try:dst.unlink(missing_ok=True)
  except Exception:pass
  print("R42A7B_RAW_FILE",json.dumps(rec),flush=True)
 return rec,dst if rec["downloaded"] else None
def extract_header(source,path):
 result={"role":source["role"],"status":"RAW_METADATA_UNAVAILABLE"}
 try:
  with h5py.File(path,"r") as f:
   keys=list(f.keys())
   ds=f["data"] if "data" in f else None
   if ds is None:raise ValueError("HDF5_DATASET_DATA_MISSING")
   attrs={str(k):value(v) for k,v in f.attrs.items()}
   d_attrs={str(k):value(v) for k,v in ds.attrs.items()}
   # Not a substitute for original TCS pointing logs.
   result.update({"status":"HDF5_READABLE","root_keys":keys,"shape":list(ds.shape),
     "dtype":str(ds.dtype),"root_attrs":attrs,"data_attrs":d_attrs,
     "native_tcs_pointing_independently_corrob":False})
   lower={k.lower():v for k,v in {**attrs,**d_attrs}.items()}
   src_name=str(lower.get("source_name",lower.get("src_name",lower.get("source","")))).strip()
   result["source_name_native"]=src_name
   for q in ["tstart","tsamp","fch1","foff","nchans","nifs","src_raj","src_dej","ra","dec","telescope_id"]:
    result["native_"+q]=lower.get(q)
   if src_name and src_name!=source["target"]:
    result["status"]="HOLD_SOURCE_NAME_DISAGREES_WITH_ARCHIVE"
   if len(ds.shape)!=3 or ds.shape[0]<16 or ds.shape[-1]<128:
    result["status"]="HOLD_UNEXPECTED_SPECTROGRAM_DIMENSIONS"
   if not isinstance(result.get("native_tstart"),(int,float)) or abs(float(result["native_tstart"])-source["mjd"])>.01:
    result["status"]="HOLD_NATIVE_START_MJD_MISMATCH_OR_MISSING"
   if result["status"]=="HDF5_READABLE":
    ntime=int(ds.shape[0]);nchan=int(ds.shape[-1])
    times=np.linspace(0,ntime-1,16,dtype=int)
    mid=nchan//2;cut=np.zeros((16,128),dtype=np.float32)
    for i,t in enumerate(times):
     v=np.asarray(ds[int(t),0,mid-64:mid+64],dtype=np.float32)
     if v.shape!=(128,) or not np.isfinite(v).all():
      raise ValueError("BAD_FROZEN_CENTER_SPECTRAL_TILE")
     cut[i]=v
    result["center_tile"]={"time_indices":times.tolist(),"frequency_slice_channels":[mid-64,mid+64],
      "array_sha256":sha(cut.tobytes()),"shape":[16,128],"selected_without_signal_search":True}
    np.save(OUT/("center_raw_"+source["role"].lower()+".npy"),cut,allow_pickle=False)
 except Exception as e:result["status"]="HOLD_HDF5_UNREADABLE";result["error"]=str(e)[:400]
 print("R42A7B_RAW_HEADER",json.dumps({"role":result["role"],"status":result["status"],
   "shape":result.get("shape"),"name":result.get("source_name_native"),
   "mjd":result.get("native_tstart"),"native_ra":result.get("native_src_raj"),
   "native_dec":result.get("native_src_dej")}),flush=True)
 return result
def main():
 results=[];headers=[]
 for input_spec in P["inputs"]:
  rec,dst=download(input_spec)
  results.append(rec)
  if dst:
   head=extract_header(input_spec,dst);headers.append(head)
   # Original archive bytes are for validation only; do not export massive raw files in artifact.
   dst.unlink(missing_ok=True)
 save("raw_download_receipt.json",{"sources":results,"original_file_count_verified":sum(x["downloaded"] for x in results)})
 save("native_header_audit.json",{"headers":headers,"native_tcs_corroborated":False,
   "pointing_note":"BL HDF5 SRC_NAME/SRC_RAJ/SRC_DEJ are original file header metadata, but may reflect source catalog rather than signed telescope TCS/encoder pointing. Need independent logs."})
 ok=len(results)==2 and all(x["downloaded"] for x in results)
 native=ok and len(headers)==2 and all(x["status"]=="HDF5_READABLE" for x in headers)
 result={"stage":"R42A7B2","status":"PASS_RAW_ARCHIVE_INTEGRITY__HOLD_INDEPENDENT_POINTING_AND_CUSTODIAN" if native else
     "PARTIAL_OR_FAILED_RAW_ARCHIVE_ADMISSION__HOLD",
     "source_urls_frozen_before_download":True,
     "official_archive_md5_verified":sum(x.get("md5_match",False) for x in results),
     "original_spectral_files_readable":sum(x["status"]=="HDF5_READABLE" for x in headers),
     "bounded_original_spectral_tiles_extracted":sum(bool(x.get("center_tile")) for x in headers),
     "native_header_source_and_time_supported":native,
     "independently_corrob_telescope_pointing":False,
     "external_blind_custodian":False,
     "independent_exam_admitted":False,
     "detector_score_comparison_run":False,
     "onoff_separation_from_archive_listing_degrees":P["listed_spatial_separation_degrees"],
     "production_modifications":False,"actual_false_alert_reduction_claim":False,
     "next":"Independent original TCS/encoder pointing and non-overlap audit, then blinded full-native-spectra equal-information evaluation. Original HDF5 data remain at frozen official archive URLs."}
 save("gate_report.json",result)
 print("R42A7B2_RESULT",json.dumps(result),flush=True)
if __name__=="__main__":main()
