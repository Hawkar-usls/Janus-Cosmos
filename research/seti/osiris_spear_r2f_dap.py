#!/usr/bin/env python3
"""OSIRIS R2F: bounded DAP public file index and optional FITS header test.

Never fetches full radio power, never attempts authorization bypass, never
turns a catalog center into actual antenna pointing. Frozen P360 file only.
"""
import hashlib,json,pathlib,requests
from urllib.parse import urlparse
P=pathlib.Path("research/seti/osiris_spear_r2f_dap_protocol.json")
S=json.loads(P.read_text())
O=pathlib.Path("runtime/osiris_spear_r2f");O.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-OSIRIS-R2F-readonly-metadata-public/1.0","Accept":"application/json"}
def sha(x):return hashlib.sha256(x).hexdigest()
def save(name,obj):(O/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def read_limited(url,limit,headers,method="GET"):
    with requests.request(method,url,headers=headers,stream=True,timeout=(8,S["limits"]["request_timeout_seconds"]),allow_redirects=True) as r:
        data=bytearray()
        if method=="GET":
            for c in r.iter_content(chunk_size=2048):
                data.extend(c)
                if len(data)>limit:break
        info={"status":r.status_code,"final_host":urlparse(r.url).hostname,
          "content_type":r.headers.get("Content-Type"),
          "content_length":r.headers.get("Content-Length"),
          "content_range":r.headers.get("Content-Range"),
          "content_disposition":r.headers.get("Content-Disposition"),
          "bytes_bounded_read":min(len(data),limit),
          "response_over_limit":len(data)>limit,
          "http_redirected":len(r.history)>0}
        return info,bytes(data[:limit])
def parse_cards(bytes_):
    # FITS standard primary header has 80 character ASCII cards; inspect only
    # the opening cards, never decode or search a radio DATA extension.
    if len(bytes_)<80 or not bytes_[:80].startswith(b"SIMPLE"):
        return {"FITS_PRIMARY_HEADER":False,"reason":"No SIMPLE card in bounded bytes"}
    fields=set(S["expected_metadata_fields"])
    fields.update({"SIMPLE","BITPIX","NAXIS","EXTEND","TELESCOP","DATE-OBS"})
    out={}
    for i in range(0,min(len(bytes_),5760)-79,80):
        a=bytes_[i:i+80]
        try:line=a.decode("ascii","replace")
        except Exception:continue
        key=line[:8].strip()
        if key=="END":break
        if key in fields:
            out[key]=line[10:80].split("/")[0].strip()[:120]
    return {"FITS_PRIMARY_HEADER":True,"cards":out}
def main():
    index=S["official_api_route"]
    state={"stage":"OSIRIS_SPEAR_R2F","collection_id":19838,
      "selected":S["selection"],"index":{"state":"UNAVAILABLE"},
      "raw_file_header":{"state":"NOT_ATTEMPTED"},
      "authentic_full_file_checksum_verified":False,
      "native_radiometric_power_read":False}
    try:
        resp,content=read_limited(index,S["limits"]["max_index_bytes"],UA)
        t=content.decode("utf-8","replace")
        fname=S["selection"]["filename"]
        match=fname in t
        state["index"].update({"state":"INDEX_READ" if resp["status"]==200 and not resp["response_over_limit"] else "UNAVAILABLE_OR_INCOMPLETE",
          "http":resp,"contains_exact_filename":match,"first_4kb_sha256":sha(content[:4096]),
          "body_bytes_read":len(content)})
        try:
            parsed=json.loads(t)
            state["index"]["top_level_data_shape"]=list(parsed)[:25] if isinstance(parsed,dict) else f"list[{len(parsed)}]" if isinstance(parsed,list) else type(parsed).__name__
        except Exception:state["index"]["parse"]="NON_JSON"
        if resp["status"]==200 and match and not resp["response_over_limit"]:
            # Public filename is an immutable lock from DAP TAP. Official API
            # single-file endpoint returns a raw download or a 401/403; do not
            # obtain signed URLs, private cookies or unapproved tokens.
            durl="https://data.csiro.au/ws/v2/collections/19838/downloadfile?filename="+fname
            hr,hb=read_limited(durl,0,{"User-Agent":UA["User-Agent"]},method="HEAD")
            state["raw_file_header"]={"state":"HEAD_DONE","head":hr}
            if hr["status"] in (200,206):
                response,data=read_limited(durl,S["limits"]["max_fits_header_bytes"],
                    {"User-Agent":UA["User-Agent"],"Range":S["limits"]["range_header"]})
                state["raw_file_header"].update({"state":"BOUNDED_PUBLIC_HEADER_READ","range_get":response,
                     "prefix_sha256":sha(data),"fits_header":parse_cards(data)})
            else:state["raw_file_header"]["state"]="NO_PUBLIC_FILE_ACCESS_VERIFIED"
    except Exception as ex:
        state["index"]["errors"]=[type(ex).__name__+":"+str(ex)[:300]]
    save("dap_file_access.json",state)
    report={"stage":"OSIRIS_SPEAR_R2F",
      "status":"PUBLIC_ORIGINAL_PSRFITS_HEADER_ACCESSED__HOLD_FULL_BYTES_TCS_ONOFF" if state["raw_file_header"].get("fits_header",{}).get("FITS_PRIMARY_HEADER") else
          "DAP_PUBLIC_FILE_ACCESS_NOT_ESTABLISHED__HOLD_ORIGINAL_SPECTRUM_AND_TCS",
      "catalog_raw_file_found_prior":True,
      "public_data_index_state":state["index"]["state"],
      "index_contains_filename":state["index"].get("contains_exact_filename"),
      "header_state":state["raw_file_header"]["state"],
      "FITS_primary_header_verified":state["raw_file_header"].get("fits_header",{}).get("FITS_PRIMARY_HEADER",False),
      "original_full_file_sha256_verified":False,
      "actual_antenna_pointing_checked":False,
      "matched_OFF_found":False,"native_Hz_resolution_absent":True,
      "raw_spectral_power_fetched":False,"source_detection_claim":False,
      "production_changed":False,
      "next":"Obtain official public download or authorized scientific file access and original antenna per-beam TCS. Coarse 3 MHz channels cannot evaluate frozen Hz-sensitive JANUS even if target lies in main beam."}
    save("r2f_gate_report.json",report)
    print("OSIRIS_R2F_RESULT",json.dumps(report),flush=True)
if __name__=="__main__":main()
