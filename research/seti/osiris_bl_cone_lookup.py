"""OSIRIS SPEAR R1C: official BL RA/Dec archive HDF5 metadata coverage.
Read-only and externally provided archive labels, not telescope pointing.
"""
import concurrent.futures,hashlib,json,math,pathlib,re,urllib.parse,requests
from bs4 import BeautifulSoup
import osiris_spear_radio_coverage as src
P=json.loads(pathlib.Path("research/seti/osiris_bl_cone_protocol.json").read_text())
O=pathlib.Path("runtime/osiris_spear_bl_r1c");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,v):(O/n).write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def one(job):
    target,scope=job
    ra,dec=target["ra_deg"],target["dec_deg"]
    half=.5/max(.2,math.cos(math.radians(dec)))
    args={"project":scope,"file_type":"HDF5","ra":f"{ra:.9f}",
      "ra_range":f"{half:.9f}","decl":f"{dec:.9f}",
      "decl_range":"0.5","perPage":"100","page":"1","search":"Search"}
    url=P["archive"]+"?"+urllib.parse.urlencode(args)
    result={"id":target["id"],"base_id":target["base_id"],"kind":target["kind"],
       "telescope":scope,"ra_deg":ra,"dec_deg":dec,"status":"UNAVAILABLE",
       "query_string_sha256":sha(url.encode()),"only_first_100_rows":True,
       "original_telescope_encoder_pointing_verified":False,"original_spectra_verified":False}
    try:
        r=requests.get(url,headers={"User-Agent":"JANUS-SPEAR-BL-RADIO-CONE/1.0"},
          timeout=(8,24))
        r.raise_for_status()
        if len(r.content)>P["response_bytes_max"]:raise ValueError("EXCEEDED_HTML_BYTE_CAP")
        soup=BeautifulSoup(r.content,"html.parser")
        rows=[];total=0;filter_violation=0
        for tr in soup.find_all("tr"):
            cells=tr.find_all("td",recursive=False)
            if len(cells)<10:continue
            vals=[c.get_text(" ",strip=True) for c in cells]
            if not re.fullmatch(r"20\d\d-\d\d-\d\d.*",vals[0]):continue
            try:
                ra0=float(vals[4]);dec0=float(vals[5])
                mjd=float(vals[1]);freq=float(vals[6])
            except (IndexError,TypeError,ValueError):continue
            total+=1
            if not (0<=ra0<360 and -90<=dec0<=90):continue
            sep=src.angular(ra,dec,ra0,dec0)
            if sep>.55:filter_violation+=1
            if sep>.50001:continue
            links=[urllib.parse.urljoin(r.url,a["href"]) for a in tr.select("a[href]") if a.get("href")]
            link=next((x for x in links if x.startswith("https://") and re.search(r"\.(h5|hdf5)(\?|$)",x,re.I)),None)
            md5=next((v.lower() for v in vals[8:] if re.fullmatch("[a-f0-9]{32}",v,re.I)),None)
            rows.append({"ra_deg_catalog":ra0,"dec_deg_catalog":dec0,"separation_deg":sep,
             "utc":vals[0],"mjd":mjd,"telescope":vals[2],"target_name":vals[3],
             "center_mhz":freq,"file_type":vals[7],"size":vals[8],
             "md5_archive":md5,"download_url":link})
        rows.sort(key=lambda v:(v["separation_deg"],v["mjd"],v["target_name"]))
        result.update({"status":"OK" if filter_violation==0 else "SUSPECT_SEARCH_NOT_CONSTRAINED",
         "http":r.status_code,"html_sha256":sha(r.content),
         "rows_examined":total,"filter_outside_0p55":filter_violation,
         "at_least_100_rows_truncation_risk":total>=100,"matches_within_0p1":sum(z["separation_deg"]<=.1 for z in rows),
         "matches_within_0p5":len(rows),"closest_sep_deg":rows[0]["separation_deg"] if rows else None,
         "verified_archive_checksum_fields":sum(bool(z["md5_archive"]) for z in rows),
         "direct_original_file_links":sum(bool(z["download_url"]) for z in rows),
         "closest_matches":rows[:10]})
    except Exception as e:result["error"]=str(e)[:200]
    print("OSIRIS_SPEAR_BL",json.dumps({k:v for k,v in result.items() if k not in ("closest_matches",)},ensure_ascii=False),flush=True)
    return result
def main():
    positions=src.positions()
    jobs=[(p,telescope) for p in positions for telescope in P["telescopes"]]
    assert len(jobs)==60
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        data=list(pool.map(one,jobs))
    save("bl_radio_coverage.json",{"official_archive":P["archive"],"metadata":data,
      "source_status":"ARCHIVE_TARGET_RA_DEC_NOT_ANTENNA_ENCODER_POSITION",
      "real_radio_spectra_downloaded":0,"JANUS_full_candidate_detector_run":False})
    original=[d for d in data if d["kind"]=="original"]
    summaries=[]
    for key in P["telescopes"]:
        sub=[d for d in original if d["telescope"]==key]
        summaries.append({"telescope":key,"target_queries":len(sub),
         "queries_succeeded":sum(d["status"]=="OK" for d in sub),
         "primary_positions_with_0p1_matches":sum(d.get("matches_within_0p1",0)>0 for d in sub),
         "matches":{d["base_id"]:{"within_0p1":d.get("matches_within_0p1"),
          "within_0p5":d.get("matches_within_0p5"),"closest_deg":d.get("closest_sep_deg"),
          "has_download_url":bool(d.get("direct_original_file_links"))}
          for d in sub}})
    ok=sum(d["status"]=="OK" for d in data)
    candidate=sum(d.get("matches_within_0p1",0)>0 for d in original)
    result={"stage":"OSIRIS_SPEAR_RADIO_R1C",
      "status":"PASS_BL_OFFICIAL_TARGET_COORDINATE_SURVEY__HOLD_NATIVE_JANUS_SPECTRAL_EXAM" if ok==60 else
       "PARTIAL_BL_CATALOG_QUERY_OR_UNCONSTRAINED__HOLD_NATIVE_JANUS_SPECTRAL_EXAM",
      "queries_attempted":len(data),"successful_constrained_queries":ok,
      "primary_positions_with_archive_0p1deg_metadata":candidate,
      "telescope_summaries":summaries,
      "real_spectral_fits_or_hdf5_downloaded":0,
      "native_actual_antenna_pointing_verified":False,
      "same_frequency_ON_OFF_cadences_qualified":0,
      "R42_radio_detector_run":False,
      "read_only":True,"historical_optical_geometry_PASSED":False,
      "no_technosignature_or_nondetection_claim":True,"production_unchanged":True,
      "next":"For any direct original target with HDF5+MD5, freeze single selected science-vs-control group before power download, verify byte hashes, native scan pointing and frequency/time, then invoke frozen JANUS native detector, retain all candidates and OFF comparisons. Mirrors/controls are new sky positions, not independent observed source confirmations."}
    save("bl_gate_report.json",result)
    print("OSIRIS_SPEAR_BL_RESULT",json.dumps(result,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
