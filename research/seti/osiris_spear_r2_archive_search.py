#!/usr/bin/env python3
"""OSIRIS SPEAR R2A: authentic BL archival radio spectral footprint and beam *preflight*.
Read-only metadata. One catalogue coordinate and central frequency never prove
actual antenna pointing or observed target flux. No detector scores computed.
"""
import concurrent.futures,hashlib,json,math,pathlib,re,urllib.parse,requests
from bs4 import BeautifulSoup
import osiris_spear_radio_coverage as r1
P=pathlib.Path("research/seti/osiris_spear_r2_protocol.json")
Q=json.loads(P.read_text())
OUT=pathlib.Path("runtime/osiris_spear_r2a");OUT.mkdir(parents=True,exist_ok=True)
HEADERS={"User-Agent":"JANUS-OSIRIS-SPEAR-R2A-catalog-beam-preflight/1.0"}
def sha(x):return hashlib.sha256(x).hexdigest()
def save(n,v):(OUT/n).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
def parse_size(s):
    try:return int(re.sub(r"[\s,]","",s))
    except Exception:return None
def beam(ang,center_mhz,telescope):
    if not math.isfinite(center_mhz) or center_mhz<=0:return {"status":"UNRESOLVED_FREQUENCY"}
    ghz=center_mhz/1000
    fwhm_arcmin=(12.6 if telescope=="GBT" else 14.3*1.42)/ghz
    radius_arcmin=fwhm_arcmin/2
    if ang*60<=radius_arcmin:cat="IN_NOMINAL_HALF_POWER_CATALOG_ONLY"
    elif ang*60<=fwhm_arcmin:cat="BETWEEN_HALF_POWER_AND_FWHM"
    else:cat="NEARBY_OUTSIDE_FWHM"
    return {"estimated_fwhm_arcmin":fwhm_arcmin,"estimated_half_power_radius_arcmin":radius_arcmin,
      "separation_arcmin":ang*60,"offset_in_FWHM":(ang*60)/fwhm_arcmin,
      "status":cat,
      "qualification":"Approximate beam model from archived center frequency, not measured receiver, feed offset, antenna encoder pointing or actual spectral band"}
def cols_match(headers, expected):
    return all(any(s in x.lower() for x in headers) for s in expected)
def parse_html(raw,url,item,file_type,page):
    soup=BeautifulSoup(raw,"html.parser")
    candidates=[]
    observed=0;outbound=0
    for tr in soup.find_all("tr"):
        cells=tr.find_all("td",recursive=False)
        if len(cells)<10:continue
        vals=[c.get_text(" ",strip=True) for c in cells]
        if not re.fullmatch(r"20\d{2}[-/]\d{2}[-/]\d{2}.*",vals[0]):continue
        try:
            ra=float(vals[4]);dec=float(vals[5]);freq=float(vals[6]);mjd=float(vals[1])
            if not(0<=ra<360 and -90<=dec<=90):continue
        except(ValueError,TypeError,IndexError):continue
        observed+=1
        dist=r1.angular(item["ra_deg"],item["dec_deg"],ra,dec)
        if dist>Q["public_archive_search"]["radius_deg"]+0.05:outbound+=1
        if dist>Q["public_archive_search"]["radius_deg"]+1e-5:continue
        links=[urllib.parse.urljoin(url,a["href"]) for a in tr.select("a[href]") if a.get("href")]
        download=next((x for x in links if x.startswith("https://") and (
            re.search(r"\.(h5|hdf5|fil|raw)(\?|$)",x,re.I) or "bldata.berkeley.edu" in x)),None)
        md5=next((x.lower() for x in vals[8:] if re.fullmatch("[0-9a-f]{32}",x,re.I)),None)
        bb=beam(dist,freq,item["telescope"])
        o={"id":item["id"],"kind":item["kind"],"base_id":item["base_id"],
          "catalog_target_name":vals[3],"catalog_ra_deg":ra,"catalog_dec_deg":dec,
          "catalog_center_frequency_mhz":freq,"catalog_utc":vals[0],"catalog_mjd":mjd,
          "telescope":item["telescope"],"reported_file_type":vals[7],
          "file_bytes":parse_size(vals[8]),"archive_md5":md5,"original_download_url":download,
          "separation_deg":dist,"beam_estimate":bb,"page":page}
        candidates.append(o)
    return candidates,observed,outbound
def one(job):
    pos,telescope,typ=job
    item={**pos,"telescope":telescope}
    r=Q["public_archive_search"]["radius_deg"]
    halfra=r/max(.12,math.cos(math.radians(item["dec_deg"])))
    params={"project":telescope,"file_type":typ,"ra":f'{item["ra_deg"]:.9f}',
      "ra_range":f"{halfra:.9f}","decl":f'{item["dec_deg"]:.9f}',
      "decl_range":f"{r:.9f}","search":"Search","perPage":"100"}
    record={"target":pos["id"],"base_id":pos["base_id"],"kind":pos["kind"],"telescope":telescope,
      "file_type":typ,"status":"UNAVAILABLE","pages":[],"rows_examined":0,
      "offcone_returned":0,"matches":[],"telescope_actual_pointing_proven":False}
    for page in (1,2):
        try:
            url=Q["source"]+"?"+urllib.parse.urlencode({**params,"page":page})
            resp=requests.get(url,headers=HEADERS,timeout=(8,Q["public_archive_search"]["timeout_seconds"]))
            resp.raise_for_status()
            if len(resp.content)>Q["public_archive_search"]["max_html_response_bytes"]:raise ValueError("MAX_HTML_BYTES")
            matches,n,out=parse_html(resp.content,resp.url,item,typ,page)
            record["pages"].append({"page":page,"html_sha256":sha(resp.content),"rows_examined":n,
              "matches_within_cone":len(matches),"out_of_cone_rows":out})
            record["rows_examined"]+=n;record["offcone_returned"]+=out
            record["matches"].extend(matches)
            if n<100:break
        except Exception as ex:
            record.setdefault("errors",[]).append({"page":page,"error":str(ex)[:300]})
            break
    record["truncated"]=bool(record["pages"] and record["pages"][-1]["rows_examined"]>=100)
    if record["pages"] and record["offcone_returned"]==0:record["status"]="OK" if not record["truncated"] else "PAGED_LIMIT_POSSIBLE_TRUNCATION"
    elif record["pages"]:record["status"]="SUSPECT_ARCHIVE_FILTER"
    record["matches"].sort(key=lambda x:(x["separation_deg"],x["file_bytes"] or 10**15,x["catalog_utc"]))
    record["match_count"]=len(record["matches"])
    record["within_half_power"]=sum(x["beam_estimate"]["status"]=="IN_NOMINAL_HALF_POWER_CATALOG_ONLY" for x in record["matches"])
    print("OSIRIS_R2A_QUERY",json.dumps({"id":pos["id"],"telescope":telescope,"type":typ,
       "status":record["status"],"pages":len(record["pages"]),"rows_examined":record["rows_examined"],
       "matches":record["match_count"],"within_nominal_halfpower":record["within_half_power"]},
       ensure_ascii=False),flush=True)
    return record
def main():
    positions=[x for x in r1.positions() if x["kind"] in ("original","control")]
    assert len(positions)==15
    jobs=[(item,telescope,typ) for item in positions for telescope in Q["public_archive_search"]["telescopes"]
          for typ in Q["public_archive_search"]["file_types"]]
    assert len(jobs)==90 and len(jobs)*2<=Q["public_archive_search"]["max_total_requests"]
    with concurrent.futures.ThreadPoolExecutor(max_workers=Q["public_archive_search"]["concurrency"]) as pool:
        result=list(pool.map(one,jobs))
    save("expanded_public_radio_inventory.json",{
      "archive":Q["source"],"protocol_sha256":sha(P.read_bytes()),
      "records":result,"catalog_ra_dec_not_encoder_pointing":True,
      "original_spectral_bytes_downloaded":0,
      "JANUS_native_R39_R42_scored":False})
    orig=[r for r in result if r["kind"]=="original"]
    rows=[v for rec in orig for v in rec["matches"]]
    # Dedupe listings for same physical spectrum appearing in different archive modes.
    unique={}
    for obj in rows:
        key=obj["original_download_url"] or (
         obj["base_id"],obj["telescope"],obj["catalog_target_name"],obj["catalog_utc"],obj["catalog_center_frequency_mhz"],obj["reported_file_type"])
        unique.setdefault(str(key),obj)
    values=list(unique.values())
    rank={"IN_NOMINAL_HALF_POWER_CATALOG_ONLY":0,"BETWEEN_HALF_POWER_AND_FWHM":1,
      "NEARBY_OUTSIDE_FWHM":2,"UNRESOLVED_FREQUENCY":3}
    values.sort(key=lambda x:(rank[x["beam_estimate"]["status"]],
       x["separation_deg"],x["file_bytes"] if x["file_bytes"] is not None else 10**15,
       x["catalog_utc"]))
    eligible=[x for x in values if x["beam_estimate"]["status"]=="IN_NOMINAL_HALF_POWER_CATALOG_ONLY"]
    save("beam_screened_candidates.json",{
      "preregistered_target_scope":Q["target_scope"],"unique_source_files":len(values),
      "provisional_half_power":eligible,
      "other_nearby": [x for x in values if x not in eligible][:80],
      "NOT_CONFIRMED_TARGET_POINTING":True,"original_spectral_bytes_verified":0,
      "beam_model_is_approximate":True})
    summaries=[]
    for tar in Q["target_scope"]:
        d=[x for x in values if x["base_id"]==tar]
        summaries.append({"target":tar,"unique_files":len(d),
          "catalog_nominal_half_power":sum(x["beam_estimate"]["status"]=="IN_NOMINAL_HALF_POWER_CATALOG_ONLY" for x in d),
          "closest_deg":min((x["separation_deg"] for x in d),default=None),
          "eligible_hdf5_with_archive_md5":sum(x["beam_estimate"]["status"]=="IN_NOMINAL_HALF_POWER_CATALOG_ONLY" and
            x["reported_file_type"]=="HDF5" and x["archive_md5"] is not None and x["original_download_url"] is not None for x in d)})
    ok=sum(x["status"]=="OK" for x in result)
    report={"stage":"OSIRIS_SPEAR_R2A","status":
      "ARCHIVE_HALF_POWER_CANDIDATES_PRESENT__HOLD_ACTUAL_POINTING_AND_NATIVE_INTEGRITY" if eligible else
      "NO_ARCHIVE_HALF_POWER_CANDIDATE_IN_BOUNDED_SAMPLE__HOLD_NATIVE_JANUS_SIGNAL_TEST",
      "queries":len(result),"successful_untruncated_and_constrained":ok,
      "unavailable_or_bad_or_truncated":len(result)-ok,
      "catalog_unique_original_file_candidates":len(values),
      "catalog_only_nominal_half_power_candidates":len(eligible),
      "targets":summaries,"native_original_radio_files_downloaded":0,
      "real_actual_antenna_pointing_verified":False,"qualified_ON_OFF_groups":0,
      "JANUS_native_detector_run":False,"OSIRIS_technosignature_detected":False,
      "scientific_nondetection_claim_permitted":False,"production_unchanged":True,
      "historical_JANUS_reference":{"false_alerts":51,"detections":318,"weak":69},
      "next":"Select a suitable catalog-only within-beam HDF5 + actual checksum and candidate OFF control if one exists. Before scoring, separately verify source file full bytes and instrument pointing; otherwise source remains NOT_TESTABLE, no SETI detection or nondetection."}
    save("r2a_gate_report.json",report)
    print("OSIRIS_SPEAR_R2A_RESULT",json.dumps(report,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
