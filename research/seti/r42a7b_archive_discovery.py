#!/usr/bin/env python3
"""R42A7B1: read-only catalog discovery of candidate *new* raw ON/OFF files.
Archive listing is NOT original telescope pointing. No performance inference.
"""
import hashlib,json,math,os,pathlib,re,urllib.parse
from datetime import datetime,timezone
import requests
from bs4 import BeautifulSoup

P=json.loads(pathlib.Path("research/seti/r42a7b_independent_onoff_protocol.json").read_text())
OUT=pathlib.Path("runtime/r42a7b")
OUT.mkdir(parents=True,exist_ok=True)
ARCHIVE=P["archive"]
HEADERS={"User-Agent":"JANUS-COSMOS-R42A7B/1.0 research-reproducibility read-only"}
def sha(s):return hashlib.sha256(s).hexdigest()
def save(name,d):(OUT/name).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
def parse_html(raw,url):
    soup=BeautifulSoup(raw,"html.parser")
    links=[]
    records=[]
    for tr in soup.find_all("tr"):
        cells=tr.find_all(["td"],recursive=False)
        if len(cells)<9:continue
        a=[" ".join(c.stripped_strings) for c in cells]
        # Expected archive result: UTC | MJD | Telescope | Target | RA | Dec | Freq | Type | Size | MD5 | Download
        date=a[0].strip()
        if not re.match(r"^20\d{2}-\d\d-\d\d(?:\s+\d\d:\d\d:\d\d)?$",date):continue
        try:
            mjd=float(a[1]);telescope=a[2].strip();target=a[3].strip()
            ra=float(a[4]);dec=float(a[5]);freq=float(a[6])
            typ=a[7].strip();size=int(re.sub(r"[^0-9]","",a[8]))
        except (IndexError,TypeError,ValueError):continue
        if not (50000<mjd<70000 and telescope in ("GBT","Parkes") and 0<=ra<360 and -90<=dec<=90 and
            100<freq<30000 and size>=P["search"]["min_bytes"] and size<=P["search"]["max_bytes"]):continue
        try:
            datetime.fromisoformat(date).replace(tzinfo=timezone.utc)
        except ValueError:continue
        md5=next((x.strip().lower() for x in a[9:] if re.fullmatch("[a-fA-F0-9]{32}",x.strip())),None)
        possible=[]
        for c in cells:
            for anchor in c.find_all("a",href=True):
                href=urllib.parse.urljoin(url,anchor["href"])
                if href.startswith("https://") and "opendatasearch?" not in href:
                    possible.append(href)
        # Avoid links to user search, javascript or unnecessary pages.
        possible=[x for x in possible if re.search(r"\.(h5|hdf5|fil|gz|raw)(?:[?#]|$)",x,re.I)]
        file_url=possible[0] if possible else None
        if not file_url:
            # Archive often serves files via javascript click buttons: preserve as unqualified.
            for c in cells:
                for tag in c.find_all(attrs={"onclick":True}):
                    links.append(str(tag.get("onclick"))[:300])
        records.append({"utc":date,"mjd":mjd,"telescope":telescope,"target":target,
            "ra_deg_catalog":ra,"dec_deg_catalog":dec,"freq_mhz":freq,"file_type":typ,
            "size_bytes":size,"md5_archive":md5,"url":file_url,
            "off_label":bool(re.search(r"_OFF[A-Z0-9]*$",target.upper())),
            "archive_page":url,
            "archive_entry_id_sha256":sha(("|".join([date,str(mjd),telescope,target,str(freq),str(size),md5 or ""])).encode())})
    return records,links
def angular(ra1,dec1,ra2,dec2):
    a,b,c,d=map(math.radians,[ra1,dec1,ra2,dec2])
    x=math.sin((b-d)/2)**2+math.cos(b)*math.cos(d)*math.sin((a-c)/2)**2
    return 2*math.degrees(math.asin(min(1,math.sqrt(max(0,x)))))
def discover():
    allrows=[];pages=[]
    session=requests.Session()
    session.headers.update(HEADERS)
    for page in range(1,P["search"]["max_pages"]+1):
        params={"file_type":"HDF5","project":"GBT","sort":"utc","order":"desc",
                "perPage":P["search"]["per_page"],"page":page,"search":"Search"}
        url=ARCHIVE+"?"+urllib.parse.urlencode(params)
        try:
            r=session.get(url,timeout=(10,P["search"]["request_timeout_seconds"]))
            r.raise_for_status()
            if len(r.content)>3000000:raise ValueError("HTML_OVER_3_MB")
            rows,onclick=parse_html(r.content,r.url)
            pages.append({"page":page,"status":"OK","bytes":len(r.content),"sha256":sha(r.content),
                "matching_records":len(rows),"unparsed_download_buttons":len(onclick),"url":r.url,
                "sample_onclick":onclick[:2]})
            allrows.extend(rows)
            print("R42A7B_ARCHIVE_PAGE",json.dumps(pages[-1]),flush=True)
        except Exception as e:
            pages.append({"page":page,"status":"UNAVAILABLE","url":url,"error":str(e)[:280]})
            print("R42A7B_ARCHIVE_PAGE",json.dumps(pages[-1]),flush=True)
    unique={}
    for x in allrows:unique.setdefault(x["archive_entry_id_sha256"],x)
    rows=list(unique.values())
    rows.sort(key=lambda x:(x["mjd"],x["target"],x["archive_entry_id_sha256"]))
    save("archive_inventory.json",{"source":ARCHIVE,"page_receipts":pages,"rows":rows,
       "observations_from_new_archive":True,"headers_read":False,"source_file_contents_read":False,
       "catalog_ra_is_not_guaranteed_actual_pointing":True})
    return rows,pages
def pairings(rows):
    # Conservative: _OFF records are potential spatial controls, never verified actual OFF by their names alone.
    offs=[x for x in rows if x["off_label"]]
    ons=[x for x in rows if not x["off_label"]]
    out=[]
    for off in offs:
        prefix=re.sub(r"_OFF[A-Z0-9]*$","",off["target"].upper())
        options=[]
        for on in ons:
            if on["target"].upper()!=prefix:continue
            if on["telescope"]!=off["telescope"] or abs(on["freq_mhz"]-off["freq_mhz"])>1e-6:continue
            gap=abs(on["mjd"]-off["mjd"])*1440
            if gap>P["pair_gate"]["max_time_minutes"]:continue
            sep=angular(on["ra_deg_catalog"],on["dec_deg_catalog"],off["ra_deg_catalog"],off["dec_deg_catalog"])
            if not (P["pair_gate"]["min_separation_deg"]<=sep<=P["pair_gate"]["max_separation_deg"]):continue
            options.append((gap,on,sep))
        if options:
            gap,on,sep=min(options,key=lambda y:(y[0],y[1]["archive_entry_id_sha256"]))
            out.append({"on":on,"off":off,"time_gap_minutes":gap,
               "listed_coordinate_separation_deg":sep,
               "status":"PROVISIONAL_FROM_ARCHIVE_LABELS_NOT_ORIGINAL_POINTING",
               "verified_original_header":False,"independent_TCS_pointing":False,
               "physically_qualified_for_blind_exam":False})
    save("provisional_onoff_pairs.json",{"pairs":out,"count":len(out),
       "note":"Archive target labels and listed RA/Dec are not independent physical per-scan pointing proof. Native file header and independently corroborated telescope logs required."})
    print("R42A7B_PROVISIONAL",json.dumps({"pairs":len(out),"examples":[
        {"target":x["on"]["target"],"utc":x["on"]["utc"],"off":x["off"]["target"],
         "gap_min":round(x["time_gap_minutes"],2),"sep_deg":round(x["listed_coordinate_separation_deg"],4)}
        for x in out[:10]]}),flush=True)
    return out
def gate(rows,pages,pairs):
    complete=sum(x["status"]=="OK" for x in pages)
    status="HOLD_NO_OFF_ROLE_PAIRS" if not pairs else "HOLD_NATIVE_SPECTRA_AND_INDEPENDENT_TCS_PENDING"
    if not complete:status="HOLD_OFFICIAL_ARCHIVE_UNAVAILABLE"
    report={"stage":"R42A7B1","status":status,
      "archive_pages_accessed":complete,"archive_pages_requested":len(pages),
      "new_archive_metadata_files":len(rows),
      "rows_with_checksum":sum(bool(x["md5_archive"]) for x in rows),
      "rows_with_direct_download_url":sum(bool(x["url"]) for x in rows),
      "archive_off_label_rows":sum(x["off_label"] for x in rows),
      "provisional_spatial_onoff_pairs":len(pairs),
      "original_file_byte_hashes_verified":0,
      "native_source_spectra_verified":0,
      "native_header_per_scan_pointing_verified":0,
      "independent_TCS_or_backend_pointing_verified":0,
      "independent_blind_custodian_confirmed":False,
      "qualified_fresh_groups":0,
      "R42A7B2_detector_allowed":False,
      "historical_audit_rescored":False,
      "reference_production_metrics_unchanged":{"false_alerts":51,"localized":318,"weak":69},
      "shadow_model_production_veto":False,
      "no_extraterrestrial_detection_claim":True,
      "next":"Choose a bounded direct-download ON/OFF pair; verify full file checksum, nonzero per-scan native coordinates/time, original observing TCS/back-end pointing separate from source catalog, disjoint observation hashes/dates and custodian-signed blind manifest. Without these hold; never fabricate ON/OFF."}
    save("gate_report.json",report)
    print("R42A7B_GATE",json.dumps(report),flush=True)
    return report
if __name__=="__main__":
    rows,pages=discover()
    pairs=pairings(rows)
    gate(rows,pages,pairs)
