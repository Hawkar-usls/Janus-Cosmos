#!/usr/bin/env python3
"""R42A7B1b: source discovery diagnostic (not new blind test).
Inspects the earliest public archive index and the raw HTML structure of download metadata.
"""
import json,requests,urllib.parse,hashlib
from bs4 import BeautifulSoup
from pathlib import Path
import r42a7b_archive_discovery as first
P=json.loads(Path("research/seti/r42a7b_independent_onoff_protocol.json").read_text())
first.P["search"]["max_bytes"]=20_000_000_000 # listing only; NEVER download these large files
out=Path("runtime/r42a7b_early");out.mkdir(parents=True,exist_ok=True)
allrows=[];pages=[];samples=[]
for file_type in ("HDF5","baseband data"):
    for page in range(1,5):
        params={"file_type":file_type,"project":"GBT","sort":"utc","order":"asc","perPage":100,"page":page,"search":"Search"}
        url=first.ARCHIVE+"?"+urllib.parse.urlencode(params)
        try:
            r=requests.get(url,timeout=(9,30),headers={"User-Agent":"JANUS-R42A7B-archive-index-study/1.0"})
            r.raise_for_status()
            assert len(r.content)<3_000_000
            rows,_=first.parse_html(r.content,r.url)
            allrows+=rows
            soup=BeautifulSoup(r.content,"html.parser")
            if not samples:
                for tr in soup.select("tr"):
                    td=tr.find_all("td",recursive=False)
                    if len(td)>9 and td[0].get_text(" ",strip=True).startswith("20"):
                        samples.append({"td_texts":[t.get_text(" ",strip=True)[:150] for t in td[:12]],
                          "urls":[a.get("href","")[:400] for a in tr.select("a[href]")][:5],
                          "html":str(td[-1])[:600]})
                        if len(samples)>=3:break
            pages.append({"type":file_type,"page":page,"status":"OK","total_rows":len(rows),
                "off":sum(x["off_label"] for x in rows),"with_url":sum(bool(x["url"]) for x in rows),
                "with_md5":sum(bool(x["md5_archive"]) for x in rows),"raw_sha256":hashlib.sha256(r.content).hexdigest()})
        except Exception as e:
            pages.append({"type":file_type,"page":page,"status":"UNAVAILABLE","error":str(e)[:240]})
        print("R42A7B_EARLY_PAGE",json.dumps(pages[-1]),flush=True)
unique={r["archive_entry_id_sha256"]:r for r in allrows}
allrows=list(unique.values())
pairs=first.pairings(allrows)
report={"stage":"R42A7B1b","status":"PUBLIC_ARCHIVE_INDEX_DISCOVERY_ONLY","page_receipts":pages,
 "rows":len(allrows),"off_labels":sum(x["off_label"] for x in allrows),
 "direct_urls":sum(bool(x["url"]) for x in allrows),"md5s":sum(bool(x["md5_archive"]) for x in allrows),
 "provisional_role_pairs":len(pairs),"sample_html":samples,
 "qualified_original_spectra":0,"qualified_real_telescope_pointing":0,"blind_exam_started":False,
 "note":"The early epoch expansion is a distinct discovery preflight, NOT a heldout scientific result. Listing coordinates alone cannot validate independent telescope pointing."}
(out/"early_archive_provenance.json").write_text(json.dumps(report,indent=2)+"\n")
(out/"early_archive_inventory.json").write_text(json.dumps({"records":allrows},indent=2)+"\n")
print("R42A7B_EARLY_SUMMARY",json.dumps({k:v for k,v in report.items() if k not in ("page_receipts","sample_html")}),flush=True)
print("R42A7B_ARCHIVE_HTML_SAMPLE",json.dumps(samples)[:2000],flush=True)
