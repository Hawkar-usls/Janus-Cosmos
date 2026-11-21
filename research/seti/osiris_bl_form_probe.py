"""Read-only interface discovery for official Breakthrough Listen RA/Dec query form."""
import hashlib,json,pathlib,requests
from bs4 import BeautifulSoup
P=json.loads(pathlib.Path("research/seti/osiris_bl_form_probe_protocol.json").read_text())
O=pathlib.Path("runtime/osiris_bl_form_probe");O.mkdir(parents=True,exist_ok=True)
results=[]
for url in P["archive_urls"]:
    item={"url":url,"status":"UNAVAILABLE"}
    try:
        r=requests.get(url,headers={"User-Agent":"JANUS-COSMOS OSIRIS open radio metadata audit/1.0"},
                       timeout=(9,20));r.raise_for_status()
        if len(r.content)>P["limits"]["max_response_bytes"]:raise ValueError("HTML_MAX_BYTES_EXCEEDED")
        soup=BeautifulSoup(r.content,"html.parser")
        item.update({"status":"HTML_OK","effective_url":r.url,
                     "html_sha256":hashlib.sha256(r.content).hexdigest(),
                     "forms":[]})
        for form in soup.find_all("form")[:8]:
            controls=[]
            for tag in form.find_all(["input","select","textarea"]):
                controls.append({"tag":tag.name,"name":tag.get("name"),"id":tag.get("id"),
                                 "placeholder":tag.get("placeholder"),"type":tag.get("type"),
                                 "value":tag.get("value"),"options":[{"value":o.get("value"),"label":o.get_text(" ",strip=True)[:75]} for o in tag.find_all("option")[:12]] if tag.name=="select" else None})
            item["forms"].append({"method":form.get("method","GET"),"action":form.get("action"),"controls":controls})
        item["table_headers"]=[[th.get_text(" ",strip=True) for th in t.find_all("th")[:15]] for t in soup.find_all("table")[:2]]
    except Exception as e:item["error"]=str(e)[:220]
    results.append(item)
    print("OSIRIS_BL_FORM",json.dumps(item,ensure_ascii=False),flush=True)
(O/"form_probe.json").write_text(json.dumps({"stage":"OSIRIS_SPEAR_RADIO_R1B","results":results,
    "original_radio_spectra_read":0,"JANUS_detector_applied":False},indent=2,ensure_ascii=False)+"\n")
