"""Read-only browser smoke on the source-derived FERMI-V7 dataset."""
from __future__ import annotations
import re
import os
import threading
from http.server import SimpleHTTPRequestHandler,ThreadingHTTPServer
from io import BytesIO
from PIL import Image
from playwright.sync_api import sync_playwright

class Handler(SimpleHTTPRequestHandler):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,directory=os.getcwd(),**kwargs)
    def log_message(self,*args):
        return

def run():
    server=ThreadingHTTPServer(("127.0.0.1",0),Handler)
    threading.Thread(target=server.serve_forever,daemon=True).start()
    try:
        with sync_playwright() as p:
            browser=p.chromium.launch(headless=True)
            page=browser.new_page(viewport={"width":1500,"height":980},device_scale_factor=1)
            errors=[];page.on("pageerror",lambda e:errors.append(str(e)))
            page.goto("http://127.0.0.1:%d/docs/fermi/"%server.server_port,
                      wait_until="domcontentloaded")
            page.wait_for_function("window.__JANUS_FERMI_V7_DIAG && window.__JANUS_FERMI_V7_DIAG().loaded",
                                   timeout=30000)
            init=page.evaluate("window.__JANUS_FERMI_V7_DIAG()")
            assert init["sourceEvents"]==175,init
            assert init["filteredEvents"]==175,init
            assert init["burstCount"]==3,init
            assert init["scientificStatus"]=="BLOCKED",init
            assert init["rawStatus"]=="RAW_FITS_NOT_RECHECKED_THIS_BUILD",init
            assert "6520"==re.sub(r"\D","",page.locator("#metric-roi").inner_text()),page.locator("#metric-roi").inner_text()
            png=page.screenshot(path="janus-fermi-v7-preview.png",full_page=True)
            pic=Image.open(BytesIO(png)).convert("RGB")
            assert pic.width>=1500
            page.locator("#burst").select_option("090510")
            after=page.evaluate("window.__JANUS_FERMI_V7_DIAG()")
            assert after["filteredEvents"]==32,after
            page.locator("#min-energy").evaluate("(el) => {el.value=10; el.dispatchEvent(new Event('input',{bubbles:true}));}")
            tight=page.evaluate("window.__JANUS_FERMI_V7_DIAG()")
            assert 0<tight["filteredEvents"]<32,tight
            page.locator("#reset").click()
            assert page.evaluate("window.__JANUS_FERMI_V7_DIAG().filteredEvents")==175
            page.locator("#ledger tr").nth(1).click()
            assert page.evaluate("window.__JANUS_FERMI_V7_DIAG().selectedEvent")
            assert "LAT Extended Pass8" in page.locator("#event-detail").inner_text()
            assert page.locator("#source-fits").get_attribute("href").startswith(
                "https://fermi.gsfc.nasa.gov/FTP/fermi/data/lat/queries/")
            assert not errors,errors
            print("PASS FERMI V7 BROWSER: archived events, 3 GRBs, source-bound filters, scientific ceilings")
            browser.close()
    finally:
        server.shutdown()

if __name__=="__main__":
    run()
