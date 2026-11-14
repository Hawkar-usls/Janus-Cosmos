"""Read-only end-to-end Galaxy Volume browser smoke; real Fermi data from local repository."""
import json
import os
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from io import BytesIO
from PIL import Image
from playwright.sync_api import sync_playwright

class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=os.getcwd(), **kwargs)
    def log_message(self, *_args):
        pass

def run():
    server=ThreadingHTTPServer(("127.0.0.1",0),Handler)
    threading.Thread(target=server.serve_forever,daemon=True).start()
    try:
        with sync_playwright() as pw:
            browser=pw.chromium.launch(headless=True,args=[
                "--enable-webgl","--use-angle=swiftshader","--enable-unsafe-swiftshader"
            ])
            page=browser.new_page(viewport={"width":1440,"height":900},device_scale_factor=1)
            errors=[]
            page.on("pageerror",lambda e:errors.append(str(e)))
            page.goto("http://127.0.0.1:%d/docs/galaxy/"%server.server_port,wait_until="domcontentloaded")
            page.wait_for_function("document.querySelector('#fermi-count').textContent.includes('38')",timeout=30000)
            engine=page.locator("#engine").inner_text()
            assert "WebGL2" in engine,engine
            assert "ILLUSTRATIVE" in page.locator("#view").inner_text()
            png=page.screenshot(path="janus-cosmos-volume-preview.png")
            pic=Image.open(BytesIO(png)).convert("RGB")
            box=pic.crop((180,210,980,745))
            bright=sum(1 for r,g,b in box.getdata() if max(r,g,b)>110 and r+g+b>165)
            assert bright>100,("No luminous galaxy viewport",bright)
            page.locator('[data-mode="fermi"]').click()
            assert "REAL DATA" in page.locator("#view").inner_text()
            assert "38" in page.locator("#fermi-count").inner_text()
            assert "log10" in page.locator("#axes").inner_text()
            assert page.locator("#points").inner_text()=="281"
            assert "/blob/6b902d4cd71971aacb4b504d680692fa798600a1/" in page.locator("#source").get_attribute("href")
            # Validate that the source-bound HI importer only grants *angular* presentation.
            sample={
                "schema":"janus.cosmos.hi.angular_samples.v1",
                "source_url":"https://lambda.gsfc.nasa.gov/product/foreground/fg_hi4pi_info.html",
                "nhi_unit":"cm-2",
                "samples":[{"l_deg":10,"b_deg":5,"nhi_cm2":1e21}]
            }
            page.locator("#hi-upload").set_input_files({
                "name":"test-observation.json",
                "mimeType":"application/json",
                "buffer":json.dumps(sample).encode("utf-8")
            })
            assert "ANGULAR SKY" in page.locator("#view").inner_text()
            assert "OBSERVED_ANGULAR_SHELL" in page.locator("#scientific-note").inner_text()
            page.locator("#reset").click()
            page.locator("#slice").evaluate("(el) => {el.value = 60; el.dispatchEvent(new Event(\u0027input\u0027, {bubbles:true}));}")
            assert page.locator("#slice-value").inner_text()=="60%"
            page.locator("#spin").click()
            assert "Включить" in page.locator("#spin").inner_text()
            assert not errors,errors
            print("PASS COSMOS GPU GALAXY + 38 REAL FERMI EVENTS + HI ANGULAR IMPORT",bright)
            browser.close()
    finally:
        server.shutdown()

if __name__=="__main__":
    run()
