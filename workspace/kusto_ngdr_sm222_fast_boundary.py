#!/usr/bin/env python3
import requests,re,json,time
BASE="https://geodataindia.gov.in/"
S=requests.Session(); S.headers.update({"User-Agent":"JANUS-KUSTO-SM222-FAST/1.0","Referer":BASE+"guestuser"})
h=S.get(BASE+"guestuser",timeout=30); h.raise_for_status()
m=re.search(r'id=["\']csrfvalue["\']\s+value=["\']([^"\']+)',h.text)
tok=m.group(1)
for ident in [5192,5837]:
    try:
        r=S.post(BASE+"guestuser/getgridBoundary_basereport?_csrf="+tok,
          data={"tablename":"baseline_data_gis_view","clmname":"baseid","grid_no":str(ident)},timeout=20)
        try: js=r.json()
        except: js=None
        print(json.dumps({"id":ident,"status":r.status_code,"json":js,"text":r.text[:1000]},ensure_ascii=False))
    except Exception as e:
        print(json.dumps({"id":ident,"error":repr(e)}))
