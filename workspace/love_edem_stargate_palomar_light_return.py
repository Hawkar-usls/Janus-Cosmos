from __future__ import annotations

import gzip, hashlib, io, json, math
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import requests

TARGETS={
 "LOVE":(204.30267916666668,-36.78240527777778),
 "EDEM":(139.22409686590188,30.26038779947318),
 "STARGATE_ABYDOS_GEOMETRY":(223.415064157,33.979315670),
}
RADII=[2,5,10,30,60,300,600]
BASE="https://raw.githubusercontent.com/jannefi/poss1-plate-slice/4005e200541b321ead3d6608f0162a14430ef1c2/results/s0-642-20260814"
STAGE_URL=f"{BASE}/stage_S0.csv.gz"
TILE_URL=f"{BASE}/tile_manifest.csv.gz"
EXPECTED_ROWS=122820
EXPECTED_SHA="2ff92f2210acb387ef9ef4b88d561595d3883e9aab27065042627272b96590f0"
OUT=Path("data/love/LOVE-EDEM-STARGATE-PALOMAR-LIGHT-RETURN-v1-LATEST-RECEIPT.json")
C=299792458.0
JULIAN_YEAR_S=365.25*86400.0
LY_M=C*JULIAN_YEAR_S

def sep_arcsec(ra0,de0,ra,de):
    r0=np.deg2rad(ra0); d0=np.deg2rad(de0); r=np.deg2rad(ra); d=np.deg2rad(de)
    x=np.sin(d0)*np.sin(d)+np.cos(d0)*np.cos(d)*np.cos(r0-r)
    return np.rad2deg(np.arccos(np.clip(x,-1,1)))*3600.0

def dl(url):
    r=requests.get(url,timeout=120); r.raise_for_status(); return r.content

def main():
    stage_gz=dl(STAGE_URL); stage_raw=gzip.decompress(stage_gz)
    sha=hashlib.sha256(stage_raw).hexdigest(); stage=pd.read_csv(io.BytesIO(stage_raw))
    if sha!=EXPECTED_SHA or len(stage)!=EXPECTED_ROWS: raise RuntimeError(f"stage identity mismatch rows={len(stage)} sha={sha}")
    tile_raw=gzip.decompress(dl(TILE_URL)); tiles=pd.read_csv(io.BytesIO(tile_raw))
    tile_cols=list(tiles.columns)
    if "tile_id" not in tiles.columns: raise RuntimeError(f"tile manifest lacks tile_id: {tile_cols}")
    tile_index=tiles.drop_duplicates("tile_id").set_index("tile_id",drop=False)
    results={}
    for name,(ra0,de0) in TARGETS.items():
        seps=sep_arcsec(ra0,de0,stage["ra"].to_numpy(float),stage["dec"].to_numpy(float))
        order=np.argsort(seps)
        nearest=[]
        for idx in order[:20]:
            row=stage.iloc[int(idx)]
            tid=str(row.get("tile_id"))
            meta={}
            if tid in tile_index.index:
                tr=tile_index.loc[tid]
                if isinstance(tr,pd.DataFrame): tr=tr.iloc[0]
                for c in tile_cols:
                    v=tr[c]
                    if pd.isna(v): meta[c]=None
                    elif isinstance(v,(np.integer,)): meta[c]=int(v)
                    elif isinstance(v,(np.floating,)): meta[c]=float(v)
                    else: meta[c]=str(v)
            nearest.append({
                "separation_arcsec":float(seps[int(idx)]),"ra_deg":float(row["ra"]),"dec_deg":float(row["dec"]),
                "src_id":str(row.get("src_id")),"tile_id":tid,"object_id":str(row.get("object_id")),"tile_manifest":meta,
            })
        results[name]={
            "radius_counts_arcsec":{str(r):int(np.sum(seps<=r)) for r in RADII},
            "nearest":nearest[0],"nearest_20":nearest,
        }
    physics={
      "speed_of_light_m_s":C,"light_year_m":LY_M,
      "one_way":{"1_ly_seconds":JULIAN_YEAR_S,"rule":"t_years=D_ly"},
      "round_trip_reflection":{"rule":"delay_years=2*D_ly","distance_ly_for_delay_years":{str(y):y/2 for y in [1,2,5,10,20,50,76,80,100,1000,9000]}},
      "palomar_era_constraint":"An Earth-origin flash observed back on Earth within years/decades requires a reflector only light-years away; kpc-scale neighbours cannot return that same flash in the 1950s."
    }
    payload={
      "schema":"janus.cosmos.love_edem_stargate.palomar_light_return.receipt.v1",
      "experiment_id":"LOVE-EDEM-STARGATE-PALOMAR-LIGHT-RETURN-v1","status":"COMPLETE","run_time_utc":datetime.now(timezone.utc).isoformat(),
      "targets":{k:{"ra_deg":v[0],"dec_deg":v[1]} for k,v in TARGETS.items()},
      "source":{"stage_rows":len(stage),"stage_uncompressed_sha256":sha,"tile_manifest_rows":len(tiles),"tile_manifest_columns":tile_cols},
      "results":results,"light_travel_physics":physics,
      "firewall":{"nearest_candidate_is_verified_transient":False,"candidate_is_artificial_mirror":False,"palomar_nuclear_association_proves_causation":False,"direction_match_alone_proves_reflection":False,"planet_claim":False,"claim_ceiling":"SPATIAL_CROSSMATCH_AND_LIGHT_TRAVEL_CAUSALITY_ONLY"}
    }
    OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(payload,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:{"nearest_arcsec":v["nearest"]["separation_arcsec"],"counts":v["radius_counts_arcsec"],"tile_manifest":v["nearest"]["tile_manifest"]} for k,v in results.items()},ensure_ascii=False,indent=2))
if __name__=="__main__": main()
