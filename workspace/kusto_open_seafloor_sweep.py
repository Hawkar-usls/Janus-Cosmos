#!/usr/bin/env python3
import json, math, os, re, sys, tempfile
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import requests
import xarray as xr

OUT = Path(os.environ.get("KUSTO_OUT", "workspace/kusto_open_seafloor_out"))
OUT.mkdir(parents=True, exist_ok=True)

TARGETS = [
    {
        "id": "CD169_TOBI_LOCUS",
        "lat": -3.865418,
        "lon": -12.14244,
        "role": "PRIMARY_OPEN_SEAFLOOR_TARGET",
        "origin": "CD169_TOPOLOGY_AND_TOBI_LINE19",
    },
    {
        "id": "KUSTO_AFRICA_LOCUS",
        "lat": -3.865418,
        "lon": 3.854924,
        "role": "LEGACY_FROZEN_KUSTO_TARGET",
        "origin": "LOVE_EDEM_BISECTOR_LEGACY_LONGITUDE",
    },
    {
        "id": "H10N1_CALIBRATION_CONTROL",
        "lat": -7.845673,
        "lon": -14.480230,
        "role": "KNOWN_CALIBRATION_CONTROL",
        "origin": "ASCENSION_HA10_CALIBRATION_LANE",
    },
]

SESSION = requests.Session()
SESSION.headers.update({"User-Agent":"JANUS-KUSTO-open-seafloor-sweep/1.0 research reproducibility"})

def haversine_m(lat1, lon1, lat2, lon2):
    R=6371008.8
    p1,p2=math.radians(lat1),math.radians(lat2)
    dp=math.radians(lat2-lat1)
    dl=math.radians(lon2-lon1)
    a=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R*math.asin(math.sqrt(a))

def request(url, params=None, binary=False, timeout=90):
    r=SESSION.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    return r.content if binary else r.text

def detect_xy_z(ds):
    # GMRT files generally expose x/y/z; keep this tolerant.
    xname = next((n for n in ["x","lon","longitude"] if n in ds.variables or n in ds.coords), None)
    yname = next((n for n in ["y","lat","latitude"] if n in ds.variables or n in ds.coords), None)
    if not xname or not yname:
        raise RuntimeError(f"Could not find coordinate axes: {list(ds.variables)}")
    candidates=[]
    for n,v in ds.data_vars.items():
        if v.ndim>=2:
            candidates.append(n)
    if not candidates:
        raise RuntimeError(f"Could not find 2D data variable: {list(ds.variables)}")
    zname = "z" if "z" in candidates else candidates[0]
    return xname,yname,zname

def normalize_grid(ds, xname, yname, zname):
    x=np.asarray(ds[xname].values, dtype=float)
    y=np.asarray(ds[yname].values, dtype=float)
    z=np.asarray(ds[zname].values, dtype=float)
    z=np.squeeze(z)
    if z.shape == (len(x),len(y)):
        z=z.T
    if z.shape != (len(y),len(x)):
        raise RuntimeError(f"grid shape mismatch z={z.shape} y={len(y)} x={len(x)}")
    return x,y,z

def local_stats(x,y,z,lat,lon):
    ix=int(np.argmin(np.abs(x-lon)))
    iy=int(np.argmin(np.abs(y-lat)))
    nearest=float(z[iy,ix]) if np.isfinite(z[iy,ix]) else None
    nearest_cell_dist=haversine_m(lat,lon,float(y[iy]),float(x[ix]))
    valid=np.isfinite(z)
    valid_fraction=float(valid.mean())
    if valid.any():
        ys,xs=np.where(valid)
        dists=np.array([haversine_m(lat,lon,float(y[j]),float(x[i])) for j,i in zip(ys,xs)])
        k=int(np.argmin(dists))
        near_valid={
            "distance_m":float(dists[k]),
            "lat":float(y[ys[k]]),
            "lon":float(x[xs[k]]),
            "value_m":float(z[ys[k],xs[k]])
        }
    else:
        near_valid=None
    stats={
        "shape":[int(z.shape[0]),int(z.shape[1])],
        "valid_fraction":valid_fraction,
        "nearest_grid_cell":{"lat":float(y[iy]),"lon":float(x[ix]),"distance_m":nearest_cell_dist,"value_m":nearest},
        "nearest_valid_cell":near_valid,
        "finite_min_m":float(np.nanmin(z)) if valid.any() else None,
        "finite_max_m":float(np.nanmax(z)) if valid.any() else None,
        "finite_median_m":float(np.nanmedian(z)) if valid.any() else None,
    }
    # morphology only around direct high-res target neighborhoods
    if nearest is not None:
        lat_scale=111320.0
        lon_scale=111320.0*math.cos(math.radians(lat))
        X=(x-lon)*lon_scale
        Y=(y-lat)*lat_scale
        xx,yy=np.meshgrid(X,Y)
        rr=np.hypot(xx,yy)
        multi={}
        for radius in [1000,2000,5000,10000]:
            m=(rr<=radius)&valid
            if m.sum()>=5:
                vals=z[m]
                multi[str(radius)] = {
                    "n":int(m.sum()),
                    "min_m":float(np.nanmin(vals)),
                    "max_m":float(np.nanmax(vals)),
                    "relief_m":float(np.nanmax(vals)-np.nanmin(vals)),
                    "median_m":float(np.nanmedian(vals)),
                    "p05_m":float(np.nanpercentile(vals,5)),
                    "p95_m":float(np.nanpercentile(vals,95))
                }
        stats["local_relief"]=multi
        # conservative gradient summary; not a target classification
        if len(x)>2 and len(y)>2:
            dx=np.nanmedian(np.diff(x))*lon_scale
            dy=np.nanmedian(np.diff(y))*lat_scale
            if abs(dx)>0 and abs(dy)>0:
                dzdy,dzdx=np.gradient(z,dy,dx)
                slope=np.degrees(np.arctan(np.hypot(dzdx,dzdy)))
                m=(rr<=5000)&np.isfinite(slope)&valid
                if m.sum()>=5:
                    stats["slope_5km_deg"]={
                        "median":float(np.nanmedian(slope[m])),
                        "p90":float(np.nanpercentile(slope[m],90)),
                        "max":float(np.nanmax(slope[m]))
                    }
    return stats

def gmrt(target):
    lat,lon=target["lat"],target["lon"]
    pad=0.12
    common={
        "west":lon-pad,"east":lon+pad,"south":lat-pad,"north":lat+pad,
        "format":"netcdf","resolution":"max"
    }
    base="https://www.gmrt.org/services/GridServer"
    md_url=base+"/metadata"
    md_params=dict(common)
    md_params.pop("format",None)
    md_params.update({"format":"netcdf","mformat":"json"})
    try:
        md_text=request(md_url,md_params)
        try: md=json.loads(md_text)
        except Exception: md={"raw":md_text[:20000]}
    except Exception as e:
        md={"error":repr(e)}
    products={}
    for layer in ["topo-mask","topo"]:
        params=dict(common); params["layer"]=layer
        fn=OUT/f'{target["id"].lower()}_gmrt_{layer.replace("-","_")}.nc'
        try:
            data=request(base,params,binary=True,timeout=180)
            fn.write_bytes(data)
            with xr.open_dataset(fn) as ds:
                xname,yname,zname=detect_xy_z(ds)
                x,y,z=normalize_grid(ds,xname,yname,zname)
                stats=local_stats(x,y,z,lat,lon)
                products[layer]={
                    "url":base+"?"+urlencode(params),
                    "bytes":len(data),
                    "coords":{"x":xname,"y":yname,"z":zname},
                    "stats":stats
                }
        except Exception as e:
            products[layer]={"error":repr(e),"url":base+"?"+urlencode(params)}
    return {"metadata":md,"products":products}

def point_to_polyline_distance_m(lat,lon,geom):
    # local equirectangular approximation, sufficient for <=~30 km search windows.
    sx=111320.0*math.cos(math.radians(lat)); sy=111320.0
    px=lon*sx; py=lat*sy
    best=None
    for path in geom.get("paths",[]):
        for a,b in zip(path,path[1:]):
            ax,ay=a[0]*sx,a[1]*sy; bx,by=b[0]*sx,b[1]*sy
            vx,vy=bx-ax,by-ay; wx,wy=px-ax,py-ay
            den=vx*vx+vy*vy
            t=0 if den==0 else max(0,min(1,(wx*vx+wy*vy)/den))
            qx,qy=ax+t*vx,ay+t*vy
            d=math.hypot(px-qx,py-qy)
            if best is None or d<best: best=d
    return best

def ncei(target):
    lat,lon=target["lat"],target["lon"]
    pad=0.25
    endpoint="https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/multibeam_dynamic/MapServer/0/query"
    params={
        "geometry":f"{lon-pad},{lat-pad},{lon+pad},{lat+pad}",
        "geometryType":"esriGeometryEnvelope",
        "inSR":"4326",
        "outSR":"4326",
        "spatialRel":"esriSpatialRelIntersects",
        "outFields":"SURVEY_ID,PLATFORM,SURVEY_YEAR,SOURCE,NGDC_ID,CHIEF_SCIENTIST,INSTRUMENT,FILE_COUNT,TRACK_LENGTH,TOTAL_TIME,BATHY_BEAMS,AMP_BEAMS,SIDESCANS,DOWNLOAD_URL,START_TIME,END_TIME",
        "returnGeometry":"true",
        "f":"json"
    }
    try:
        data=json.loads(request(endpoint,params))
    except Exception as e:
        return {"error":repr(e),"query_url":endpoint+"?"+urlencode(params)}
    feats=[]
    for f in data.get("features",[]):
        a=f.get("attributes",{})
        dist=point_to_polyline_distance_m(lat,lon,f.get("geometry",{}))
        a["min_trackline_distance_m"]=dist
        feats.append(a)
    feats.sort(key=lambda a: float("inf") if a.get("min_trackline_distance_m") is None else a["min_trackline_distance_m"])
    return {
        "query_url":endpoint+"?"+urlencode(params),
        "candidate_count":len(feats),
        "nearest_candidates":feats[:50],
        "trackline_within_1km_count":sum(1 for a in feats if a.get("min_trackline_distance_m") is not None and a["min_trackline_distance_m"]<=1000),
        "trackline_within_5km_count":sum(1 for a in feats if a.get("min_trackline_distance_m") is not None and a["min_trackline_distance_m"]<=5000),
        "warning":"TRACKLINE_PROXIMITY_IS_NOT_SWATH_COVERAGE"
    }


def ncei_footprints(target):
    lat,lon=target["lat"],target["lon"]
    endpoint="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
    exact_params={
        "geometry":f"{lon},{lat}",
        "geometryType":"esriGeometryPoint",
        "inSR":"4326",
        "outSR":"4326",
        "spatialRel":"esriSpatialRelIntersects",
        "outFields":"NCEI_ID,SURVEY_ID,PLATFORM,SOURCE,CHIEF_SCIENTIST,INSTRUMENT,START_TIME,END_TIME,SURVEY_YEAR,DOWNLOAD_URL,SURVEY_AND_VERSION",
        "returnGeometry":"false",
        "f":"json"
    }
    try:
        exact=json.loads(request(endpoint,exact_params))
        feats=[x.get("attributes",{}) for x in exact.get("features",[])]
    except Exception as e:
        feats=[]
        exact={"error":repr(e)}
    nearby_params={
        "geometry":f"{lon-0.05},{lat-0.05},{lon+0.05},{lat+0.05}",
        "geometryType":"esriGeometryEnvelope",
        "inSR":"4326",
        "outSR":"4326",
        "spatialRel":"esriSpatialRelIntersects",
        "outFields":"NCEI_ID,SURVEY_ID,PLATFORM,SOURCE,CHIEF_SCIENTIST,INSTRUMENT,START_TIME,END_TIME,SURVEY_YEAR,DOWNLOAD_URL,SURVEY_AND_VERSION",
        "returnGeometry":"false",
        "f":"json"
    }
    try:
        near=json.loads(request(endpoint,nearby_params))
        nearfeats=[x.get("attributes",{}) for x in near.get("features",[])]
    except Exception as e:
        nearfeats=[]
        near={"error":repr(e)}
    return {
        "exact_point_intersection_count":len(feats),
        "exact_point_intersections":feats,
        "nearby_0p05deg_candidate_count":len(nearfeats),
        "nearby_0p05deg_candidates":nearfeats,
        "exact_query_url":endpoint+"?"+urlencode(exact_params),
        "nearby_query_url":endpoint+"?"+urlencode(nearby_params),
        "authority":"NCEI_MULTIBEAM_FOOTPRINT_POLYGON_LAYER",
        "firewall":"FOOTPRINT_POLYGON_INTERSECTION_IS_SURVEY_COVERAGE_EVIDENCE__NOT_EXACT_PING_BEAM_IDENTITY"
    }

def parse_feature_value(fi):
    txt=fi.get("response","") if isinstance(fi,dict) else ""
    m=re.search(r"value_list\\s*=\\s*'([^']+)'",txt)
    return None if not m else m.group(1)

def gebco_feature_info(target, layer):
    lat,lon=target["lat"],target["lon"]
    pad=0.02
    endpoint="https://wms.gebco.net/mapserv"
    params={
        "SERVICE":"WMS","VERSION":"1.3.0","REQUEST":"GetFeatureInfo",
        "LAYERS":layer,"QUERY_LAYERS":layer,
        "CRS":"EPSG:4326",
        # EPSG:4326 axis order in WMS 1.3 is lat,lon.
        "BBOX":f"{lat-pad},{lon-pad},{lat+pad},{lon+pad}",
        "WIDTH":"101","HEIGHT":"101","I":"50","J":"50",
        "INFO_FORMAT":"text/plain",
        "FEATURE_COUNT":"10"
    }
    try:
        txt=request(endpoint,params)
        return {"url":endpoint+"?"+urlencode(params),"response":txt[:20000]}
    except Exception as e:
        return {"url":endpoint+"?"+urlencode(params),"error":repr(e)}

def gebco(target):
    return {
        "grid_feature_info":gebco_feature_info(target,"gebco_latest_2"),
        "tid_feature_info":gebco_feature_info(target,"gebco_latest_tid_2"),
        "measured_mask_feature_info":gebco_feature_info(target,"gebco_latest_tid"),
        "note":"GEBCO_2026 latest WMS; TID response retained raw for provenance. Do not infer exact source survey from TID alone."
    }

def classify(target, gm, nc):
    mask=gm.get("products",{}).get("topo-mask",{})
    stats=mask.get("stats") if isinstance(mask,dict) else None
    direct=bool(stats and stats.get("nearest_grid_cell",{}).get("value_m") is not None)
    near_valid=None if not stats else stats.get("nearest_valid_cell")
    track1=nc.get("trackline_within_1km_count") if isinstance(nc,dict) else None
    track5=nc.get("trackline_within_5km_count") if isinstance(nc,dict) else None
    if direct:
        status="GMRT_HIGH_RES_AT_TARGET__MORPHOLOGY_ALLOWED_WITH_SOURCE_CAVEATS"
    elif near_valid and near_valid.get("distance_m") is not None and near_valid["distance_m"]<=1000:
        status="GMRT_HIGH_RES_WITHIN_1KM__TARGET_DIRECT_COVERAGE_UNCONFIRMED"
    else:
        status="NO_GMRT_HIGH_RES_TARGET_CELL__GLOBAL_GRID_ONLY_UNLESS_OTHER_SOURCE_FOUND"
    return {
        "coverage_status":status,
        "gmrt_high_res_direct":direct,
        "nearest_gmrt_high_res_distance_m":None if not near_valid else near_valid.get("distance_m"),
        "ncei_trackline_within_1km_count":track1,
        "ncei_trackline_within_5km_count":track5,
        "firewall":"NCEI_TRACKLINE_NE_SWATH__GMRT_MASK_NE_SOURCE_IDENTITY__GEBCO_CELL_NE_MEASURED_SURVEY_FOOTPRINT"
    }

def main():
    out={
        "artifact_id":"JANUS-KUSTO-OPEN-SEAFLOOR-SWEEP-2026-09-21-v1.0",
        "status":"RUNNING",
        "targets":[],
        "hard_rules":[
            "COVERAGE_FIRST_MORPHOLOGY_SECOND",
            "TRACKLINE_NE_SWATH",
            "GLOBAL_GRID_CELL_NE_MEASURED_CELL",
            "GMRT_HIGH_RES_MASK_NE_EXACT_SOURCE_IDENTITY",
            "NO_ARTIFICIALITY_CLASSIFICATION_FROM_TOPOGRAPHY_ALONE",
            "NO_RETARGETING_FROM_INTERESTING_SHAPE",
            "NEGATIVE_COVERAGE_RESULTS_PRESERVED"
        ]
    }
    for t in TARGETS:
        gm=gmrt(t)
        nc=ncei(t)
        gb=gebco(t)
        fp=ncei_footprints(t)
        out["targets"].append({
            **t,
            "gmrt":gm,
            "ncei":nc,
            "ncei_multibeam_footprints":fp,
            "gebco_2026_wms":gb,
            "classification":classify(t,gm,nc)
        })
    out["status"]="COMPLETE"
    outp=OUT/"JANUS-KUSTO-OPEN-SEAFLOOR-SWEEP-2026-09-21-v1.0.json"
    outp.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(outp)
    for t in out["targets"]:
        print("\n",t["id"],json.dumps(t["classification"],ensure_ascii=False))
        m=t["gmrt"].get("products",{}).get("topo-mask",{})
        if "stats" in m:
            print("GMRT mask:",json.dumps(m["stats"],ensure_ascii=False)[:3000])
        print("NCEI nearest:",json.dumps(t["ncei"].get("nearest_candidates",[])[:5],ensure_ascii=False)[:5000])
        print("NCEI exact footprints:",json.dumps(t["ncei_multibeam_footprints"],ensure_ascii=False)[:8000])
        print("GEBCO elevation:",json.dumps(t["gebco_2026_wms"]["grid_feature_info"],ensure_ascii=False)[:2000])
        print("GEBCO TID:",json.dumps(t["gebco_2026_wms"]["tid_feature_info"],ensure_ascii=False)[:2000])
        print("GEBCO parsed:", {"elevation_value":parse_feature_value(t["gebco_2026_wms"]["grid_feature_info"]), "tid_value":parse_feature_value(t["gebco_2026_wms"]["tid_feature_info"])})

if __name__=="__main__":
    main()
