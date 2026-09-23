#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"
OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-KUSTO-USGS-natural-controls-source-gate/1.2"}

DATASETS=[
    {
        "id":"USGS_MIAMI_KEY_BISCAYNE_POCKMARKS",
        "sciencebase_id":"5a96f5cee4b06990606c4ff0",
        "doi":"10.5066/F72J6B4Z",
        "known_child_item_ids":[]
    },
    {
        "id":"USGS_LAKE_CRESCENT_2016",
        "sciencebase_id":"586d3165e4b0f5ce109faa51",
        "doi":"10.5066/F7B56GW5",
        "known_child_item_ids":[
            "5eea74cb82ce3bd58d8572af",
            "5eea74d882ce3bd58d8572dd"
        ]
    }
]

def getj(url,params=None):
    r=requests.get(url,params=params,headers=UA,timeout=90)
    r.raise_for_status()
    return r.json()

def file_role(name):
    n=name.lower()
    roles=[]
    if any(x in n for x in ["bathy","bathym","dem","depth","elevation"]):
        roles.append("BATHYMETRY")
    if any(x in n for x in ["backscatter","back_scatter","reflectivity","reflect"]):
        roles.append("BACKSCATTER")
    if any(x in n for x in [".xml",".txt","metadata","readme","fgdc"]):
        roles.append("METADATA_OR_TEXT")
    if any(x in n for x in [".tif",".tiff",".asc",".grd",".xyz",".nc",".zip"]):
        roles.append("NUMERIC_RASTER_OR_ARCHIVE_PRODUCT")
    return roles

def norm_file(f,item_id,item_title):
    name=f.get("name") or ""
    return {
        "item_id":item_id,
        "item_title":item_title,
        "name":name,
        "url":f.get("url"),
        "size":f.get("size"),
        "contentType":f.get("contentType"),
        "checksum":f.get("checksum"),
        "roles":file_role(name)
    }

def children(parent_id):
    attempts=[
        ("parentId_param",{"parentId":parent_id,"format":"json","max":1000}),
        ("parentId_filter",{"filter":f"parentId={parent_id}","format":"json","max":1000})
    ]
    audit=[]
    for mode,params in attempts:
        try:
            j=getj("https://www.sciencebase.gov/catalog/items",params=params)
            items=j.get("items") or []
            exact=[x for x in items if str(x.get("parentId") or "")==str(parent_id)]
            audit.append({"mode":mode,"returned":len(items),"exact_children":len(exact)})
            if exact:
                return exact,audit
        except Exception as e:
            audit.append({"mode":mode,"error":type(e).__name__+": "+str(e)})
    return [],audit

def crawl(root_id,max_depth=3):
    queue=[(root_id,0)]
    seen=set()
    items=[]
    query_audit=[]
    while queue:
        iid,depth=queue.pop(0)
        if iid in seen or depth>max_depth:
            continue
        seen.add(iid)
        j=getj(f"https://www.sciencebase.gov/catalog/item/{iid}?format=json")
        items.append(j)
        if j.get("hasChildren") or not (j.get("files") or []):
            kids,audit=children(iid)
            query_audit.append({"parent_id":iid,"audit":audit,"exact_children":len(kids)})
            for k in kids:
                kid=k.get("id")
                if kid and kid not in seen:
                    queue.append((kid,depth+1))
    return items,query_audit

outsets=[]
for ds in DATASETS:
    root_url=f"https://www.sciencebase.gov/catalog/item/{ds['sciencebase_id']}?format=json"
    try:
        items,query_audit=crawl(ds["sciencebase_id"])
        root=items[0]
        recovered_ids={str(x.get("id")) for x in items}
        explicit_children=[]
        for kid in ds.get("known_child_item_ids") or []:
            if str(kid) in recovered_ids:
                continue
            cj=getj(f"https://www.sciencebase.gov/catalog/item/{kid}?format=json")
            items.append(cj)
            recovered_ids.add(str(kid))
            explicit_children.append(str(kid))

        files=[]
        item_summaries=[]
        for item in items:
            iid=item.get("id")
            title=item.get("title")
            fs=[norm_file(f,iid,title) for f in (item.get("files") or [])]
            files.extend(fs)
            item_summaries.append({
                "id":iid,
                "title":title,
                "parentId":item.get("parentId"),
                "hasChildren":item.get("hasChildren"),
                "files":len(fs)
            })

        if explicit_children:
            status="BOUND_WITH_EXPLICIT_LEGACY_CHILD_BINDING"
        elif len(items)>1:
            status="BOUND_WITH_CHILD_TRAVERSAL"
        elif files:
            status="BOUND"
        else:
            status="BOUND_NO_FILES_LISTED"

        outsets.append({
            **ds,
            "sciencebase_api":root_url,
            "title":root.get("title"),
            "root_hasChildren":root.get("hasChildren"),
            "explicit_child_ids_from_independent_usgs_metadata":explicit_children,
            "items_traversed":len(items),
            "child_query_audit":query_audit,
            "item_summaries":item_summaries,
            "files_total":len(files),
            "bathymetry_candidates":[f for f in files if "BATHYMETRY" in f["roles"]],
            "backscatter_candidates":[f for f in files if "BACKSCATTER" in f["roles"]],
            "all_files":files,
            "status":status
        })
    except Exception as e:
        outsets.append({
            **ds,
            "sciencebase_api":root_url,
            "status":"ERROR",
            "error":type(e).__name__+": "+str(e)
        })

out={
    "artifact_id":"JANUS-KUSTO-USGS-NATURAL-CONTROLS-SOURCE-GATE-RUN-2026-09-23-v1.2",
    "supersedes":"JANUS-KUSTO-USGS-NATURAL-CONTROLS-SOURCE-GATE-RUN-2026-09-23-v1.1",
    "stage":"G0_SOURCE_AND_PROVENANCE_INVENTORY",
    "depth_values_read":False,
    "detector_run":False,
    "datasets":outsets,
    "hard_rules":[
        "ScienceBase metadata/file inventory only; no raster values read.",
        "Parent item with zero files is not treated as data absence.",
        "Legacy child IDs may be bound from independent USGS metadata explicitly citing the same DOI; such bindings are recorded separately.",
        "File-name role classification is routing metadata, not scientific classification.",
        "Processed products may calibrate detection but cannot independently validate their own processing lineage."
    ],
    "next_gate":"FREEZE_NUMERIC_PRODUCTS_AND_NATIVE_RESOLUTION_BEFORE_BLIND_DETECTION",
    "claim_ceiling":"SOURCE_BINDING_ONLY"
}

p=OUT/"JANUS-KUSTO-USGS-NATURAL-CONTROLS-SOURCE-GATE-RUN-2026-09-23-v1.2.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
    "artifact_id":out["artifact_id"],
    "datasets":[{
        "id":d["id"],
        "status":d["status"],
        "items_traversed":d.get("items_traversed"),
        "files_total":d.get("files_total"),
        "bathymetry_candidates":len(d.get("bathymetry_candidates") or []),
        "backscatter_candidates":len(d.get("backscatter_candidates") or [])
    } for d in outsets],
    "depth_values_read":False
},indent=2))
