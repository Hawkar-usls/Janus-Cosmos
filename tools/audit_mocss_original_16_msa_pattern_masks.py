#!/usr/bin/env python3
"""Audit ORIGINAL NIRSpec GO-3290 MSA shutter topology by exact (MSAMETID,PATT_NUM).

For each of the 16 independently SHA-256-verified archival RATE files, read
its *actual* MSAMETID and PATT_NUM FITS headers. Select the original matching
MSA SHUTTER_INFO rows, respecting DITHER_POINT_INDEX, OPEN shutter state and
unique (quadrant,row,column) sites. Output counts and exact binary 365x171
open-shutter raster for each of 4 quadrants and each of 16 dither exposures.

This constitutes a PHYSICAL-MSA-PLANE shutter mask, NOT the illumination map
on the dispersed detector. Every open shutter fans into many detector pixels
according to wavelength-dependent CRDS/WCS; no 2D detector weight can be
claimed until that projection, leakage, overlaps and throughput are handled.
No original FITS bytes are changed.
"""
import argparse, concurrent.futures, hashlib, json
from pathlib import Path
from urllib.parse import quote
import numpy as np
from astropy.io import fits
from fetch_jwst_mocss_mast_subset import BASE,ROOT_URI,read_and_store

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"domains/space/MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json"

def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for v in iter(lambda:f.read(4*1024*1024),b""):h.update(v)
    return h.hexdigest()

def fetch_verified(args):
    name,hash_expected,root=args
    dest=root/name
    read_and_store(BASE+quote(ROOT_URI+name,safe=""),dest)
    got=sha256(dest)
    if got!=hash_expected:raise RuntimeError(f"SHA-256 archive mismatch {name}")
    return str(dest)

def normalized(s):
    return str(s).strip().upper()

def slice_input(row,raw,msa):
    with fits.open(raw,memmap=False) as f:
        hdr=f[0].header
        # Header not assumed to have same MSA id as some unwrapped receipt:
        # original FITS SHUTTER_INFO uses signed 16-bit.
        header_id=int(hdr["MSAMETID"])
        patt=int(hdr["PATT_NUM"])
        assert hdr["DETECTOR"].lower()==row["expected_science_detector_from_WCS"]
    signed=int(row["source_msametid_signed16"])
    assert header_id % 65536==signed % 65536,(header_id,signed)
    with fits.open(msa,memmap=False) as f:
        tab=f["SHUTTER_INFO"].data
        names=set(tab.names)
        assert {"msa_metadata_id","dither_point_index","shutter_quadrant","shutter_row",
                "shutter_column","shutter_state","background","primary_source"}.issubset(names)
        by_id=tab[tab["msa_metadata_id"]==signed]
        chosen=by_id[by_id["dither_point_index"]==patt]
        if not len(chosen):
            raise ValueError(f"No MSA rows for FITS metadata_id={signed} pattern={patt}, exposure={row['exposure_index']}")
        n_msid=len(by_id);n_chosen=len(chosen)
        rasters={}
        cfg={}
        shutter=row["detectors"][row["expected_science_detector_from_WCS"]]["selected_virtual_shutter"]
        expected_coordinate=[int(shutter["quadrant"]),int(shutter["row"]),int(shutter["column"])]
        all_unique=set()
        for quad in range(1,5):
            subset=chosen[np.asarray(chosen["shutter_quadrant"],dtype=int)==quad]
            raster=np.zeros((365,171),dtype=np.uint8)
            mask_states={}
            duplicate_sites=0
            for t in subset:
                r,c=int(t["shutter_row"]),int(t["shutter_column"])
                assert 1<=r<=365 and 1<=c<=171
                flag=normalized(t["shutter_state"])
                mask_states[flag]=mask_states.get(flag,0)+1
                key=(quad,r,c)
                if key in all_unique:duplicate_sites+=1
                all_unique.add(key)
                if flag=="OPEN":raster[r-1,c-1]=1
            row_counts=raster.sum(axis=1)
            col_counts=raster.sum(axis=0)
            rasters[str(quad)]=raster
            cfg[str(quad)]={
                "rows_in_original_SHUTTER_INFO":len(subset),
                "unique_sites":len(set((int(t["shutter_row"]),int(t["shutter_column"])) for t in subset)),
                "duplicate_sites":duplicate_sites,
                "states":mask_states,
                "unique_OPEN_sites":int(raster.sum()),
                "open_fraction_of_365x171":float(raster.mean()),
                "fraction_rows_open_median":float(np.median(row_counts/171)),
                "fraction_rows_open_05_95":np.quantile(row_counts/171,[.05,.95]).tolist(),
                "fraction_columns_open_median":float(np.median(col_counts/365)),
                "shutter_plane_raster_sha256":hashlib.sha256(raster.tobytes()).hexdigest()
            }
        q,r,c=expected_coordinate
        raster=rasters[str(q)]
        source_line={"quadrant":q,"row":r,"column":c,
                     "center_is_marked_open":bool(raster[r-1,c-1]),
                     "three_dispersion_rows":[
                         {"row":j,"marked_open":bool(raster[j-1,c-1])}
                         for j in range(max(1,r-1),min(365,r+1)+1)],
                     "row_open_fraction_near_frozen_source":float(raster[r-1].mean()),
                     "row_open_count":int(raster[r-1].sum())}
        return {"exposure_index":row["exposure_index"],"group":row["group"],
                "detector":row["expected_science_detector_from_WCS"],
                "archive_rate_file":Path(raw).name,
                "original_RATE_sha256":row["detectors"][row["expected_science_detector_from_WCS"]]["sha256"],
                "original_MSA_file":Path(msa).name,
                "actual_RATE_FITS_MSAMETID":header_id,
                "matching_original_MSA_signed_metadata_ID":signed,
                "actual_RATE_FITS_PATT_NUM":patt,
                "MSA_all_rows_for_metadata_id_before_pattern_filter":n_msid,
                "MSA_rows_for_metadata_id_AND_dither_pattern":n_chosen,
                "original_shutter_record_state_by_quadrant":cfg,
                "frozen_source_nearest_shutter":source_line,
                "original_MSA_contains_heliocentric_emission_line_source":None,
                "map_at_detector_pixel":False},rasters

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--data-dir",type=Path,default=Path("/tmp/janus_mocss_msa_original_topology"))
    p.add_argument("--output",type=Path,default=Path("/tmp/janus_mocss_16_original_MSA_pattern_topology.json"))
    p.add_argument("--raster",type=Path,default=Path("/tmp/janus_mocss_16_original_MSA_pattern_topology.npz"))
    p.add_argument("--workers",type=int,default=4)
    a=p.parse_args()
    rec=json.loads(SOURCE.read_text())
    assert len(rec["exposures"])==16
    a.data_dir.mkdir(parents=True,exist_ok=True)
    msa={k:v["sha256"] for k,v in rec["two_full_msa_reference_files"].items()}
    targets={}
    for e in rec["exposures"]:
        det=e["expected_science_detector_from_WCS"];v=e["detectors"][det]
        targets[v["rate_file"]]=v["sha256"]
    targets.update(msa)
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as ex:
        list(ex.map(fetch_verified,[(k,v,a.data_dir) for k,v in sorted(targets.items())]))
    output={"schema":"janus-cosmos-nirspec-16-original-exposure-matched-MSA-plane-raster-v1",
            "original_source_receipt":"domains/space/MOCSS_16_SKY_REGISTERED_RATE_RECEIPT_2026-10-16.json",
            "original_RATE_SHA_verified_count":16,
            "original_MSA_SHA_verified_count":2,
            "matching_msa_criterion":"SHUTTER_INFO.MSA_METADATA_ID=original RATE.MSAMETID AND DITHER_POINT_INDEX=original RATE.PATT_NUM",
            "raster_shape_per_quadrant":[365,171],
            "raster_exposure_count":16,"original_FITS_mutated":False,
            "science":"ORIGINAL_MSA_PLANE_BINARY_APERTURE_MAP_ONLY",
            "original_MSA_config_full_detector_illumination_projected":False,
            "crds_NIRSpec_wavelength_detector_projection":False,
            "source_spectral_flux_or_snr_calibrated":False,
            "source_class":None,"redshift":None,
            "exposures":[]}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.raster.parent.mkdir(parents=True,exist_ok=True)
    matrices={}
    for e in rec["exposures"]:
        name=("jw03290001001_05_msa.fits" if e["orientation"]==1
              else "jw03290002001_01_msa.fits")
        det=e["expected_science_detector_from_WCS"]
        source=e["detectors"][det]["rate_file"]
        info,rasters=slice_input(e,a.data_dir/source,a.data_dir/name)
        for q,mask in rasters.items():
            matrices[f"exposure_{e['exposure_index']:02d}_quad_{q}"]=mask
        output["exposures"].append(info)
        a.output.write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
        print("EXPOSURE",e["exposure_index"],"group",e["group"],
              "MSAMETID",info["actual_RATE_FITS_MSAMETID"],
              "PATT_NUM",info["actual_RATE_FITS_PATT_NUM"],
              "center_open",info["frozen_source_nearest_shutter"]["center_is_marked_open"],flush=True)
    np.savez_compressed(a.raster,**matrices)
    output["original_64_quadrant_rasters_npz_sha256"]=sha256(a.raster)
    a.output.write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
    print("PASS: 16 original FITS exact MSA-id+PATT_NUM configurations, 64 physical aperture rasters; DETECTOR_PROJECTION=NOT_PERFORMED",flush=True)
if __name__=="__main__":main()
