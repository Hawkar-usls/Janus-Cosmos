#!/usr/bin/env python3
"""Independent geometric spectral-alias audit from ORIGINAL JWST CRDS WCS receipts.

The inputs are four GitHub Actions ZIP artifacts, each containing four official
JWST v1.19.1 CRDS jwst_1413.pmap detector-frame WCS footprints for source-
adjacent OPEN original NIRSpec MSA shutter coordinates. ZIP SHA256s are pinned.

Use centers of real shutter-coordinate, VIRTUAL BACKGROUND-source slit WCS
traces to estimate the wavelength in an adjacent slit whose detector X
coordinate equals the nominal target slit at 16004 A. This is a local linear
interpolation/extrapolation from official 15850.7933 A and 16004 A samples,
NOT an exact inverse WCS solution, calibrated spectrum, proof of leakage or
identification of any individual real neighboring galaxy.

No source flux, pixel DQ, PSF convolution, or true MSA source-target spectrum
is assessed. It is physically possible for two different slit/wavelength
coordinates to map to one detector coordinate; whether light is present is
entirely unmeasured here.
"""
import argparse
import hashlib
import json
import math
import statistics
import zipfile
from pathlib import Path

PRIMARY = 16004.0
SECONDARY = 15850.7932988
ARTIFACTS = {
    "A": {"artifact_id": 12186247363,
          "zip_sha256": "b5a1d60dcc5680ad889f7027cbebd13739e33c67b8367b86bbe1460ade8fe680",
          "workflow_run_id": 39251950005},
    "B": {"artifact_id": 12184714901,
          "zip_sha256": "1fecd6dca98b556c78be86d0564a145673b07823cc45419c155d14903ca7cb25",
          "workflow_run_id": 39251950005},
    "C": {"artifact_id": 12185693370,
          "zip_sha256": "b3f8bcd48b7410a41d87f69b34a40657fc81850ad4882c3efb75e63bfdbd8aa3",
          "workflow_run_id": 39251950005},
    "D": {"artifact_id": 12186316819,
          "zip_sha256": "bd0279d59364d7317318c54d2365de458d4fcbc1f3ae4e2176c767d112417257",
          "workflow_run_id": 39251950005},
}


def read_artifact(group, path):
    data=path.read_bytes()
    sha=hashlib.sha256(data).hexdigest()
    assert sha==ARTIFACTS[group]["zip_sha256"], (group,"ZIP_SHA256_CHANGED",sha)
    with zipfile.ZipFile(path) as bundle:
        name=f"janus_mocss_{group}_four_dither_gwcs_neighbor_map.json"
        assert bundle.namelist()==[name], bundle.namelist()
        receipt=json.loads(bundle.read(name))
    assert receipt["jwst_version"]=="1.19.1"
    assert receipt["crds_context"]=="jwst_1413.pmap"
    assert receipt["test_MSA_inputs_are_synthetic_background_source_zero"] is True
    assert receipt["real_MSA_shutter_coordinates_and_open_states"] is True
    assert receipt["exact_full_detector_illumination_for_all_open_shutters"] is False
    assert receipt["STPSF_weighted_detector_pixel_throughput"] is False
    assert abs(receipt["wavelengths_A"]["primary"]-PRIMARY)<1e-6
    assert abs(receipt["wavelengths_A"]["OIII_4960"]-SECONDARY)<1e-6
    assert len(receipt["records"])==4
    return receipt


def near_line(slit, label):
    line=slit["lines"][label]
    assert line["status"]=="FINITE_GWCS_SAMPLE" and line["on_detector_center"]
    return float(line["center_x_px"]),float(line["center_y_px"])


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--zip-dir",type=Path,default=Path("/tmp"))
    p.add_argument("--output",type=Path,default=Path("/tmp/janus_mocss_gwcs_aliases.json"))
    a=p.parse_args()
    exposures=[]
    source_sha=set()
    for g,spec in ARTIFACTS.items():
        artifact=a.zip_dir/f"janus_mocss_neighbor_{g.lower()}.zip"
        receipt=read_artifact(g,artifact)
        for exp in receipt["records"]:
            assert exp["group"]==g and exp["three_rows_open_in_original_MSA"]
            assert exp["status"]=="OFFICIAL_CRDS_WCS_SOURCE_NEIGHBOR_ON_DETECTOR"
            assert len(exp["samples_by_shutter"])==3
            assert all(s["wcs_slit_exists"] for s in exp["samples_by_shutter"])
            source_sha.add(exp["original_rate_sha256"])
            sorted_s=sorted(exp["samples_by_shutter"],key=lambda s:s["row"])
            assert all(s["quadrant"]==sorted_s[0]["quadrant"] and
                       s["column"]==sorted_s[0]["column"] for s in sorted_s)
            assert sorted_s[0]["row"]+1==sorted_s[1]["row"]==sorted_s[2]["row"]-1
            central=sorted_s[1]
            central_x,central_y=near_line(central,"primary")
            row_aliases=[]
            for relative,neighbor in [(-1,sorted_s[0]),(1,sorted_s[2])]:
                at_primary=near_line(neighbor,"primary")
                at_other=near_line(neighbor,"OIII_4960")
                dx_dA=(at_primary[0]-at_other[0])/(PRIMARY-SECONDARY)
                dy_dA=(at_primary[1]-at_other[1])/(PRIMARY-SECONDARY)
                assert 0.05<dx_dA<1 and math.isfinite(dy_dA)
                predicted_wavelength=PRIMARY+(central_x-at_primary[0])/dx_dA
                predicted_y=at_primary[1]+(predicted_wavelength-PRIMARY)*dy_dA
                yerr=predicted_y-central_y
                # A box at λ16004 is approx ~2 detector pixels in x; pure
                # closest-edge gap quantifies separation at the SAME λ.
                box_a=central["lines"]["primary"]["detector_box"]
                box_b=neighbor["lines"]["primary"]["detector_box"]
                gap=max(0.0,float(box_a["x_min"])-float(box_b["x_max"]),
                        float(box_b["x_min"])-float(box_a["x_max"]))
                row_aliases.append({
                    "neighbor_relative_row":relative,
                    "neighbor_row":neighbor["row"],
                    "same_line_center_x_gap_pixels":at_primary[0]-central_x,
                    "same_line_footprint_x_gap_pixels":gap,
                    "local_neighbor_dispersion_dx_pixels_per_angstrom":dx_dA,
                    "first_order_neighbor_true_lambda_for_primary_x_A":predicted_wavelength,
                    "first_order_neighbor_wavelength_alias_shift_A":predicted_wavelength-PRIMARY,
                    "first_order_residual_y_pixels":yerr,
                    "linear_approximation_not_exact_inverse_GWCS":True,
                    "actual_neighbor_photon_flux_measured":False,
                })
            exposures.append({
                "group":g,
                "index":int(exp["exposure_index"]),
                "original_rate_sha256":exp["original_rate_sha256"],
                "original_msa_sha256":exp["original_msa_sha256"],
                "main_shutter_row":central["row"],
                "main_shutter_column":central["column"],
                "target_line_center_pixel_xy":[central_x,central_y],
                "adjacent_open_shutter_aliases":row_aliases,
            })
    assert [e["index"] for e in exposures]==list(range(1,17))
    assert len(source_sha)==16,"RATE SHA indicates duplicate exposures"
    plus=[r for e in exposures for r in e["adjacent_open_shutter_aliases"] if r["neighbor_relative_row"]==1]
    minus=[r for e in exposures for r in e["adjacent_open_shutter_aliases"] if r["neighbor_relative_row"]==-1]

    # Authors observed two dispersion orientations roughly 89deg apart and
    # already required an identical SKY host plus <=6A line consistency.
    # Wavelength-only consistency does NOT substitute for that sky-host test.
    cross_orientation=[]
    for shift in (-1,1):
        first=[r["first_order_neighbor_true_lambda_for_primary_x_A"]
               for e in exposures if e["group"] in ("A","B")
               for r in e["adjacent_open_shutter_aliases"]
               if r["neighbor_relative_row"]==shift]
        second=[r["first_order_neighbor_true_lambda_for_primary_x_A"]
                for e in exposures if e["group"] in ("C","D")
                for r in e["adjacent_open_shutter_aliases"]
                if r["neighbor_relative_row"]==shift]
        assert len(first)==len(second)==8
        cross_orientation.append({
            "neighbor_relative_row_in_each_instrument_orientation":shift,
            "orientation_1_AB_wavelength_median_A":statistics.median(first),
            "orientation_2_CD_wavelength_median_A":statistics.median(second),
            "absolute_wavelength_median_difference_A":abs(statistics.median(first)-statistics.median(second)),
            "wavelength_only_within_paper_6A_tolerance":abs(statistics.median(first)-statistics.median(second))<=6,
            "same_physical_sky_location_between_orientations_proven":False,
            "paper_already_tested_host_shutter_and_wavelength_cross_orientation":True,
        })

    def summarise(r):
        return {key:{"min":min(z[key] for z in r),"median":statistics.median(z[key] for z in r),
                     "max":max(z[key] for z in r)} for key in (
            "same_line_center_x_gap_pixels", "same_line_footprint_x_gap_pixels",
            "local_neighbor_dispersion_dx_pixels_per_angstrom",
            "first_order_neighbor_true_lambda_for_primary_x_A",
            "first_order_neighbor_wavelength_alias_shift_A",
            "first_order_residual_y_pixels")}

    out={
        "schema":"janus-mocss-original-crds-16rate-adjacent-shutter-spectral-alias-v1",
        "status":"GEOMETRIC_NEIGHBORING_OPEN_SHUTTER_WAVELENGTH_ALIAS_POSSIBLE_NOT_PHOTOMETRICALLY_DEMONSTRATED",
        "primary_assumed_measured_wavelength_A":PRIMARY,
        "secondary_crds_sample_wavelength_A":SECONDARY,
        "official_wcs":"jwst 1.19.1 + jwst_1413.pmap",
        "input_provenance":ARTIFACTS,
        "read_inputs_all_four_zip_original_digests_verified":True,
        "original_rate_exposures_distinct_sha256":len(source_sha),
        "exposure_count":len(exposures),
        "official_msa_original_three_source_adjacent_rows_open_every_exposure":True,
        "virtual_background_source_wcs_warning":"The official pipeline WCS was calculated on a synthetic 9-slit MSA with real OPEN shutter geometry and source_id=0. It is not the author's complete science spec2 WCS / cross-contamination solution.",
        "method":"In each adjacent OPEN row, linearly estimate lambda solving x_neighbor(lambda)=x_target(16004A), using two official detector WCS x(lambda) evaluations at 15850.7933 and 16004 A. Then evaluate y residual using same local slope. A geometric degeneracy is possible; no photon flux or astrophysical identification follows.",
        "summary_nearest_open_row_minus_1":summarise(minus),
        "summary_nearest_open_row_plus_1":summarise(plus),
        "cross_orientation_wavelength_only_not_host_solution":cross_orientation,
        "per_exposure":exposures,
        "science_cautions":{
            "actual_open_shutter_flux_overlap":None,
            "full_MSAs_all_open_illumination_model":None,
            "real_sky_sources_projected_onto_neighbors":None,
            "source_calibrated_spectral_line_identification":None,
            "wavelength_solution_is_exact_inverse_GWCS":False,
            "author_cross_orientation_host_test_replicated":False,
            "line_16004_must_be_neighbor_contamination":False,
            "redshift":None,"physical_propagation_vector":None,"trail_size":None,
            "interpretation":"Original paper Sec. III.2 explicitly anticipated slitless wavelength/sky degeneracy and used two near-perpendicular orientations plus common sky-shutter and <=6A matching. Our official-WCS audit confirms the per-orientation geometric degeneracy but does NOT reverse or replicate the author's two-orientation source matching. Both sky and wavelength must be matched; first-order same-sign 1-row aliases match by wavelength alone, which is insufficient. No actual neighboring photons or motion measured."
        },
        "next_test":"Reproduce author Sec. III.2 two-orientation common-SKY-host-and-<6A gate independently with all open-shutter WCS and measured line centroids, before entertaining adjacent-shutter alias alternatives; then fit original calibrated science detector data with STPSF leakage and sky DQ/background."
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status":out["status"],"num_records":len(exposures),
        "minus_row_alias_A":out["summary_nearest_open_row_minus_1"]["first_order_neighbor_true_lambda_for_primary_x_A"],
        "plus_row_alias_A":out["summary_nearest_open_row_plus_1"]["first_order_neighbor_true_lambda_for_primary_x_A"],
        "max_abs_y_difference_px":max(abs(z["first_order_residual_y_pixels"]) for z in minus+plus),
        "output":str(a.output)},indent=2))

if __name__=="__main__":
    main()
