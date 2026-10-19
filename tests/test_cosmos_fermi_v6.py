"""V6 offline OGIP RSP2 decoding, on-off likelihood and science nonpromotion."""

import copy
import json
import math
import struct
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workspace"))
sys.path.insert(0, str(ROOT / "scripts"))

from janus_cosmos_fermi_raw_v4 import HDU, sha256
from janus_cosmos_fermi_response_v6 import (
    AuditError, decode_drm, weighted_drm, profile_background,
    poisson_deviance, photon_bin_flux, folded_rate, cstat_onoff,
    fit_segment, SCHEMA as GBM_SCHEMA,
)
from janus_cosmos_fermi_lat_diagnostics_v6 import (
    exposure_proxy, SCHEMA as LAT_SCHEMA,
)
from janus_cosmos_fermi_joint_v6 import run as joint_run, SCHEMA as JOINT_SCHEMA
from janus_cosmos_simulation_audit import digest
from fetch_fermi_v6_inputs import safe_archive_url

BASE = ROOT / "data/cosmos"
GBM = BASE / "COSMOS-ALI-V6-GBM-RSP2-ONOFF-SPECTRAL-2026-10-19.json"
LAT = BASE / "COSMOS-ALI-V6-LAT-SOURCE-CONTROL-EXPOSURE-2026-10-19.json"
JOINT = BASE / "COSMOS-ALI-V6-JOINT-OBSERVATION-READINESS-2026-10-19.json"
V4 = BASE / "COSMOS-ALI-V4-FERMI-RAW-EVENT-DIAGNOSTICS-2026-10-18.json"
V5 = BASE / "COSMOS-ALI-V5-CALDB-GBM-LAT-INSTRUMENT-PULSE-2026-10-19.json"
MANIFEST = BASE / "COSMOS-ALI-V5-FERMI-LAT-GBM-CALDB-MANIFEST-2026-10-19.json"


def example_rsp2(correction: float = 1., start: float = 0., stop: float = 10.) -> HDU:
    rows, heap = [], bytearray()
    for i in range(140):
        base = len(heap)
        heap.extend(struct.pack(">h", 1))  # 1-based PHA F_CHAN
        heap.extend(struct.pack(">h", 128))
        heap.extend(struct.pack(">128f", *([correction] * 128)))
        rows.append(struct.pack(">ffhiiiiii", 5. + i, 6. + i, 1,
                                1, base, 1, base + 2, 128, base + 4))
    h = {"EXTNAME": "SPECRESP MATRIX", "NAXIS1": 34, "NAXIS2": 140,
         "DETCHANS": 128, "TFORM1": "1E", "TFORM2": "1E", "TFORM3": "1I",
         "TFORM4": "PI(128)", "TFORM5": "PI(128)", "TFORM6": "PE(128)",
         "TSTART": start, "TSTOP": stop}
    return HDU(h, b"".join(rows) + bytes(heap))


class TestV6OGIP(unittest.TestCase):
    def test_ogip_heap_matrix_decodes_and_one_based_channel_maps_to_zero_based(self):
        matrix = decode_drm(example_rsp2())
        self.assertEqual(matrix["channels"], 128)
        self.assertEqual(matrix["true_energy_bins"], 140)
        self.assertEqual(matrix["response_cm2"][0][0], 1.)
        self.assertEqual(matrix["response_cm2"][139][127], 1.)
        self.assertEqual(matrix["true_energy_edges_keV"][0], (5., 6.))

    def test_corrupted_fits_heap_pointer_fails_closed(self):
        h = example_rsp2()
        data = bytearray(h.content)
        struct.pack_into(">ii", data, 26, 128, 9999999)
        with self.assertRaisesRegex(AuditError, "heap"):
            decode_drm(HDU(h.header, bytes(data)))

    def test_weighted_attitude_rsp2_exact_integral(self):
        one = example_rsp2(2., 0., 2.)
        two = example_rsp2(6., 2., 10.)
        energy, drm, info = weighted_drm([one, two], 1., 5.)
        self.assertEqual(len(energy), 140)
        self.assertAlmostEqual(drm[0][0], (2. * 1. + 6. * 3.) / 4.)
        self.assertEqual([i["matrix_number"] for i in info["rsp2_matrices_used"]], [1, 2])
        self.assertEqual(info["rsp2_total_weight"], 1.)
        with self.assertRaisesRegex(AuditError, "not entirely covered"):
            weighted_drm([one, two], 1., 12.)

    def test_poisson_onoff_profile_is_finite_and_null_exact(self):
        self.assertAlmostEqual(profile_background(20, 30, 0., 2.), 50 / 3.)
        self.assertEqual(profile_background(0, 0, 0., 2.), 0.)
        self.assertAlmostEqual(poisson_deviance(10, 10.), 0.)
        self.assertGreater(poisson_deviance(10, 15.), 0.)
        self.assertAlmostEqual(cstat_onoff([20] + [0] * 127, [30] + [0] * 127,
                                          [0.] * 128, 0., 1., 2., [True] + [False] * 127),
                               poisson_deviance(20, 50 / 3) + poisson_deviance(30, 100 / 3))

    def test_injected_photon_spectrum_forward_folds_to_channel(self):
        flux = photon_bin_flux([(100., 200.), (200., 300.)], -1., 500.)
        expected = [0.] * 128
        expected[0] = 2 * flux[0]
        expected[127] = 3 * flux[1]
        response = [[2.] + [0.] * 127, [0.] * 127 + [3.]]
        calculated = folded_rate([(100., 200.), (200., 300.)], response, -1., 500.)
        self.assertAlmostEqual(calculated[0], expected[0])
        self.assertAlmostEqual(calculated[127], expected[127])

    def test_pilot_three_detector_source_improves_null_only_as_diagnostic(self):
        energy = [(10. + i * 100., 20. + i * 100.) for i in range(140)]
        response = [[1.] + [0.] * 127 for _ in energy]
        detectors = [{"name": label, "energy": energy, "response": response,
                      "on": [100] + [0] * 127, "off": [10] + [0] * 127,
                      "duration": 1., "bg_duration": 15., "allowed": [True] + [False] * 127}
                     for label in ("n6", "n7", "b1")]
        fit = fit_segment(detectors)
        self.assertGreater(fit["source_vs_no_source_diagnostic_delta_deviance"], 0)
        self.assertTrue(fit["not_a_discovery_significance_or_calibrated_p_value"])


class TestV6RealReceipts(unittest.TestCase):
    def test_unchanged_manifest_links_to_historical_v5_sources(self):
        m = json.loads(MANIFEST.read_text())
        self.assertEqual(m["lat_caldb"]["irf_name"], "P8R3_SOURCE_V3")
        self.assertEqual(len(m["lat_caldb"]["caldb_files"]), 9)
        self.assertEqual(sum(len(x["detectors"]) * 3 for x in m["gbm_bursts"]), 27)
        self.assertEqual(m["v4_frozen_observed_json_sha256"], sha256(V4))

    def test_archive_fetcher_refuses_unallowlisted_urls(self):
        self.assertTrue(safe_archive_url("https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2009/bn090510016/current/a.fit", True))
        self.assertFalse(safe_archive_url("https://heasarc.gsfc.nasa.gov.evil.net/FTP/fermi/data/gbm/triggers/2009/a.fit", True))
        self.assertFalse(safe_archive_url("http://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2009/a.fit", True))
        self.assertFalse(safe_archive_url("https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2009/a.fit?token=a", True))

    def test_actual_gbm_rsp2_fit_nonpromotion_and_hot_burst_residual(self):
        x = json.loads(GBM.read_text())
        self.assertEqual(x["schema"], GBM_SCHEMA)
        self.assertEqual(sum(len(b["time_resolved_spectral_fit_segments"]) for b in x["gbm_bursts"]), 17)
        stress = {b["grb"]: sum(row["onoff_spectral_grid_fit"]["spectral_mismatch_stress_gate"] == "HIGH_RESIDUAL_MODEL_REVISE" for row in b["time_resolved_spectral_fit_segments"]) for b in x["gbm_bursts"]}
        self.assertEqual(stress, {"090510": 0, "090902B": 1, "130427A": 2})
        self.assertIsNone(x["observational_joint_p_value"])
        self.assertEqual(x["scientific_gates"]["GBM_MODEL_COMPLETENESS"], "NOT_ESTABLISHED_SINGLE_CPL_NOT_BAND_PLUS_PL")

    def test_actual_lat_roi_and_control_exposure_not_pvalue(self):
        x = json.loads(LAT.read_text())
        self.assertEqual(x["schema"], LAT_SCHEMA)
        self.assertEqual(x["counts"]["v4_lat_source_region_events_recovered_exactly"], 175)
        self.assertEqual([b["five_degree_source_region_candidate_total"] for b in x["GRB_sources"]], [32, 67, 76])
        self.assertEqual([next(w for w in b["three_time_windows"] if w["time_window"] == "prompt")["source_circle_5deg_count"] for b in x["GRB_sources"]], [19, 42, 47])
        for b in x["GRB_sources"]:
            self.assertEqual(b["event_type_3_partitions_decoded_count"], b["five_degree_source_region_candidate_total"])
            for window in b["three_time_windows"]:
                self.assertFalse(window["is_real_calibrated_lat_exposure"])
        self.assertIsNone(x["real_joint_observed_liv_p_value"])

    def test_exposure_proxy_uses_gti_and_livetime_nonnegative(self):
        grid = {"ENERG_LO": [1000.], "ENERG_HI": [10000.],
                "CTHETA_LO": [.2], "CTHETA_HI": [1.], "EFFAREA": [.5]}
        spacecraft = [{"START": 0., "STOP": 10., "LIVETIME": 5., "DATA_QUAL": 1,
                       "LAT_CONFIG": 1, "IN_SAA": False, "RA_SCZ": 0., "DEC_SCZ": 0.}]
        res = exposure_proxy(spacecraft, [{"START": 0., "STOP": 10.}], (2., 8.), 0., 0.,
                             {"FRONT": grid, "BACK": grid})
        self.assertAlmostEqual(res["nearest_FB_aeff_times_livetime_cm2_s_by_energy_MeV"]["1000"], 30000.)

    def test_joint_gate_cannot_turn_toy_mimic_into_real_liv_pvalue(self):
        b = json.loads(GBM.read_text()); l = json.loads(LAT.read_text())
        p = json.loads(V5.read_text()); r = json.loads(V4.read_text())
        actual = json.loads(JOINT.read_text())
        self.assertEqual(joint_run(b, l, p, r), actual)
        self.assertEqual(actual["schema"], JOINT_SCHEMA)
        self.assertIsNone(actual["observational_liv_p_value"])
        self.assertIsNone(actual["conditional_hypothetical_GeV_lag_prior"])
        self.assertEqual(actual["conditional_source_mimic_negative_control_SYNTHETIC_ONLY"]["unconstrained_eta"]["twice_delta_log_likelihood"], 0.)
        self.assertEqual(actual["scientific_gates"]["OBSERVED_MULTIBURST_LIV_LRT"], "BLOCKED_NO_FREQUENTIST_COVERAGE_OR_VALID_SOURCE_MODEL")
        for value in (b, l, actual):
            copy_value = dict(value)
            claimed = copy_value.pop("result_sha256")
            self.assertEqual(claimed, digest(copy_value))
        broken = copy.deepcopy(b)
        broken["scientific_gates"]["OBSERVED_LAT_GBM_LIV_HYPOTHESIS_TEST"] = "DISCOVERED"
        with self.assertRaises(AuditError):
            joint_run(broken, l, p, r)


if __name__ == "__main__":
    unittest.main()
