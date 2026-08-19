from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _f(d: dict[str, Any] | None, key: str) -> float | None:
    if not d:
        return None
    try:
        v = d.get(key)
        return None if v is None else float(v)
    except Exception:
        return None


def _catalog_selected(group: dict[str, Any], key: str) -> dict[str, Any] | None:
    return group.get("catalog_queries", {}).get(key, {}).get("selected")


def _sdss_class(group: dict[str, Any]) -> int | None:
    s = _catalog_selected(group, "SDSS_DR16")
    if not s:
        return None
    try:
        v = s.get("features", {}).get("class")
        return None if v is None else int(v)
    except Exception:
        return None


def _nature(name: str, group: dict[str, Any]) -> dict[str, Any]:
    c = group.get("classification", {})
    gd = c.get("gaia_derived") or {}
    sdss_class = _sdss_class(group)
    gaia = _catalog_selected(group, "GAIA_DR3")
    wise = _catalog_selected(group, "ALLWISE")
    ps1 = _catalog_selected(group, "PANSTARRS_DR1")
    tm = _catalog_selected(group, "2MASS_PSC")

    # SDSS PhotoType convention: 3 = galaxy, 6 = star.
    if sdss_class == 3 and gaia is None:
        return {
            "nature": "GALAXY",
            "confidence": "MEDIUM_HIGH",
            "plain_ru": "Это очень вероятно галактика: SDSS морфологически помечает источник как galaxy (class=3), а Pan-STARRS независимо видит объект в той же позиции.",
            "evidence": {
                "sdss_photo_class": 3,
                "sdss_mapping": "3=GALAXY",
                "panstarrs_counterpart": ps1 is not None,
                "gaia_counterpart_within_2arcsec": False,
                "allwise_counterpart_within_2arcsec": wise is not None,
                "twomass_counterpart_within_2arcsec": tm is not None,
            },
            "ceiling": "MORPHOLOGICAL_GALAXY_CLASSIFICATION_NOT_SPECTROSCOPIC_REDSHIFT"
        }

    if c.get("class") == "GALACTIC_STAR":
        teff = None
        logg = None
        feh = None
        if gaia:
            gf = gaia.get("features", {})
            try: teff = None if gf.get("Teff") is None else float(gf.get("Teff"))
            except Exception: pass
            try: logg = None if gf.get("logg") is None else float(gf.get("logg"))
            except Exception: pass
            try: feh = None if gf.get("[Fe/H]") is None else float(gf.get("[Fe/H]"))
            except Exception: pass

        subtype = "GALACTIC_STAR"
        plain = "Это обычная галактическая звезда, подтверждённая астрометрией Gaia."
        if teff is not None and logg is not None and 3900 <= teff <= 5200 and logg >= 4.0:
            subtype = "COOL_K_TYPE_DWARF_STAR"
            plain = "Это холодный карлик K-типа: Gaia видит значимое движение/параллакс, а температура и высокая поверхностная гравитация соответствуют карликовой звезде."
        elif _f(gd, "bp_rp") is not None and _f(gd, "bp_rp") >= 1.5:
            subtype = "RED_GALACTIC_STAR"
            plain = "Это красная галактическая звезда: Gaia фиксирует очень значимое собственное движение; точный спектральный подтип пока не заморожен."

        return {
            "nature": subtype,
            "confidence": "HIGH",
            "plain_ru": plain,
            "evidence": {
                "gaia_source_id": None if gaia is None else gaia.get("source_id"),
                "parallax_mas": _f(gd, "parallax_mas"),
                "parallax_snr": _f(gd, "parallax_snr"),
                "proper_motion_vector_snr": _f(gd, "proper_motion_vector_snr"),
                "ruwe": _f(gd, "ruwe"),
                "bp_rp": _f(gd, "bp_rp"),
                "g_mag": _f(gd, "g_mag"),
                "naive_distance_pc_if_high_snr": _f(gd, "naive_inverse_parallax_distance_pc_if_snr_ge_5"),
                "teff_K": teff,
                "logg": logg,
                "feh": feh,
                "sdss_photo_class": sdss_class,
                "sdss_mapping_if_present": "6=STAR" if sdss_class == 6 else None,
            },
            "ceiling": "GALACTIC_STELLAR_SOURCE_CLASSIFICATION"
        }

    if c.get("class") == "FAINT_INFRARED_SOURCE_UNCLASSIFIED":
        wf = {} if wise is None else wise.get("features", {})
        return {
            "nature": "FAINT_INFRARED_POINT_SOURCE_UNCLASSIFIED",
            "confidence": "LOW_MEDIUM",
            "plain_ru": "Это реальный слабый инфракрасный точечный источник в AllWISE, но без Gaia/оптического/2MASS совпадения в нашем жёстком радиусе его нельзя честно назвать звездой или галактикой.",
            "evidence": {
                "allwise_source_id": None if wise is None else wise.get("source_id"),
                "W1mag": wf.get("W1mag"),
                "W2mag": wf.get("W2mag"),
                "photometric_quality": wf.get("qph"),
                "contamination_flags": wf.get("ccf"),
                "extended_flag": wf.get("ex"),
                "gaia_counterpart_within_2arcsec": gaia is not None,
                "panstarrs_counterpart_within_2arcsec": ps1 is not None,
                "twomass_counterpart_within_2arcsec": tm is not None,
            },
            "ceiling": "INFRARED_CATALOG_SOURCE_ONLY"
        }

    return {
        "nature": "UNRESOLVED_SOURCE_CLASS",
        "confidence": "LOW",
        "plain_ru": "Источник реален в каталогах, но текущих данных недостаточно для более узкой физической классификации.",
        "evidence": {"upstream_class": c.get("class")},
        "ceiling": "UNRESOLVED"
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    src = json.loads(args.input.read_text(encoding="utf-8"))
    groups = src.get("groups", {})
    out_groups = {name: {
        "target": g.get("target"),
        "offset_from_frozen_target_center_arcsec": g.get("offset_from_frozen_target_center_arcsec"),
        "representative_icrs": g.get("representative_icrs"),
        "nature_result": _nature(name, g),
    } for name, g in groups.items()}

    payload = {
        "schema": "janus.cosmos.love_edem.source_nature.receipt.v1",
        "experiment_id": "LOVE-EDEM-SOURCE-NATURE-v1",
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETE",
        "input_receipt": str(args.input),
        "input_status": src.get("status"),
        "input_service_failures": src.get("service_failures", []),
        "anomaly_scoring_used": False,
        "groups": out_groups,
        "sdss_photo_type_mapping": {
            "3": "GALAXY",
            "6": "STAR",
            "reference": "SDSS SkyServer PhotoObj/PhotoType convention"
        },
        "summary": {
            "high_confidence_galactic_stars": [n for n, g in out_groups.items() if g["nature_result"]["confidence"] == "HIGH" and "STAR" in g["nature_result"]["nature"]],
            "galaxy_classifications": [n for n, g in out_groups.items() if g["nature_result"]["nature"] == "GALAXY"],
            "infrared_unclassified": [n for n, g in out_groups.items() if g["nature_result"]["nature"] == "FAINT_INFRARED_POINT_SOURCE_UNCLASSIFIED"],
        },
        "firewall": {
            "near_LOVE_source_is_identity_LOVE": False,
            "near_EDEM_source_is_identity_EDEM": False,
            "source_classification_establishes_planet": False,
            "source_classification_establishes_anomaly": False,
            "source_classification_establishes_physical_LOVE_EDEM_ORION_link": False,
            "claim_ceiling": "NATURE_OF_CATALOG_SOURCE_GROUPS_AROUND_FROZEN_COORDINATES_ONLY"
        }
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload["summary"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
