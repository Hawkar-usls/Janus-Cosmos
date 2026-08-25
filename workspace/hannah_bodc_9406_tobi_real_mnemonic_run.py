#!/usr/bin/env python3
"""Run Cousteau Synesthetic Memory Core on real Hannah/BODC CD169 TOBI bytes.

The adapter was frozen before sensory scoring. It streams the full target raw
file once, hashes it, summarizes predetermined time windows, and emits mnemonic
passports. Raw sonar bytes are never written to output or uploaded.
"""
from __future__ import annotations

import argparse
import ftplib
import hashlib
import json
import math
import posixpath
import statistics
import struct
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cousteau_synesthetic_memory_core import build_passport, compare_passports  # noqa:E402

HOST = "livftp.noc.ac.uk"
FILE = "/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11285/TOBI.DAT"
REL = "sd11285/TOBI.DAT"
BLOCK = 40960
TARGET = datetime(2005, 2, 28, 1, 7, 25, tzinfo=timezone.utc)
TARGET_TOWFISH_FROZEN = (-3.8654180644718967, -12.142441475)
SCALES = (60, 300, 1800, 7200)
ARRAYS = {
    "port_sidescan": 0x0240,
    "stbd_sidescan": 0x2180,
    "profiler": 0x40C0,
    "port_swath": 0x6000,
    "stbd_swath": 0x7F40,
}


def ftp() -> ftplib.FTP:
    f = ftplib.FTP(timeout=90)
    f.connect(HOST, 21)
    f.login("anonymous", "janus-probe@example.invalid")
    f.voidcmd("TYPE I")
    return f


def dos_dt(d: int, t: int) -> datetime | None:
    try:
        return datetime(
            1980 + ((d >> 9) & 127), (d >> 5) & 15, d & 31,
            (t >> 11) & 31, (t >> 5) & 63, (t & 31) * 2,
            tzinfo=timezone.utc,
        )
    except Exception:
        return None


def coord(deg: int, mins: float) -> float | None:
    if not math.isfinite(mins) or abs(mins) > 60.5:
        return None
    sign = -1.0 if deg < 0 or mins < 0 else 1.0
    return sign * (abs(float(deg)) + abs(float(mins)) / 60.0)


def circular_mean_deg(vals: list[float]) -> float | None:
    if not vals:
        return None
    s = sum(math.sin(math.radians(v % 360.0)) for v in vals)
    c = sum(math.cos(math.radians(v % 360.0)) for v in vals)
    if abs(s) < 1e-15 and abs(c) < 1e-15:
        return None
    return math.degrees(math.atan2(s, c)) % 360.0


def median(xs: list[float]) -> float | None:
    return float(statistics.median(xs)) if xs else None


def mad(xs: list[float]) -> float | None:
    if not xs:
        return None
    m = statistics.median(xs)
    return float(statistics.median(abs(x - m) for x in xs))


def haversine_km(a: tuple[float, float], b: tuple[float, float]) -> float:
    r = 6371.0088
    lat1, lon1 = map(math.radians, a)
    lat2, lon2 = map(math.radians, b)
    dlat, dlon = lat2 - lat1, lon2 - lon1
    h = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return 2 * r * math.asin(min(1.0, math.sqrt(h)))


@dataclass
class ArrayAgg:
    count: int = 0
    total: int = 0
    total_sq: int = 0
    hist: list[int] = field(default_factory=lambda: [0] * 256)

    def add_block(self, values: tuple[int, ...]) -> None:
        self.count += len(values)
        self.total += sum(values)
        self.total_sq += sum(v * v for v in values)
        for v in values:
            self.hist[((v + 32768) >> 8) & 255] += 1

    def merge_stats(self, other: "ArrayAgg") -> None:
        self.count += other.count
        self.total += other.total
        self.total_sq += other.total_sq
        for i, n in enumerate(other.hist):
            self.hist[i] += n

    def summary(self) -> dict[str, float | None]:
        if not self.count:
            return {"mean": None, "std": None, "entropy8": None}
        mean = self.total / self.count
        var = max(0.0, self.total_sq / self.count - mean * mean)
        ent = 0.0
        for n in self.hist:
            if n:
                p = n / self.count
                ent -= p * math.log2(p)
        return {"mean": mean, "std": math.sqrt(var), "entropy8": ent}


def block_array_aggs(block: bytes) -> dict[str, ArrayAgg]:
    out: dict[str, ArrayAgg] = {}
    for name, off in ARRAYS.items():
        vals = struct.unpack_from("<4000h", block, off)
        agg = ArrayAgg()
        agg.add_block(vals)
        out[name] = agg
    return out


def parse_telemetry(block: bytes) -> dict[str, Any] | None:
    try:
        tm, dat, alt = struct.unpack_from("<HHH", block, 0x32)
        dt = dos_dt(dat, tm)
        if dt is None or not (2004 <= dt.year <= 2006):
            return None
        lon_deg, lat_deg, lon_min, lat_min = struct.unpack_from("<hhff", block, 0x38)
        magx = list(struct.unpack_from("<8i", block, 0x44))
        magy = list(struct.unpack_from("<8i", block, 0x64))
        magz = list(struct.unpack_from("<8i", block, 0x84))
        roll = list(struct.unpack_from("<8h", block, 0xA4))
        pitch = list(struct.unpack_from("<8h", block, 0xB4))
        gyro = list(struct.unpack_from("<8h", block, 0xD4))
        press = list(struct.unpack_from("<8H", block, 0xE4))
        temp = list(struct.unpack_from("<8H", block, 0xF4))
        cond = list(struct.unpack_from("<8H", block, 0x104))
        water, wire = struct.unpack_from("<hh", block, 0x114)
        lss = list(struct.unpack_from("<8i", block, 0x118))
        return {
            "dt": dt,
            "altitude": float(alt),
            "latitude": coord(lat_deg, lat_min),
            "longitude": coord(lon_deg, lon_min),
            "gyro_deg": [(x / 10.0 - 10.1) % 360.0 for x in gyro],
            "roll_deg": [x / 6.4 for x in roll],
            "pitch_raw": [float(x) for x in pitch],
            "pressure_dbar": [x / 10.0 - 5.0 for x in press],
            "temperature_C": [x / 2000.0 - 2.0 for x in temp],
            "conductivity_mmho": [x / 1000.0 - 2.0 for x in cond],
            "magx_raw": [float(x) for x in magx],
            "magy_raw": [float(x) for x in magy],
            "magz_raw": [float(x) for x in magz],
            "lss_volts": [x * 5.0 / 524288.0 for x in lss],
            "water_path_ms": float(water),
            "wire_out_m": float(wire),
        }
    except Exception:
        return None


@dataclass
class WindowAcc:
    name: str
    nominal_seconds: int
    direction: str
    scale_label: str
    blocks: int = 0
    raw_sha256: Any = field(default_factory=hashlib.sha256)
    raw_blake2: Any = field(default_factory=lambda: hashlib.blake2b(digest_size=32))
    timestamps: list[float] = field(default_factory=list)
    lat: list[float] = field(default_factory=list)
    lon: list[float] = field(default_factory=list)
    alt: list[float] = field(default_factory=list)
    wire: list[float] = field(default_factory=list)
    water: list[float] = field(default_factory=list)
    gyro: list[float] = field(default_factory=list)
    roll: list[float] = field(default_factory=list)
    pitch: list[float] = field(default_factory=list)
    press: list[float] = field(default_factory=list)
    temp: list[float] = field(default_factory=list)
    cond: list[float] = field(default_factory=list)
    magx: list[float] = field(default_factory=list)
    magy: list[float] = field(default_factory=list)
    magz: list[float] = field(default_factory=list)
    lss: list[float] = field(default_factory=list)
    arrays: dict[str, ArrayAgg] = field(default_factory=lambda: {k: ArrayAgg() for k in ARRAYS})

    def add(self, block: bytes, t: dict[str, Any], aa: dict[str, ArrayAgg]) -> None:
        self.blocks += 1
        self.raw_sha256.update(block)
        self.raw_blake2.update(block)
        self.timestamps.append(t["dt"].timestamp())
        if t["latitude"] is not None: self.lat.append(t["latitude"])
        if t["longitude"] is not None: self.lon.append(t["longitude"])
        self.alt.append(t["altitude"]); self.wire.append(t["wire_out_m"]); self.water.append(t["water_path_ms"])
        self.gyro.extend(t["gyro_deg"]); self.roll.extend(t["roll_deg"]); self.pitch.extend(t["pitch_raw"])
        self.press.extend(t["pressure_dbar"]); self.temp.extend(t["temperature_C"]); self.cond.extend(t["conductivity_mmho"])
        self.magx.extend(t["magx_raw"]); self.magy.extend(t["magy_raw"]); self.magz.extend(t["magz_raw"]); self.lss.extend(t["lss_volts"])
        for k in ARRAYS: self.arrays[k].merge_stats(aa[k])

    def summarize(self) -> tuple[dict[str, Any], dict[str, Any]]:
        diffs = [b - a for a, b in zip(self.timestamps, self.timestamps[1:]) if b > a]
        expected = max(1, round(self.nominal_seconds / 4)) if self.nominal_seconds else 1
        missing = max(0.0, min(1.0, 1.0 - self.blocks / expected)) if self.nominal_seconds else 0.0
        arr = {k: self.arrays[k].summary() for k in ARRAYS}
        pm, sm = arr["port_sidescan"]["mean"], arr["stbd_sidescan"]["mean"]
        psm, ssm = arr["port_swath"]["mean"], arr["stbd_swath"]["mean"]
        asym_ss = None if pm is None or sm is None else (pm-sm) / max(abs(pm)+abs(sm), 1e-12)
        asym_sw = None if psm is None or ssm is None else (psm-ssm) / max(abs(psm)+abs(ssm), 1e-12)
        heading = circular_mean_deg(self.gyro)
        physical = {
            "latitude": median(self.lat), "longitude": median(self.lon), "heading_deg": heading,
            "altitude_m": median(self.alt), "wire_out_m": median(self.wire), "water_path_ms": median(self.water),
            "pressure_dbar": median(self.press), "temperature_C": median(self.temp), "conductivity_mmho": median(self.cond),
            "roll_deg": median(self.roll), "pitch_raw": median(self.pitch),
            "magx_raw": median(self.magx), "magy_raw": median(self.magy), "magz_raw": median(self.magz), "lss_volts": median(self.lss),
            "timestamp_cadence_s": median(diffs), "cadence_jitter_mad_s": mad(diffs), "missing_fraction": missing,
            "array_summaries": arr, "sidescan_asymmetry": asym_ss, "swath_asymmetry": asym_sw,
        }
        payload: dict[str, Any] = {
            "latitude": physical["latitude"], "longitude": physical["longitude"], "heading": heading,
            "timestamp_cadence": physical["timestamp_cadence_s"], "cadence_jitter": physical["cadence_jitter_mad_s"],
            "missing_fraction": missing,
            "acquisition_tobi_altitude_scaled": None if physical["altitude_m"] is None else physical["altitude_m"] / 1000.0,
            "acquisition_tobi_wire_out_scaled": None if physical["wire_out_m"] is None else physical["wire_out_m"] / 10000.0,
            "acquisition_tobi_water_path_scaled": None if physical["water_path_ms"] is None else physical["water_path_ms"] / 10000.0,
            "acquisition_tobi_pressure_scaled": None if physical["pressure_dbar"] is None else physical["pressure_dbar"] / 5000.0,
            "acquisition_tobi_temperature_scaled": None if physical["temperature_C"] is None else physical["temperature_C"] / 10.0,
            "acquisition_tobi_conductivity_scaled": None if physical["conductivity_mmho"] is None else physical["conductivity_mmho"] / 50.0,
            "acquisition_tobi_roll_scaled": None if physical["roll_deg"] is None else physical["roll_deg"] / 30.0,
            "acquisition_tobi_pitch_raw_scaled": None if physical["pitch_raw"] is None else physical["pitch_raw"] / 100.0,
            "acquisition_tobi_magx_raw_scaled": None if physical["magx_raw"] is None else physical["magx_raw"] / 200000.0,
            "acquisition_tobi_magy_raw_scaled": None if physical["magy_raw"] is None else physical["magy_raw"] / 200000.0,
            "acquisition_tobi_magz_raw_scaled": None if physical["magz_raw"] is None else physical["magz_raw"] / 200000.0,
            "acquisition_tobi_lss_scaled": None if physical["lss_volts"] is None else physical["lss_volts"] / 10.0,
            "acquisition_tobi_sidescan_asymmetry": asym_ss,
            "acquisition_tobi_swath_asymmetry": asym_sw,
        }
        for k, s in arr.items():
            payload[f"acquisition_tobi_{k}_mean_scaled"] = None if s["mean"] is None else s["mean"] / 32768.0
            payload[f"acquisition_tobi_{k}_std_scaled"] = None if s["std"] is None else s["std"] / 32768.0
            payload[f"acquisition_tobi_{k}_entropy8_scaled"] = None if s["entropy8"] is None else s["entropy8"] / 8.0
        payload = {k:v for k,v in payload.items() if v is not None}
        meta = {
            "window": self.name, "direction": self.direction, "scale": self.scale_label,
            "block_count": self.blocks, "nominal_seconds": self.nominal_seconds,
            "start_utc": datetime.fromtimestamp(min(self.timestamps), timezone.utc).isoformat() if self.timestamps else None,
            "end_utc": datetime.fromtimestamp(max(self.timestamps), timezone.utc).isoformat() if self.timestamps else None,
            "raw_window_sha256": self.raw_sha256.hexdigest(), "raw_window_blake2b_256": self.raw_blake2.hexdigest(),
        }
        return payload, {"meta": meta, "physical_summary": physical}


def passport_for(acc: WindowAcc) -> dict[str, Any]:
    payload, info = acc.summarize()
    # The core needs raw bytes for a true source hash. We intentionally do not retain
    # the window bytes; use the already-computed raw hash as provenance and canonical
    # payload source identity. Exact source hash remains in info.meta.
    pp = build_passport(
        payload, direction=acc.direction, scale=acc.scale_label,
        provenance={"source": "BODCREQ-9406/CD169_TOBI/" + REL, **info["meta"]},
        expected_feature_count=33,
    )
    pp["source_identity"]["raw_window_sha256_authoritative"] = info["meta"]["raw_window_sha256"]
    pp["source_identity"]["raw_window_blake2b_256_authoritative"] = info["meta"]["raw_window_blake2b_256"]
    pp["tobi_physical_summary"] = info["physical_summary"]
    pp["tobi_semantic_warning"] = "TOBI_SONAR_TEXTURE_AND_PRESSURE_CUES_ARE_MNEMONIC_ONLY__NOT_SEAFLOOR_GEOMETRY_OR_DEPTH_VERDICTS"
    return pp


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--output", required=True, type=Path); args = ap.parse_args()
    result: dict[str, Any] = {
        "artifact_id": "JANUS-HANNAH-BODC-MEASUREMENT-RECEIPT-001-REAL-CD169-TOBI-SENSORY-PASSPORTS",
        "schema": "janus.cosmos.cousteau.hannah_bodc.real_tobi_mnemonic_run.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "STARTED", "scientific_convergence_claim": False,
        "adapter_contract": "data/cousteau/JANUS-HANNAH-CD169-TOBI-SYNESTHETIC-ADAPTER-CONTRACT-2026-08-25-v1.0.json",
        "adapter_contract_commit": "f33acaa3336d07b84caa06de223b6d8426f81907",
        "locator_run": 32893244352,
        "target_utc": TARGET.isoformat(),
        "target_towfish_coordinate_frozen": {"latitude": TARGET_TOWFISH_FROZEN[0], "longitude": TARGET_TOWFISH_FROZEN[1], "role": "CONTEXT_ONLY_NOT_FINGERPRINT"},
        "source": {"host": HOST, "relative_path": REL, "raw_bytes_redistributed": False},
    }
    windows: dict[str, WindowAcc] = {}
    for s in SCALES:
        label = f"{s}s"
        windows[f"CENTER_{label}"] = WindowAcc(f"CENTER_{label}", s, "CENTER", label)
        windows[f"PRE_{label}"] = WindowAcc(f"PRE_{label}", s, "HEAD_FORWARD", label)
        windows[f"POST_{label}"] = WindowAcc(f"POST_{label}", s, "TAIL_REVERSE", label)
    file_sha = hashlib.sha256(); file_b2 = hashlib.blake2b(digest_size=32); file_bytes=0; file_blocks=0
    nearest: tuple[float, bytes, dict[str, Any], dict[str, ArrayAgg], int] | None = None
    chronology_deltas: list[float] = []
    prev_dt: datetime | None = None
    try:
        f = ftp()
        try:
            size = f.size(FILE)
            result["source"]["size_bytes"] = size
            sock = f.transfercmd("RETR " + FILE)
            buf = bytearray(); idx = 0
            try:
                while True:
                    chunk = sock.recv(1024 * 1024)
                    if not chunk: break
                    file_sha.update(chunk); file_b2.update(chunk); file_bytes += len(chunk); buf.extend(chunk)
                    while len(buf) >= BLOCK:
                        block = bytes(buf[:BLOCK]); del buf[:BLOCK]; file_blocks += 1
                        t = parse_telemetry(block); idx_now = idx; idx += 1
                        if t is None: continue
                        if prev_dt is not None: chronology_deltas.append((t["dt"] - prev_dt).total_seconds())
                        prev_dt = t["dt"]
                        delta = (t["dt"] - TARGET).total_seconds()
                        if abs(delta) > 7200: continue
                        aa = block_array_aggs(block)
                        if nearest is None or abs(delta) < abs(nearest[0]): nearest = (delta, block, t, aa, idx_now)
                        for s in SCALES:
                            if abs(delta) <= s / 2.0: windows[f"CENTER_{s}s"].add(block, t, aa)
                            if -s <= delta < 0: windows[f"PRE_{s}s"].add(block, t, aa)
                            if 0 < delta <= s: windows[f"POST_{s}s"].add(block, t, aa)
            finally:
                try: sock.close()
                except Exception: pass
        finally:
            try: f.close()
            except Exception: pass

        result["raw_file_integrity"] = {
            "bytes_streamed": file_bytes, "block_count": file_blocks,
            "sha256": file_sha.hexdigest(), "blake2b_256": file_b2.hexdigest(),
            "size_matches_server": (result["source"].get("size_bytes") == file_bytes),
            "size_mod_block": file_bytes % BLOCK,
        }
        result["whole_file_cadence"] = {"median_seconds": median(chronology_deltas), "mad_seconds": mad(chronology_deltas)}
        if nearest is None: raise RuntimeError("no valid block within +/-7200s of frozen target")
        delta, raw_exact, t_exact, aa_exact, idx_exact = nearest
        exact = WindowAcc("EXACT_NEAREST_BLOCK", 0, "CENTER", "custom"); exact.add(raw_exact, t_exact, aa_exact)
        exact_pp = passport_for(exact)
        result["exact_nearest_block"] = {
            "record_index": idx_exact, "delta_seconds_from_frozen_target": delta,
            "datetime_utc": t_exact["dt"].isoformat(), "raw_block_sha256": hashlib.sha256(raw_exact).hexdigest(),
            "passport": exact_pp,
        }
        ship = (exact_pp["tobi_physical_summary"]["latitude"], exact_pp["tobi_physical_summary"]["longitude"])
        result["ship_vs_frozen_towfish_context"] = {
            "ship_coordinate_from_raw_header": {"latitude": ship[0], "longitude": ship[1]},
            "frozen_towfish_coordinate": {"latitude": TARGET_TOWFISH_FROZEN[0], "longitude": TARGET_TOWFISH_FROZEN[1]},
            "great_circle_separation_km": haversine_km(ship, TARGET_TOWFISH_FROZEN),
            "rule": "SEPARATION_IS_EXPECTED_TOWED_INSTRUMENT_CONTEXT__NOT_A_NAVIGATION_ERROR_BY_ITSELF",
        }
        result["windows"] = {name: passport_for(acc) for name, acc in windows.items()}
        result["local_mirror_comparisons"] = {}
        for s in SCALES:
            result["local_mirror_comparisons"][f"{s}s"] = compare_passports(result["windows"][f"PRE_{s}s"], result["windows"][f"POST_{s}s"])
        result["exact_vs_center"] = {f"{s}s": compare_passports(exact_pp, result["windows"][f"CENTER_{s}s"]) for s in SCALES}

        # Frozen anti-leakage and direction-isolation tests on the same real measurement.
        exact_payload, _ = exact.summarize()
        base_a = build_passport(exact_payload, direction="CENTER", scale="custom")
        base_b = build_passport(exact_payload, direction="CENTER", scale="custom")
        polluted = dict(exact_payload); polluted.update({"H1": 1, "pyramid": 999, "verdict": 1, "prediction": 42})
        forbidden = build_passport(polluted, direction="CENTER", scale="custom")
        dh = build_passport(exact_payload, direction="HEAD_FORWARD", scale="custom")
        dtail = build_passport(exact_payload, direction="TAIL_REVERSE", scale="custom")
        result["real_data_regression_gates"] = {
            "same_input_same_passport": base_a["passport_sha256"] == base_b["passport_sha256"],
            "forbidden_label_invariance": base_a["measurement_fingerprint"]["sha256"] == forbidden["measurement_fingerprint"]["sha256"],
            "direction_does_not_change_measurement_fingerprint": dh["measurement_fingerprint"]["sha256"] == dtail["measurement_fingerprint"]["sha256"],
            "direction_does_change_context_passport": dh["passport_sha256"] != dtail["passport_sha256"],
            "raw_target_block_hash_present": bool(result["exact_nearest_block"]["raw_block_sha256"]),
        }
        result["all_real_data_regression_gates_pass"] = all(result["real_data_regression_gates"].values())
        result["status"] = "REAL_TOBI_SENSORY_PASSPORTS_READY"
    except Exception as exc:
        result["status"] = "REAL_TOBI_MNEMONIC_RUN_FAILED"; result["error_type"] = type(exc).__name__; result["error"] = str(exc)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "raw_file_sha256": (result.get("raw_file_integrity") or {}).get("sha256"),
        "exact": {k:v for k,v in (result.get("exact_nearest_block") or {}).items() if k != "passport"},
        "mirror_similarity": {k:v.get("common_measurement_similarity") for k,v in (result.get("local_mirror_comparisons") or {}).items()},
        "real_data_gates_pass": result.get("all_real_data_regression_gates_pass"),
        "scientific_convergence_claim": False,
    }, indent=2))
    return 0 if result["status"] == "REAL_TOBI_SENSORY_PASSPORTS_READY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
