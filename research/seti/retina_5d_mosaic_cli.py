"""Opt-in read-only native radio or TOPA frame-manifest entrypoint."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from retina_5d_mosaic import analyze_spectrogram, inspect_retina_video_tiles


def from_topa_frame_manifest(filename: Path) -> dict:
    """Verify TOPA Retina sampled frame hashes; do not decode the original video."""
    manifest = json.loads(filename.read_text(encoding="utf-8"))
    if manifest.get("schema") != "janus.topa.retina.frame-manifest.v1":
        raise ValueError("not an authenticated-format TOPA frame manifest")
    if len(str(manifest.get("source_sha256",""))) != 64:
        raise ValueError("TOPA source video digest missing")
    records = manifest.get("frames")
    if not isinstance(records, list) or not 2 <= len(records) <= 256:
        raise ValueError("require 2..256 sampled frames")
    try:
        from PIL import Image
    except ImportError as exc:
        raise RuntimeError("Pillow required only for TOPA JPEG decoding") from exc
    frames=[]
    for entry in records:
        name=entry.get("file")
        if not isinstance(name,str) or Path(name).name != name or name.startswith("."):
            raise ValueError("unsafe or missing TOPA filename")
        file=filename.parent / "frames" / name
        if hashlib.sha256(file.read_bytes()).hexdigest() != entry.get("sha256"):
            raise ValueError("TOPA sampled frame integrity mismatch")
        if entry.get("pts_time_s") is None:
            raise ValueError("TOPA sampled frame PTS missing")
        with Image.open(file) as image:
            frames.append(np.asarray(image.convert("L"), dtype=np.float64))
    if any(x.shape != frames[0].shape for x in frames):
        raise ValueError("TOPA samples are not same-sized; no silent resampling")
    result=inspect_retina_video_tiles(np.stack(frames))
    result["topa_video_sha256"]=manifest["source_sha256"]
    result["topa_frame_sha256"]=[entry["sha256"] for entry in records]
    result["topa_frame_pts_time_s"]=[entry["pts_time_s"] for entry in records]
    result["original_video_rescanned"]=False
    return result


def main() -> None:
    ap=argparse.ArgumentParser(description="Shadow-only 5D mosaic on native COSMOS power or TOPA Retina QC")
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--native-on",type=Path,help="read-only 2D .npy native time/frequency power")
    mode.add_argument("--retina-manifest",type=Path,help="TOPA frame_manifest.json with local verified frames/")
    ap.add_argument("--native-off",type=Path)
    ap.add_argument("--source-id")
    ap.add_argument("--reference-id")
    ap.add_argument("--native-grid-verified",action="store_true",help="caller attestation; not proof of instrument metadata")
    ap.add_argument("--baseline-alert",action="store_true",help="original detector alert; never altered")
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()
    if args.retina_manifest:
        if args.native_off:
            ap.error("radio OFF cannot be used on video mode")
        result=from_topa_frame_manifest(args.retina_manifest)
    else:
        on=np.load(args.native_on,allow_pickle=False)
        off=np.load(args.native_off,allow_pickle=False) if args.native_off else None
        result=analyze_spectrogram(on,off_power=off,source_id=args.source_id,
                                   reference_id=args.reference_id,
                                   native_grid_verified=args.native_grid_verified,
                                   baseline_alert=args.baseline_alert)
    payload=json.dumps(result,indent=2,ensure_ascii=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(payload,encoding="utf-8")
    else:
        print(payload,end="")


if __name__=="__main__":
    main()
