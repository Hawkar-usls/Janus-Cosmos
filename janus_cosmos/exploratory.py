from __future__ import annotations

import math
import shutil
from collections import Counter

import numpy as np


def _entropy(text: str) -> float:
    if not text:
        return 0.0
    counts = Counter(text)
    n = len(text)
    return float(-sum((c / n) * math.log2(c / n) for c in counts.values()))


def run_exploratory(image: np.ndarray) -> dict:
    """Post-gate exploratory enrichment that cannot change the blind decision."""
    out = {
        "affects_blind_gate": False,
        "ocr": {"status": "CAPABILITY_UNAVAILABLE"},
        "face_like_detection": {"status": "CAPABILITY_UNAVAILABLE", "identity_search": False},
        "semantic_analysis": {"status": "CAPABILITY_UNAVAILABLE", "model_download": False},
        "cipher_scan": {"status": "NO_OCR_TEXT"},
        "post_hoc_tuning": {"status": "FORBIDDEN_ON_EVALUATION_DATA"},
    }

    try:
        import pytesseract
        from PIL import Image

        if shutil.which("tesseract"):
            x = np.asarray(image, dtype=np.float32)
            lo, hi = np.percentile(x, [1, 99])
            vis = np.clip((x - lo) / max(float(hi - lo), 1e-9), 0, 1)
            u8 = (vis * 255).astype(np.uint8)
            text = pytesseract.image_to_string(Image.fromarray(u8), config="--psm 11").strip()
            out["ocr"] = {"status": "OK", "text": text, "character_count": len(text)}
            cleaned = "".join(ch for ch in text if ch.isprintable()).strip()
            if cleaned:
                tokens = [t for t in cleaned.split() if t]
                repetitions = {k: v for k, v in Counter(tokens).items() if v > 1}
                out["cipher_scan"] = {
                    "status": "DESCRIPTIVE_ONLY",
                    "character_entropy_bits": _entropy(cleaned),
                    "repeated_tokens": repetitions,
                    "claim_ceiling": "Pattern statistics only; not evidence of a message or cipher.",
                }
    except Exception as exc:
        out["ocr"] = {"status": "CAPABILITY_UNAVAILABLE", "error": f"{type(exc).__name__}: {exc}"}

    try:
        import cv2
        x = np.asarray(image, dtype=np.float32)
        lo, hi = np.percentile(x, [1, 99])
        u8 = (np.clip((x - lo) / max(float(hi - lo), 1e-9), 0, 1) * 255).astype(np.uint8)
        cascade_path = cv2.data.haarcascades + "haarcascade_frontalface_default.xml"
        cascade = cv2.CascadeClassifier(cascade_path)
        boxes = cascade.detectMultiScale(u8, scaleFactor=1.1, minNeighbors=5, minSize=(24, 24))
        out["face_like_detection"] = {
            "status": "OK",
            "identity_search": False,
            "count": int(len(boxes)),
            "boxes_xywh": [[int(v) for v in box] for box in boxes],
            "claim_ceiling": "Pareidolia-prone detections only; no identity or astrophysical interpretation.",
        }
    except Exception as exc:
        out["face_like_detection"] = {
            "status": "CAPABILITY_UNAVAILABLE",
            "identity_search": False,
            "error": f"{type(exc).__name__}: {exc}",
        }

    return out
