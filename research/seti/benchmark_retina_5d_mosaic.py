"""Fully synthetic 5D mosaic diagnostic; no astronomical FPR or deployment claim."""
import json
import numpy as np
from retina_5d_mosaic import analyze_spectrogram, _normalize_rows, SLOPES


def matched_mean(x):
    z = _normalize_rows(x)
    n, f = z.shape
    best = float("-inf")
    for slope in SLOPES:
        drift = np.rint(slope * (np.arange(n) - (n - 1) / 2.)).astype(int)
        for start in range(2-int(drift.min()), f-2-int(drift.max())):
            p = start + drift
            t = np.arange(n)
            val = np.mean(z[t,p] - (z[t,p-2] + z[t,p+2])*.5)
            best = max(best, val)
    return float(best)


def sample(rng, kind):
    on = rng.normal(10, 1, (32,64))
    off = rng.normal(10, 1, (32,64))
    if kind in ("weak", "strong", "shared"):
        amp = {"weak":1.0, "strong":3.0, "shared":3.0}[kind]
        on[:,30] += amp
        if kind == "shared":
            off[:,30] += amp
    elif kind == "impulse":
        for _ in range(2):
            on[rng.integers(0,32),rng.integers(3,61)] += 18
    return on, off


def scores(x, y):
    v = analyze_spectrogram(x, off_power=y, source_id="on",
                            reference_id="off", native_grid_verified=True)
    return (matched_mean(x), v["evidence"]["native_research_score"],
            v["evidence"]["shared_pattern_advisory_only"])


def main():
    rng = np.random.default_rng(2026112301)
    cal = np.array([scores(*sample(rng,"noise"))[:2] for _ in range(80)])
    threshold = np.quantile(cal, 0.95, axis=0)
    rng = np.random.default_rng(2026112302)
    rows={}
    for kind in ("noise","impulse","weak","strong","shared"):
        current=np.array([scores(*sample(rng,kind)) for _ in range(40)])
        rows[kind]={
            "count":40,
            "direct_exceed":int(sum(current[:,0]>threshold[0])),
            "mosaic_exceed":int(sum(current[:,1]>threshold[1])),
            "shared_advisory":int(sum(current[:,2]>0))}
    result={
        "schema":"janus.cosmos.retina_5d_mosaic.synthetic_diagnostic.v1",
        "seed_calibration":2026112301,"seed_evaluation":2026112302,
        "calibration_noise_only":80,"evaluation_cases":200,
        "threshold_basis":"synthetic in-sample noise p95 for each method; NOT SCIENTIFIC FPR",
        "direct_95th":float(threshold[0]),"mosaic_95th":float(threshold[1]),
        "metrics":rows,"source_kind":"SYNTHETIC_GAUSSIAN_ONLY",
        "production_adoption":False,
        "scientific_status":"HOLD_FRESH_REAL_INDEPENDENT_REPLAY"}
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
