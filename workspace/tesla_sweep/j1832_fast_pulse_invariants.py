#!/usr/bin/env python3
"""Deterministic invariants for the published FAST 0.2 s J1832-0911 pulse.

This script does not load or reconstruct the FAST waveform. It only derives
geometry/timescale invariants from published scalar measurements, preserving
an explicit claim ceiling for Tesla Sweep / TOPA.
"""

from __future__ import annotations

import json
import math

C_KM_S = 299_792.458
PERIOD_S = 2656.23
PULSE_WIDTH_S = 0.2
SAMPLE_S = 196e-6


def compute() -> dict:
    fraction = PULSE_WIDTH_S / PERIOD_S
    light_cylinder_km = C_KM_S * PERIOD_S / (2.0 * math.pi)
    causal_km = C_KM_S * PULSE_WIDTH_S
    return {
        "period_s": PERIOD_S,
        "pulse_width_s": PULSE_WIDTH_S,
        "sample_time_s": SAMPLE_S,
        "pulse_fraction_of_rotation": fraction,
        "pulse_percent_of_rotation": 100.0 * fraction,
        "rotational_phase_width_deg_if_geometric": 360.0 * fraction,
        "samples_across_pulse": PULSE_WIDTH_S / SAMPLE_S,
        "causal_light_crossing_upper_scale_km": causal_km,
        "light_cylinder_radius_km": light_cylinder_km,
        "causal_scale_fraction_of_light_cylinder": causal_km / light_cylinder_km,
        "boundary": (
            "Derived from published scalars only. Pulse width need not equal physical "
            "source size because beaming, propagation, scattering, relativistic motion "
            "and geometric sweep can modify observed duration."
        ),
        "raw_fast_waveform_processed": False,
    }


if __name__ == "__main__":
    print(json.dumps(compute(), indent=2, sort_keys=True))
