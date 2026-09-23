#!/usr/bin/env python3
import importlib.util, random
from pathlib import Path

P=Path(__file__).with_name("kusto_synthetic_multibeam_nuisance_v2.py")
spec=importlib.util.spec_from_file_location("sampler",P)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

assert m.bin_for_across(0)=="INNER"
assert m.bin_for_across(2999.9)=="INNER"
assert m.bin_for_across(3000)=="MID"
assert m.bin_for_across(4500)=="OUTER"
assert m.bin_for_across(5500)=="EXTREME"
assert m.bin_for_across(-6100)=="EXTREME"

c={"support_cell_count":4,"signed_residual_m":[-12.0,7.0,125.0,-160.0]}
rng=random.Random(1)
for _ in range(100):
    x,src=m.draw_one(rng,c,"CD169_REPLAY")
    assert x in c["signed_residual_m"]
    assert src=="SIGNED_EMPIRICAL_BOOTSTRAP"

rng=random.Random(2)
allowed={12.0,7.0,125.0,160.0}
for _ in range(100):
    x,src=m.draw_one(rng,c,"GENERIC_UNSEEN_CALIBRATED")
    assert abs(x) in allowed
    assert src=="EMPIRICAL_MAGNITUDE_BOOTSTRAP_RANDOMIZED_SIGN"

rng=random.Random(3)
for _ in range(100):
    x,src=m.draw_one(rng,c,"ADVERSARIAL_REAL_TAILS")
    assert abs(x) in {125.0,160.0}
    assert src=="EMPIRICAL_GE100_TAIL_STRESS_RANDOMIZED_SIGN"

try:
    m.draw_one(random.Random(4),{"support_cell_count":2,"signed_residual_m":[-12.0,20.0]},"ADVERSARIAL_REAL_TAILS")
    raise AssertionError("expected missing-tail ValueError")
except ValueError:
    pass

print("PASS: empirical multibeam nuisance sampler self-test")
