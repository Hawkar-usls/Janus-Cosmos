#!/usr/bin/env python3
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
direct = root / "experiments" / "direct"
legacy = direct / "s_phallus_h_gate_2_bounded_k_scaling_holdout_budget_guard.py"
spiral = root / "workspace" / "osiris_spiral_runtime_strict.py"

args = list(sys.argv[1:])
legacy_mode = "--legacy" in args
if legacy_mode:
    args.remove("--legacy")
    target = legacy
    cwd = direct
else:
    target = spiral
    cwd = root

if not target.exists():
    print("OSIRIS source snapshot is not initialized.", file=sys.stderr)
    print("Run: git submodule update --init --recursive", file=sys.stderr)
    raise SystemExit(2)

# No arguments means: exercise the current canonical runtime.
if not args:
    args = ["--self-test"]

cmd = [sys.executable, str(target), *args]
raise SystemExit(subprocess.call(cmd, cwd=cwd))
