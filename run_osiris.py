#!/usr/bin/env python3
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
direct = root / "experiments" / "direct"
target = direct / "s_phallus_h_gate_2_bounded_k_scaling_holdout_budget_guard.py"
if not target.exists():
    print("OSIRIS source snapshot is not initialized.", file=sys.stderr)
    print("Run: git submodule update --init --recursive", file=sys.stderr)
    raise SystemExit(2)
cmd = [sys.executable, str(target), "--self-test", *sys.argv[1:]]
raise SystemExit(subprocess.call(cmd, cwd=direct))
