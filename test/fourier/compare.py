#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Compares the waves of calcStreamFunctionWave with those of the independent
# implementation fenton.py, for each input type (uType), shallow to deep water
# and low to steep waves. Run with the OpenFOAM environment sourced and
# streamFunctionWave built. Usage: compare.py <workDir>

import os
import re
import subprocess
import sys

import fenton

work = sys.argv[1]
os.makedirs(work, exist_ok=True)

# name: entries of the input dictionary (uType 0 and 3 take lambda)
waves = {
    "example": dict(H=5, d=50, T=10, uType=2, u=0),
    "shallowEulerian": dict(H=0.5, d=1, T=6, uType=1, u=0.2),
    "lab": dict(H=0.1, d=0.6, T=1.2, uType=2, u=0.1),
    "lengthOnly": dict(H=0.2, d=1, **{"lambda": 2}, uType=0),
    "lengthAndPeriod": dict(H=0.2, d=1, **{"lambda": 2}, T=1.2, uType=3),
    "steepDeep": dict(H=2, d=100, **{"lambda": 20}, uType=0),
}
N, nHSteps = 35, 20
quantities = ["k", "c", "uBar", "u1", "u2", "Q", "R", "h_crest", "h_trough"]
# Both solve the same discrete equations; they agree to round-off, except for
# steep waves in deep water, where the difference Jacobians of both limit the
# agreement to about 1e-8
tolerance = 1e-7

fails = 0
print(f"{'wave':16s} {'max relative difference':>24s}  quantity")
for name, w in waves.items():
    inFile, outFile = os.path.join(work, name), os.path.join(work, name + ".out")
    with open(inFile, "w") as f:
        f.write("FoamFile { format ascii; class dictionary; object waveProperties; }\n")
        for key, value in dict(w, N=N, nHSteps=nHSteps).items():
            f.write(f"{key} {value};\n")
    if os.path.exists(outFile):
        os.remove(outFile)
    log = subprocess.run(["calcStreamFunctionWave", inFile, outFile],
                         capture_output=True, text=True)
    if log.returncode != 0:
        print(f"{name:16s} calcStreamFunctionWave failed:\n{log.stdout[-2000:]}{log.stderr[-2000:]}")
        fails += 1
        continue
    text = open(outFile).read()
    foam = {q: float(re.search(rf"^{q}\s+([^;]+);", text, re.M).group(1)) for q in quantities}

    args = dict(H=w["H"], d=w["d"], N=N, nSteps=nHSteps)
    if w["uType"] in (0, 3):
        args["lam"] = w["lambda"]
        if w["uType"] == 3:
            args["T"] = w["T"]
    else:
        args["T"] = w["T"]
        args["u1" if w["uType"] == 1 else "u2"] = w["u"]
    ref = fenton.solve(**args)

    # Relative to the size of each quantity's natural scale
    scale = {"k": ref["k"], "Q": abs(ref["Q"]) + w["d"]*(9.81*w["d"])**0.5,
             "R": ref["R"], "h_crest": w["H"], "h_trough": w["H"]}
    diff = {q: abs(foam[q] - ref[q])/scale.get(q, (9.81*w["d"])**0.5) for q in quantities}
    worst = max(diff, key=diff.get)
    ok = diff[worst] < tolerance
    fails += not ok
    print(f"{name:16s} {diff[worst]:24.1e}  {worst} {'ok' if ok else 'FAIL'}")

print("ALL WAVES AGREE" if not fails else f"{fails} WAVES DIFFER")
raise SystemExit(fails)
