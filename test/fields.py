#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Checks the initial velocity that streamFunFoam writes against the stream
# function wave, Fenton (1999), eqs. (3.24)-(3.25), evaluated here from the
# coefficients of calcStreamFunctionWave. In every cell full of water (alpha 1)
# the written velocity must be the wave velocity at the cell centre. Run after
# streamFunFoam and postProcess -func writeCellCentres at time 0, with ascii
# fields. Usage: fields.py <case> <coefficient file>

import re
import sys

import numpy as np

case, outFile = sys.argv[1:3]


def entry(text, key):
    return re.search(r"^\s*" + key + r"\s+([^;]*);", text, re.M).group(1)


def internalField(name):
    """Internal field of an ascii volField at time 0, as an (n, k) array"""
    text = open(f"{case}/0/{name}").read()
    m = re.search(r"internalField\s+nonuniform\s+List<\w+>\s*(\d+)\s*\(", text)
    n = int(m.group(1))
    values = re.findall(r"[-+0-9.eE]+", text[m.end():].split("\n)", 1)[0])
    return np.array(values, dtype=float).reshape(n, -1)


coeffs = open(outFile).read()
d, k, omega, uBar = (float(entry(coeffs, key)) for key in ("d", "k", "omega", "uBar"))
B = np.array(re.search(r"^\s*B\s+\d+\s*\(([^)]*)\)", coeffs, re.M).group(1).split(), dtype=float)
j = np.arange(1, B.size + 1)[:, None]

wave = open(f"{case}/constant/waveProperties").read()
kDir = np.array(entry(wave, "direction").strip("() ").split(), dtype=float)
kDir /= np.linalg.norm(kDir)
gVec = np.array(entry(open(f"{case}/constant/g").read(), "value").strip("() ").split(), dtype=float)
g = np.linalg.norm(gVec)
up = -gVec/g
log = open(f"{case}/log.streamFunFoam").read()
phase = float(re.search(r"Wave phase: (\S+)", log).group(1))
xWaveMaker = float(re.search(r"Wave maker position: (\S+)", log).group(1))
seabed = float(re.search(r"Seabed height: (\S+)", log).group(1))

C = internalField("C")
U = internalField("U")
alpha = internalField("alpha.water")[:, 0]

x, y = C @ kDir - xWaveMaker, C @ up - seabed
theta = k*x + 2*np.pi*phase
a = j*B[:, None]*np.sqrt(g/k)/np.cosh(j*k*d)
u = np.sum(a*np.cosh(j*k*y)*np.cos(j*theta), axis=0) + omega/k - uBar
v = np.sum(a*np.sinh(j*k*y)*np.sin(j*theta), axis=0)

water = alpha == 1
uRef = np.abs(u[water]).max()
err = max(np.abs(U[water] @ kDir - u[water]).max(), np.abs(U[water] @ up - v[water]).max())/uRef
# The fields are written with 6 significant digits
ok = water.sum() > 0 and err < 1e-5
print(f"{case.rstrip('/').split('/')[-1]}: max velocity error in {water.sum()} water cells "
      f"{err:.1e} of the largest velocity {'ok' if ok else 'FAIL'}")
raise SystemExit(not ok)
