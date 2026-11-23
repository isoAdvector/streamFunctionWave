#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# The figures of the README, from waves calculated by calcStreamFunctionWave
# at three relative depths d/lambda (shallow, intermediate, deep):
#
#   waveScan.png   the surface of waves from 0.1 to 1 times the highest wave
#                  of the table, H_top, at each depth
#   accuracy.png   how well the waves satisfy the free-surface conditions,
#                  from 0.05 H_top to H_top, for tolerances 1e-6 to 1e-10,
#                  and the number of Fourier modes used
#
# The velocity potential of a stream function wave is a sum of terms that
# each satisfy the Laplace equation and the condition of no flow through the
# bed exactly, so the accuracy of a wave is how well it satisfies the two
# conditions at the surface, in the frame moving with the wave:
#
#   kinematic  psi(x, eta) = -Q            (the surface is a streamline)
#   dynamic    (U^2 + V^2)/2 + g eta = R   (constant pressure, Bernoulli)
#
# They are evaluated here, independently of the C++ code, from the Fourier
# coefficients written by calcStreamFunctionWave, at 4001 points over half a
# wave length, mostly between the points where the equations were solved.
# The residuals are scaled to compare with the tolerance: the kinematic one
# by uBar H, giving about the displacement of the surface relative to H, and
# the dynamic one by uBar^2, giving about the error in the speed at the
# surface relative to the mean speed uBar.
#
# A wave fails if a residual exceeds its tolerance, or 3e-9 for the waves
# within 2 % of H_top, which are determined only to about 1e-9.
#
# Run with the OpenFOAM environment sourced and streamFunctionWave built.
# Usage: figures.py <workDir> <figureDir>

import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D

work, figDir = sys.argv[1], sys.argv[2]
os.makedirs(work, exist_ok=True)
os.makedirs(figDir, exist_ok=True)

# d/lambda of the three depths, from the bottom of the scan figure up
depths = [("shallow", 0.05), ("intermediate", 0.15), ("deep", 1.0)]
scanSigmas = [0.1*i for i in range(1, 11)]
sigmas = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 1.0]
tolerances = [1e-6, 1e-8, 1e-10]
lam = 1.0

# Colours: the blue ordinal ramp for the heights (light to dark), and the
# first three categorical slots for the tolerances
heightColours = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
                 "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
tolColours = ["#2a78d6", "#eb6834", "#1baf7a"]
ink, inkMuted, grid = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": inkMuted, "axes.labelcolor": ink,
    "xtick.color": inkMuted, "ytick.color": inkMuted, "text.color": ink,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
})


def parse(text):
    """Scalars and lists of a calcStreamFunctionWave output dictionary"""
    out = {}
    for m in re.finditer(r"^(\w+)\s+([-+0-9.eE]+)\s*;", text, re.M):
        out[m.group(1)] = float(m.group(2))
    for m in re.finditer(r"^(\w+)\s+(\d+)\s*\(([^)]*)\)", text, re.M):
        out[m.group(1)] = np.array(m.group(3).split(), dtype=float)
    return out


def calc(name, **entries):
    """The wave of the given input entries, calculated by
    calcStreamFunctionWave, or None with the log printed if it failed"""
    inFile, outFile = os.path.join(work, name), os.path.join(work, name + ".out")
    with open(inFile, "w") as f:
        f.write("FoamFile { format ascii; class dictionary; object waveProperties; }\n")
        for key, value in entries.items():
            f.write(f"{key} {value!r};\n")
    if os.path.exists(outFile):
        os.remove(outFile)
    log = subprocess.run(["calcStreamFunctionWave", inFile, outFile],
                         capture_output=True, text=True)
    if log.returncode != 0:
        print(f"{name}: calcStreamFunctionWave failed:\n{log.stdout[-1500:]}{log.stderr[-1500:]}")
        return None
    return parse(open(outFile).read())


def ratio(a, b, sign):
    """sinh(a)/cosh(b) (sign -1) or cosh(a)/cosh(b) (sign 1) for 0 <= a,
    without overflow for large b"""
    return (np.exp(a - b) + sign*np.exp(-a - b))/(1 + np.exp(-2*b))


def surface(w, theta):
    """Surface elevation above the bed at phase theta"""
    m = np.arange(len(w["E"]))
    return 2*(w["E"][:, None]*np.cos(m[:, None]*theta)).sum(axis=0)


def residuals(w, n=4001):
    """Scaled kinematic and dynamic residuals, and the relative error of the
    height, of wave w"""
    k, d, g, H = w["k"], w["d"], w["g"], w["H"]
    uBar, Q, R, B = w["uBar"], w["Q"], w["R"], w["B"]
    theta = np.linspace(0, np.pi, n)
    eta = surface(w, theta)
    j = np.arange(1, len(B) + 1)[:, None]
    sh = ratio(j*k*eta, j*k*d, -1)
    ch = ratio(j*k*eta, j*k*d, 1)
    Bj = B[:, None]
    psi = -uBar*eta + np.sqrt(g/k**3)*(Bj*sh*np.cos(j*theta)).sum(axis=0)
    U = -uBar + np.sqrt(g/k)*(j*Bj*ch*np.cos(j*theta)).sum(axis=0)
    V = np.sqrt(g/k)*(j*Bj*sh*np.sin(j*theta)).sum(axis=0)
    kinematic = np.abs(psi + Q).max()/(uBar*H)
    dynamic = np.abs(0.5*(U**2 + V**2) + g*eta - R).max()/uBar**2
    height = abs(eta[0] - eta[-1] - H)/H
    return kinematic, dynamic, height


# The highest wave of the table at each depth, from the H/H_top of a low wave
Htop = {}
for label, dl in depths:
    w = calc(f"probe_{label}", H=0.002*lam, d=dl*lam, **{"lambda": lam}, uType=0)
    if w is None:
        raise SystemExit(1)
    Htop[label] = 0.002*lam/w["HbyHtop"]

# All waves, calculated in parallel
jobs = {}
for label, dl in depths:
    for s in scanSigmas:
        jobs[("scan", label, round(s, 3), 1e-8)] = (label, dl, s, 1e-8)
    for s in sigmas:
        for tol in tolerances:
            jobs[("acc", label, s, tol)] = (label, dl, s, tol)


def run(key):
    label, dl, s, tol = jobs[key]
    name = f"{key[0]}_{label}_{s:.3f}_{tol:.0e}"
    return key, calc(name, H=s*Htop[label], d=dl*lam, **{"lambda": lam},
                     uType=0, tolerance=tol)


with ThreadPoolExecutor(max_workers=min(16, os.cpu_count() or 1)) as pool:
    waves = dict(pool.map(run, jobs))

fails = sum(w is None for w in waves.values())

# ----------------------------------------------------------------------------
# Scan of heights at the three depths, stacked with equal scales in x and y

fig, ax = plt.subplots(figsize=(7.2, 2.7))
x = np.linspace(-0.5, 0.5, 801)
theta = 2*np.pi*x/lam
offset, top = 0.0, None
gap = 0.025*lam
for label, dl in depths:
    curves = [surface(waves[("scan", label, round(s, 3), 1e-8)], theta) - dl*lam
              for s in scanSigmas if waves[("scan", label, round(s, 3), 1e-8)] is not None]
    if top is not None:
        offset = top + gap - min(c.min() for c in curves)/lam
    ax.axhline(offset, color=grid, lw=1, zorder=0)
    for c, colour in zip(curves, heightColours):
        ax.plot(x, c/lam + offset, color=colour, lw=1.2)
    top = offset + max(c.max() for c in curves)/lam
    ax.text(-0.52, offset, f"{label}\nd/λ = {dl:g}\nH_top/λ = {Htop[label]/lam:.4f}",
            ha="right", va="center", color=inkMuted, fontsize=8)
ax.set_aspect("equal")
ax.set_xlim(-0.5, 0.5)
ax.set_xlabel("x/λ")
ax.set_yticks([])
ax.spines["left"].set_visible(False)
sm = plt.cm.ScalarMappable(
    cmap=ListedColormap(heightColours),
    norm=BoundaryNorm(np.linspace(0.05, 1.05, 11), 10))
cb = fig.colorbar(sm, ax=ax, ticks=scanSigmas, fraction=0.025, pad=0.02,
                  aspect=18, shrink=0.95)
cb.ax.set_yticklabels([f"{s:.1f}" for s in scanSigmas])
cb.set_label("H/H_top")
cb.outline.set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(figDir, "waveScan.png"), dpi=160)
plt.close(fig)

# ----------------------------------------------------------------------------
# Accuracy: residuals of the surface conditions, and the number of modes

fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.6), sharex=True,
                        gridspec_kw=dict(height_ratios=[2, 1]))
print(f"{'depth':13s} {'H/H_top':>7s} {'tolerance':>9s} {'N':>4s} {'kinematic':>10s} "
      f"{'dynamic':>10s} {'height':>10s}")
for col, (label, dl) in enumerate(depths):
    ax, axN = axs[0, col], axs[1, col]
    for tol, colour in zip(tolerances, tolColours):
        ss, kin, dyn, N = [], [], [], []
        for s in sigmas:
            w = waves[("acc", label, s, tol)]
            if w is None:
                continue
            r = residuals(w)
            bound = tol if s <= 0.98 else max(tol, 3e-9)
            ok = max(r) <= bound
            fails += not ok
            print(f"{label:13s} {s:7.3f} {tol:9.0e} {int(w['N']):4d} {r[0]:10.2e} "
                  f"{r[1]:10.2e} {r[2]:10.2e} {'ok' if ok else 'FAIL'}")
            ss.append(s)
            kin.append(r[0])
            dyn.append(r[1])
            N.append(w["N"])
        ax.axhline(tol, color=colour, lw=0.8, ls=":")
        ax.semilogy(ss, kin, color=colour, lw=1.5, marker="o", ms=3.5)
        ax.semilogy(ss, dyn, color=colour, lw=1.5, ls="--", marker="s", ms=3.5,
                    mfc="white")
        axN.plot(ss, N, color=colour, lw=1.5, marker="o", ms=3.5)
    ax.set_title(f"{label}, d/λ = {dl:g}", fontsize=9)
    ax.set_ylim(1e-15, 1e-5)
    ax.set_yticks([1e-14, 1e-12, 1e-10, 1e-8, 1e-6])
    ax.grid(axis="y", color=grid, lw=0.6)
    axN.set_yscale("log")
    axN.set_ylim(4, 400)
    axN.set_yticks([5, 10, 20, 50, 100, 200])
    axN.set_yticklabels(["5", "10", "20", "50", "100", "200"])
    axN.grid(axis="y", color=grid, lw=0.6)
    axN.set_xlabel("H/H_top")
    axN.set_xlim(0, 1.02)
    if col == 0:
        ax.set_ylabel("residual of the surface conditions")
        axN.set_ylabel("Fourier modes N")
    else:
        ax.set_yticklabels([])
        axN.set_yticklabels([])
handles = [Line2D([], [], color=c, lw=1.5, label=f"tolerance {t:.0e}")
           for t, c in zip(tolerances, tolColours)]
handles += [Line2D([], [], color=inkMuted, lw=1.5, marker="o", ms=3.5,
                   label="kinematic, max |ψ(x, η) + Q|/(ū H)"),
            Line2D([], [], color=inkMuted, lw=1.5, ls="--", marker="s", ms=3.5,
                   mfc="white", label="dynamic, max |(U² + V²)/2 + gη − R|/ū²"),
            Line2D([], [], color=inkMuted, lw=0.8, ls=":", label="the tolerance")]
fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
           fontsize=7.5, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.89))
fig.savefig(os.path.join(figDir, "accuracy.png"), dpi=160)
plt.close(fig)

print("ALL WAVES ACCURATE" if not fails else f"{fails} WAVES FAILED")
raise SystemExit(fails)
