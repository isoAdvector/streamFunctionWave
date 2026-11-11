#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# An independent implementation of the Fourier approximation method of
# section 3.1.2 of Fenton (1999), Numerical methods for nonlinear waves,
# Advances in Coastal and Ocean Engineering 5, 241-324, written from the paper
# for testing calcStreamFunctionWave. It solves the same discrete equations
# (N Fourier modes, N + 1 points over half a wave length) with scipy, in units
# of the depth d and gravity g, stepping up the wave height from linear theory.
#
# Frame moving with the wave, X from a crest, Y up from the bed:
#   psi(X, Y) = -uBar Y + sum_j B_j sinh(j k Y)/cosh(j k d) cos(j k X)
# Surface eta_m = eta(X_m), X_m = m pi/(N k), m = 0..N:
#   kinematic:  psi(X_m, eta_m) = -Q
#   dynamic:    (u^2 + v^2)/2 + g eta_m = R
#   mean depth: trapezoidal mean of eta_m = d
#   height:     eta_0 - eta_N = H
# and two of: k lambda = 2 pi; k c T = 2 pi; c - uBar = u1 (Eulerian current);
# c - Q/d = u2 (mass transport velocity); c = 0 (wave seen in its own frame).

import numpy as np
from scipy.optimize import brentq


def ratio(a, b, f):
    """sinh(a)/cosh(b) (f = -1) or cosh(a)/cosh(b) (f = 1) for 0 <= a <~ b,
    without overflow"""
    return (np.exp(a - b) + f*np.exp(-a - b))/(1 + np.exp(-2*b))


def residual(z, N, H, spec):
    k, c, uBar, Q, R = z[:5]
    eta, B = z[5:N + 6], z[N + 6:]
    j = np.arange(1, N + 1)[None, :]
    theta = (np.pi*np.arange(N + 1)/N)[:, None]
    a, b = j*k*eta[:, None], j*k
    sh, ch = ratio(a, b, -1), ratio(a, b, 1)
    psi = -uBar*eta + (B*sh*np.cos(j*theta)).sum(1)
    u = -uBar + (j*k*B*ch*np.cos(j*theta)).sum(1)
    v = (j*k*B*sh*np.sin(j*theta)).sum(1)
    mean = (eta[0] + eta[-1] + 2*eta[1:-1].sum())/(2*N)
    eqs = [mean - 1, eta[0] - eta[-1] - H]
    for name, value in spec:
        eqs.append({
            "lambda": lambda: k*value - 2*np.pi,
            "T": lambda: k*c*value - 2*np.pi,
            "u1": lambda: c - uBar - value,
            "u2": lambda: c - Q - value,
            "c": lambda: c - value,
        }[name]())
    return np.concatenate([eqs, psi + Q, 0.5*(u**2 + v**2) + eta - R])


def newton(f, z, args, tol=1e-13, maxIter=50):
    """Newton iteration with a forward difference Jacobian. The step is
    absolute for small unknowns: the high Fourier coefficients are near
    round-off, where a relative step (as in scipy's fsolve) gives no
    derivative."""
    z = z.copy()
    last = np.inf
    for _ in range(maxIter):
        r = f(z, *args)
        err = np.abs(r).max()
        # Converged, or at the round-off floor of the difference Jacobian
        if err < tol or (err < 1e-9 and err > 0.5*last):
            return z
        last = err
        h = 1e-7*np.maximum(np.abs(z), 1e-3)
        J = np.empty((len(r), len(z)))
        for i in range(len(z)):
            zi = z.copy()
            zi[i] += h[i]
            J[:, i] = (f(zi, *args) - r)/h[i]
        z -= np.linalg.solve(J, r)
    raise RuntimeError(f"no convergence, residual {np.abs(f(z, *args)).max():.1e}")


def linear(N, H, spec):
    """Linear wave as the first guess (dimensionless)"""
    s = dict(spec)
    if "lambda" in s:
        k = 2*np.pi/s["lambda"]
        c = s.get("c", 2*np.pi/(k*s["T"]) if "T" in s else 0)
        uBar = np.sqrt(np.tanh(k)/k)
    else:
        omega = 2*np.pi/s["T"]
        u = s.get("u1", s.get("u2", 0))
        # Doppler shifted dispersion relation, the root with omega - k u > 0
        dispersion = lambda k: omega - k*u - np.sqrt(k*np.tanh(k))
        kMax = 1.0
        while dispersion(kMax) > 0:
            kMax *= 2
        k = brentq(dispersion, 1e-12, kMax)
        c = omega/k
        uBar = c - u
    theta = np.pi*np.arange(N + 1)/N
    B = np.zeros(N)
    B[0] = uBar*0.5*H/np.tanh(k)
    return np.concatenate([[k, c, uBar, uBar, 0.5*uBar**2 + 1],
                           1 + 0.5*H*np.cos(theta), B])


def solve(H, d, N=35, nSteps=20, g=9.81, T=None, lam=None, u1=None, u2=None):
    """Wave of height H on depth d; give T or lambda (or both), and with T
    alone u1 or u2. Returns the dimensional quantities"""
    tScale, lScale, uScale = np.sqrt(d/g), d, np.sqrt(g*d)
    spec = []
    if lam is not None:
        spec.append(("lambda", lam/lScale))
        spec.append(("T", T/tScale) if T is not None else ("c", 0.0))
    else:
        spec.append(("T", T/tScale))
        spec.append(("u1", u1/uScale) if u1 is not None else ("u2", (u2 or 0)/uScale))
    Hd = H/d
    z = linear(N, Hd/nSteps, spec)
    previous = None
    for i in range(1, nSteps + 1):
        # Linear extrapolation from the last two solved heights
        guess = z if previous is None else 2*z - previous
        sol = newton(residual, guess, (N, Hd*i/nSteps, spec))
        previous, z = (z if i > 1 else None), sol
    k, c, uBar, Q, R = z[:5]
    eta = z[5:N + 6]
    return {
        "k": k/lScale, "lambda": 2*np.pi*lScale/k, "c": c*uScale,
        "uBar": uBar*uScale, "u1": (c - uBar)*uScale, "u2": (c - Q)*uScale,
        "Q": Q*uScale*lScale, "R": R*g*d,
        "h_crest": (eta[0] - 1)*lScale, "h_trough": (eta[-1] - 1)*lScale,
        "maxResidual": np.abs(residual(z, N, Hd, spec)).max(),
    }


if __name__ == "__main__":
    # Checks of this implementation against classical results
    fails = 0

    # Linear limit: dispersion relation omega^2 = g k tanh(k d)
    for kd in (0.3, 1.0, 3.0):
        d, g = 10.0, 9.81
        T = 2*np.pi/np.sqrt(g*kd/d*np.tanh(kd))
        w = solve(1e-6*d, d, N=10, nSteps=1, T=T, u2=0)
        err = abs(w["k"]*d/kd - 1)
        ok = err < 1e-9
        fails += not ok
        print(f"linear limit kd = {kd}: relative error of k {err:.1e} {'ok' if ok else 'FAIL'}")

    # Deep water, Stokes: c = c0 (1 + (k a)^2/2), error O((k a)^4)
    d, g, lam = 1000.0, 9.81, 100.0
    k0 = 2*np.pi/lam
    for ka in (0.02, 0.05, 0.1):
        w = solve(2*ka/k0, d, N=20, nSteps=5, lam=lam)
        c0 = np.sqrt(g/k0)
        stokes = c0*(1 + 0.5*ka**2)
        err = abs(w["uBar"]/stokes - 1)
        ok = err < 2*ka**4
        fails += not ok
        print(f"Stokes deep water ka = {ka}: relative difference of c {err:.1e} "
              f"(bound 2 (ka)^4 = {2*ka**4:.1e}) {'ok' if ok else 'FAIL'}")
    raise SystemExit(fails)
