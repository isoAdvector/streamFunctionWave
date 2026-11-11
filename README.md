<!--
  Copyright (C) 2012-2017 DHI A/S
  Copyright (C) 2017-2026 Johan Roenby
  SPDX-License-Identifier: GPL-3.0-or-later
-->

# streamFunctionWave

Developed and maintained by Johan Roenby, [STROMNING](https://stromning.com).

Steady nonlinear water waves of any height up to near breaking, on any depth,
for OpenFOAM.com: a wave generating boundary condition, a utility that sets the
fields of a case to the exact wave at every time step, and a utility that
calculates the wave. The waves are stream function waves, calculated with the
Fourier approximation method of section 3.1.2 of

> Fenton, J. D. (1999). Numerical methods for nonlinear waves. *Advances in
> Coastal and Ocean Engineering*, 5, 241–324.

whose equation numbers the code refers to. The method assumes two-dimensional,
inviscid and irrotational flow over a horizontal bed, in a frame in which the
wave is steady, but no small wave height, unlike linear (Airy) or Stokes
waves.

streamFunctionWave is the successor of
[fentonWave](https://github.com/roenby/fentonWave), developed at DHI; see
[History](#history).

## Contents

| Directory | What |
|---|---|
| [src/streamFunctionWave](src/streamFunctionWave) | The class `streamFunctionWave`: calculates the wave and gives the surface elevation, velocity and pressure at any point and time |
| [src/streamFunctionWaveFvPatchField](src/streamFunctionWaveFvPatchField) | The boundary condition `streamFunctionWave` for the volume fraction and the velocity, in the library `libstreamFunctionWaveFvPatchField` |
| [applications/streamFunFoam](applications/streamFunFoam) | Sets `alpha`, `U`, `phi`, `p` and `p_rgh` of an interFoam-type case to the wave at every time step: initial conditions, or a reference solution sampled as a simulation would be |
| [applications/calcStreamFunctionWave](applications/calcStreamFunctionWave) | Calculates a wave from an input dictionary and writes its parameters and Fourier coefficients |
| [run](run) | Example cases |
| [test](test) | Tests |

## Installation

With OpenFOAM.com sourced (v2306 to v2606):

```sh
git clone https://github.com/isoAdvector/streamFunctionWave.git
cd streamFunctionWave
./Allwmake
```

This builds `libstreamFunctionWaveFvPatchField` into `$FOAM_USER_LIBBIN` and
`streamFunFoam` and `calcStreamFunctionWave` into `$FOAM_USER_APPBIN`. Build it
again with each OpenFOAM version used.

## Usage

### The wave

The wave is described in a dictionary, by default `constant/waveProperties`
(see [applications/calcStreamFunctionWave/waveProperties](applications/calcStreamFunctionWave/waveProperties)):

| Entry | Meaning | Default |
|---|---|---|
| `H` | wave height (crest to trough) | required |
| `d` | mean water depth | required |
| `uType` | what gives the length: `0` wave length `lambda`, the wave seen in its own frame; `1` period `T` and Eulerian mean current `u`; `2` period `T` and mass transport velocity `u`; `3` wave length `lambda` and period `T` | `2` |
| `T` | wave period (`uType` 1, 2, 3) | |
| `lambda` | wave length (`uType` 0, 3) | |
| `u` | current (`uType` 1, 2) | `0` |
| `N` | number of Fourier modes | `20` |
| `nHSteps` | number of steps in which the wave height is raised to `H` | `20` |
| `g` | gravity | `9.81` |
| `direction` | direction of propagation (boundary condition, `streamFunFoam`) | required |
| `seabedHeight` | vertical position of the bed | lowest point of the patch or mesh |
| `waveMakerPos` | origin of the horizontal coordinate | first point of the patch or mesh along `direction` |
| `phase` | phase, in wave lengths | the surface at the mean level at the origin at time 0 |

`streamFunFoam` also reads `uAir` (velocity of the air, default zero) and
`comoving` and `c` (to give the water velocity in a frame moving with
celerity `c`). Steep waves need more modes and height steps: the tests use
`N 35` and `nHSteps 20` up to H/λ = 0.1, about 70 % of the breaking height.

### The boundary condition

In `0/U` and `0/alpha.water` (the patch also works for `alpha1`):

```
inlet
{
    type            streamFunctionWave;
    value           uniform (0 0 0);    // uniform 0 for alpha.water
    waveDictName    waveProperties;     // default: waveProperties
    waveDictPath    "constant";         // default: "constant"
}
```

and in `system/controlDict`:

```
libs            (streamFunctionWaveFvPatchField);
```

A `fixedFluxPressure` condition suits `p_rgh` on the patch. The wave is
calculated when the patch is first constructed and stored in its entry, with
its direction, phase and origin, so `decomposePar`, restarts and
`reconstructPar` read it without recalculating it. Restore `0/` from `0.orig/`
to recalculate the wave after changing the dictionary.

### streamFunFoam and calcStreamFunctionWave

`streamFunFoam`, run in a case like a solver, writes the wave fields at every
write time; its time step follows `maxCo`, `maxAlphaCo` and `maxDeltaT` as in
the VoF solvers. With `endTime 0` it writes the initial fields of a
simulation. `calcStreamFunctionWave <inputDict> <outputDict>` writes the wave
length, celerity, currents, crest and trough heights and the Fourier
coefficients of a wave, and recalculates only when the input describes another
wave than the output.

### Example cases

| Case | What |
|---|---|
| [run/waveFlume](run/waveFlume) | A flume with the boundary condition at the inlet, run with `interIsoFoam` |
| [run/streamFunWaveFlume](run/streamFunWaveFlume) | The same flume with the fields set by `streamFunFoam` |
| [run/periodicWaveTank](run/periodicWaveTank) | A periodic tank of one wave length, initialised by `streamFunFoam` and run with `interIsoFoam` |

Each case has an `Allrun` and an `Allclean`. The wave is that of H = 0.24 m
and T = 1.8 s on 1 m of water, for 30 s.

## Tests

With OpenFOAM sourced and streamFunctionWave built, `test/Allrun` (about half
a minute) checks:

- **[test/fourier/fenton.py](test/fourier/fenton.py)**, an independent
  implementation of Fenton's method in Python (numpy, scipy), against the
  linear dispersion relation for small waves and against the celerity of
  Stokes waves in deep water.
- **`calcStreamFunctionWave` against fenton.py**
  ([test/fourier/compare.py](test/fourier/compare.py)) for each `uType`, from
  shallow water (kd = 0.28) to deep water (kd = 31), up to H/λ = 0.1. The wave
  numbers, celerities, currents, volume fluxes, Bernoulli constants and crest
  and trough heights agree to 1e-13 relative for moderate waves, to 7e-10 at
  H/λ = 0.1 and kd = π, and to 1e-8 at H/λ = 0.1 in deep water, where the
  difference Jacobians of both limit the agreement.
- **The example cases** in [run](run) run for 0.2 s.

`test/Allclean` removes the work directory.

## History

Johan Roenby wrote **fentonWave**
([github.com/roenby/fentonWave](https://github.com/roenby/fentonWave)), a stream
function wave boundary condition for interFoam, at **DHI Water & Environment**
in 2012 and released it as free software. Its successor, streamFunctionWave,
was developed at DHI from 2014 to 2017: the wave calculation moved into a class
of its own without the dependency on GSL, with the currents and wave length
input of Fenton's method, the coefficients of a patch were calculated once and
shared between its fields, and `calcStreamFunctionWave` and `streamFunFoam`
were added. We gratefully acknowledge DHI, where the code was first developed.

Since 2017 Johan Roenby has developed it further: ports to the current OpenFOAM.com
releases, the wave stored in the boundary condition for restarts and parallel
runs, the input of both wave length and period (`uType` 3), `streamFunFoam`
setting the fields at every time step, and in this release a more robust
calculation (steep waves, and deep water, where it overflowed for kd > 20) and
the tests. streamFunctionWave generates the waves of the test cases of
[surfacePosition](https://github.com/isoAdvector/surfacePosition) and of
fluxinator, a VoF solver to be released shortly.

## License

Copyright (C) 2012-2017 DHI A/S, Copyright (C) 2017-2026 Johan Roenby.
streamFunctionWave is an extension of OpenFOAM and, like OpenFOAM, free
software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version. It is
distributed in the hope that it will be useful, but without any warranty; see
[LICENSE](LICENSE) for details.

## Author

Johan Roenby, [STROMNING](https://stromning.com). STROMNING offers consulting,
custom development, training and support for interfacial and multiphase flow
simulation with OpenFOAM. Bug reports and questions are welcome as GitHub issues;
for anything else, contact johan@stromning.com.
