/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Johan Roenby
-------------------------------------------------------------------------------
License
    This file is an extension of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    Test-streamFunctionWaveRequests

Description
    Random requests to streamFunctionWave over the range of the table:
    relative depths d/lambda from 0.03 to 3, heights up to the top of the
    table (half of them in its upper 10 %), tolerances 1e-6 to 1e-10 and all
    input types. For period input (uType 1, 2) the period is that of the
    wave of the drawn length with a drawn current of up to 20 % of the wave
    speed, so the wave length is known.

    Each wave is compared with a request for the same wave with a tolerance
    1000 times smaller (at least 1e-11): the surface elevation relative to
    H, the velocities at the surface relative to the mean fluid speed in the
    frame of the wave, and the wave length. A request passes if the
    difference is below its tolerance, or below the round-off level of the
    highest waves (3e-9) where the reference warned that it reached it.

Usage
    Test-streamFunctionWaveRequests [-n 200] [-seed 1]

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Random.H"
#include "clockTime.H"
#include "mathematicalConstants.H"
#include "streamFunctionWave.H"
#include "streamFunctionWaveTable.H"

using namespace Foam;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// The wave of dictionary d, calculated quietly. Returns false on failure.
bool calc(const dictionary& d, autoPtr<streamFunctionWave>& w, scalar& t)
{
    w.reset(new streamFunctionWave(d));
    const bool throwing = FatalError.throwing(true);
    clockTime timer;
    bool ok = true;
    try
    {
        w->calcCoeffs();
    }
    catch (const Foam::error& err)
    {
        Info<< "    FAILED: " << err.message().c_str() << nl;
        ok = false;
    }
    t = timer.elapsedTime();
    FatalError.throwing(throwing);
    return ok;
}


// Largest difference of wave a to the reference b
scalar difference(const streamFunctionWave& a, const streamFunctionWave& b)
{
    const scalar ka = a.k(), kb = b.k();
    const scalar H = b.H();
    const scalar uBar = b.uBar();

    const label n = 400;
    scalarField x(n);
    forAll(x, i)
    {
        x[i] = 2*pi/kb*i/n;
    }
    const scalarField eta(b.getEta(x, 0, 0));

    scalar diff = mag(ka - kb)/kb;
    diff = max(diff, max(mag(a.getEta(x, 0, 0) - eta))/H);
    diff = max(diff, max(mag(a.getU(x, eta, 0, 0) - b.getU(x, eta, 0, 0)))/uBar);
    diff = max(diff, max(mag(a.getV(x, eta, 0, 0) - b.getV(x, eta, 0, 0)))/uBar);
    return diff;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("n", "n", "Number of requests (200)");
    argList::addOption("seed", "seed", "Random seed (1)");
    argList args(argc, argv);

    // Values put into dictionaries go through streams of this precision
    IOstream::defaultPrecision(17);

    const label n = args.getOrDefault<label>("n", 200);
    Random rnd(args.getOrDefault<label>("seed", 1));
    const streamFunctionWaveTable& table = streamFunctionWaveTable::builtIn();
    const scalar d = 1, g = 9.81;
    const scalar roundOff = 3e-9;

    label nFail = 0, nOk = 0;
    scalarList times;
    scalarList timesTop;

    for (label i = 0; i < n; ++i)
    {
        const scalar dl = 0.03*std::pow(100.0, rnd.sample01<scalar>());
        const scalar sigma =
            rnd.sample01<scalar>() < 0.5
          ? rnd.sample01<scalar>()
          : 0.9 + 0.1*rnd.sample01<scalar>();
        const scalar tolerance = std::pow(10.0, -6 - 2*label(3*rnd.sample01<scalar>()));
        const label uType = label(4*rnd.sample01<scalar>());
        const scalar lambda = d/dl;
        const scalar H = sigma*table.HtopByLambda(dl)*lambda;

        dictionary in;
        in.add("H", H);
        in.add("d", d);
        in.add("g", g);
        in.add("tolerance", tolerance);
        in.add("uType", uType);

        string what;
        if (uType == 0 || uType == 3)
        {
            in.add("lambda", lambda);
            if (uType == 3)
            {
                in.add("T", lambda/(std::sqrt(g*std::tanh(2*pi*dl)/(2*pi/lambda))));
            }
        }
        else
        {
            // The period of the wave of length lambda with the drawn current
            dictionary l(in);
            l.set("uType", 0);
            l.set("lambda", lambda);
            l.set("tolerance", 1e-12);
            autoPtr<streamFunctionWave> w;
            scalar t;
            if (!calc(l, w, t))
            {
                ++nFail;
                continue;
            }
            const scalar uBar = w->uBar();
            const scalar u = 0.2*(2*rnd.sample01<scalar>() - 1)*uBar;
            const scalar c = uType == 1 ? u + uBar : u + w->Q()/d;
            in.add("T", lambda/c);
            in.add("u", u);
        }

        Info<< "d/lambda " << dl << ", H/H_top " << sigma << ", uType "
            << uType << ", tolerance " << tolerance << nl;

        autoPtr<streamFunctionWave> w, ref;
        scalar t, tRef;
        if (!calc(in, w, t))
        {
            ++nFail;
            continue;
        }
        dictionary inRef(in);
        inRef.set("tolerance", max(1e-3*tolerance, 1e-11));
        if (!calc(inRef, ref, tRef))
        {
            ++nFail;
            continue;
        }

        const scalar refErr = ref->error();
        const scalar lambdaW = 2*pi/w->k();
        const scalar diff = difference(*w, *ref);
        const bool lengthOk =
            uType == 0 || uType == 3
         || mag(lambdaW - lambda) < max(tolerance, roundOff)*lambda;
        const bool ok =
            lengthOk
         && (diff <= tolerance || (refErr > 0.01*tolerance && diff <= roundOff));

        Info<< "    tolerance " << tolerance << ", lambda " << lambdaW
            << " (rel. error " << mag(lambdaW - lambda)/lambda << "), N "
            << w->N() << ", estimated error " << w->error()
            << ", difference to reference " << diff
            << " (reference N " << ref->N() << ", error "
            << refErr << "), " << t << " s: " << (ok ? "ok" : "BAD") << nl;

        if (ok) ++nOk; else ++nFail;
        (sigma > 0.98 ? timesTop : times).append(t);
    }

    Foam::sort(times);
    Foam::sort(timesTop);
    Info<< nl << nOk << " of " << n << " requests ok, " << nFail << " failed"
        << nl;
    if (times.size())
    {
        Info<< "time H/H_top <= 0.98: median " << times[times.size()/2]
            << " s, max " << times.last() << " s" << nl;
    }
    if (timesTop.size())
    {
        Info<< "time H/H_top > 0.98: median " << timesTop[timesTop.size()/2]
            << " s, max " << timesTop.last() << " s" << nl;
    }

    return nFail ? 1 : 0;
}


// ************************************************************************* //
