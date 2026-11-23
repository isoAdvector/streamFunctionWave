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
    generateStreamFunctionWaveTable

Group
    grpPreProcessingUtilities

Description
    Generates the table of stream function waves that streamFunctionWave
    starts a requested wave from (src/streamFunctionWave/table).

    For each relative depth d/lambda (rows logarithmically spaced between
    -dlMin and -dlMax), the wave is continued from a small wave in the crest
    parameter omega = 1 - U_crest^2/U_trough^2 towards the highest wave. A
    step is accepted when the solution agrees to -tol with the solution of
    about 25 % more modes; the number of modes grows as needed up to -Nmax.
    The highest wave accepted is H_top(d/lambda), the top of the table.

    Then, at the relative heights sigma = H/H_top = 1 - (1 - i/nSigma)^2,
    i = 1..nSigma, the wave is solved with an increasing number of modes;
    the node records the number of modes from which on the solutions agree
    to each tolerance level, and the solution with the fewest modes that
    meets the tightest level reached. (With more modes than needed, the
    high Fourier coefficients are not unique.)

Usage
    \b generateStreamFunctionWaveTable [OPTIONS] \<outputFile\>

    Options:
      - \par -nRows \<n\>   Number of relative depths (21)
      - \par -dlMin \<d/lambda\>  Shallowest relative depth (0.03)
      - \par -dlMax \<d/lambda\>  Deepest relative depth (3)
      - \par -nSigma \<n\>  Number of relative heights per row (16)
      - \par -tol \<tol\>   Tolerance of the top of the table (1e-10)
      - \par -Nmax \<N\>    Largest number of modes (197)
      - \par -rows \<list\> Only these rows (for testing)
      - \par -verbose      Report every solve of the nodes

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OFstream.H"
#include "clockTime.H"
#include "mathematicalConstants.H"
#include "streamFunctionWaveSolver.H"
#include "streamFunctionWaveTable.H"

using namespace Foam;
using namespace Foam::constant::mathematical;

typedef streamFunctionWaveSolver solver;

static bool verbose = false;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// A wave accepted by the continuation
struct contWave
{
    scalar omega;
    scalar kH;
    label N;
    solver::modal md;
};


// Continuation in omega from a small wave to the top for k d = kd. Returns
// the accepted waves, the last of which is the top.
List<contWave> continuation
(
    const scalar kd,
    const scalar tol,
    const label Nmax
)
{
    solver::parameters ph;
    ph.kd = kd;

    // A small wave of given height, which selects the fundamental (a start
    // in omega can land on the wave of length lambda/n), then the switch to
    // its crest parameter
    label N = 8;
    ph.h = 1e-3*sqr(min(scalar(1), kd));
    solver s0(N, ph);
    scalarField z(s0.linearWave(kd, ph.h));
    s0.solve(z, 1e-14);

    solver::parameters pc(ph);
    pc.heightCondition = false;
    auto crest = [&](const scalar omega)
    {
        solver::parameters p(pc);
        p.omega = omega;
        return p;
    };

    scalar omega = s0.crestParameter(z);
    solver s(N, crest(omega));
    s.solve(z, 1e-14);

    List<contWave> waves;
    scalar sv = -std::log(1 - omega);
    scalar ds = 0.2;
    scalar sp = 0;
    scalarField zp;
    bool havePrev = false;

    while (true)
    {
        const scalar sn = sv + ds;
        const scalar wn = 1 - std::exp(-sn);

        // Guess extrapolated in s = -log(1 - omega)
        autoPtr<solver> s1(new solver(N, crest(wn)));
        scalarField z1(z);
        if (havePrev)
        {
            z1 += (z - zp)*(sn - sv)/(sv - sp);
        }
        solver::result r1 = s1->solve(z1, tol);

        // Check against grow(N) modes; grow N while that is within Nmax
        // (N itself never above Nmax, so that a step that fails at Nmax is
        // retried with a smaller step at Nmax)
        bool ok = false;
        while (r1.residual < GREAT)
        {
            const label N2 = solver::grow(N);
            autoPtr<solver> s2(new solver(N2, crest(wn)));
            scalarField z2(s2->fromModal(s1->toModal(z1)));
            const solver::result r2 = s2->solve(z2, tol);

            const scalar err =
                (r1.converged && r2.converged)
              ? max(s1->difference(z1, *s2, z2), max(r1.step, r2.step))
              : GREAT;

            if (err <= tol && r1.residual < tol)
            {
                ok = true;
                break;
            }
            if (N2 > Nmax)
            {
                break;
            }

            // Continue with more modes, from the last accepted wave
            solver sOld(N2, crest(1 - std::exp(-sv)));
            z = sOld.fromModal(s.toModal(z));
            s = sOld;
            havePrev = false;
            s1 = std::move(s2);
            z1 = z2;
            r1 = r2;
            N = N2;
        }

        if (!ok || !s1->regular(z1))
        {
            ds *= 0.5;
            if (ds < 1e-4)
            {
                break;
            }
            continue;
        }

        const scalar kH = s1->kH(z1);
        if (!waves.empty() && kH <= waves.last().kH)
        {
            // Past the highest wave
            break;
        }

        zp = z;
        sp = sv;
        havePrev = true;
        z = z1;
        sv = sn;
        s = *s1;
        waves.append(contWave{wn, kH, N, s.toModal(z)});

        ds = min(1.3*ds, scalar(0.25));
    }

    return waves;
}


// Node at relative height sigma of the row with continuation waves
streamFunctionWaveTable::nodeData makeNode
(
    const scalar kd,
    const scalar sigma,
    const List<contWave>& waves,
    const scalarList& levels,
    const label Nmax
)
{
    const scalar kHtop = waves.last().kH;
    const scalar kH = sigma*kHtop;
    solver::parameters p;
    p.kd = kd;
    p.h = kH/kd;

    // Start from the continuation wave nearest in height
    label nearest = 0;
    forAll(waves, i)
    {
        if (mag(waves[i].kH - kH) < mag(waves[nearest].kH - kH))
        {
            nearest = i;
        }
    }
    solver::modal best = waves[nearest].md;

    // Solutions with an increasing number of modes, and the difference of
    // each to the next
    struct sol
    {
        label N;
        autoPtr<solver> s;
        scalarField z;
        bool ok;
        scalar step;
    };
    DynamicList<sol*> sols;
    auto get = [&](const label n) -> sol&
    {
        for (sol* x : sols)
        {
            if (x->N == n) return *x;
        }
        sol* x = new sol{n, autoPtr<solver>(new solver(n, p)), scalarField(), false, 0};
        x->z = x->s->fromModal(best);
        x->s->fitB(x->z);
        const solver::result r = x->s->solve(x->z, 1e-14);
        x->ok = r.converged && x->s->regular(x->z);
        x->step = r.step;
        if (verbose)
        {
            Info<< "        N " << n << ": " << r.iterations
                << " iterations, residual " << r.residual << ", converged "
                << r.converged << ", regular " << x->s->regular(x->z) << nl;
        }
        if (x->ok)
        {
            best = x->s->toModal(x->z);
        }
        sols.append(x);
        return *x;
    };

    const label Ncheck = solver::grow(Nmax);
    DynamicList<Tuple2<label, scalar>> chain;
    label Nc = 8;
    while (true)
    {
        sol& a = get(Nc);
        sol& b = get(solver::grow(Nc));
        // The difference, or the uncertainty of a solution whose iteration
        // stalled at round-off if larger
        const scalar e =
            (a.ok && b.ok)
          ? max(a.s->difference(a.z, *b.s, b.z), max(a.step, b.step))
          : GREAT;
        chain.append(Tuple2<label, scalar>(Nc, e));
        if (verbose)
        {
            Info<< "      sigma " << sigma << ": N " << Nc << " against "
                << solver::grow(Nc) << ": " << e << nl;
        }
        if (e < 1e-13 || solver::grow(solver::grow(Nc)) > Ncheck)
        {
            break;
        }
        Nc = solver::grow(Nc);
    }

    streamFunctionWaveTable::nodeData nd;
    nd.sigma = sigma;
    nd.kH = kH;
    nd.N.resize(levels.size());
    forAll(levels, l)
    {
        // The smallest N from which on every N in the chain meets the level
        label n = -1;
        for (label i = chain.size() - 1; i >= 0; --i)
        {
            if (chain[i].second() <= levels[l]) n = chain[i].first();
            else break;
        }
        nd.N[l] = n;
    }

    // The solution with the fewest modes that meets the tightest level
    // reached
    label Ns = Nc;
    for (label l = levels.size() - 1; l >= 0; --l)
    {
        if (nd.N[l] > 0)
        {
            Ns = nd.N[l];
            break;
        }
    }
    sol& x = get(Ns);
    nd.md = x.s->toModal(x.z);

    // Drop the trailing surface coefficients that do not matter
    label nE = nd.md.E.size();
    while (nE > 1 && mag(nd.md.E[nE - 1]) < 1e-13*kH)
    {
        --nE;
    }
    nd.md.E.resize(nE);
    nd.md.B.clear();

    for (sol* y : sols)
    {
        delete y;
    }

    return nd;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate the table of stream function waves that streamFunctionWave"
        " starts a requested wave from"
    );
    argList::noBanner();
    argList::noParallel();
    argList::addArgument("outputFile");
    argList::addOption("nRows", "n", "Number of relative depths (21)");
    argList::addOption("dlMin", "d/lambda", "Shallowest depth (0.03)");
    argList::addOption("dlMax", "d/lambda", "Deepest depth (3)");
    argList::addOption("nSigma", "n", "Relative heights per row (16)");
    argList::addOption("tol", "tol", "Tolerance of the top (1e-10)");
    argList::addOption("Nmax", "N", "Largest number of modes (197)");
    argList::addOption("rows", "list", "Only these rows (for testing)");
    argList::addBoolOption("verbose", "Report every solve of the nodes");

    argList args(argc, argv);

    const label nRows = args.getOrDefault<label>("nRows", 21);
    const scalar dlMin = args.getOrDefault<scalar>("dlMin", 0.03);
    const scalar dlMax = args.getOrDefault<scalar>("dlMax", 3);
    const label nSigma = args.getOrDefault<label>("nSigma", 16);
    const scalar tol = args.getOrDefault<scalar>("tol", 1e-10);
    const label Nmax = args.getOrDefault<label>("Nmax", 197);
    verbose = args.found("verbose");
    labelList rows(identity(nRows));
    args.readIfPresent("rows", rows);

    const scalarList levels({1e-4, 1e-6, 1e-8, 1e-10, 1e-12});

    clockTime total;
    List<streamFunctionWaveTable::rowData> table(rows.size());
    forAll(rows, ri)
    {
        clockTime t;
        const label i = rows[ri];
        const scalar dl =
            nRows > 1
          ? dlMin*std::pow(dlMax/dlMin, scalar(i)/(nRows - 1))
          : dlMin;
        const scalar kd = 2*pi*dl;

        const List<contWave> waves(continuation(kd, tol, Nmax));
        const scalar tCont = t.timeIncrement();

        streamFunctionWaveTable::rowData& row = table[ri];
        row.dl = dl;
        row.kd = kd;
        row.kHtop = waves.last().kH;
        row.nodes.resize(nSigma);
        forAll(row.nodes, k)
        {
            const scalar sigma = 1 - sqr(1 - scalar(k + 1)/nSigma);
            row.nodes[k] = makeNode(kd, sigma, waves, levels, Nmax);
        }

        Info<< "d/lambda " << dl
            << ": H_top/lambda " << row.kHtop/(2*pi)
            << ", omega " << waves.last().omega
            << ", N " << waves.last().N
            << " (" << tCont << " + " << t.timeIncrement() << " s)" << nl;
        for (const auto& nd : row.nodes)
        {
            Info<< "    sigma " << nd.sigma << ": N per level " << nd.N
                << ", nE " << nd.md.E.size() << nl;
        }
        Info<< flush;
    }

    OFstream os(args.get<fileName>(1));
    os  << "// Stream function wave table written by"
        << " generateStreamFunctionWaveTable" << nl
        << "// nRows " << nRows << ", dlMin " << dlMin << ", dlMax " << dlMax
        << ", nSigma " << nSigma << ", tol " << tol << ", Nmax " << Nmax
        << nl;
    streamFunctionWaveTable::write(os, levels, tol, Nmax, table);

    Info<< "Wrote " << args.get<fileName>(1) << " in " << total.elapsedTime()
        << " s" << endl;

    return 0;
}


// ************************************************************************* //
