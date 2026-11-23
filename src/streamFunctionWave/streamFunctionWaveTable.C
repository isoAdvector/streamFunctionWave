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

\*---------------------------------------------------------------------------*/

#include "streamFunctionWaveTable.H"
#include "mathematicalConstants.H"
#include "Ostream.H"
#include "error.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Index of the interval of the sorted xs that holds x, clipped to the
// intervals that exist
Foam::label interval(const Foam::UList<Foam::scalar>& xs, const Foam::scalar x)
{
    const Foam::label n = xs.size();
    Foam::label i = 0;
    while (i < n - 2 && x >= xs[i + 1])
    {
        ++i;
    }
    return i;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::streamFunctionWaveTable::streamFunctionWaveTable
(
    const double* data,
    const label size
)
:
    data_(data),
    genTol_(0),
    genNmax_(0)
{
    label pos = 0;
    auto next = [&]() -> double
    {
        if (pos >= size)
        {
            FatalErrorInFunction
                << "The stream function wave table ends early"
                << abort(FatalError);
        }
        return data_[pos++];
    };

    const label version = label(next());
    if (version != 1)
    {
        FatalErrorInFunction
            << "Unknown stream function wave table version " << version
            << abort(FatalError);
    }

    levels_.resize(label(next()));
    for (scalar& l : levels_)
    {
        l = next();
    }
    genTol_ = next();
    genNmax_ = label(next());

    const label nRows = label(next());
    x_.resize(nRows);
    kd_.resize(nRows);
    HtopL_.resize(nRows);
    node_.resize(nRows);

    const label L = levels_.size();
    for (label i = 0; i < nRows; ++i)
    {
        x_[i] = std::log(next());
        kd_[i] = next();
        HtopL_[i] = next()/(2*pi);

        labelList& nodes = node_[i];
        nodes.resize(label(next()));
        scalarList t(nodes.size() + 1, Zero);
        forAll(nodes, k)
        {
            nodes[k] = pos;
            t[k + 1] = 1 - std::sqrt(1 - data_[pos]);
            pos += 2 + L + 3;
            const label nE = label(next());
            pos += nE;
        }

        if (i == 0)
        {
            t_ = t;
        }
        else if (t.size() != t_.size() || max(mag(t - t_)) > 1e-12)
        {
            FatalErrorInFunction
                << "Row " << i << " of the stream function wave table has"
                << " other nodes than the first"
                << abort(FatalError);
        }
    }

    if (pos != size)
    {
        FatalErrorInFunction
            << "The stream function wave table has " << size - pos
            << " numbers more than its rows"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::streamFunctionWaveTable::weights
(
    const UList<scalar>& xs,
    const scalar x,
    labelList& idx,
    scalarList& w
)
{
    const label n = xs.size();
    const label order = min(3, n - 1);
    const label i = interval(xs, x);
    const label lo = max(0, min(i - (order - 1)/2, n - order - 1));

    idx.resize(order + 1);
    w.resize(order + 1);
    forAll(idx, a)
    {
        idx[a] = lo + a;
    }
    forAll(idx, a)
    {
        scalar p = 1;
        forAll(idx, b)
        {
            if (b != a)
            {
                p *= (x - xs[idx[b]])/(xs[idx[a]] - xs[idx[b]]);
            }
        }
        w[a] = p;
    }
}


void Foam::streamFunctionWaveTable::node
(
    const label i,
    const label k,
    const UList<scalar>& theta,
    scalarField& ke,
    scalar& u,
    scalar& qr,
    scalar& rr
) const
{
    ke.resize(theta.size());

    if (k == 0)
    {
        // Still water
        ke = Zero;
        u = std::sqrt(std::tanh(kd_[i]));
        qr = 0;
        rr = 0.5*u*u;
        return;
    }

    const label L = levels_.size();
    const double* nd = data_ + node_[i][k - 1];
    u = nd[2 + L];
    qr = nd[3 + L] - u*kd_[i];
    rr = nd[4 + L] - kd_[i];

    const label nE = min(label(nd[5 + L]), theta.size());
    const double* E = nd + 6 + L;
    forAll(theta, p)
    {
        scalar s = -kd_[i];
        for (label m = 0; m < nE; ++m)
        {
            s += E[m]*std::cos(m*theta[p]);
        }
        ke[p] = s;
    }
}


Foam::label Foam::streamFunctionWaveTable::nodeN
(
    const label i,
    const label k,
    const label l
) const
{
    if (k == 0)
    {
        return 1;
    }
    return label(data_[node_[i][k - 1] + 2 + l]);
}


void Foam::streamFunctionWaveTable::coords
(
    const scalar dl,
    const scalar sigma,
    scalar& x,
    scalar& t
) const
{
    x = min(max(std::log(dl), x_.first()), x_.last());
    t = 1 - std::sqrt(1 - min(max(sigma, scalar(0)), scalar(1)));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::streamFunctionWaveTable::HtopByLambda(const scalar dl) const
{
    const scalar x = min(max(std::log(dl), x_.first()), x_.last());
    labelList idx;
    scalarList w;
    weights(x_, x, idx, w);
    scalar H = 0;
    forAll(idx, a)
    {
        H += w[a]*HtopL_[idx[a]];
    }
    return H;
}


Foam::scalar Foam::streamFunctionWaveTable::sigma
(
    const scalar dl,
    const scalar h
) const
{
    if (dl < dlMin())
    {
        return h*dlMin()/HtopL_.first();
    }
    return h*dl/HtopByLambda(dl);
}


Foam::label Foam::streamFunctionWaveTable::guessN
(
    const scalar dl,
    const scalar sigma
) const
{
    scalar x, t;
    coords(dl, sigma, x, t);

    // The tightest level that is no tighter than the guess level
    label l = 0;
    forAll(levels_, li)
    {
        if (levels_[li] >= guessLevel*(1 - 1e-6))
        {
            l = li;
        }
    }

    // The corners of the cell
    const label ix = interval(x_, x);
    const label it = interval(t_, t);
    label N = 4;
    for (label i = ix; i <= min(ix + 1, x_.size() - 1); ++i)
    {
        for (label k = it; k <= min(it + 1, t_.size() - 1); ++k)
        {
            const label n = nodeN(i, k, l);
            N = max(N, n < 0 ? genNmax_ : n);
        }
    }
    return N;
}


Foam::tmp<Foam::scalarField> Foam::streamFunctionWaveTable::guess
(
    const scalar dl,
    const scalar sigma,
    const streamFunctionWaveSolver& s
) const
{
    const label M = s.M();
    const label N = s.N();
    const scalarField& theta = s.theta();
    const scalar kd = 2*pi*dl;

    scalar x, t;
    coords(dl, sigma, x, t);
    labelList ix, it;
    scalarList wx, wt;
    weights(x_, x, ix, wx);
    weights(t_, t, it, wt);

    // Interpolated surface, u, q - u k d and r - k d
    scalarField ke(M, Zero), keNode;
    scalar u = 0, qr = 0, rr = 0;
    forAll(ix, a)
    {
        forAll(it, b)
        {
            const scalar w = wx[a]*wt[b];
            scalar un, qn, rn;
            node(ix[a], it[b], theta, keNode, un, qn, rn);
            ke += w*keNode;
            u += w*un;
            qr += w*qn;
            rr += w*rn;
        }
    }
    ke += kd;
    const scalar q = qr + u*kd;
    const scalar r = rr + kd;

    auto tz = tmp<scalarField>::New(s.nUnknowns(), Zero);
    auto& z = tz.ref();
    for (label p = 0; p < M; ++p)
    {
        z[p] = ke[p];
    }
    z[s.iu()] = u;
    z[s.iq()] = q;
    z[s.ir()] = r;
    if (s.params().freeKd)
    {
        z[M + N + 3] = kd;
    }

    // The B_j from the kinematic condition, at about the accuracy of the
    // guess
    s.fitB(z, 1e-8);

    return tz;
}


void Foam::streamFunctionWaveTable::uq
(
    const scalar kd,
    const scalar h,
    scalar& u,
    scalar& q
) const
{
    const scalar dl = kd/(2*pi);
    scalar x, t;
    coords(dl, sigma(dl, h), x, t);
    labelList ix, it;
    scalarList wx, wt;
    weights(x_, x, ix, wx);
    weights(t_, t, it, wt);

    const scalarList noPoints;
    scalarField ke;
    u = 0;
    scalar qr = 0;
    forAll(ix, a)
    {
        forAll(it, b)
        {
            const scalar w = wx[a]*wt[b];
            scalar un, qn, rn;
            node(ix[a], it[b], noPoints, ke, un, qn, rn);
            u += w*un;
            qr += w*qn;
        }
    }
    q = qr + u*kd;
}


Foam::scalar Foam::streamFunctionWaveTable::kdForPeriod
(
    const streamFunctionWaveSolver::parameters& p,
    scalar& kdOther
) const
{
    auto f = [&](const scalar kd)
    {
        scalar u, q;
        uq(kd, p.h, u, q);
        return streamFunctionWaveSolver::currentEquation(p, u, q, kd);
    };

    // Bracket the roots on a logarithmic grid of d/lambda from 0.01 to 10
    const label n = 60;
    scalarList kd(n), fs(n);
    forAll(kd, i)
    {
        kd[i] = 2*pi*0.01*std::pow(1000.0, scalar(i)/(n - 1));
        fs[i] = f(kd[i]);
    }
    labelList roots;
    for (label i = 0; i < n - 1; ++i)
    {
        if ((fs[i] < 0) != (fs[i + 1] < 0))
        {
            roots.append(i);
        }
    }

    kdOther = -1;
    if (roots.empty())
    {
        return -1;
    }
    if (roots.size() > 1)
    {
        kdOther = 0.5*(kd[roots[1]] + kd[roots[1] + 1]);
    }

    // The first root, by regula falsi with the Illinois modification
    scalar a = kd[roots[0]], b = kd[roots[0] + 1];
    scalar fa = fs[roots[0]], fb = fs[roots[0] + 1];
    scalar m = b;
    for (label iter = 0; iter < 100; ++iter)
    {
        m = b - fb*(b - a)/(fb - fa);
        const scalar fm = f(m);
        if (mag(b - a) < 1e-14*b || fm == 0)
        {
            break;
        }
        if ((fm < 0) == (fb < 0))
        {
            fa *= 0.5;
        }
        else
        {
            a = b;
            fa = fb;
        }
        b = m;
        fb = fm;
    }
    return m;
}


void Foam::streamFunctionWaveTable::write
(
    Ostream& os,
    const UList<scalar>& levels,
    const scalar genTol,
    const label genNmax,
    const UList<rowData>& rows
)
{
    const int oldPrecision = os.precision();
    os.precision(17);

    label count = 0;
    auto put = [&](const scalar v)
    {
        os  << v << ',';
        if (++count % 4 == 0)
        {
            os  << nl;
        }
        else
        {
            os  << ' ';
        }
    };
    auto newLine = [&]()
    {
        if (count % 4 != 0)
        {
            os  << nl;
            count = 0;
        }
    };

    put(1);
    put(levels.size());
    for (const scalar l : levels)
    {
        put(l);
    }
    put(genTol);
    put(genNmax);
    put(rows.size());
    newLine();

    for (const rowData& row : rows)
    {
        os  << "// d/lambda " << row.dl << nl;
        put(row.dl);
        put(row.kd);
        put(row.kHtop);
        put(row.nodes.size());
        newLine();
        for (const nodeData& nd : row.nodes)
        {
            put(nd.sigma);
            put(nd.kH);
            for (const label N : nd.N)
            {
                put(N);
            }
            put(nd.md.u);
            put(nd.md.q);
            put(nd.md.r);
            put(nd.md.E.size());
            for (const scalar E : nd.md.E)
            {
                put(E);
            }
            newLine();
        }
    }

    os.precision(oldPrecision);
}


// ************************************************************************* //
