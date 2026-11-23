/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 DHI A/S
    Copyright (C) 2017-2026 Johan Roenby
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

#include "streamFunctionWave.H"
#include "streamFunctionWaveSolver.H"
#include "streamFunctionWaveTable.H"
#include "mathematicalConstants.H"
#include "OFstream.H"
#include "clockTime.H"
#include "SubField.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// sinh(a)/cosh(b) and cosh(a)/cosh(b) for 0 <= a, written so that they do
// not overflow for the large arguments of deep water (j k d > 710)

Foam::tmp<Foam::scalarField> sinhByCosh
(
    const Foam::scalarField& a,
    const Foam::scalar b
)
{
    return (Foam::exp(a - b) - Foam::exp(-a - b))/(1 + Foam::exp(-2*b));
}

Foam::tmp<Foam::scalarField> coshByCosh
(
    const Foam::scalarField& a,
    const Foam::scalar b
)
{
    return (Foam::exp(a - b) + Foam::exp(-a - b))/(1 + Foam::exp(-2*b));
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::streamFunctionWave::streamFunctionWave(const dictionary& waveDict)
{
    readParmsFromDict(waveDict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::streamFunctionWave::readParmsFromDict(const dictionary& waveDict)
{
    // The results of an earlier calculation are cleared
    *this = streamFunctionWave();

    // Once per run, not for every field read by a post-processing utility
    static bool warnedNModes = false;
    if (!warnedNModes && (waveDict.found("N") || waveDict.found("nHSteps")))
    {
        warnedNModes = true;
        WarningInFunction
            << "N and nHSteps are no longer used and are ignored: the"
            << " number of Fourier modes follows from the tolerance"
            << " (default 1e-8), and the wave is started from a table"
            << " of precomputed waves" << endl;
    }

    H_ = waveDict.getOrDefault<scalar>("H", 0);
    d_ = waveDict.getOrDefault<scalar>("d", 0);
    uType_ = waveDict.getOrDefault<label>("uType", 2);
    tolerance_ = waveDict.getOrDefault<scalar>("tolerance", 1e-8);
    g_ = waveDict.getOrDefault<scalar>("g", 9.81);

    if (tolerance_ <= 0 || tolerance_ > 1e-3)
    {
        FatalIOErrorInFunction(waveDict)
            << "tolerance has to be in (0, 1e-3], but was " << tolerance_
            << exit(FatalIOError);
    }

    if (uType_ == 0)
    {
        const scalar lambda = waveDict.get<scalar>("lambda");
        k_ = 2*pi/lambda;
        omega_ = 0;     // Co-moving frame
        // u is not used for anything in this case
    }
    else if (uType_ == 3)
    {
        const scalar lambda = waveDict.get<scalar>("lambda");
        k_ = 2*pi/lambda;
        const scalar tau = waveDict.get<scalar>("T");
        omega_ = 2*pi/tau;
        // The current follows from lambda and T, so u is not used
    }
    else if ((uType_ == 1) || (uType_ == 2))
    {
        // k_ is calculated by calcCoeffs()
        const scalar tau = waveDict.get<scalar>("T");
        omega_ = 2*pi/tau;
        u_ = waveDict.getOrDefault<scalar>("u", 0);
    }
    else
    {
        FatalIOErrorInFunction(waveDict)
            << "uType has to be 0, 1, 2 or 3, but was " << uType_
            << exit(FatalIOError);
    }
}


bool Foam::streamFunctionWave::sameInput(const streamFunctionWave& w) const
{
    // k_ is the input for wave length input only: for period input
    // calcCoeffs() calculates it
    return
        H_ == w.H_ && d_ == w.d_ && uType_ == w.uType_
     && tolerance_ == w.tolerance_ && g_ == w.g_ && omega_ == w.omega_
     && u_ == w.u_ && (!lengthInput() || k_ == w.k_);
}


void Foam::streamFunctionWave::calcCoeffs()
{
    typedef streamFunctionWaveSolver solver;

    clockTime timer;
    const scalar g = mag(g_);
    const streamFunctionWaveTable& table = streamFunctionWaveTable::builtIn();

    if (table.empty())
    {
        FatalErrorInFunction
            << "The stream function wave table is empty. Generate it with"
            << " generateStreamFunctionWaveTable and rebuild the library."
            << exit(FatalError);
    }
    if (!(H_ > 0) || !(d_ > 0))
    {
        FatalErrorInFunction
            << "The wave height H and the depth d have to be positive, but"
            << " were " << H_ << " and " << d_
            << exit(FatalError);
    }

    const scalar h = H_/d_;
    solver::parameters p;
    p.h = h;

    if (lengthInput())
    {
        // k_ set to 2*pi/lambda in readParmsFromDict(). omega_ is zero for
        // uType 0, giving the co-moving frame, and 2*pi/T for uType 3,
        // giving the frame in which the wave moves with celerity lambda/T.
        p.kd = k_*d_;
    }
    else
    {
        // k d is an unknown, which the table gives a first value of
        p.freeKd = true;
        p.uType = uType_;
        p.tau = 2*pi/omega_*Foam::sqrt(g/d_);
        p.uStar = u_/Foam::sqrt(g*d_);

        scalar kdOther;
        p.kd = table.kdForPeriod(p, kdOther);
        if (p.kd < 0)
        {
            FatalErrorInFunction
                << "There is no wave of period T = " << 2*pi/omega_
                << " with the current u = " << u_ << " (uType " << uType_
                << ") in depth d = " << d_ << ". An opposing current may"
                << " be too strong for the wave to propagate."
                << exit(FatalError);
        }
        if (kdOther > 0)
        {
            Info<< "Note: a short wave of about lambda = "
                << 2*pi*d_/kdOther << " has the same period and current;"
                << " the longest wave, lambda about " << 2*pi*d_/p.kd
                << ", is calculated" << endl;
        }
    }

    const scalar dl = p.kd/(2*pi);
    sigma_ = table.sigma(dl, h);

    if (sigma_ > 1)
    {
        const scalar HtopL =
            dl < table.dlMin()
          ? table.HtopByLambda(table.dlMin())*dl/table.dlMin()
          : table.HtopByLambda(dl);
        FatalErrorInFunction
            << "The wave of height H = " << H_ << " (H/lambda = "
            << H_*p.kd/(2*pi*d_) << ", H/d = " << h << ") is above the"
            << " highest wave of the table at this depth, H/lambda = "
            << HtopL << " (H/d = " << HtopL/dl << "), by "
            << 100*(sigma_ - 1) << " %." << nl
            << "The highest possible wave is about 1 % (deep water) to 9 %"
            << " (shallow water) higher than that, but too close to"
            << " breaking to be calculated reliably."
            << exit(FatalError);
    }
    if (dl < table.dlMin())
    {
        WarningInFunction
            << "The relative depth d/lambda = " << dl << " is below the"
            << " table (" << table.dlMin() << "). The wave is started from"
            << " the shallowest waves of the table, which may be slow or"
            << " fail." << endl;
    }

    // The round-off level of the error measure, which is relative to the
    // wave height, so higher for small waves
    const scalar tolMin = 1e-12 + 1e-14/h;
    scalar tol = tolerance_;
    if (tol < tolMin)
    {
        WarningInFunction
            << "The tolerance " << tolerance_ << " is below the round-off"
            << " level " << tolMin << " of this wave and is raised to it"
            << endl;
        tol = tolMin;
    }

    // The table's number of modes for the accuracy of its guess, then more
    // until two successive numbers of modes agree to half the tolerance
    label N = table.guessN(dl, sigma_);
    autoPtr<solver> s(new solver(N, p));
    scalarField z(table.guess(dl, sigma_, *s));
    solver::result r = s->solve(z, tol);
    label nIter = r.iterations;

    const label Nmax = 2*table.genNmax();
    DynamicList<scalar> errs;
    label roundOff = 0;

    // The coarser solution of the pair with the smallest error estimate
    autoPtr<solver> best;
    scalarField bestZ;
    scalar bestErr = GREAT;
    while (true)
    {
        const label N2 = solver::grow(N);
        if (N2 > Nmax && best.valid() && bestErr < 1e-9)
        {
            // At round-off without having detected it: the best level
            WarningInFunction
                << "The stream function wave reached round-off at an"
                << " estimated error of " << bestErr << ", above the"
                << " tolerance " << tol << endl;
            s = std::move(best);
            z = bestZ;
            N = s->N();
            errs.append(bestErr);
            break;
        }
        if (N2 > Nmax)
        {
            FatalErrorInFunction
                << "The stream function wave did not converge to the"
                << " tolerance " << tol << " with " << N << " modes (the"
                << " last difference between successive numbers of modes"
                << " was " << (errs.empty() ? GREAT : errs.last()) << ")."
                << " The wave is probably too close to breaking for this"
                << " tolerance."
                << exit(FatalError);
        }

        autoPtr<solver> s2(new solver(N2, p));
        scalarField z2(s2->fromModal(s->toModal(z)));
        const solver::result r2 = s2->solve(z2, tol);
        nIter += r2.iterations;

        // The difference, or the uncertainty of a solution whose iteration
        // stalled at round-off if larger
        const bool ok =
            r.converged && r2.converged && s->regular(z) && s2->regular(z2);
        const scalar diff = ok ? s->difference(z, *s2, z2) : GREAT;
        const scalar stall = max(r.step, r2.step);
        const scalar err = max(diff, stall);
        errs.append(err);
        if (err < bestErr)
        {
            bestErr = err;
            best.reset(new solver(N, p));
            bestZ = z;
        }

        if (err <= 0.5*tol)
        {
            break;
        }

        // Round-off reached before the tolerance, at two successive numbers
        // of modes: the error estimate is the uncertainty of the stalled
        // iteration, not the truncation, or it is below 1e-9 and has stopped
        // falling. Near the highest wave the round-off is about 1e-9, since
        // the height there hardly determines the wave, and for steep waves
        // generally about 1e-11.
        const label n = errs.size();
        const bool flat =
            n >= 2 && err < 1e-9 && errs[n - 1] > 0.5*errs[n - 2];
        roundOff = (ok && (stall >= diff || flat)) ? roundOff + 1 : 0;
        if (roundOff >= 2)
        {
            WarningInFunction
                << "The stream function wave reached round-off at an"
                << " estimated error of " << err << ", above the tolerance "
                << tol << ". Close to the highest wave of the table the"
                << " prescribed height determines the wave only to about"
                << " 1e-9." << endl;
            break;
        }

        s = std::move(s2);
        z = z2;
        r = r2;
        N = N2;
    }

    // The coarser of the two: its high Fourier coefficients are determined,
    // and the wave is cheaper to evaluate. Its error is about the difference.
    error_ = errs.last();
    N_ = N;
    const label M = s->M();
    if (!lengthInput())
    {
        k_ = s->kdOf(z)/d_;
    }

    B_ = SubField<scalar>(z, N, M);

    // E_j of eta = sum_j 2 E_j cos(j theta), without the trailing ones that
    // together are below 1 % of the tolerance
    const solver::modal md(s->toModal(z));
    label nE = md.E.size();
    scalar tail = 0;
    while (nE > 2 && tail + mag(md.E[nE - 1]) < 0.01*tol*s->kH(z))
    {
        tail += mag(md.E[--nE]);
    }
    E_.resize(nE);
    forAll(E_, j)
    {
        E_[j] = 0.5*md.E[j]/k_;
    }

    const scalar u = z[s->iu()];
    const scalar q = z[s->iq()];
    const scalar rr = z[s->ir()];
    uBar_ = u/Foam::sqrt(k_/g);
    Q_ = q/Foam::sqrt(k_*k_*k_/g);
    R_ = rr*g/k_;

    Info<< "Stream function wave: lambda = " << 2*pi/k_
        << ", d/lambda = " << d_*k_/(2*pi)
        << ", H/lambda = " << H_*k_/(2*pi)
        << ", H/H_top = " << sigma_ << nl
        << "    " << N_ << " Fourier modes, estimated error " << error_
        << " (tolerance " << tol << "), " << nIter << " iterations, "
        << timer.elapsedTime() << " s" << endl;
}


Foam::scalar Foam::streamFunctionWave::zeroPointPhase() const
{
    const scalar lambda = 2.0*pi/k_;
    const scalar tol = 1e-6;

    scalarField x0(1, Zero);
    scalarField x1(1, lambda/2.0);
    scalarField xMean(1, Zero);

    scalar err = GREAT;

    while (mag(err) > tol)
    {
        xMean[0] = 0.5*(x0[0] + x1[0]);

        const scalarField etatmp(getEta(xMean, 0.0, 0.0));
        err = (etatmp[0] - d_)/d_;

        if (err > 0)
        {
            x0[0] = xMean[0];
        }
        else
        {
            x1[0] = xMean[0];
        }
    }

    return xMean[0]/lambda;
}


Foam::scalarField Foam::streamFunctionWave::getEta
(
    const scalarField& x,
    const scalar t,
    const scalar phi
) const
{
    const scalarField theta(k_*x - omega_*t + 2*pi*phi);

    scalarField eta(x.size(), Zero);
    forAll(E_, jj)
    {
        // Fenton (1999), eq. (3.27)
        eta += 2*E_[jj]*cos(jj*theta);
    }

    return eta;
}


Foam::scalarField Foam::streamFunctionWave::getU
(
    const scalarField& x,
    const scalarField& y,
    const scalar t,
    const scalar phi
) const
{
    const scalarField theta(k_*x - omega_*t + 2*pi*phi);

    scalarField u(y.size(), Zero);
    forAll(B_, ii)
    {
        const label jj = ii + 1;

        // Fenton (1999), eq. (3.24)
        u += jj*B_[ii]*coshByCosh(jj*k_*y, jj*k_*d_)*cos(jj*theta);
    }
    u *= sqrt(g_/k_);
    u += omega_/k_ - uBar_;

    return u;
}


Foam::scalarField Foam::streamFunctionWave::getV
(
    const scalarField& x,
    const scalarField& y,
    const scalar t,
    const scalar phi
) const
{
    const scalarField theta(k_*x - omega_*t + 2*pi*phi);

    scalarField v(y.size(), Zero);
    forAll(B_, ii)
    {
        const label jj = ii + 1;

        // Fenton (1999), eq. (3.25)
        v += jj*B_[ii]*sinhByCosh(jj*k_*y, jj*k_*d_)*sin(jj*theta);
    }
    v *= sqrt(g_/k_);

    return v;
}


Foam::scalarField Foam::streamFunctionWave::getPbyRho
(
    const scalarField& x,
    const scalarField& y,
    const scalar t,
    const scalar phi
) const
{
    const scalarField u(getU(x, y, t, phi));
    const scalarField v(getV(x, y, t, phi));
    const scalar c = omega_/(k_ + SMALL);

    // Fenton (1999), eq. (3.28)
    return R_ - g_*y - 0.5*(sqr(u - c) + sqr(v));
}


void Foam::streamFunctionWave::writeDict(const fileName& outDictName) const
{
    OFstream os(outDictName);
    os.precision(15);

    os.writeEntry("H", H_);
    os.writeEntry("d", d_);
    os.writeEntry("g", g_);
    if (omega_ > 0)     // Period is finite
    {
        os.writeEntry("T", 2*pi/omega_);
    }
    os.writeEntry("uType", uType_);
    os.writeEntry("u", u_);
    os.writeEntry("tolerance", tolerance_);
    os.writeEntry("N", N_);
    os.writeEntry("error", error_);
    os.writeEntry("HbyHtop", sigma_);

    const scalar lambda = 2*pi/k_;
    os.writeEntry("lambda", lambda);
    os.writeEntry("k", k_);
    os.writeEntry("omega", omega_);
    os.writeEntry("c", omega_/k_);

    // Eulerian mean current and mass transport velocity in the frame in
    // which the wave has celerity c
    os.writeEntry("u1", omega_/k_ - uBar_);
    os.writeEntry("u2", omega_/k_ - Q_/d_);
    os.writeEntry("R", R_);
    os.writeEntry("uBar", uBar_);
    os.writeEntry("B", B_);
    os.writeEntry("E", E_);
    os.writeEntry("Q", Q_);
    os.writeEntry("T0", lambda*d_/Q_);

    scalarField x(2);
    x[0] = 0;
    x[1] = lambda/2.0;

    const scalarField eta(getEta(x, 0, 0));
    const scalarField u(getU(x, eta, 0, 0));

    os.writeEntry("h_crest", eta[0] - d_);
    os.writeEntry("h_trough", eta[1] - d_);
    os.writeEntry("u_crest", u[0]);
    os.writeEntry("u_trough", u[1]);
    os.writeEntry("steepness", H_/lambda);
}


void Foam::streamFunctionWave::write(Ostream& os) const
{
    os.writeEntry("H", H_);
    os.writeEntry("d", d_);
    if (omega_ > 0)     // Period is finite
    {
        os.writeEntry("T", 2*pi/omega_);
    }
    os.writeEntry("lambda", 2*pi/k_);
    os.writeEntry("uType", uType_);
    os.writeEntry("u", u_);
    os.writeEntry("tolerance", tolerance_);
    os.writeEntry("N", N_);
    os.writeEntry("error", error_);
    os.writeEntry("HbyHtop", sigma_);
}


// ************************************************************************* //
