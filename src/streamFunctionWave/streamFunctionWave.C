/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2026 Johan Roenby
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
#include "mathematicalConstants.H"
#include "OFstream.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::streamFunctionWave::streamFunctionWave(const dictionary& waveDict)
{
    readParmsFromDict(waveDict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::streamFunctionWave::readParmsFromDict(const dictionary& waveDict)
{
    N_ = waveDict.getOrDefault<label>("N", 20);

    const scalarField Btmp(N_, Zero);
    const scalarField Etmp(N_ + 1, Zero);

    E_ = waveDict.getOrDefault<scalarField>("E", Etmp);
    B_ = waveDict.getOrDefault<scalarField>("B", Btmp);
    H_ = waveDict.getOrDefault<scalar>("H", 0);
    d_ = waveDict.getOrDefault<scalar>("d", 0);
    uType_ = waveDict.getOrDefault<label>("uType", 2);
    nHSteps_ = waveDict.getOrDefault<label>("nHSteps", 20);
    g_ = waveDict.getOrDefault<scalar>("g", 9.81);
    k_ = waveDict.getOrDefault<scalar>("k", 0);
    R_ = waveDict.getOrDefault<scalar>("R", 0);
    uBar_ = waveDict.getOrDefault<scalar>("uBar", 0);
    Q_ = waveDict.getOrDefault<scalar>("Q", 0);

    if (uType_ == 0)
    {
        const scalar lambda = waveDict.get<scalar>("lambda");
        k_ = 2*pi/lambda;
        omega_ = 0;     // Co-moving frame
        // u is not used for anything in this case
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
            << "uType has to be 0, 1 or 2, but was " << uType_
            << exit(FatalIOError);
    }
}


bool Foam::streamFunctionWave::sameWave(const dictionary& waveDict) const
{
    // A dictionary written by an older version has no nHSteps entry, so it
    // cannot be matched and the wave is recalculated
    label nHSteps(0);
    if (!waveDict.readIfPresent("nHSteps", nHSteps))
    {
        return false;
    }

    const scalar H(waveDict.get<scalar>("H"));
    const scalar d(waveDict.get<scalar>("d"));
    const label N(waveDict.get<label>("N"));
    const label uType(waveDict.get<label>("uType"));

    bool sameWave(true);

    const scalar prec(Foam::pow(10.0, 1 - scalar(IOstream::defaultPrecision())));

    sameWave = sameWave && (N == N_);
    sameWave = sameWave && (nHSteps == nHSteps_);
    sameWave = sameWave && (uType == uType_);
    sameWave = sameWave && (mag(H/d - H_/d_) < prec);

    const scalar k(waveDict.get<scalar>("k"));
    const scalar omega(waveDict.get<scalar>("omega"));
    const scalar u(waveDict.get<scalar>("u"));

    if (uType_ == 0)
    {
        sameWave = sameWave && (mag(k/d - k_/d_) < prec);
    }
    else
    {
        sameWave = sameWave && (mag(omega - omega_) < prec);
        sameWave = sameWave && (mag(u - u_) < prec);
    }

    return sameWave;
}


void Foam::streamFunctionWave::calcCoeffs()
{
    scalar lort;
    const scalar g = mag(g_);

    if (uType_ == 0)
    {
        // k_ set to 2*pi/lambda in readParmsFromDict()
        const scalar lambda = 2*pi/k_;
        lort = lambda/d_;   // Non-dimensionalised wave length
        omega_ = 0;
        // u_ not used in this case
    }
    else
    {
        const scalar tau = 2*pi/omega_;

        // Non-dimensionalised wave period, T*sqrt(g/d)
        lort = tau*Foam::sqrt(g/d_);

        k_ = sqr(omega_)
           /(
                g
               *Foam::pow
                (
                    Foam::tanh(Foam::pow(omega_*Foam::sqrt(d_/g), 3.0/2.0)),
                    2.0/3.0
                )
            );

        u_ /= Foam::sqrt(g*d_);
    }

    // Initial guess for the variables from 1st order theory
    // (Fenton 1999, p. 14)

    const scalar kd = k_*d_;    // Mean depth normalised by wave number, k*d

    // Non-dimensionalised surface elevation above sea level, keta = k*eta
    scalarField keta(N_ + 1);
    forAll(keta, ii)
    {
        keta[ii] = kd + 0.5*k_*(H_/nHSteps_)*Foam::cos(ii*pi/N_);
    }

    // Non-dimensionalised mean horizontal velocity in the co-moving frame,
    // Ubar*sqrt(k/g)
    scalar u = Foam::sqrt(Foam::tanh(kd));

    B_ = 0;     // Fourier coefficients
    B_[0] = 0.5*k_*(H_/nHSteps_)/u;

    // Non-dimensionalised volume flux in the co-moving frame, q = Q*sqrt(k^3/g)
    scalar q = kd*u;

    // Non-dimensionalised constant in the Bernoulli equation, R*k/g
    scalar r = kd + 0.5*sqr(u);

    // Putting the guessed values into the x vector
    scalarField x(2*N_ + 5);
    forAll(keta, ii)
    {
        x[ii] = keta[ii];
    }
    forAll(B_, ii)
    {
        x[ii + N_ + 1] = B_[ii];
    }
    x[2*N_ + 1] = u;
    x[2*N_ + 2] = kd;
    x[2*N_ + 3] = q;
    x[2*N_ + 4] = r;

    // Finding root in parameter space
    for (label ii = 1; ii <= nHSteps_; ++ii)
    {
        findRoot(x, ii*H_/nHSteps_/d_, lort);

        Info<< "Iteration " << ii << " of " << nHSteps_
            << " with H/d = " << ii*H_/nHSteps_/d_
            << " gave kd = " << x[2*N_ + 2] << endl;
    }

    // Converting the x vector into named variables
    forAll(keta, ii)
    {
        keta[ii] = x[ii];
    }
    forAll(B_, ii)
    {
        B_[ii] = x[ii + N_ + 1];
    }
    u = x[2*N_ + 1];
    q = x[2*N_ + 3];

    k_ = x[2*N_ + 2]/d_;
    R_ = x[2*N_ + 4]*g/k_;

    // Calculating E_ using eq. (3.26) (with j running from 0 to N) and
    // eq. (3.27) (with the RHS divided by N to correct the equation) in
    // Fenton 1999. With these corrections A_[j] = 2/N*E_{j+1} for
    // j = [0,N-2], and A_[N-1] = E_{N}/N.
    scalarField etapp = keta/k_;
    etapp[0] *= 0.5;
    etapp[N_] *= 0.5;

    E_ = 0;
    forAll(E_, jj)
    {
        forAll(etapp, mm)
        {
            E_[jj] += etapp[mm]*Foam::cos(jj*mm*pi/N_);
        }
    }
    E_[0] *= 0.5;
    E_[N_] *= 0.5;
    E_ /= N_;

    uBar_ = u/Foam::sqrt(k_/g);
    Q_ = q/Foam::sqrt(k_*k_*k_/g);

    Info<< "Done calculating stream function wave" << endl;
}


void Foam::streamFunctionWave::findRoot
(
    scalarField& X,
    scalar h,
    scalar lort
)
{
    const scalar tol(1e-5);     // HARDCODED TOLERANCE
    const label maxIter(50);    // HARDCODED ITERATION LIMIT

    scalar err(2*tol);
    label nIter(0);

    while (err > tol)
    {
        scalarField F0(fun(X, h, lort));
        scalarSquareMatrix M(-jac(X, h, lort));

        LUsolve(M, F0);         // LUsolve puts dX into F0

        X += F0;
        err = max(mag(F0));
        ++nIter;

        Info<< "Root iteration: " << nIter << " Error = " << err << endl;

        if (nIter >= maxIter && err > tol)
        {
            FatalErrorInFunction
                << "Stream function wave did not converge in " << maxIter
                << " iterations. Error = " << err
                << " with tolerance " << tol << nl
                << "The wave of height H/d = " << h
                << " is probably too steep. Try increasing nHSteps or N."
                << exit(FatalError);
        }
    }
}


Foam::scalarField Foam::streamFunctionWave::fun
(
    const scalarField& x,
    scalar h,
    scalar lort
) const
{
    // Converting from the x vector into named variables
    scalarField keta(N_ + 1);
    forAll(keta, ii)
    {
        keta[ii] = x[ii];
    }

    scalarField B(N_);
    forAll(B, ii)
    {
        B[ii] = x[ii + N_ + 1];
    }

    const scalar u = x[2*N_ + 1];
    const scalar kd = x[2*N_ + 2];
    const scalar q = x[2*N_ + 3];
    const scalar r = x[2*N_ + 4];

    // Kinematic boundary condition, (3.7) in Fenton 1999
    scalarField kinBC(N_ + 1);

    // Dynamic boundary condition, (3.8) in Fenton 1999
    scalarField dynBC(N_ + 1);

    // Mean depth equation, (3.9) in Fenton 1999
    scalar meanDepthEq = -(0.5/N_)*(keta[0] + keta[N_]) - kd;

    forAll(kinBC, ii)
    {
        const scalar kx = ii*pi/N_;
        const scalar ke = keta[ii];

        kinBC[ii] = -u*ke + q;
        dynBC[ii] = ke - r;

        scalar U = -u;
        scalar V = 0;

        forAll(B, jj)
        {
            const label j = jj + 1;

            kinBC[ii] +=
                B[jj]*(Foam::sinh(j*ke)/Foam::cosh(j*kd))*Foam::cos(j*kx);
            U += j*B[jj]*(Foam::cosh(j*ke)/Foam::cosh(j*kd))*Foam::cos(j*kx);
            V += j*B[jj]*(Foam::sinh(j*ke)/Foam::cosh(j*kd))*Foam::sin(j*kx);
        }

        dynBC[ii] += 0.5*sqr(U) + 0.5*sqr(V);
        meanDepthEq += keta[ii]/N_;
    }

    // Wave height equation, (3.10) in Fenton 1999
    const scalar waveHeightEq = keta[0] - keta[N_] - kd*h;

    scalar lastEq = 0;
    if (uType_ == 0)
    {
        // Relative wave length equation, (3.11) in Fenton 1999
        lastEq = kd - 2.0*pi/lort;
    }
    else if (uType_ == 1)
    {
        // Equation for the mean current in the laboratory frame,
        // (3.15) in Fenton 1999
        lastEq = Foam::sqrt(kd)*u + kd*u_ - 2.0*pi/lort;
    }
    else if (uType_ == 2)
    {
        // Equation for the mass transport velocity, (3.16) in Fenton 1999
        lastEq = q/Foam::sqrt(kd) + kd*u_ - 2.0*pi/lort;
    }

    // Writing the wave conditions into the solution vector F
    scalarField F(2*N_ + 5);
    forAll(kinBC, ii)
    {
        F[ii] = kinBC[ii];
        F[ii + N_ + 1] = dynBC[ii];
    }
    F[2*N_ + 2] = meanDepthEq;
    F[2*N_ + 3] = waveHeightEq;
    F[2*N_ + 4] = lastEq;

    return F;
}


Foam::scalarSquareMatrix Foam::streamFunctionWave::jac
(
    const scalarField& x,
    scalar h,
    scalar lort
) const
{
    const scalar d(1e-10);      // HARDCODED STEP SIZE

    // Constant across the columns, so evaluated once
    const scalarField F0(fun(x, h, lort));

    scalarSquareMatrix J(x.size());
    forAll(x, ii)
    {
        scalarField dx(x.size(), Zero);
        dx[ii] = d;

        const scalarField dF((fun(x + dx, h, lort) - F0)/d);

        forAll(dF, jj)
        {
            J[jj][ii] = dF[jj];
        }
    }

    return J;
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
        u += jj*B_[ii]*cosh(jj*k_*y)/cosh(jj*k_*d_)*cos(jj*theta);
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
        v += jj*B_[ii]*sinh(jj*k_*y)/cosh(jj*k_*d_)*sin(jj*theta);
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
    Info<< "Enter writeDict(const fileName& outDictName)" << endl;

    IOstream::defaultPrecision(14);

    dictionary outDict(outDictName);

    outDict.set("H", H_);
    outDict.set("d", d_);
    if (omega_ > 0)     // Period is finite
    {
        outDict.set("T", 2*pi/omega_);
    }
    outDict.set("uType", uType_);
    outDict.set("u", u_);
    outDict.set("N", N_);
    outDict.set("nHSteps", nHSteps_);

    const scalar lambda = 2*pi/k_;
    outDict.set("lambda", lambda);
    outDict.set("k", k_);
    outDict.set("omega", omega_);
    outDict.set("c", omega_/k_);
    outDict.set("R", R_);
    outDict.set("uBar", uBar_);
    outDict.set("B", B_);
    outDict.set("E", E_);
    outDict.set("Q", Q_);
    outDict.set("T0", lambda*d_/Q_);

    scalarField x(2);
    x[0] = 0;
    x[1] = lambda/2.0;

    const scalarField eta(getEta(x, 0, 0));
    const scalarField u(getU(x, eta, 0, 0));

    outDict.set("h_crest", eta[0] - d_);
    outDict.set("h_trough", eta[1] - d_);
    outDict.set("u_crest", u[0]);
    outDict.set("u_trough", u[1]);
    outDict.set("steepness", H_/lambda);

    outDict.write(OFstream(outDictName)(), false);
}


// ************************************************************************* //
