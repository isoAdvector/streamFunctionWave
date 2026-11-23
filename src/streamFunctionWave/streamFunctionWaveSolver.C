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

#include "streamFunctionWaveSolver.H"
#include "mathematicalConstants.H"
#include "SubField.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

using Foam::scalar;
using Foam::label;

// sinh(a)/cosh(b) and cosh(a)/cosh(b) for 0 <= a, written so that they do
// not overflow for the large arguments of deep water (j k d > 710)

inline scalar sinhByCosh(const scalar a, const scalar b)
{
    return (std::exp(a - b) - std::exp(-a - b))/(1 + std::exp(-2*b));
}

inline scalar coshByCosh(const scalar a, const scalar b)
{
    return (std::exp(a - b) + std::exp(-a - b))/(1 + std::exp(-2*b));
}


bool allFinite(const scalar* v, const label n)
{
    for (label i = 0; i < n; ++i)
    {
        if (!std::isfinite(v[i]))
        {
            return false;
        }
    }
    return true;
}

bool allFinite(const Foam::UList<scalar>& f)
{
    return allFinite(f.cdata(), f.size());
}

bool allFinite(const Foam::scalarRectangularMatrix& A)
{
    return allFinite(A.cdata(), A.size());
}


} // End anonymous namespace


// Least-squares solution of J x = b with the structure of the Jacobian:
// the surface unknowns k eta_m (columns 0..M-1) appear only in rows m and
// M + m (the two surface conditions at point m) and in the global rows 2M..
// They are eliminated exactly by orthogonal transformations, which leave a
// small problem in the remaining unknowns. J and b are destroyed.
Foam::tmp<Foam::scalarField> Foam::streamFunctionWaveSolver::blockLeastSquares
(
    scalarRectangularMatrix& J,
    scalarField& b,
    const label M,
    const scalar rcond
)
{
    const label nEq = J.m();
    const label n = J.n();
    const label nG = nEq - 2*M;

    // 1. A Givens rotation of rows k and M + k removes k eta_k from row
    //    M + k. Both rows have no other surface unknown.
    for (label k = 0; k < M; ++k)
    {
        const scalar a = J(k, k);
        const scalar c = J(M + k, k);
        const scalar r = std::hypot(a, c);
        if (r == 0)
        {
            continue;
        }
        const scalar cs = a/r;
        const scalar sn = c/r;

        auto rotate = [&](const label col)
        {
            const scalar x = J(k, col);
            const scalar y = J(M + k, col);
            J(k, col) = cs*x + sn*y;
            J(M + k, col) = -sn*x + cs*y;
        };
        rotate(k);
        for (label col = M; col < n; ++col)
        {
            rotate(col);
        }
        J(M + k, k) = 0;

        const scalar x = b[k];
        const scalar y = b[M + k];
        b[k] = cs*x + sn*y;
        b[M + k] = -sn*x + cs*y;
    }

    // 2. A Householder reflection of row k and the global rows removes
    //    k eta_k from the global rows. Row k fills in towards the later
    //    surface unknowns, which leaves it upper triangular.
    Foam::scalarField v(1 + nG);
    for (label k = 0; k < M; ++k)
    {
        v[0] = J(k, k);
        for (label g = 0; g < nG; ++g)
        {
            v[1 + g] = J(2*M + g, k);
        }
        scalar nrm = 0;
        for (const scalar vi : v)
        {
            nrm = std::hypot(nrm, vi);
        }
        if (nrm == 0)
        {
            continue;
        }
        const scalar alpha = (v[0] >= 0 ? -nrm : nrm);
        v[0] -= alpha;
        scalar vv = 0;
        for (const scalar vi : v)
        {
            vv += vi*vi;
        }
        if (vv > 0)
        {
            const scalar f = 2/vv;
            auto reflect = [&](auto get)
            {
                scalar s = v[0]*get(k);
                for (label g = 0; g < nG; ++g)
                {
                    s += v[1 + g]*get(2*M + g);
                }
                s *= f;
                get(k) -= s*v[0];
                for (label g = 0; g < nG; ++g)
                {
                    get(2*M + g) -= s*v[1 + g];
                }
            };
            for (label col = k + 1; col < n; ++col)
            {
                reflect([&](const label row) -> scalar& { return J(row, col); });
            }
            reflect([&](const label row) -> scalar& { return b[row]; });
        }
        J(k, k) = alpha;
        for (label g = 0; g < nG; ++g)
        {
            J(2*M + g, k) = 0;
        }
    }

    // 3. The reduced problem: rows M..2M-1 and the global rows, in the
    //    columns M..n-1. Stored with the columns as rows for the QR.
    const label nRed = M + nG;
    const label nCol = n - M;
    Foam::scalarRectangularMatrix AT(nCol, nRed);
    Foam::scalarField bRed(nRed);
    for (label i = 0; i < nRed; ++i)
    {
        const label row = M + i;    // M..2M-1, then 2M..
        for (label c = 0; c < nCol; ++c)
        {
            AT(c, i) = J(row, M + c);
        }
        bRed[i] = b[row];
    }
    const Foam::scalarField xRed
    (
        Foam::streamFunctionWaveSolver::leastSquares(AT, bRed, rcond)
    );

    // 4. Back substitution for the surface unknowns
    auto tx = Foam::tmp<Foam::scalarField>::New(n, Foam::Zero);
    auto& x = tx.ref();
    for (label c = 0; c < nCol; ++c)
    {
        x[M + c] = xRed[c];
    }
    for (label k = M - 1; k >= 0; --k)
    {
        scalar s = b[k];
        for (label col = k + 1; col < n; ++col)
        {
            s -= J(k, col)*x[col];
        }
        x[k] = s/J(k, k);
    }

    return tx;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::streamFunctionWaveSolver::streamFunctionWaveSolver
(
    const label N,
    const parameters& p
)
:
    N_(N),
    M_(2*N),
    p_(p),
    theta_(M_),
    w_(M_, 1.0/(M_ - 1)),
    C_(M_, N_),
    S_(M_, N_)
{
    w_[0] = w_[M_ - 1] = 0.5/(M_ - 1);

    for (label m = 0; m < M_; ++m)
    {
        theta_[m] = pi*m/(M_ - 1);
        for (label jj = 0; jj < N_; ++jj)
        {
            C_(m, jj) = std::cos((jj + 1)*theta_[m]);
            S_(m, jj) = std::sin((jj + 1)*theta_[m]);
        }
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::streamFunctionWaveSolver::basis
(
    const scalarField& z,
    scalarRectangularMatrix& Sh,
    scalarRectangularMatrix& Ch
) const
{
    const scalar kD = kd(z);
    Sh.resize(M_, N_);
    Ch.resize(M_, N_);
    for (label m = 0; m < M_; ++m)
    {
        for (label jj = 0; jj < N_; ++jj)
        {
            const label j = jj + 1;
            Sh(m, jj) = sinhByCosh(j*z[m], j*kD);
            Ch(m, jj) = coshByCosh(j*z[m], j*kD);
        }
    }
}


void Foam::streamFunctionWaveSolver::surfaceVelocity
(
    const scalarField& z,
    const UList<scalar>& theta,
    const UList<scalar>& ke,
    scalarField& U,
    scalarField& V
) const
{
    const scalar kD = kd(z);
    const scalar u = z[iu()];
    U.resize(theta.size());
    V.resize(theta.size());
    forAll(theta, i)
    {
        U[i] = -u;
        V[i] = 0;
        for (label jj = 0; jj < N_; ++jj)
        {
            const label j = jj + 1;
            const scalar Bj = z[M_ + jj];
            U[i] += j*Bj*coshByCosh(j*ke[i], j*kD)*std::cos(j*theta[i]);
            V[i] += j*Bj*sinhByCosh(j*ke[i], j*kD)*std::sin(j*theta[i]);
        }
    }
}


Foam::scalar Foam::streamFunctionWaveSolver::compare
(
    const scalarField& z,
    const scalarField& ke,
    const streamFunctionWaveSolver& s2,
    const scalarField& z2
) const
{
    const label M2 = s2.M();
    const scalarField ke2(SubField<scalar>(z2, M2));

    scalarField U, V, U2, V2;
    surfaceVelocity(z, s2.theta(), ke, U, V);
    s2.surfaceVelocity(z2, s2.theta(), ke2, U2, V2);

    // In depth units: elevations by d, velocities by sqrt(g d), q by
    // sqrt(g d^3), r by g d
    const scalar kd1 = kd(z);
    const scalar kd2 = s2.kd(z2);
    const scalar v1 = 1/std::sqrt(kd1);
    const scalar v2 = 1/std::sqrt(kd2);
    const scalar H = s2.kH(z2)/kd2;
    const scalar u = z2[s2.iu()]*v2;

    scalar diff = 0;
    for (label i = 0; i < M2; ++i)
    {
        diff = max(diff, mag(ke[i]/kd1 - ke2[i]/kd2)/H);
        diff = max(diff, mag(U[i]*v1 - U2[i]*v2)/u);
        diff = max(diff, mag(V[i]*v1 - V2[i]*v2)/u);
    }
    diff = max(diff, mag(z[iu()]*v1 - z2[s2.iu()]*v2)/u);
    const scalar q2 = z2[s2.iq()]/(kd2*std::sqrt(kd2));
    diff = max(diff, mag(z[iq()]/(kd1*std::sqrt(kd1)) - q2)/mag(q2));
    const scalar r2 = z2[s2.ir()]/kd2;
    diff = max(diff, mag(z[ir()]/kd1 - r2)/mag(r2));
    diff = max(diff, mag(kd1 - kd2)/kd2);

    return diff;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::streamFunctionWaveSolver::currentEquation
(
    const parameters& p,
    const scalar u,
    const scalar q,
    const scalar kd
)
{
    if (p.uType == 1)
    {
        // Mean current in the laboratory frame, (3.15) in Fenton 1999
        return std::sqrt(kd)*u + kd*p.uStar - 2*pi/p.tau;
    }

    // Mass transport velocity, (3.16) in Fenton 1999
    return q/std::sqrt(kd) + kd*p.uStar - 2*pi/p.tau;
}


Foam::tmp<Foam::scalarField> Foam::streamFunctionWaveSolver::residual
(
    const scalarField& z
) const
{
    const scalar kD = kd(z);
    const scalar u = z[iu()];
    const scalar q = z[iq()];
    const scalar r = z[ir()];

    scalarRectangularMatrix Sh, Ch;
    basis(z, Sh, Ch);

    auto tF = tmp<scalarField>::New(nEquations(), Zero);
    auto& F = tF.ref();

    scalar U0 = 0, UM = 0, mean = -kD;
    for (label m = 0; m < M_; ++m)
    {
        const scalar ke = z[m];

        // Kinematic and dynamic conditions, (3.7) and (3.8) in Fenton 1999
        scalar kin = -u*ke + q;
        scalar U = -u;
        scalar V = 0;
        for (label jj = 0; jj < N_; ++jj)
        {
            const label j = jj + 1;
            const scalar Bj = z[M_ + jj];
            kin += Bj*Sh(m, jj)*C_(m, jj);
            U += j*Bj*Ch(m, jj)*C_(m, jj);
            V += j*Bj*Sh(m, jj)*S_(m, jj);
        }
        F[m] = kin;
        F[M_ + m] = ke - r + 0.5*(U*U + V*V);
        mean += w_[m]*ke;

        if (m == 0) U0 = U;
        if (m == M_ - 1) UM = U;
    }

    // Mean depth, (3.9) in Fenton 1999
    F[2*M_] = W_*mean;

    // Wave height, (3.10) in Fenton 1999, or the crest parameter
    F[2*M_ + 1] =
        p_.heightCondition
      ? W_*(z[0] - z[M_ - 1] - kD*p_.h)
      : W_*(U0*U0 - (1 - p_.omega)*UM*UM);

    if (p_.freeKd)
    {
        F[2*M_ + 2] = W_*currentEquation(p_, u, q, kD);
    }

    return tF;
}


void Foam::streamFunctionWaveSolver::jacobian
(
    const scalarField& z,
    scalarRectangularMatrix& J
) const
{
    const scalar kD = kd(z);
    const scalar u = z[iu()];
    const scalar q = z[iq()];
    const label ikd = M_ + N_ + 3;

    scalarRectangularMatrix Sh, Ch;
    basis(z, Sh, Ch);

    J.resize(nEquations(), nUnknowns());
    J = Zero;

    // Derivatives of U at the crest and the trough, for the crest condition
    scalarField dU0(nUnknowns(), Zero), dUM(nUnknowns(), Zero);
    scalar U0 = 0, UM = 0;

    for (label m = 0; m < M_; ++m)
    {
        const scalar ke = z[m];
        scalar U = -u, V = 0;
        scalar dUdke = 0, dVdke = 0, dkindkd = 0, dUdkd = 0, dVdkd = 0;
        for (label jj = 0; jj < N_; ++jj)
        {
            const label j = jj + 1;
            const scalar Bj = z[M_ + jj];
            const scalar sc = Sh(m, jj)*C_(m, jj);
            const scalar cc = Ch(m, jj)*C_(m, jj);
            const scalar ss = Sh(m, jj)*S_(m, jj);
            const scalar cs = Ch(m, jj)*S_(m, jj);
            const scalar tj = j*std::tanh(j*kD);
            U += j*Bj*cc;
            V += j*Bj*ss;
            dUdke += j*j*Bj*sc;
            dVdke += j*j*Bj*cs;
            dkindkd -= Bj*sc*tj;
            dUdkd -= j*Bj*cc*tj;
            dVdkd -= j*Bj*ss*tj;
        }

        // Kinematic condition
        J(m, m) = U;
        for (label jj = 0; jj < N_; ++jj)
        {
            J(m, M_ + jj) = Sh(m, jj)*C_(m, jj);
        }
        J(m, iu()) = -ke;
        J(m, iq()) = 1;

        // Dynamic condition
        J(M_ + m, m) = 1 + U*dUdke + V*dVdke;
        for (label jj = 0; jj < N_; ++jj)
        {
            const label j = jj + 1;
            J(M_ + m, M_ + jj) =
                U*j*Ch(m, jj)*C_(m, jj) + V*j*Sh(m, jj)*S_(m, jj);
        }
        J(M_ + m, iu()) = -U;
        J(M_ + m, ir()) = -1;

        if (p_.freeKd)
        {
            J(m, ikd) = dkindkd;
            J(M_ + m, ikd) = U*dUdkd + V*dVdkd;
        }

        if (m == 0 || m == M_ - 1)
        {
            scalarField& g = (m == 0 ? dU0 : dUM);
            (m == 0 ? U0 : UM) = U;
            g[m] = dUdke;
            for (label jj = 0; jj < N_; ++jj)
            {
                g[M_ + jj] = (jj + 1)*Ch(m, jj)*C_(m, jj);
            }
            g[iu()] = -1;
            if (p_.freeKd)
            {
                g[ikd] = dUdkd;
            }
        }

        // Mean depth
        J(2*M_, m) = W_*w_[m];
    }

    if (p_.freeKd)
    {
        J(2*M_, ikd) = -W_;
    }

    // Height or crest condition
    if (p_.heightCondition)
    {
        J(2*M_ + 1, 0) = W_;
        J(2*M_ + 1, M_ - 1) = -W_;
        if (p_.freeKd)
        {
            J(2*M_ + 1, ikd) = -W_*p_.h;
        }
    }
    else
    {
        for (label c = 0; c < nUnknowns(); ++c)
        {
            J(2*M_ + 1, c) =
                W_*(2*U0*dU0[c] - (1 - p_.omega)*2*UM*dUM[c]);
        }
    }

    // Current
    if (p_.freeKd)
    {
        const label row = 2*M_ + 2;
        if (p_.uType == 1)
        {
            J(row, iu()) = W_*std::sqrt(kD);
            J(row, ikd) = W_*(0.5*u/std::sqrt(kD) + p_.uStar);
        }
        else
        {
            J(row, iq()) = W_/std::sqrt(kD);
            J(row, ikd) = W_*(-0.5*q/(kD*std::sqrt(kD)) + p_.uStar);
        }
    }
}


Foam::streamFunctionWaveSolver::result Foam::streamFunctionWaveSolver::solve
(
    scalarField& z,
    const scalar tol,
    label maxIter
) const
{
    const label n = nUnknowns();
    result res;
    scalarRectangularMatrix J;
    scalar last = GREAT;

    for (label it = 1; it <= maxIter; ++it)
    {
        res.iterations = it;

        const scalarField F(residual(z));
        jacobian(z, J);
        if (!allFinite(F) || !allFinite(J))
        {
            res.residual = GREAT;
            res.converged = false;
            return res;
        }

        // Columns scaled to their largest element
        scalarField D(n);
        for (label c = 0; c < n; ++c)
        {
            scalar mx = 0;
            for (label i = 0; i < J.m(); ++i)
            {
                mx = max(mx, mag(J(i, c)));
            }
            D[c] = (mx > 0 ? 1/mx : 1);
            for (label i = 0; i < J.m(); ++i)
            {
                J(i, c) *= D[c];
            }
        }

        scalarField b(-F);
        scalarField dz(blockLeastSquares(J, b, M_, 1e-13));
        dz *= D;
        if (!allFinite(dz))
        {
            res.residual = GREAT;
            res.converged = false;
            return res;
        }

        // Backtracking: halve the update until the residual does not grow
        // (much), so that a poor guess cannot run away
        const scalar f0 = max(mag(F));
        scalar lambda = 1;
        scalar f1 = GREAT;
        for (label i = 0; i < 30; ++i)
        {
            f1 = max(mag(residual(z + lambda*dz)));
            if (std::isfinite(f1) && f1 <= max(1.1*f0, 1e-13))
            {
                break;
            }
            lambda *= 0.5;
        }

        const scalarField z0(z);
        z += lambda*dz;

        // Size of the update in the measure of difference(): the change of
        // the surface and its velocities on the points, not of the
        // individual B_j, which are not unique with more modes than needed
        const scalar step =
            compare(z0, scalarField(SubField<scalar>(z0, M_)), *this, z);

        if (!std::isfinite(step))
        {
            res.residual = GREAT;
            res.converged = false;
            return res;
        }
        res.residual = f1;
        res.step = step;
        if (step < max(0.01*tol, 1e-14))
        {
            res.converged = true;
            return res;
        }
        // (A shortened update is a stall too once it is at round-off: the
        // residual is then noise, which the backtracking cannot reduce.)
        if (it >= 3 && step > 0.5*last && (lambda == 1 || step < 1e-8))
        {
            // Stalled at round-off. Near the highest wave this can be well
            // above machine precision: prescribing the height there leaves a
            // nearly undetermined direction (the height hardly changes with
            // the crest parameter), along which the update wanders by about
            // 1e-9. The solution is then known to about the update, which
            // callers add to their error estimate.
            res.converged = step < 1e-6;
            return res;
        }
        last = step;
    }

    res.converged = false;
    return res;
}


void Foam::streamFunctionWaveSolver::fitB
(
    scalarField& z,
    const scalar rcond
) const
{
    const scalar kD = kd(z);
    const scalar u = z[iu()];
    const scalar q = z[iq()];

    // sum_j B_j sinh(j k eta)/cosh(j k d) cos(j theta) = u k eta - q, with
    // the columns scaled to their largest element
    scalarRectangularMatrix AT(N_, M_);
    scalarField rhs(M_), D(N_);
    for (label jj = 0; jj < N_; ++jj)
    {
        const label j = jj + 1;
        scalar mx = 0;
        for (label m = 0; m < M_; ++m)
        {
            AT(jj, m) = sinhByCosh(j*z[m], j*kD)*C_(m, jj);
            mx = max(mx, mag(AT(jj, m)));
        }
        D[jj] = (mx > 0 ? 1/mx : 1);
        for (label m = 0; m < M_; ++m)
        {
            AT(jj, m) *= D[jj];
        }
    }
    for (label m = 0; m < M_; ++m)
    {
        rhs[m] = u*z[m] - q;
    }
    const scalarField B(D*leastSquares(AT, rhs, rcond));
    for (label jj = 0; jj < N_; ++jj)
    {
        z[M_ + jj] = B[jj];
    }
}


Foam::streamFunctionWaveSolver::modal Foam::streamFunctionWaveSolver::toModal
(
    const scalarField& z
) const
{
    // Cosine coefficients of k eta on the points: exact interpolation by the
    // discrete cosine transform (DCT-I)
    modal md;
    md.E.resize(M_);
    for (label m = 0; m < M_; ++m)
    {
        scalar s = 0;
        for (label i = 0; i < M_; ++i)
        {
            const scalar wi = (i == 0 || i == M_ - 1 ? 0.5 : 1.0);
            s += wi*z[i]*std::cos(m*theta_[i]);
        }
        md.E[m] = 2*s/(M_ - 1);
    }
    md.E[0] *= 0.5;
    md.E[M_ - 1] *= 0.5;

    md.B = SubField<scalar>(z, N_, M_);
    md.u = z[iu()];
    md.q = z[iq()];
    md.r = z[ir()];
    md.kd = kd(z);
    return md;
}


Foam::tmp<Foam::scalarField> Foam::streamFunctionWaveSolver::fromModal
(
    const modal& md
) const
{
    auto tz = tmp<scalarField>::New(nUnknowns(), Zero);
    auto& z = tz.ref();

    const label nE = min(md.E.size(), M_);
    for (label i = 0; i < M_; ++i)
    {
        scalar s = 0;
        for (label m = 0; m < nE; ++m)
        {
            s += md.E[m]*std::cos(m*theta_[i]);
        }
        z[i] = s;
    }
    const label nB = min(md.B.size(), N_);
    for (label jj = 0; jj < nB; ++jj)
    {
        z[M_ + jj] = md.B[jj];
    }
    z[iu()] = md.u;
    z[iq()] = md.q;
    z[ir()] = md.r;
    if (p_.freeKd)
    {
        z[M_ + N_ + 3] = md.kd;
    }
    return tz;
}


Foam::scalar Foam::streamFunctionWaveSolver::difference
(
    const scalarField& z,
    const streamFunctionWaveSolver& s2,
    const scalarField& z2
) const
{
    const scalarField ke(SubField<scalar>(s2.fromModal(toModal(z))(), s2.M()));
    return compare(z, ke, s2, z2);
}


bool Foam::streamFunctionWaveSolver::regular(const scalarField& z) const
{
    const scalar tol = 1e-3*mag(kH(z));
    for (label m = 0; m < M_ - 1; ++m)
    {
        if (z[m + 1] - z[m] > tol)
        {
            return false;
        }
    }
    return true;
}


Foam::scalar Foam::streamFunctionWaveSolver::crestParameter
(
    const scalarField& z
) const
{
    const scalarList theta({theta_[0], theta_[M_ - 1]});
    const scalarList ke({z[0], z[M_ - 1]});
    scalarField U, V;
    surfaceVelocity(z, theta, ke, U, V);
    return 1 - sqr(U[0])/sqr(U[1]);
}


Foam::tmp<Foam::scalarField> Foam::streamFunctionWaveSolver::leastSquares
(
    scalarRectangularMatrix& AT,
    scalarField& b,
    const scalar rcond
)
{
    const label n = AT.m();     // columns of A
    const label m = AT.n();     // rows of A
    const label sz = min(m, n);

    labelList p(identity(n));

    // Squared column norms, downdated, and recomputed when cancellation
    // makes the downdated value unreliable (as in LAPACK's xGEQP3)
    scalarField norms(n, Zero), normsRef(n);
    for (label c = 0; c < n; ++c)
    {
        for (label i = 0; i < m; ++i)
        {
            norms[c] += sqr(AT(c, i));
        }
    }
    normsRef = norms;

    scalarField diag(sz, Zero);
    label rank = sz;
    scalar R00 = 0;

    for (label k = 0; k < sz; ++k)
    {
        // Pivot: the remaining column of largest norm
        label j = k;
        for (label c = k + 1; c < n; ++c)
        {
            if (norms[c] > norms[j]) j = c;
        }
        if (j != k)
        {
            for (label i = 0; i < m; ++i)
            {
                std::swap(AT(k, i), AT(j, i));
            }
            std::swap(p[k], p[j]);
            std::swap(norms[k], norms[j]);
            std::swap(normsRef[k], normsRef[j]);
        }

        scalar nrm = 0;
        for (label i = k; i < m; ++i)
        {
            nrm = std::hypot(nrm, AT(k, i));
        }
        if (k == 0)
        {
            R00 = nrm;
        }
        if (nrm <= rcond*R00 || nrm == 0)
        {
            rank = k;
            break;
        }

        // Householder reflection I - 2 v v^T/(v^T v), v = x - alpha e_k
        const scalar alpha = (AT(k, k) >= 0 ? -nrm : nrm);
        AT(k, k) -= alpha;
        scalar vv = 0;
        for (label i = k; i < m; ++i)
        {
            vv += sqr(AT(k, i));
        }
        const scalar f = (vv > 0 ? 2/vv : 0);

        for (label c = k + 1; c < n; ++c)
        {
            scalar s = 0;
            for (label i = k; i < m; ++i)
            {
                s += AT(k, i)*AT(c, i);
            }
            s *= f;
            for (label i = k; i < m; ++i)
            {
                AT(c, i) -= s*AT(k, i);
            }

            // Downdate the norm of the rest of column c
            norms[c] -= sqr(AT(c, k));
            if (norms[c] <= 1e-6*normsRef[c])
            {
                norms[c] = 0;
                for (label i = k + 1; i < m; ++i)
                {
                    norms[c] += sqr(AT(c, i));
                }
                normsRef[c] = norms[c];
            }
        }
        {
            scalar s = 0;
            for (label i = k; i < m; ++i)
            {
                s += AT(k, i)*b[i];
            }
            s *= f;
            for (label i = k; i < m; ++i)
            {
                b[i] -= s*AT(k, i);
            }
        }
        diag[k] = alpha;
    }

    // The truncated R = [R11 R12] (rank x n, in pivoted order) as a dense
    // matrix. Element (i, c) of R is AT(c, i) for c > i.
    scalarRectangularMatrix R(rank, n, Zero);
    for (label i = 0; i < rank; ++i)
    {
        R(i, i) = diag[i];
        for (label c = i + 1; c < n; ++c)
        {
            R(i, c) = AT(c, i);
        }
    }

    // The solution of least norm (as an SVD would give), not the basic
    // solution with the truncated unknowns zero: with the basic solution,
    // the retained, nearly dependent columns take up round-off noise, which
    // near the highest wave stalls the iteration at about 1e-10. Householder
    // reflections H_i from the right turn [R11 R12] into [T 0] (as LAPACK's
    // xTZRZF), T upper triangular: [R11 R12] H_rank-1 ... H_0 = [T 0], so
    // x = H_rank-1 ... H_0 [T^-1 c, 0].
    const label nTail = n - rank;
    List<scalarField> vs(rank);
    scalarList vvs(rank, Zero);
    if (nTail > 0)
    {
        for (label i = rank - 1; i >= 0; --i)
        {
            scalarField& v = vs[i];
            v.resize(1 + nTail);
            v[0] = R(i, i);
            scalar nrm = mag(v[0]);
            for (label t = 0; t < nTail; ++t)
            {
                v[1 + t] = R(i, rank + t);
                nrm = std::hypot(nrm, v[1 + t]);
            }
            const scalar alpha = (v[0] >= 0 ? -nrm : nrm);
            v[0] -= alpha;
            scalar vv = 0;
            for (const scalar vi : v)
            {
                vv += vi*vi;
            }
            vvs[i] = vv;
            if (vv > 0)
            {
                const scalar f = 2/vv;
                for (label k = 0; k < i; ++k)
                {
                    scalar s = v[0]*R(k, i);
                    for (label t = 0; t < nTail; ++t)
                    {
                        s += v[1 + t]*R(k, rank + t);
                    }
                    s *= f;
                    R(k, i) -= s*v[0];
                    for (label t = 0; t < nTail; ++t)
                    {
                        R(k, rank + t) -= s*v[1 + t];
                    }
                }
            }
            R(i, i) = alpha;
            for (label t = 0; t < nTail; ++t)
            {
                R(i, rank + t) = 0;
            }
        }
    }

    // T w = c
    scalarField y(n, Zero);
    for (label i = rank - 1; i >= 0; --i)
    {
        scalar s = b[i];
        for (label c = i + 1; c < rank; ++c)
        {
            s -= R(i, c)*y[c];
        }
        y[i] = s/R(i, i);
    }

    // y = H_rank-1 ... H_0 [w, 0]: H_0 first
    if (nTail > 0)
    {
        for (label i = 0; i < rank; ++i)
        {
            if (vvs[i] == 0) continue;
            const scalarField& v = vs[i];
            scalar s = v[0]*y[i];
            for (label t = 0; t < nTail; ++t)
            {
                s += v[1 + t]*y[rank + t];
            }
            s *= 2/vvs[i];
            y[i] -= s*v[0];
            for (label t = 0; t < nTail; ++t)
            {
                y[rank + t] -= s*v[1 + t];
            }
        }
    }

    auto tx = tmp<scalarField>::New(n, Zero);
    for (label i = 0; i < n; ++i)
    {
        tx.ref()[p[i]] = y[i];
    }
    return tx;
}


Foam::tmp<Foam::scalarField> Foam::streamFunctionWaveSolver::linearWave
(
    const scalar kD,
    const scalar h
) const
{
    auto tz = tmp<scalarField>::New(nUnknowns(), Zero);
    auto& z = tz.ref();
    const scalar u = std::sqrt(std::tanh(kD));
    for (label m = 0; m < M_; ++m)
    {
        z[m] = kD + 0.5*kD*h*std::cos(theta_[m]);
    }
    z[M_] = 0.5*kD*h/u;
    z[iu()] = u;
    z[iq()] = kD*u;
    z[ir()] = kD + 0.5*u*u;
    if (p_.freeKd)
    {
        z[M_ + N_ + 3] = kD;
    }
    return tz;
}


// ************************************************************************* //
