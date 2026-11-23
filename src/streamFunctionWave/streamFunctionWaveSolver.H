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

Class
    Foam::streamFunctionWaveSolver

Description
    Least-squares version of Fenton's Fourier approximation method for a
    stream function wave with N Fourier modes (Fenton 1999, section 3.1.2).

    The surface is collocated at M = 2N points over half a wave length,
    theta_m = pi m/(M - 1), so there are more equations (2M + 2, or 2M + 3
    with period input) than unknowns (M + N + 3, or M + N + 4). They are
    solved in the least-squares sense by Gauss-Newton iteration. Unlike
    Fenton's square system, adding modes then does not spoil the solution,
    and the method reaches machine precision also for steep waves.

    Everything is non-dimensionalised with the wave number k and g as in
    Fenton (1999). The unknowns are

        z = [k eta_0..M-1, B_1..N, u, q, r(, k d)]

    with u = ubar sqrt(k/g), q = Q sqrt(k^3/g) and r = R k/g. k d is an
    unknown for period input (uType 1 or 2) and fixed for wave length input.
    The size of the wave is fixed by its height H/d, or by the crest
    parameter omega = 1 - U_crest^2/U_trough^2 (velocities in the frame of
    the wave), which is used to continue towards the highest wave.

    The linear least-squares problem of each iteration exploits that k eta_m
    appears only in the two surface conditions at point m and in the few
    global equations: these unknowns are eliminated exactly by orthogonal
    transformations, and the remaining problem in B, u, q, r(, k d) is
    solved by a QR decomposition with column pivoting, truncated where the
    diagonal falls below rcond times its first element.

SourceFiles
    streamFunctionWaveSolver.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_streamFunctionWaveSolver_H
#define Foam_streamFunctionWaveSolver_H

#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class streamFunctionWaveSolver Declaration
\*---------------------------------------------------------------------------*/

class streamFunctionWaveSolver
{
public:

    // Public Data Types

        //- The parameters of the wave, apart from the number of modes
        struct parameters
        {
            //- k d fixed (wave length input) or an unknown (period input)
            bool freeKd = false;

            //- k d for wave length input, else not used
            scalar kd = 0;

            //- The wave height H/d fixes the size of the wave, or else the
            //- crest parameter omega
            bool heightCondition = true;

            //- Wave height divided by the depth, H/d
            scalar h = 0;

            //- Crest parameter omega = 1 - U_crest^2/U_trough^2
            scalar omega = 0;

            //- 1: Eulerian mean current, 2: mass transport velocity
            label uType = 2;

            //- Wave period, T sqrt(g/d)
            scalar tau = 0;

            //- Current, u/sqrt(g d)
            scalar uStar = 0;
        };

        //- A solution in a form independent of N: the cosine coefficients
        //- of k eta (k eta = sum_m E_m cos(m theta), E_0 = k d), B and u,
        //- q, r, k d
        struct modal
        {
            scalarField E;
            scalarField B;
            scalar u = 0;
            scalar q = 0;
            scalar r = 0;
            scalar kd = 0;
        };

        //- Outcome of solve()
        struct result
        {
            label iterations = 0;
            scalar residual = 0;
            bool converged = false;

            //- Size of the last update (in the measure of difference()),
            //- the uncertainty of a solution whose iteration stalled
            scalar step = 0;
        };


private:

    // Private Data

        //- Number of Fourier modes
        label N_;

        //- Number of collocation points over half a wave length
        label M_;

        //- Wave parameters
        parameters p_;

        //- Weight of the global equations (mean depth, height, current)
        static constexpr scalar W_ = 10;

        //- Collocation points theta_m and trapezoidal weights of the mean
        scalarField theta_;
        scalarField w_;

        //- cos(j theta_m) and sin(j theta_m), M x N
        scalarRectangularMatrix C_;
        scalarRectangularMatrix S_;


    // Private Member Functions

        //- k d of solution z
        scalar kd(const scalarField& z) const
        {
            return p_.freeKd ? z[M_ + N_ + 3] : p_.kd;
        }

        //- sinh(j k eta)/cosh(j k d) and cosh(j k eta)/cosh(j k d) at the
        //- surface points (M x N)
        void basis
        (
            const scalarField& z,
            scalarRectangularMatrix& Sh,
            scalarRectangularMatrix& Ch
        ) const;

        //- Horizontal and vertical velocity in the frame of the wave at
        //- the points theta with surface elevation ke, of solution z
        void surfaceVelocity
        (
            const scalarField& z,
            const UList<scalar>& theta,
            const UList<scalar>& ke,
            scalarField& U,
            scalarField& V
        ) const;

        //- Largest difference between solution z of this problem and z2
        //- of problem s2, at the points of s2 where this surface is ke
        scalar compare
        (
            const scalarField& z,
            const scalarField& ke,
            const streamFunctionWaveSolver& s2,
            const scalarField& z2
        ) const;


public:

    // Constructors

        //- Construct for N modes
        streamFunctionWaveSolver(const label N, const parameters& p);


    // Member Functions

        // Access

            label N() const noexcept { return N_; }
            label M() const noexcept { return M_; }
            const parameters& params() const noexcept { return p_; }
            const scalarField& theta() const noexcept { return theta_; }

            //- Number of unknowns and of equations
            label nUnknowns() const noexcept
            {
                return M_ + N_ + 3 + (p_.freeKd ? 1 : 0);
            }
            label nEquations() const noexcept
            {
                return 2*M_ + 2 + (p_.freeKd ? 1 : 0);
            }

            //- Positions of u, q, r and k d in z
            label iu() const noexcept { return M_ + N_; }
            label iq() const noexcept { return M_ + N_ + 1; }
            label ir() const noexcept { return M_ + N_ + 2; }

            //- k d of solution z
            scalar kdOf(const scalarField& z) const { return kd(z); }

            //- k H of solution z
            scalar kH(const scalarField& z) const
            {
                return z[0] - z[M_ - 1];
            }


        // Equations

            //- Residual of the equations
            tmp<scalarField> residual(const scalarField& z) const;

            //- Analytic Jacobian of residual()
            void jacobian(const scalarField& z, scalarRectangularMatrix& J)
                const;

            //- Last equation for period input, without the weight
            static scalar currentEquation
            (
                const parameters& p,
                const scalar u,
                const scalar q,
                const scalar kd
            );


        // Solution

            //- Gauss-Newton iteration from z until the update of the surface
            //- and its velocities is below 0.01 tol, or stalls at round-off
            result solve(scalarField& z, const scalar tol, label maxIter = 40)
                const;

            //- Replace the B_j of z by the least-squares fit of the
            //- kinematic condition to its surface, u and q (linear in the
            //- B_j), truncated at rcond: directions that the surface does
            //- not determine to rcond are left empty. Directions below the
            //- truncation of solve() (1e-13) are not corrected by it, so a
            //- start that carries noise there can stall it above the
            //- truncation error; this gives a clean start.
            void fitB(scalarField& z, const scalar rcond = 1e-8) const;

            //- Modal form of solution z
            modal toModal(const scalarField& z) const;

            //- Solution vector of this problem from a modal form, truncated
            //- or padded with zeros
            tmp<scalarField> fromModal(const modal& m) const;

            //- Largest difference between solution z and solution z2 of
            //- problem s2, on the points of s2: surface elevation relative to
            //- H, surface velocities and mean speed relative to the mean
            //- speed, q, r and k d relative, all in depth units
            scalar difference
            (
                const scalarField& z,
                const streamFunctionWaveSolver& s2,
                const scalarField& z2
            ) const;

            //- True if the surface falls from crest to trough (one crest per
            //- wave length; trough wiggles well below H allowed)
            bool regular(const scalarField& z) const;

            //- Crest parameter omega of solution z
            scalar crestParameter(const scalarField& z) const;


        // Linear algebra

            //- Least-squares solution of A x = b by a QR decomposition with
            //- column pivoting, truncated where the column norm falls below
            //- rcond times the first, and of least norm (a complete
            //- orthogonal decomposition, as LAPACK's xGELSY). A is m x n,
            //- with the columns of A as the rows of AT, which is destroyed,
            //- as is b.
            static tmp<scalarField> leastSquares
            (
                scalarRectangularMatrix& AT,
                scalarField& b,
                const scalar rcond
            );


            //- Least-squares solution of the linearised equations J x = b,
            //- exploiting that k eta_m appears only in rows m, M + m and the
            //- global rows: these unknowns are eliminated exactly by
            //- orthogonal transformations, the rest by leastSquares(). J and
            //- b are destroyed.
            static tmp<scalarField> blockLeastSquares
            (
                scalarRectangularMatrix& J,
                scalarField& b,
                const label M,
                const scalar rcond
            );


        // Helpers

            //- Number of modes in the next level of refinement
            static label grow(const label N)
            {
                return label(1.25*N) + 4;
            }

            //- Linear (small amplitude) wave of height H/d = h and depth kd
            tmp<scalarField> linearWave(const scalar kd, const scalar h)
                const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
