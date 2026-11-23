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
    Foam::streamFunctionWaveTable

Description
    Table of precomputed stream function waves over the relative depth
    d/lambda and the relative height sigma = H/H_top(d/lambda), from which
    a requested wave is started.

    H_top(d/lambda) is the highest wave that the table generator
    (generateStreamFunctionWaveTable) could compute to its tolerance with
    its largest number of modes, a little below the highest possible wave.
    Each row of the table is a relative depth, with nodes at fixed sigma,
    uniform in t = 1 - sqrt(1 - sigma) (dense towards the top). Each node
    holds the surface as cosine coefficients of k eta, u, q, r and the
    number of modes needed for each of a set of tolerances.

    A guess for a wave is interpolated, cubically in log(d/lambda) and t,
    from the surfaces, u, q - u k d and r - k d of the nodes; the B_j follow
    from the kinematic condition, which is linear in them. Interpolating the
    B_j themselves fails for steep waves, where they decay like exp(-a j)
    with a that changes with the height.

    The table is compiled into the library (streamFunctionWaveTableData.C
    includes table/streamFunctionWaveTable.dat). Only an index of it is
    built, when the first wave is calculated; the node data are read in
    place, so only the pages of the nodes used are ever loaded.

    Layout of the data, a flat list of numbers:
    \verbatim
        version (1), nLevels, levels..., generator tolerance, generator
        largest N, nRows,
        per row:  d/lambda, k d, k H_top, nNodes,
        per node: sigma, k H, N per level (-1: not reached), u, q, r, nE,
                  E_0..E_nE-1
    \endverbatim
    with E_0 = k d (k eta = sum_m E_m cos(m theta)).

SourceFiles
    streamFunctionWaveTable.C
    streamFunctionWaveTableData.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_streamFunctionWaveTable_H
#define Foam_streamFunctionWaveTable_H

#include "streamFunctionWaveSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Ostream;

/*---------------------------------------------------------------------------*\
                   Class streamFunctionWaveTable Declaration
\*---------------------------------------------------------------------------*/

class streamFunctionWaveTable
{
public:

    // Public Data Types

        //- A node, as written by the generator
        struct nodeData
        {
            scalar sigma = 0;
            scalar kH = 0;
            labelList N;
            streamFunctionWaveSolver::modal md;
        };

        //- A row, as written by the generator
        struct rowData
        {
            scalar dl = 0;
            scalar kd = 0;
            scalar kHtop = 0;
            List<nodeData> nodes;
        };


private:

    // Private Data

        //- The data
        const double* data_;

        //- Tolerance levels for which the nodes give N
        scalarList levels_;

        //- Tolerance and largest number of modes of the generator
        scalar genTol_;
        label genNmax_;

        //- Per row: log(d/lambda), k d and H_top/lambda
        scalarList x_;
        scalarList kd_;
        scalarList HtopL_;

        //- Coordinate t of the nodes, with t = 0 (still water) first
        scalarList t_;

        //- Offset of each node in the data, per row
        List<labelList> node_;


    // Private Member Functions

        //- Lagrange weights of the 4 (or fewer) nodes of xs around x
        static void weights
        (
            const UList<scalar>& xs,
            const scalar x,
            labelList& idx,
            scalarList& w
        );

        //- Node k of row i (k = 0 is still water): u, q - u k d, r - k d and
        //- the surface relative to k d at the points theta
        void node
        (
            const label i,
            const label k,
            const UList<scalar>& theta,
            scalarField& ke,
            scalar& u,
            scalar& qr,
            scalar& rr
        ) const;

        //- Modes node k of row i needs for level l, or -1
        label nodeN(const label i, const label k, const label l) const;

        //- Interpolation coordinates of a wave, clipped to the table
        void coords
        (
            const scalar dl,
            const scalar sigma,
            scalar& x,
            scalar& t
        ) const;


public:

    // Static Data

        //- Tolerance level of the guess, which sets the first number of
        //- modes of a request
        static constexpr scalar guessLevel = 1e-6;


    // Constructors

        //- Construct the index of the given data
        streamFunctionWaveTable(const double* data, const label size);

        //- The table compiled into the library
        static const streamFunctionWaveTable& builtIn();


    // Member Functions

        //- True if the table has no rows
        bool empty() const noexcept { return x_.empty(); }

        //- Range of relative depths
        scalar dlMin() const { return std::exp(x_.first()); }
        scalar dlMax() const { return std::exp(x_.last()); }

        //- Largest number of modes of the generator
        label genNmax() const noexcept { return genNmax_; }

        //- H_top/lambda at relative depth d/lambda (the deepest row's value
        //- in deeper water)
        scalar HtopByLambda(const scalar dl) const;

        //- sigma = H/H_top of the wave of relative depth d/lambda and
        //- H/d = h. Below the table it is H/d against H_top/d of the
        //- shallowest row.
        scalar sigma(const scalar dl, const scalar h) const;

        //- Number of modes that suits the accuracy of the guess
        label guessN(const scalar dl, const scalar sigma) const;

        //- The guess for problem s, of relative depth d/lambda (k d =
        //- 2 pi d/lambda) and relative height sigma
        tmp<scalarField> guess
        (
            const scalar dl,
            const scalar sigma,
            const streamFunctionWaveSolver& s
        ) const;

        //- u and q of the table wave of depth k d and H/d = h
        void uq(const scalar kd, const scalar h, scalar& u, scalar& q) const;

        //- k d of the table wave with the period and current of p. Against
        //- a current there may also be a short wave of the same period;
        //- the longest is returned and kdOther set to the next one (or -1).
        //- Returns -1 if there is none.
        scalar kdForPeriod
        (
            const streamFunctionWaveSolver::parameters& p,
            scalar& kdOther
        ) const;

        //- Write a table in the layout above, as a C initialiser list
        static void write
        (
            Ostream& os,
            const UList<scalar>& levels,
            const scalar genTol,
            const label genNmax,
            const UList<rowData>& rows
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
