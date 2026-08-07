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

Class
    Foam::streamFunctionWave

Description
    Calculates Fourier coefficients for a stream function theory wave.

    The calculation follows the outline in section 3.1.2 of:

    \verbatim
        Fenton, J. D. (1999).
        Numerical methods for nonlinear waves.
        Advances in Coastal and Ocean Engineering, 5, 241-324.
    \endverbatim

    Equation numbering in the source files refers to this paper.

SourceFiles
    streamFunctionWave.C

Author
    Johan Roenby, DHI

\*---------------------------------------------------------------------------*/

#ifndef Foam_streamFunctionWave_H
#define Foam_streamFunctionWave_H

#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class streamFunctionWave Declaration
\*---------------------------------------------------------------------------*/

class streamFunctionWave
{
protected:

    // Protected Data

        //- Number of Fourier modes
        label N_;

        //- Fourier coefficients for velocity field and pressure
        scalarField B_;

        //- Fourier coefficients for surface elevation
        scalarField E_;

        //- Wave height
        scalar H_;

        //- Mean water level
        scalar d_;

        //- Current type: 0 for wave length input, 1 for wave period and
        //- "Eulerian" current, 2 for wave period and "Stokes" current
        label uType_;

        //- Number of steps to reach the final wave height
        label nHsteps_;

        //- Gravitational acceleration
        scalar g_;

        //- Current - either u1 or u2 - type determined by uType_
        scalar u_;

        //- Wave number, k = 2*pi/lambda
        scalar k_;

        //- Angular wave frequency
        scalar omega_;

        //- Mean fluid speed in the co-moving frame
        scalar uBar_;

        //- Bernoulli constant
        scalar R_;

        //- Volume flow rate per unit span underneath the waves in the
        //- co-moving frame
        scalar Q_;


    // Protected Member Functions

        //- Newton-Raphson iteration towards the wave of height h
        void findRoot(scalarField& X, scalar h, scalar lort);

        //- The set of equations to be satisfied by the wave
        scalarField fun(const scalarField& X, scalar h, scalar lort) const;

        //- Numerical Jacobian of fun()
        scalarSquareMatrix jac(const scalarField& X, scalar h, scalar lort) const;


public:

    // Constructors

        //- Default construct
        streamFunctionWave();

        //- Construct from dictionary
        explicit streamFunctionWave(const dictionary& waveDict);


    // Member Functions

        //- Read the wave parameters from the given dictionary
        void readParmsFromDict(const dictionary& waveDict);

        //- True if the given dictionary describes the same wave as this one
        bool sameWave(const dictionary& waveDict) const;

        //- Calculate the Fourier coefficients of the wave
        void calcCoeffs();

        //- Surface elevation above the seabed
        scalarField getEta
        (
            const scalarField& x,
            const scalar t,
            const scalar phi
        ) const;

        //- Velocity component along the direction of wave propagation
        scalarField getU
        (
            const scalarField& x,
            const scalarField& y,
            const scalar t,
            const scalar phi
        ) const;

        //- Vertical velocity component
        scalarField getV
        (
            const scalarField& x,
            const scalarField& y,
            const scalar t,
            const scalar phi
        ) const;

        //- Pressure divided by density
        scalarField getPbyRho
        (
            const scalarField& x,
            const scalarField& y,
            const scalar t,
            const scalar phi
        ) const;

        //- Find the x/lambda value where eta = d
        scalar zeroPointPhase() const;

        //- Write the calculated wave to the named dictionary file
        void writeDict(const fileName& outDictName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
