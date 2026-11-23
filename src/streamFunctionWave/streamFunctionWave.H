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

Class
    Foam::streamFunctionWave

Description
    Calculates Fourier coefficients for a stream function theory wave.

    The calculation follows Fenton's Fourier approximation method, section
    3.1.2 of:

    \verbatim
        Fenton, J. D. (1999).
        Numerical methods for nonlinear waves.
        Advances in Coastal and Ocean Engineering, 5, 241-324.
    \endverbatim

    Equation numbering in the source files refers to this paper.

    The wave is specified by its height H, the mean water depth d and one of
    the following, selected with uType:
      - 0: the wave length lambda. The wave is given in the frame moving with
           the wave, so it is stationary.
      - 1: the wave period T and the Eulerian mean current u (u1 in Fenton).
      - 2: the wave period T and the mass transport velocity u (u2 in
           Fenton).
      - 3: the wave length lambda and the wave period T. The wave is
           calculated as for uType 0 and then given in the frame in which it
           moves with celerity lambda/T. The Eulerian mean current and the
           mass transport velocity follow from this.

    The number of Fourier modes follows from the tolerance (default 1e-8):
    the largest error allowed in the surface elevation relative to H, and in
    the velocities at the surface and the mean fluid speed relative to the
    mean fluid speed. The wave is accepted when the solution agrees to half
    the tolerance with the solution of about 25 % more modes.

    The equations are solved in the least-squares sense on more surface
    points than modes (streamFunctionWaveSolver), from a guess interpolated
    in a table of precomputed waves (streamFunctionWaveTable) over the
    relative depth d/lambda and the relative height H/H_top. H_top is the
    highest wave of the table, about 1 % below the highest possible wave in
    deep water and 2-9 % below it in shallow water; higher waves are an
    error. The table covers d/lambda from 0.03 to 3; deeper water is as
    deep water, shallower water is tried with a warning.

SourceFiles
    streamFunctionWave.C

Author
    Johan Roenby, DHI and STROMNING

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
        label N_ = 0;

        //- Tolerance of the wave
        scalar tolerance_ = 1e-8;

        //- Estimated error of the calculated wave
        scalar error_ = 0;

        //- Relative height H/H_top of the calculated wave
        scalar sigma_ = 0;

        //- Fourier coefficients for velocity field and pressure
        scalarField B_;

        //- Fourier coefficients for surface elevation
        scalarField E_;

        //- Wave height
        scalar H_ = 0;

        //- Mean water level
        scalar d_ = 0;

        //- Input type: 0 for wave length input, 1 for wave period and
        //- "Eulerian" current, 2 for wave period and "Stokes" current,
        //- 3 for wave length and wave period
        label uType_ = 0;

        //- Gravitational acceleration
        scalar g_ = 0;

        //- Current - either u1 or u2 - type determined by uType_
        scalar u_ = 0;

        //- Wave number, k = 2*pi/lambda
        scalar k_ = 0;

        //- Angular wave frequency
        scalar omega_ = 0;

        //- Mean fluid speed in the co-moving frame
        scalar uBar_ = 0;

        //- Bernoulli constant
        scalar R_ = 0;

        //- Volume flow rate per unit span underneath the waves in the
        //- co-moving frame
        scalar Q_ = 0;


    // Protected Member Functions

        //- True if the wave length is given, so the wave period does not
        //- enter the calculation of the Fourier coefficients
        bool lengthInput() const noexcept
        {
            return uType_ == 0 || uType_ == 3;
        }


public:

    // Constructors

        //- Default construct
        streamFunctionWave() = default;

        //- Construct from dictionary
        explicit streamFunctionWave(const dictionary& waveDict);


    // Member Functions

        //- Read the wave parameters from the given dictionary
        void readParmsFromDict(const dictionary& waveDict);

        //- True if the input parameters of this wave and of w are equal.
        //- The comparison is exact, for parameters read from the same text.
        bool sameInput(const streamFunctionWave& w) const;

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

        // Access

            //- Gravitational acceleration
            scalar g() const noexcept { return g_; }

            //- Wave height
            scalar H() const noexcept { return H_; }

            //- Wave number, calculated by calcCoeffs() for period input
            scalar k() const noexcept { return k_; }

            //- Mean fluid speed in the frame of the wave
            scalar uBar() const noexcept { return uBar_; }

            //- Volume flow rate per unit span in the frame of the wave
            scalar Q() const noexcept { return Q_; }

            //- Number of Fourier modes of the calculated wave
            label N() const noexcept { return N_; }

            //- Estimated error of the calculated wave
            scalar error() const noexcept { return error_; }

        //- Find the x/lambda value where eta = d
        scalar zeroPointPhase() const;

        //- Write the calculated wave, with its Fourier coefficients, to the
        //- named dictionary file
        void writeDict(const fileName& outDictName) const;

        //- Write the input and the main results of the calculated wave as
        //- dictionary entries, for information
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
