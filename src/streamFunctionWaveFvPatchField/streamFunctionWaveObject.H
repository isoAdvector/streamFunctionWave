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
    Foam::streamFunctionWaveObject

Description
    A calculated stream function wave held by an object registry (the mesh),
    so that the fields of a patch, and the fields read again at other times
    by post-processing utilities, share one calculation of the wave, as
    OpenFOAM's waveModel::lookupOrCreate does for its wave models.

    The wave is recalculated if the input differs from that of the stored
    wave, e.g. after the wave dictionary was edited. It is not written.

SourceFiles
    streamFunctionWaveObject.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_streamFunctionWaveObject_H
#define Foam_streamFunctionWaveObject_H

#include "regIOobject.H"
#include "streamFunctionWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class streamFunctionWaveObject Declaration
\*---------------------------------------------------------------------------*/

class streamFunctionWaveObject
:
    public regIOobject
{
    // Private Data

        //- The input, before the calculation
        streamFunctionWave input_;

        //- The calculated wave
        streamFunctionWave wave_;


public:

    //- Runtime type information
    TypeName("streamFunctionWaveObject");


    // Constructors

        //- Construct from the input, calculating the wave
        streamFunctionWaveObject
        (
            const IOobject& io,
            const streamFunctionWave& input
        );


    // Member Functions

        //- The wave of the given input, from the registry if it holds the
        //- wave of this input under the given name, else calculated and
        //- stored there
        static const streamFunctionWave& lookupOrCreate
        (
            const objectRegistry& db,
            const word& name,
            const streamFunctionWave& input
        );

        //- The calculated wave
        const streamFunctionWave& wave() const noexcept { return wave_; }

        //- Not written
        virtual bool writeData(Ostream&) const { return true; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
