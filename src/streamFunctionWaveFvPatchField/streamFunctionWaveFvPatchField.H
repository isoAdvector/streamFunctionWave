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
    Foam::streamFunctionWaveFvPatchField

Description
    Wave generating boundary condition based on the Fourier approximation
    method for stream function waves described in section 3.1.2 of:

    \verbatim
        Fenton, J. D. (1999).
        Numerical methods for nonlinear waves.
        Advances in Coastal and Ocean Engineering, 5, 241-324.
    \endverbatim

    Equation numbers in the code refer to this paper. At the time of writing
    it could be downloaded from

    http://www.johndfenton.com/Papers/Fenton99Liu-Numerical-methods-for-nonlinear-waves.pdf

    The patch works with interFoam type solvers.

    The wave is calculated from the wave dictionary when the patch is first
    constructed. The patch then writes the wave parameters and Fourier
    coefficients to its entry in a \c wave sub-dictionary, together with the
    wave direction, phase and origin. When the patch is constructed from an
    entry with a \c wave sub-dictionary, e.g. by decomposePar, a restart or
    reconstructPar, the wave is read from it without being recalculated and
    the wave dictionary is not read. Remove the \c wave sub-dictionary (e.g.
    by restoring 0/ from 0.orig/) to recalculate the wave after changing the
    wave dictionary.

Usage
    The wave generator patch velocity is specified in 0/U as

    \verbatim
    <patchName>
    {
        type            streamFunctionWave;
        value           uniform (0 0 0);
        waveDictName    waveProperties;     // default: waveProperties
        waveDictPath    "constant";         // default: "constant"
    }
    \endverbatim

    and similarly for the phase fraction field in 0/alpha.water

    \verbatim
    <patchName>
    {
        type            streamFunctionWave;
        value           uniform 0;
    }
    \endverbatim

    The file <case>/system/controlDict must contain

    \verbatim
    libs            (streamFunctionWaveFvPatchField);
    \endverbatim

SourceFiles
    streamFunctionWaveFvPatchField.C

Author
    Johan Roenby, DHI

\*---------------------------------------------------------------------------*/

#ifndef Foam_streamFunctionWaveFvPatchField_H
#define Foam_streamFunctionWaveFvPatchField_H

#include "mixedFvPatchFields.H"
#include "uniformDimensionedFields.H"
#include "streamFunctionWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class streamFunctionWaveFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class streamFunctionWaveFvPatchField
:
    public mixedFvPatchField<Type>
{
protected:

    // Protected Data

        //- The stream function wave
        streamFunctionWave wave_;

        //- Wave vector
        vector K_;

        //- Unit vector pointing upwards
        vector up_;

        //- Wave phase in radians
        scalar phi_;

        //- Vertical coordinate of the seabed at the patch
        scalar seabedHeight_;

        //- Horizontal position of the wave maker patch
        scalar xWaveMaker_;

        //- Dictionary from which to read the input data
        word waveDictName_;

        //- Path to the dictionary from which to read the input data
        fileName waveDictPath_;

        //- Position of the lowest point of each face of the patch
        scalarField ymin_;

        //- Position of the highest point of each face of the patch
        scalarField ymax_;


public:

    //- Runtime type information
    TypeName("streamFunctionWave");


    // Constructors

        //- Construct from patch and internal field
        streamFunctionWaveFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        streamFunctionWaveFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given streamFunctionWaveFvPatchField
        //- onto a new patch
        streamFunctionWaveFvPatchField
        (
            const streamFunctionWaveFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        streamFunctionWaveFvPatchField
        (
            const streamFunctionWaveFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return fvPatchField<Type>::Clone(*this);
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return fvPatchField<Type>::Clone(*this, iF);
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Set the patch field values, dispatched on the field Type
        virtual void setField
        (
            const Field<Type>&,
            const scalarField& alpha,
            const scalarField& u,
            const scalarField& v
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "streamFunctionWaveFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
