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

Application
    calcStreamFunctionWave

Group
    grpPreProcessingUtilities

Description
    Calculates Fourier coefficients for a stream function theory wave.

    The calculation follows the outline in section 3.1.2 of:

    \verbatim
        Fenton, J. D. (1999).
        Numerical methods for nonlinear waves.
        Advances in Coastal and Ocean Engineering, 5, 241-324.
    \endverbatim

    Equation numbering in the source files refers to this paper.

Usage
    \b calcStreamFunctionWave \<inputDict\> \<outputDict\>

    where \<inputDict\> holds the parameters needed to calculate the Fourier
    coefficients - see the waveProperties file next to this source for an
    example. If \<outputDict\> already describes the same wave, it is left
    alone and nothing is recalculated.

Author
    Johan Roenby, DHI and STROMNING

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "dictionary.H"
#include "streamFunctionWave.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Read the specified stream function wave input dictionary, check"
        " whether it describes the same wave as the output dictionary"
        " (if that exists) and if not, calculate the wave Fourier"
        " coefficients and update the output dictionary."
    );

    argList::noBanner();
    argList::noParallel();
    argList::addArgument("inputDict");
    argList::addArgument("outputDict");

    argList args(argc, argv);

    const fileName inDictName(args[1]);
    const fileName outDictName(args[2]);

    const dictionary inDict(IFstream(inDictName)(), true);
    streamFunctionWave inWave(inDict);

    bool calcWave(true);

    IFstream fileStream(outDictName);
    if (fileStream.good())
    {
        const dictionary oldOutDict(fileStream(), true);

        if (inWave.sameWave(oldOutDict))
        {
            calcWave = false;

            Info<< "Same wave in " << inDictName << " and " << outDictName << nl
                << "Not recalculating" << endl;
        }
    }

    if (calcWave)
    {
        Info<< "Calculating stream function wave..." << endl;
        inWave.calcCoeffs();

        Info<< "Writing wave to file " << outDictName << endl;
        inWave.writeDict(outDictName);
    }

    return 0;
}


// ************************************************************************* //
