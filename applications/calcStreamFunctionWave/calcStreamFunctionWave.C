/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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
    calcStreamFuncionWave

Description
    Calculates Fourier coefficients for a stream function theory wave.
	The calculation follows the outline in section 3.1.2 of:
	
	Fenton, John D. "Numerical methods for nonlinear waves." Advances in 
	coastal and ocean engineering 5 (1999): 241-324.
	
	Equation numbering in the source files refer to this paper.

Usage
    - calcStreamFunctinoWave inputDict outputDict
	
	where inputDict contains the necessary parameters to calculate the 
	Fourier coefficients. See waveProperties for an example.

Author
    Johan Roenby, DHI, all rights reserved

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "dictionary.H"
#include "includeEntry.H"
#include "streamFunctionWave.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{	
    argList::addNote
    (
        "Read the specified streamFunctionWaveDict input file, \n"
		"check if identical to streamFunctionWaveDict output file \n"
		"(if this exists) and if not calculating wave Fourier \n"
		"coefficients and updating output file."
    );

    argList::noBanner();
    argList::noParallel();
    argList::validArgs.append("inputDict");
	argList::validArgs.append("outputDict");
    argList args(argc, argv);

    const string inDictName = args[1];
    dictionary inDict(IFstream(inDictName)(), true);
	streamFunctionWave inWave(inDict);

    const string outDictName = args[2];	
	
	bool calcWave(true);

	if ( IFstream(outDictName) )
	{
		dictionary oldOutDict(IFstream(outDictName)(), true);
		if ( inWave.sameWave(oldOutDict) )
		{
			calcWave = false;
			Info << "Same wave in " << inDictName << " and " << outDictName << endl;
			Info << "Not recalculating" << endl;
		}
	}

	if ( calcWave )
	{
		Info << "Calculating stream function wave..." << endl;
		inWave.calcCoeffs();
		Info << "Writing wave to file " << outDictName << endl;
		inWave.writeDict(outDictName);
	}
	
    return 0;
}

// ************************************************************************* //