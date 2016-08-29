/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Class
    Foam::streamFunctionWave

Description
    Calculates Fourier coefficients for a stream function theory wave.
	The calculation follows the outline in section 3.1.2 of:
	
	Fenton, John D. "Numerical methods for nonlinear waves." Advances in 
	coastal and ocean engineering 5 (1999): 241-324.
	
	Equation numbering in the source files refer to this paper.
	
SourceFiles
    streamFunctionWave.C
	
Author
	Johan Roenby, DHI

\*---------------------------------------------------------------------------*/

#ifndef streamFunctionWave_H
#define streamFunctionWave_H

#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class streamFunctionWave
{
protected:
	//- Number of Fourier modes
	scalar N_;

	//- Fourier coefficients for velocity field and pressure
	scalarField B_;

	//- Fourier coefficients for surface elevation
	scalarField E_;

	//- Wave Height
	scalar H_;

	//- Mean water level
	scalar d_;
			
	//- Current type
	label uType_; //0 for wavelength input, 1 for wave period and "Eulerian" current, 2 for wave period and "Stokes" current
	
	//- Number of steps to reach final wave height
	label nHsteps_;
				
	//- 3.1415...
	scalar pi_;
	
	//- Gravitational acceleration
	scalar g_;

	//- Current - either u1 or u2 - type determined by uType_
	scalar u_;

	//- Wave number k = 2*pi/lambda
	scalar k_;

	//- Angular wave frequency
	scalar omega_;
	
	//- Mean fluid speed in co-moving frame
	scalar uBar_;
	
	//- Bernoulli constant
	scalar R_;

	//- Volume flow rate per unit span underneath the waves in the co-mvoing frame
	scalar Q_;

	void findRoot(scalarField& X, scalar h, scalar lort);

	scalarField fun(const scalarField& X, scalar h, scalar lort);
	
	scalarSquareMatrix jac(const scalarField& X, scalar h, scalar lort);

public:

	streamFunctionWave();

	streamFunctionWave(const dictionary&);
	
	void readParmsFromDict(const dictionary&);

	bool sameWave(const dictionary&);

	void calcCoeffs();

	scalarField getEta
	(
		const scalarField& x,
		const scalar& time,
		const scalar& phi
	);

	scalarField getU
	(
		const scalarField& x,
		const scalarField& y,
		const scalar& t,
		const scalar& phi
	);

	scalarField getV
	(
		const scalarField& x,
		const scalarField& y,
		const scalar& t,
		const scalar& phi
	);

	//- Find x/lambda value where eta = d
	scalar zeroPointPhase();

	void writeDict(fileName);	

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
