/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is an extension of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    streamFunctionWaveFvPatchField

Description
	This is an implementation of the Fourier approximation method for finding
	stream function waves described in section 3.1.2 of the paper by J.D.
	Fenton, "Numerical Methods for Nonlinear Waves" (1999), Numerical Methods 
	for Nonlinear Waves, Advances in Coastal and Ocean Engineering, Vol. 5,
	pp241-324. Equation numbers in the code refer to this paper. At the time
	of writing the paper could be downloaded here:

	http://www.johndfenton.com/Papers/Fenton99Liu-Numerical-methods-for-nonlinear-waves.pdf	

	Usage: 
	
	The wave generator patch velocity should be specified in 0/U like

		patchName            
		{
			type            streamFunctionWave;
			value			uniform (0 0 0);
		}

	Similarly the dynamic pressure in 0/pd and volumen fraction field in 
	0/alpha1 must look like:

		patchName            
		{
			type            streamFunctionWave;
			value			uniform 0;
		}

	In the <case>/constant folder there should be a file named waveProperties
	containing the wave parameters for the patch in a subdictionary named after the patch:

	patchName
	{
	//	Name			Value ex.	Description
		H	          0.1;			//Wave heigth
		T   	       20;			//Wave period (only used if uType = 1 or 2)
		length			2;			//Wave length (only used if uType = 0)
		uType           2; 			//Input type. 0: wave lenght specified. 1 or 2: wave period and current specified (default = 2).
		u1or2           0; 			//Current, Eulerian if uType = 1, Stokes if uType = 2 (default = 0).

		d       	    0.3;		//Water depth
		direction       ( 1 0 0 ); 	//Direction of wave vector
		seabedHeight	-0.3; 		//Vertical coordinate of seabed on patch (default = lowest point on patch)
		waveMakerPos	-3.0; 		//Horizontal position of wave maker (default = backmost point on patch)
		phase			.1;			//Phase of wave, theta = k*x+omega*t + 2*pi*phase (default: the value that gives eta = d)

		N          40; 		//Number of Fourier modes in expansion (default = 30)
		nHSteps    20; 		//Number of height steps in search for wave solution (default = 20)
	}
	
	This folder must also contain a file named g containing the gravitation 
	vector, e.g.:
	
		dimensions      [0 1 -2 0 0 0 0];
		value           ( 0 -9.81 0 );

	The file <case>/system/controlDict must end with the line:
	
		libs ( "libstreamFunctionWaveFvPatchField.so" "libOpenFOAM.so" );	
		
SourceFiles
    streamFunctionWaveFvPatchField.C

Author
	Johan Roenby, DHI.

\*---------------------------------------------------------------------------*/

#ifndef streamFunctionWaveFvPatchField_H
#define streamFunctionWaveFvPatchField_H

#include "mixedFvPatchFields.H"
#include "uniformDimensionedFields.H" //gravity type of field
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
    public mixedFvPatchField<Type>,
	public streamFunctionWave
{
protected:
		
		//- Wave vector
        vector K_;
		
		//- Unit vector pointing upwards
        vector up_;

		//- Wave phase in radians
        scalar phi_;
		
		//- Vertical coordinate of seabed at patch
		scalar seabedHeight_;

		//- Horizontal position of wave maker patch
		scalar xWaveMaker_;

		//- Water density
		scalar rho1_;

		//- Air density
		scalar rho2_;
	
    // Protected member functions
				
		//Initialize member data from other patch (avoid constructor code duplication)
		void copyInit(const streamFunctionWaveFvPatchField<Type>&);
								
		void writeOutMembers();

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
        //  onto a new patch
        streamFunctionWaveFvPatchField
        (
            const streamFunctionWaveFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new streamFunctionWaveFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        streamFunctionWaveFvPatchField
        (
            const streamFunctionWaveFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new streamFunctionWaveFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions
				
		//- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

		//- Update the coefficients associated with the patch field
        virtual void evaluate();

        //- Set field values on patch
		virtual void setField
		(
			const Field<Type>&, 
			const scalarField&, 
			const scalarField&, 
			const scalarField&
		);
		
		//- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "streamFunctionWaveFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
