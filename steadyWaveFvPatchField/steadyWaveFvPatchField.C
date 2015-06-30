/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "steadyWaveFvPatchField.H"
#include "surfaceFields.H" //using lookupPatchField to get flux
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
steadyWaveFvPatchField<Type>::steadyWaveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF)
{
	Info << "Construct field " << this->dimensionedInternalField().name() << " on patch " 
		 << this->patch().name() << " from patch and internal field." << endl;
	Info<< "STEADYWAVEFVPATCHFIELD WARNING: THIS CONSTRUCTOR CURRENTLY DOES NOT INITIALISE PATCH PROPERLY!!!" << endl;
//	writeOutMembers();

	//Initialising to fixed value boundary with zero as the value
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 1.0;
}

template<class Type>
steadyWaveFvPatchField<Type>::steadyWaveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF)
{
	Info << "Construct field " << this->dimensionedInternalField().name() << " on patch " 
		 << this->patch().name() << " from patch, internal field and dictionary." << endl;
	
//	if ( dict.found("readFromWaveDict") && !dict.lookup("readFromWaveDict") ) //reading from constant/waveProperties
	{
		dictionary patchWaveProperties;
		{
			IOdictionary waveProp
			(
				IOobject
				(
					"waveProperties",
					this->db().time().constant(),
					this->db(),
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				)
			);
			patchWaveProperties = waveProp.subDict(this->patch().name());
		}
		//This is supposed to set N_, H_, d_, uType, nHSteps_, pi_, g_ by reading from subDict but not sure it does:
//		steadyWave sw(steadyWave::steadyWave(patchWaveProperties)); 
		readParmsFromDict(patchWaveProperties);

		K_ = vector(patchWaveProperties.lookupOrDefault<vector>("direction",vector::zero));
		phi_ = patchWaveProperties.lookupOrDefault<scalar>("phase",0.0);
	}

	vector gVector = 
	uniformDimensionedVectorField
	(
		IOobject
		(
			"g",
			this->db().time().constant(),
			this->db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	).value();
	g_ = mag(gVector);
	up_ = -gVector/g_;

	//Here we should first check if any other field on this patch has already calculated its 
	calcCoeffs();
	phi_ = 1.0-zeroPointPhase();
	
	if ( mag(up_ & K_) > SMALL )
	{
		Info << "Steady wave warning: Wave vector has non-zero vertical component which will be ignored." << endl;
		K_ -= (up_ & K_)*up_;;
	}

	if ( mag(K_) < SMALL )
	{
		Info << "Direction vector has zero length. Using face normals for wave direction." << endl; //STILL NOT IMPLEMENTED!!!!!!!
	}
	else
	{
		K_ = K_/mag(K_);
	}
	
	//Finding lowest and the "backmost" mesh point on the patch
    const fvMesh & mesh = this->dimensionedInternalField().mesh();
	const pointField& pts = mesh.points();
    const label patchID = mesh.boundaryMesh().findPatchID(this->patch().name());
	const polyPatch& patch = mesh.boundaryMesh()[patchID];
	const labelList& labels(patch.meshPoints());

	scalar ymin(1e20), xmin(1e20);
	forAll(labels, pi)
	{
		scalar yi = (up_ & (pts[labels[pi]]));
		scalar xi = (K_ & (pts[labels[pi]]));
		if (yi < ymin) ymin = yi;
		if (xi < xmin) xmin = xi;
	}
	reduce(ymin, minOp<scalar>());    
	reduce(xmin, minOp<scalar>());    
	
	//Vertical position of seabed at the boundary in coordinates used to generate mesh 
	seabedHeight_ = dict.lookupOrDefault<scalar>("seabedHeight",ymin);

	//Horizontal wave maker position
	xWaveMaker_ = dict.lookupOrDefault<scalar>("waveMakerPos",xmin);
	
	IOdictionary tp
	(
		IOobject
		(
			"transportProperties",
			this->db().time().constant(),
			this->db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	wordList defaultPhaseNames(2);
	defaultPhaseNames[0] = "phase1";
	defaultPhaseNames[1] = "phase2";
	wordList phases(wordList(tp.lookupOrDefault("phases",defaultPhaseNames))); //In OF-2.3.0 the phases keyword is required in transportProperties

	dictionary sD1(tp.subDict(phases[0]));
	rho1_ = (dimensionedScalar(sD1.lookup("rho"))).value();	
	dictionary sD2(tp.subDict(phases[1]));
	rho2_ = (dimensionedScalar(sD2.lookup("rho"))).value();	

    if (dict.found("value")) //Not sure why this is needed but it exists in many BC's derived from mixed, so better keep it (JRO).
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }

    this->refValue() = *this;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;

	evaluate();
	writeOutMembers();
}

template<class Type>
steadyWaveFvPatchField<Type>::steadyWaveFvPatchField
(
    const steadyWaveFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper)
{
	Info << "Construct field " << this->dimensionedInternalField().name() << " on patch " 
		 << this->patch().name() << " by mapping given steadyWaveFvPatchField onto a new patch." << endl;
	copyInit(ptf);
}


template<class Type>
steadyWaveFvPatchField<Type>::steadyWaveFvPatchField
(
    const steadyWaveFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF)
{
	Info << "Construct field " << this->dimensionedInternalField().name() << " on patch " 
		 << this->patch().name() << " as copy setting internal field reference." << endl;
	copyInit(ptf);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Evaluate the field on the patch
template<class Type>
void steadyWaveFvPatchField<Type>::evaluate()
{
//	Info<< "Enter steadyWaveFvPatchField<Type>::evaluate()" << endl;
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    
    Field<Type>::operator=
    (
        this->valueFraction()*this->refValue()
      +
        (1.0 - this->valueFraction())*
        (
            this->patchInternalField()
          + this->refGrad()/this->patch().deltaCoeffs()
        )
    );

    mixedFvPatchField<Type>::evaluate();
}

// Update the coefficients associated with the patch field
template<class Type>
void steadyWaveFvPatchField<Type>::updateCoeffs()
{
//	Info<< "Enter steadyWaveFvPatchField<Type>::updateCoeffs()" << endl;

	if (this->updated())
    {
        return;
    }
	
	scalarField x = (K_/mag(K_)) & (this->patch().Cf());
	x -= xWaveMaker_;
	scalarField y = up_ & (this->patch().Cf());
	y -= seabedHeight_; //Vertical coordinate as measured from seabed
	
	scalar t = this->db().time().value();
	
	scalarField eta = getEta(x,t,phi_);
	scalarField u = getU(x,y,t,phi_);
	scalarField v = getV(x,y,t,phi_);
/*		
	scalarField theta = k_*x - omega_*t + 2*pi_*phi_;


	//Calculating surface elevation above sealevel
	scalarField eta(x.size(),0.0);
	forAll(E_,jj)
	{
		eta += 2*E_[jj]*cos(jj*theta); //Fenton (1999), eq. (3.27)
	}
	
	//Calculating velocity field components
	scalarField u(y.size(),0.0);
	scalarField v(y.size(),0.0);
	int	jj;
	forAll(B_,ii)
	{
		jj = ii + 1;
		u += jj*B_[ii]*cosh(jj*k_*y)/cosh(jj*k_*d_)*cos(jj*theta); //Fenton (1999), eq. (3.24)
		v += jj*B_[ii]*sinh(jj*k_*y)/cosh(jj*k_*d_)*sin(jj*theta); //Fenton (1999), eq. (3.25)
	}
	u *= sqrt(g_/k_);
	v *= sqrt(g_/k_);
	u += omega_/k_ - uBar_;
*/
/*
	//Attenuating fields during ramp up time - CURRENTLY NOT PROPERLY IMPLEMENTED!!!!!!!!!!
	scalar fac(1.0);
	if (t < rampUpTime_)
	{
		fac = sin(pi_/2.0*t/rampUpTime_);
		u *= fac;
		v *= fac;
		eta = fac*(eta-d_) + d_;
	}
*/

	//Calculating volume fraction (assuming parallelogram shaped - e.g. rectangular - faces with horizontal sides)
	const faceList& fs = this->dimensionedInternalField().mesh().faces(); //Handle to list of all faces
	const pointField& p = this->dimensionedInternalField().mesh().points(); //Handle to list of all grid points
	const label start = this->patch().patch().start(); //index in mesh face list of first face belonging to this patch
	scalarField alpha(y.size(),0.0);
	
	scalar umax(0.0), vmax(0.0), yumax(0.0);
	forAll(y,facei)
    {
        const labelList& f = fs[facei + start];
        label nPoints = f.size();
		scalar ymin(1e20), ymax(-1e20);
		for (label pi = 0; pi < nPoints; pi++)
        {
			scalar yi = (up_ & (p[f[pi]])) - seabedHeight_;
			if (yi < ymin) ymin = yi;
			if (yi > ymax) ymax = yi;
        }
		alpha[facei] = min(max((eta[facei]-ymin)/(ymax-ymin),0.0),1.0);
		u[facei] *= pos(alpha[facei]-rho2_/rho1_);
		v[facei] *= pos(alpha[facei]-rho2_/rho1_);
		if (mag(u[facei]) > mag(umax)) 
		{
			umax = u[facei];
			vmax = v[facei];
			yumax = ymax;// + alpha[facei]*(ymax-ymin);
		}
    }
	reduce(umax, maxOp<scalar>());
	reduce(vmax, maxOp<scalar>());    
	reduce(yumax, maxOp<scalar>());    

    //Setting velocity field in air phase to decrease linearly from -(umax,vmax) just above free surface to (0,0) at roof.
    scalar yRoof(gMax(y));
	scalar magU = sqrt(umax*umax+vmax*vmax);
    forAll(y,facei)
    {
		scalar yf = y[facei];
		if (yf > yumax)
		{
				u[facei] = 0.0;//-mag(umax)*(yf-yRoof)/(yumax-yRoof);
				v[facei] = sign(vmax)*magU*(yf-yRoof)/(yumax-yRoof);
//				u[facei] = umax*(yf-yRoof)/(yumax-yRoof);
//				v[facei] = vmax*(yf-yRoof)/(yumax-yRoof);
		}
	}

	//Setting field depending on its type
	setField(this->refValue(), alpha, u, v);

    mixedFvPatchField<Type>::updateCoeffs(); //This line simply runs fvPatchField::updateCoeffs() which sets its private boolean updated_ = true;

}

//In the following we use template specialization to set field depending on its Type: 
//Inspired by http://www.parashift.com/c++-faq-lite/templates.html#faq-35.7

// Setting velocity field
template<> 
void steadyWaveFvPatchField<vector>::setField
(
	Field<vector> const& refVal, 
	const scalarField& alpha, 
	const scalarField& u, 
	const scalarField& v
)
{			
//	Info << "Enter steadyWaveFvPatchField<vector>::setField" << endl;
	this->refValue() = u*K_/mag(K_) + v*up_;
	this->refGrad() = vector::zero;
//	this->valueFraction() = 1.0;
	this->valueFraction() = pos(alpha-SMALL);
} 

// Setting volume fraction (alpha1) and pressure (pd etc) field
template<> 
void steadyWaveFvPatchField<scalar>::setField
(
	Field<scalar> const& iF, 
	const scalarField& alpha, 
	const scalarField& u, 
	const scalarField& v
)
{	
//	Info << "Enter steadyWaveFvPatchField<scalar>::setField" << endl;

	const word& fieldName = this->dimensionedInternalField().name();

	if (fieldName == "alpha1" || fieldName == "gamma")
	{			
		this->refValue() = alpha;
		this->refGrad() = 0.0;
		this->valueFraction() = 1.0;
	}
	else if (fieldName == "pd" || fieldName == "p_rgh" || fieldName == "ph_rgh")
	{
//		if (this->db().time().timeIndex() != 0)
		if (this->db().time().value() != this->db().time().startTime().value())
		{
//			Failed attempts with fixedValue for pd:
//			scalarField rho(alpha*rho1_+(1-alpha)*rho2_);
//			scalar c = omega_/k_;
//			this->refValue() =  rho*( R_ - 0.5*pow(u-c,2) - 0.5*pow(v,2) ); //Fenton (1999) eq. (3.28)
//			this->valueFraction() = 1.0; //fixedValue
//			this->valueFraction() = (1.0-alpha); //Fixed value in water and fixed Gradient in air - Runs way too slow with this

//			Taken from buoyantPressure. Corresponds to a fixed wall with no slip. Not optimal but works.
			const fvPatchField<scalar>& rho =
				patch().lookupPatchField<volScalarField, scalar>("rho");
			this->refGrad() = -rho.snGrad()*((-g_*up_) & this->patch().Cf()); 
			this->valueFraction() = 0.0; //fixedGradient as in buoyantPressure
		}
		else
		{
			this->refValue() =  0.0;
			this->refGrad() = 0.0;
			this->valueFraction() = 0.0;
		}
	}	
}

template<class Type>
void steadyWaveFvPatchField<Type>::setField
(
	Field<Type> const& refVal, 
	const scalarField& alpha, 
	const scalarField& u, 
	const scalarField& v
)
{
//	Info << "Enter steadyWaveFvPatchField<Type>::setField" << endl;
//   dummy code executed if setField is called with other Type than scalar or vector
}   

// Initialising 
template<class Type>
void steadyWaveFvPatchField<Type>::copyInit(const steadyWaveFvPatchField<Type>& ptf)
{
//	Info << "Enter steadyWaveFvPatchField<Type>::copyInit(const steadyWaveFvPatchField<Type>& ptf)" << endl;
	N_ = ptf.N_;	
	H_ = ptf.H_;
	d_ = ptf.d_;
	k_ = ptf.k_;
	omega_ = ptf.omega_;
	uBar_ = ptf.uBar_;
	E_ = ptf.E_;
	B_ = ptf.B_;
	K_ = ptf.K_;
	phi_ = ptf.phi_;

	pi_ = ptf.pi_;
	g_ = ptf.g_;

	seabedHeight_ = ptf.seabedHeight_;
	xWaveMaker_ = ptf.xWaveMaker_;
	rho1_ = ptf.rho1_;
	rho2_ = ptf.rho2_;
	
	uType_ = ptf.uType_;
	nHsteps_ = ptf.nHsteps_;
	u1or2_ = ptf.u1or2_;
	R_ = ptf.R_;
	up_ = ptf.up_;	

//	writeOutMembers();
}

// Write
template<class Type>
void steadyWaveFvPatchField<Type>::write(Ostream& os) const
{		
//	Info << "Entering write function and doing fvPatchField<Type>::write(os);" << endl;

    fvPatchField<Type>::write(os);
//	os.writeKeyword("d") << d_ << token::END_STATEMENT << nl;
//	os.writeKeyword("k") << k_ << token::END_STATEMENT << nl;
//	os.writeKeyword("omega") << omega_ << token::END_STATEMENT << nl;
//	os.writeKeyword("uBar") << uBar_ << token::END_STATEMENT << nl;
//	os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
//	os.writeKeyword("B") << B_ << token::END_STATEMENT << nl;
//	os.writeKeyword("direction") << K_ << token::END_STATEMENT << nl;
//	os.writeKeyword("phase") << phi_ << token::END_STATEMENT << nl;
//	os.writeKeyword("rampUpTime") << rampUpTime_ << token::END_STATEMENT << nl;

/*    os.writeKeyword("depth")
        << d_ << token::END_STATEMENT << nl;
    os.writeKeyword("period")
        << 2*pi_/omega_ << token::END_STATEMENT << nl;
    os.writeKeyword("length")
        << 2*pi_/k_ << token::END_STATEMENT << nl;
 */
//	Info << "Doing this->writeEntry('value', os);" << endl;
    this->writeEntry("value", os);
}

template<class Type>
void steadyWaveFvPatchField<Type>::writeOutMembers()
{
//	Info << "Enter  steadyWaveFvPatchField<Type>::writeOutMembers()" << endl;
	Info << "N_ " << N_ << endl;
	Info << "H_ " << H_ << endl;
	Info << "d_ " << d_ << endl;
	Info << "uType_ " << uType_ << endl;
	Info << "nHsteps_ " << nHsteps_ << endl;
	Info << "pi_ " << pi_ << endl;
	Info << "g_ " << g_ << endl;
	Info << "up_ " << up_ << endl;
	Info << "K_ " << K_ << endl;
	Info << "seabedHeight_ " << seabedHeight_ << endl;
	Info << "xWaveMaker_ " << xWaveMaker_ << endl;
	Info << "phi_ " << phi_ << endl;

	Info << "rho1_ " << rho1_ << endl;
	Info << "rho2_ " << rho2_ << endl;
	Info << "k_ " << k_ << endl;
	Info << "omega_ " << omega_ << endl;

	Info << "u1or2_ " << u1or2_ << endl;
	Info << "uBar_ " << uBar_ << endl;
	Info << "R_ " << R_ << endl;
	Info << "E_ " << E_ << endl;
	Info << "B_ " << B_ << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
