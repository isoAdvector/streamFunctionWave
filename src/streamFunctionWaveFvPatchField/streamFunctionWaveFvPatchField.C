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

\*---------------------------------------------------------------------------*/

#include "streamFunctionWaveFvPatchField.H"
#include "surfaceFields.H" //using lookupPatchField to get flux
#include "volFields.H"
#include "IFstream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
streamFunctionWaveFvPatchField<Type>::streamFunctionWaveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF)
{
//	Info << "Construct field " << this->dimensionedInternalField().name() << " on patch " 
//		 << this->patch().name() << " from patch and internal field." << endl;
//	Info<< "Warning: This streamFunctionWave constructor currently just sets refValue and refGrad to zero and valueFraction to 1!" << endl;
//	writeOutMembers();

	//Initialising to fixed value boundary with zero as the value
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 1.0;
}

template<class Type>
streamFunctionWaveFvPatchField<Type>::streamFunctionWaveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    waveDictName_(dict.lookupOrDefault<word>("waveDictName", "waveProperties")),
    waveDictPath_(dict.lookupOrDefault<fileName>("waveDictPath", "constant")),
    ymin_(this->patch().Cf().size()),
    ymax_(this->patch().Cf().size())
{

//	Info << "Construct field " << this->dimensionedInternalField().name() << " on patch " 
//		 << this->patch().name() << " from patch, internal field and dictionary." << endl;

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


    fileName waveDictPath;
    if (Pstream::parRun())
    {
        waveDictPath = this->db().time().path()/".."/waveDictPath_;
    }
    else
    {
        waveDictPath = this->db().time().path()/waveDictPath_;
    }
    
	IOdictionary waveProp
	(
		IOobject
		(
			waveDictName_,
			waveDictPath,
			this->db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);
    dictionary patchInWave = waveProp;
    if (waveProp.found(this->patch().name()))
    {
        patchInWave = waveProp.subDict(this->patch().name());
    }
    
	word filnam("streamFunctionWaveOutputFor" + this->patch().name());
	fileName outDictName(this->db().time().caseConstant()/filnam);

	bool calcWave(true);
	streamFunctionWave inWave(patchInWave);	
	IFstream fileStream(outDictName);
	if ( fileStream.good() )
	{
		dictionary oldOutDict(IFstream(outDictName)(), true);
		if ( inWave.sameWave(oldOutDict) )
		{
			calcWave = false;
			Info << "Wave already calculated on " << this->patch().name() << "... reusing old calculation" << endl;
			readParmsFromDict(oldOutDict);
		}
	}

	if (calcWave)
	{
		Info << "Calculating stream function wave..." << endl;
		readParmsFromDict(patchInWave);
		calcCoeffs();
		writeDict(outDictName);
	}
	
	K_ = vector(patchInWave.lookupOrDefault<vector>("direction",vector::zero));
	if ( patchInWave.found("phase") )
	{
		phi_ = patchInWave.lookupOrDefault<scalar>("phase",0.0);
	}
	else
	{
		phi_ = 1.0-zeroPointPhase();
	}
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
	
	//Setting seabedHeight_ and xWaveMaker_ with lowest and "backmost" mesh 
    //point on the patch as defaults
    {
        scalar ymin(1e20), xmin(1e20);

        const fvMesh & mesh = this->internalField().mesh();
        const pointField& pts = mesh.points();
        const label patchID = mesh.boundaryMesh().findPatchID(this->patch().name());
        const polyPatch& patch = mesh.boundaryMesh()[patchID];
        const labelList& labels(patch.meshPoints());

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
    }

    //Reading fluid densities
    {
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
        rho1_ = sD1.get<scalar>("rho");	
        dictionary sD2(tp.subDict(phases[1]));
        rho2_ = sD2.get<scalar>("rho");
    }

    //Find minimum and maximum vertical point for each face and populate ymin_and ymax_
    {
        const faceList& fs = this->internalField().mesh().faces(); //Handle to list of all faces
        const pointField& p = this->internalField().mesh().points(); //Handle to list of all grid points
        const label start = this->patch().patch().start(); //index in mesh face list of first face belonging to this patch

        forAll(ymin_,facei)
        {
            const labelList& f = fs[facei + start];
            label nPoints = f.size();
            ymin_[facei] = 1e20;
            ymax_[facei] = -1e20;
            for (label pi = 0; pi < nPoints; pi++)
            {
                scalar yi = (up_ & (p[f[pi]])) - seabedHeight_;
                if (yi < ymin_[facei]) ymin_[facei] = yi;
                if (yi > ymax_[facei]) ymax_[facei] = yi;
            }
        }
    }

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

//	writeOutMembers();
	evaluate();
}

template<class Type>
streamFunctionWaveFvPatchField<Type>::streamFunctionWaveFvPatchField
(
    const streamFunctionWaveFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper)
{
//	These lines made code crash when run in parallel:
//	Info << "Construct field " << this->dimensionedInternalField().name() << " on patch " 
//		 << this->patch().name() << " by mapping given streamFunctionWaveFvPatchField onto a new patch." << endl;
	copyInit(ptf);
}


template<class Type>
streamFunctionWaveFvPatchField<Type>::streamFunctionWaveFvPatchField
(
    const streamFunctionWaveFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF)
{
//	Info << "Construct field " << this->dimensionedInternalField().name() << " on patch " 
//		 << this->patch().name() << " as copy setting internal field reference." << endl;
	copyInit(ptf);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Evaluate the field on the patch
template<class Type>
void streamFunctionWaveFvPatchField<Type>::evaluate()
{
//	Info<< "Enter streamFunctionWaveFvPatchField<Type>::evaluate()" << endl;
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
void streamFunctionWaveFvPatchField<Type>::updateCoeffs()
{
//	Info<< "Enter streamFunctionWaveFvPatchField<Type>::updateCoeffs()" << endl;
	
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

	//Setting air velocity 
	scalarField umax = getU(x,eta,t,phi_);
	scalarField vmax = getV(x,eta,t,phi_);

	//Calculating volume fraction (assuming parallelogram shaped - e.g. 
    //rectangular - faces with horizontal sides)
	scalarField alpha(y.size(),0.0);
	
	forAll(y,facei)
    {
		alpha[facei] = min(max((eta[facei]-ymin_[facei])/(ymax_[facei]-ymin_[facei]),0.0),1.0);

        if (ymin_[facei] > eta[facei])
        {
            u[facei] = 0;
            v[facei] = 0;
        }
        else if (ymin_[facei] <= eta[facei] && eta[facei] <= ymax_[facei])
        {
            u[facei] = umax[facei];
            v[facei] = vmax[facei];
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
void streamFunctionWaveFvPatchField<vector>::setField
(
	Field<vector> const& refVal, 
	const scalarField& alpha, 
	const scalarField& u, 
	const scalarField& v
)
{			
//	Info << "Enter streamFunctionWaveFvPatchField<vector>::setField" << endl;
	this->refValue() = (u*K_/mag(K_) + v*up_);
	this->refGrad() = vector::zero;
	this->valueFraction() = 1.0;
//	this->valueFraction() = max(alpha-SMALL);
} 

// Setting volume fraction (alpha1) and pressure (pd etc) field
template<> 
void streamFunctionWaveFvPatchField<scalar>::setField
(
	Field<scalar> const& iF, 
	const scalarField& alpha, 
	const scalarField& u, 
	const scalarField& v
)
{	
//	Info << "Enter streamFunctionWaveFvPatchField<scalar>::setField" << endl;

	const word& fieldName = this->internalField().name();
	if (fieldName == "alpha1" || fieldName == "gamma" || fieldName == "alpha.water")
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
void streamFunctionWaveFvPatchField<Type>::setField
(
	Field<Type> const& refVal, 
	const scalarField& alpha, 
	const scalarField& u, 
	const scalarField& v
)
{
//	Info << "Enter streamFunctionWaveFvPatchField<Type>::setField" << endl;
//   dummy code executed if setField is called with other Type than scalar or vector
}   

// Initialising 
template<class Type>
void streamFunctionWaveFvPatchField<Type>::copyInit(const streamFunctionWaveFvPatchField<Type>& ptf)
{
//	Info << "Enter streamFunctionWaveFvPatchField<Type>::copyInit(const streamFunctionWaveFvPatchField<Type>& ptf)" << endl;
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
	u_ = ptf.u_;
	R_ = ptf.R_;
	up_ = ptf.up_;	
    waveDictName_ = ptf.waveDictName_;
    waveDictPath_ = ptf.waveDictPath_;
    ymin_ = ptf.ymin_;
    ymax_ = ptf.ymax_;
    
//	Info << "Out of streamFunctionWaveFvPatchField<Type>::copyInit(const streamFunctionWaveFvPatchField<Type>& ptf)" << endl;
//	writeOutMembers();
}

// Write
template<class Type>
void streamFunctionWaveFvPatchField<Type>::write(Ostream& os) const
{		
//	Info << "Entering write function and doing fvPatchField<Type>::write(os);" << endl;

    fvPatchField<Type>::write(os);
	os.writeKeyword("H") << H_ << token::END_STATEMENT << nl;
	os.writeKeyword("d") << d_ << token::END_STATEMENT << nl;
	if (omega_ > 0) //Period is finite 
	{
		os.writeKeyword("T") << (2*pi_/omega_) << token::END_STATEMENT << nl;
	}
	os.writeKeyword("uType") << uType_ << token::END_STATEMENT << nl;
	os.writeKeyword("u") << u_ << token::END_STATEMENT << nl;
	os.writeKeyword("N") << N_ << token::END_STATEMENT << nl;
	os.writeKeyword("nHsteps") << nHsteps_ << token::END_STATEMENT << nl;
	os.writeKeyword("pi") << pi_ << token::END_STATEMENT << nl;
	os.writeKeyword("g") << g_ << token::END_STATEMENT << nl;
	os.writeKeyword("k") << k_ << token::END_STATEMENT << nl;
	os.writeKeyword("R") << R_ << token::END_STATEMENT << nl;
	os.writeKeyword("uBar") << uBar_ << token::END_STATEMENT << nl;
	os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
	os.writeKeyword("B") << B_ << token::END_STATEMENT << nl;
	os.writeKeyword("direction") << K_ << token::END_STATEMENT << nl;
	os.writeKeyword("phase") << phi_ << token::END_STATEMENT << nl;
	os.writeKeyword("waveDictName") << waveDictName_ << token::END_STATEMENT << nl;
	os.writeKeyword("waveDictPath") << waveDictPath_ << token::END_STATEMENT << nl;

    this->writeEntry("value", os);
}

template<class Type>
void streamFunctionWaveFvPatchField<Type>::writeOutMembers()
{
//	Info << "Enter  streamFunctionWaveFvPatchField<Type>::writeOutMembers()" << endl;
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

	Info << "u_ " << u_ << endl;
	Info << "uBar_ " << uBar_ << endl;
	Info << "R_ " << R_ << endl;
	Info << "E_ " << E_ << endl;
	Info << "B_ " << B_ << endl;
	Info << "waveDictName_ " << waveDictName_ << endl;
	Info << "waveDictPath_ " << waveDictPath_ << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
