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

\*---------------------------------------------------------------------------*/

#include "streamFunctionWaveFvPatchField.H"
#include "volFields.H"
#include "mathematicalConstants.H"

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
    // Initialising to a fixed value boundary with zero as the value
    this->refValue() = Zero;
    this->refGrad() = Zero;
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
    waveDictName_(dict.getOrDefault<word>("waveDictName", "waveProperties")),
    waveDictPath_(dict.getOrDefault<fileName>("waveDictPath", "constant")),
    ymin_(p.size()),
    ymax_(p.size())
{
    const vector gVector
    (
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
        ).value()
    );

    up_ = -gVector/mag(gVector);

    if
    (
        dict.isDict("wave")
     && streamFunctionWave::hasCoeffs(dict.subDict("wave"))
    )
    {
        // The wave was written by this patch, e.g. before decomposePar,
        // reconstructPar or a restart, so it is not recalculated
        wave_.readParmsFromDict(dict.subDict("wave"));

        K_ = dict.get<vector>("direction");
        phi_ = dict.get<scalar>("phase");
    }
    else
    {
        const fileName waveDictPath
        (
            UPstream::parRun()
          ? this->db().time().path()/".."/waveDictPath_
          : this->db().time().path()/waveDictPath_
        );

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

        // Wave input either at the top level or in a sub-dictionary named
        // after the patch
        const dictionary& patchInWave =
            waveProp.optionalSubDict(this->patch().name());

        // The gravitational acceleration is taken from constant/g unless
        // given in the wave dictionary
        dictionary waveDict(patchInWave);
        if (!waveDict.found("g"))
        {
            waveDict.add("g", mag(gVector));
        }

        Info<< "Calculating stream function wave on patch "
            << this->patch().name() << endl;

        wave_.readParmsFromDict(waveDict);
        wave_.calcCoeffs();

        K_ = patchInWave.get<vector>("direction");

        if (!patchInWave.readIfPresent("phase", phi_))
        {
            phi_ = 1.0 - wave_.zeroPointPhase();
        }
    }

    if (mag(up_ & K_) > SMALL)
    {
        WarningInFunction
            << "Wave vector on patch " << this->patch().name()
            << " has a non-zero vertical component which will be ignored."
            << endl;

        K_ -= (up_ & K_)*up_;
    }

    if (mag(K_) < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "The wave direction on patch " << this->patch().name()
            << " has zero length in the horizontal plane." << nl
            << "Deriving the wave direction from the face normals is not"
            << " implemented, so a non-zero direction must be given."
            << exit(FatalIOError);
    }

    K_ /= mag(K_);

    // Setting seabedHeight_ and xWaveMaker_ with the lowest and "backmost"
    // mesh point on the patch as defaults
    {
        const pointField& pts = this->internalField().mesh().points();
        const labelList& meshPts = this->patch().patch().meshPoints();

        scalar ymin(1e20);
        scalar xmin(1e20);

        for (const label pointi : meshPts)
        {
            ymin = min(ymin, up_ & pts[pointi]);
            xmin = min(xmin, K_ & pts[pointi]);
        }

        reduce(ymin, minOp<scalar>());
        reduce(xmin, minOp<scalar>());

        // Vertical position of the seabed at the boundary, in the
        // coordinates used to generate the mesh
        seabedHeight_ = dict.getOrDefault<scalar>("seabedHeight", ymin);

        // Horizontal wave maker position
        xWaveMaker_ = dict.getOrDefault<scalar>("waveMakerPos", xmin);
    }

    // Minimum and maximum vertical position of each patch face
    {
        const faceList& fs = this->internalField().mesh().faces();
        const pointField& pts = this->internalField().mesh().points();
        const label start = this->patch().patch().start();

        forAll(ymin_, facei)
        {
            const face& f = fs[facei + start];

            ymin_[facei] = 1e20;
            ymax_[facei] = -1e20;

            for (const label pointi : f)
            {
                const scalar yi = (up_ & pts[pointi]) - seabedHeight_;

                ymin_[facei] = min(ymin_[facei], yi);
                ymax_[facei] = max(ymax_[facei], yi);
            }
        }
    }

    if (!this->readValueEntry(dict))
    {
        fvPatchField<Type>::extrapolateInternal();
    }

    this->refValue() = *this;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;

    this->evaluate();
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
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    wave_(ptf.wave_),
    K_(ptf.K_),
    up_(ptf.up_),
    phi_(ptf.phi_),
    seabedHeight_(ptf.seabedHeight_),
    xWaveMaker_(ptf.xWaveMaker_),
    waveDictName_(ptf.waveDictName_),
    waveDictPath_(ptf.waveDictPath_),
    ymin_(ptf.ymin_),
    ymax_(ptf.ymax_)
{}


template<class Type>
streamFunctionWaveFvPatchField<Type>::streamFunctionWaveFvPatchField
(
    const streamFunctionWaveFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    wave_(ptf.wave_),
    K_(ptf.K_),
    up_(ptf.up_),
    phi_(ptf.phi_),
    seabedHeight_(ptf.seabedHeight_),
    xWaveMaker_(ptf.xWaveMaker_),
    waveDictName_(ptf.waveDictName_),
    waveDictPath_(ptf.waveDictPath_),
    ymin_(ptf.ymin_),
    ymax_(ptf.ymax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void streamFunctionWaveFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    scalarField x((K_/mag(K_)) & this->patch().Cf());
    x -= xWaveMaker_;

    // Vertical coordinate as measured from the seabed
    scalarField y(up_ & this->patch().Cf());
    y -= seabedHeight_;

    const scalar t = this->db().time().value();

    const scalarField eta(wave_.getEta(x, t, phi_));
    scalarField u(wave_.getU(x, y, t, phi_));
    scalarField v(wave_.getV(x, y, t, phi_));

    // Velocity at the free surface, used above it
    const scalarField umax(wave_.getU(x, eta, t, phi_));
    const scalarField vmax(wave_.getV(x, eta, t, phi_));

    // Calculating the volume fraction, assuming parallelogram shaped -
    // e.g. rectangular - faces with horizontal sides
    scalarField alpha(y.size(), Zero);

    forAll(y, facei)
    {
        alpha[facei] =
            min
            (
                max
                (
                    (eta[facei] - ymin_[facei])/(ymax_[facei] - ymin_[facei]),
                    0.0
                ),
                1.0
            );

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

    // Setting the field depending on its type
    setField(this->refValue(), alpha, u, v);

    // Sets the updated_ flag of fvPatchField
    mixedFvPatchField<Type>::updateCoeffs();
}


// In the following, template specialisation is used to set the field
// depending on its Type

//- Setting the velocity field
template<>
void streamFunctionWaveFvPatchField<vector>::setField
(
    const Field<vector>& refVal,
    const scalarField& alpha,
    const scalarField& u,
    const scalarField& v
)
{
    this->refValue() = (u*K_/mag(K_) + v*up_);
    this->refGrad() = Zero;
    this->valueFraction() = 1.0;
}


//- Setting the volume fraction (alpha1) and the pressure (p_rgh etc) field
template<>
void streamFunctionWaveFvPatchField<scalar>::setField
(
    const Field<scalar>& iF,
    const scalarField& alpha,
    const scalarField& u,
    const scalarField& v
)
{
    const word& fieldName = this->internalField().name();

    if (fieldName == "alpha1" || fieldName == "gamma" || fieldName == "alpha.water")
    {
        this->refValue() = alpha;
        this->refGrad() = 0.0;
        this->valueFraction() = 1.0;
    }
    else if (fieldName == "pd" || fieldName == "p_rgh" || fieldName == "ph_rgh")
    {
        if (this->db().time().value() != this->db().time().startTime().value())
        {
            // Taken from buoyantPressure. Corresponds to a fixed wall with
            // no slip. Not optimal, but works.
            const fvPatchField<scalar>& rho =
                patch().lookupPatchField<volScalarField>("rho");

            this->refGrad() = -rho.snGrad()*((-wave_.g()*up_) & this->patch().Cf());
            this->valueFraction() = 0.0;    // fixedGradient as in buoyantPressure
        }
        else
        {
            this->refValue() = 0.0;
            this->refGrad() = 0.0;
            this->valueFraction() = 0.0;
        }
    }
}


//- Dummy code, executed if setField is called with a Type other than
//- scalar or vector
template<class Type>
void streamFunctionWaveFvPatchField<Type>::setField
(
    const Field<Type>& refVal,
    const scalarField& alpha,
    const scalarField& u,
    const scalarField& v
)
{}


template<class Type>
void streamFunctionWaveFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    os.beginBlock("wave");
    wave_.write(os);
    os.endBlock();

    os.writeEntry("direction", K_);
    os.writeEntry("phase", phi_);
    os.writeEntry("seabedHeight", seabedHeight_);
    os.writeEntry("waveMakerPos", xWaveMaker_);
    os.writeEntry("waveDictName", waveDictName_);
    os.writeEntry("waveDictPath", waveDictPath_);

    this->writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
