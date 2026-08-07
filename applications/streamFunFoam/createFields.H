#include "readGravitationalAcceleration.H"

Info<< "Reading waveProperties\n" << endl;

IOdictionary waveProperties
(
    IOobject
    (
        "waveProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Calculating stream function wave coefficients" << endl;

streamFunctionWave wave(waveProperties);
wave.calcCoeffs();

// Unit vector pointing up and unit vector in the wave propagation direction
const vector up(normalised(-g.value()));
const vector kDir(normalised(waveProperties.get<vector>("direction")));

// Origin of the wave coordinate system. By default the wave starts at the
// backmost mesh point and the water depth is measured from the lowest one.
const scalar xWaveMaker
(
    waveProperties.getOrDefault<scalar>
    (
        "waveMakerPos",
        gMin(kDir & mesh.points())
    )
);

const scalar seabedHeight
(
    waveProperties.getOrDefault<scalar>
    (
        "seabedHeight",
        gMin(up & mesh.points())
    )
);

// Wave phase in wave lengths. Default puts the still water level at x = 0
// for t = 0.
const scalar phase
(
    waveProperties.getOrDefault<scalar>("phase", 1.0 - wave.zeroPointPhase())
);

Info<< "Wave direction: " << kDir << nl
    << "Wave maker position: " << xWaveMaker << nl
    << "Seabed height: " << seabedHeight << nl
    << "Wave phase: " << phase << " wave lengths" << nl << endl;

// Wave coordinates of the mesh points and of the cell centres: xp/xc measured
// along the wave direction from the wave maker, yp/yc measured up from the
// seabed.
const scalarField xp((kDir & mesh.points()) - xWaveMaker);
const scalarField yp((up & mesh.points()) - seabedHeight);

const scalarField xc((kDir & mesh.C().primitiveField()) - xWaveMaker);
const scalarField yc((up & mesh.C().primitiveField()) - seabedHeight);

// Wave coordinates of the boundary face centres
List<scalarField> xb(mesh.boundary().size());
List<scalarField> yb(mesh.boundary().size());
forAll(mesh.boundary(), patchi)
{
    const vectorField Cf(mesh.boundary()[patchi].Cf());
    xb[patchi] = (kDir & Cf) - xWaveMaker;
    yb[patchi] = (up & Cf) - seabedHeight;
}

// Velocity assigned to the air phase
const vector uAir(waveProperties.getOrDefault<vector>("uAir", vector::zero));

// Velocity of the frame of reference. With comoving enabled the wave is
// stationary in the generated fields.
const scalar uFrame
(
    waveProperties.getOrDefault("comoving", false)
  ? waveProperties.get<scalar>("c")
  : 0
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);

Info<< "Reading transportProperties\n" << endl;
incompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);

// Gravitational potential of the cell centres
const volScalarField gh("gh", g & mesh.C());

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);
