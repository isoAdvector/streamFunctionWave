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

    #include "createPhi.H"


    Info<< "Reading transportProperties\n" << endl;
    twoPhaseMixture twoPhaseProperties(U, phi);

    volScalarField& alpha1(twoPhaseProperties.alpha1());

    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();


    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + (scalar(1) - alpha1)*rho2,
        alpha1.boundaryField().types()
    );
    rho.oldTime();


    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi
    );


    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);


    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );

    #include "readGravitationalAcceleration.H"

    /*
    dimensionedVector g0(g);

    // Read the data file and initialise the interpolation table
    interpolationTable<vector> timeSeriesAcceleration
    (
        runTime.path()/runTime.caseConstant()/"acceleration.dat"
    );
    */

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    fv::IOoptionList fvOptions(mesh);

    IOdictionary waveProps
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

    vector up = -g.value();
    up /= mag(up);
    vector K = vector(waveProps.lookupOrDefault<vector>("direction",vector::zero));

    scalar ymin(1e20), xmin(1e20);
    const pointField& pts = mesh.points();

    ymin = gMin(up & pts);
    xmin = gMin(K & pts);
    
	scalar xWaveMaker = waveProps.lookupOrDefault<scalar>("waveMakerPos",xmin);
	scalar seabedHeight = waveProps.lookupOrDefault<scalar>("seabedHeight",ymin);

    scalarField x = (K/mag(K)) & mesh.C();
	x -= xWaveMaker;
	scalarField y = up & mesh.C();
	y -= seabedHeight; //Vertical coordinate as measured from seabed

    streamFunctionWave wave(waveProps);
    wave.calcCoeffs();
    scalar phase = 0;
    if ( waveProps.found("phase") )
	{
		phase = waveProps.lookupOrDefault<scalar>("phase",0.0);
	}
	else
	{
		phase = 1.0 - wave.zeroPointPhase();
	}

    vector uAir = vector(waveProps.lookupOrDefault<vector>("uAir",vector::zero));

    //Setting stresam function alpha field in whole domain
	scalar t = runTime.time().value();
    scalarField xp = ((K/mag(K)) & pts);
    xp -= xWaveMaker;
	scalarField yp = (up & pts);
    yp -= seabedHeight;
	scalarField etap = wave.getEta(xp, t, phase) - yp;
    {
        isoCutCell icc(mesh, etap);
        icc.VolumeOfFluid(alpha1, 0);
    }
    alpha1.write(); //Writing volScalarField alpha1

	scalarField eta = wave.getEta(x,t,phase);
	scalarField u = wave.getU(x,y,t,phase);
	scalarField v = wave.getV(x,y,t,phase);
    scalarField umax = wave.getU(x,eta,t,phase);
    scalarField vmax = wave.getV(x,eta,t,phase);

    forAll(U,ci)
    {
        u[ci] = sign(u[ci])*min(mag(u[ci]),mag(umax[ci]));
        v[ci] = sign(v[ci])*min(mag(v[ci]),mag(vmax[ci]));
        U[ci] = pos(alpha1[ci] - 1e-3)*(u[ci]*K/mag(K) + v[ci]*up)
                + neg(alpha1[ci] - 1e-3)*uAir;
    }

    U.write(); //Writing volScalarField alpha1

    /*    
    //Find the stream function wave patch
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    label sfpLabel = -1;
    forAll(patches, patchI)
    {
        if(waveProps.found(patches[patchI].name()))
        {
            sfpLabel = patchI;
            break;
        }
    }
    
    if (isA<streamFunctionWave>(U.boundaryField()[sfpLabel]))
    {
        //Set the alpha field in the whole domain
        scalarField eta = U.boundaryField()[sfpLabel].getEta(mesh.C(),runTime.time().value(),0)
        //Set the velocity field in the whole domain
    }
 */
