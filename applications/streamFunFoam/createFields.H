    #include "readGravitationalAcceleration.H"
    vector up = -g.value();
    up /= mag(up);

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

    //Setting stresam function alpha field in whole domain
	scalar t = runTime.time().value();
    scalarField xp = ((K/mag(K)) & pts);
    xp -= xWaveMaker;
	scalarField yp = (up & pts);
    yp -= seabedHeight;
	scalarField etap = wave.getEta(xp, t, phase) - yp;
    
    Info<< "Reading field U\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    isoCutCell icc(mesh, etap);
    icc.VolumeOfFluid(alpha1, 0);
    alpha1.write(); //Writing volScalarField alpha1
    
    // Calculate rho
    Info<< "Reading fluid densities\n" << endl;
    IOdictionary tp
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    wordList defaultPhaseNames(2);
    defaultPhaseNames[0] = "phase1";
    defaultPhaseNames[1] = "phase2";
    //From OF-2.3.0 the phases keyword is required in transportProperties
    wordList phases(wordList(tp.lookupOrDefault("phases", defaultPhaseNames)));

    dictionary sD1(tp.subDict(phases[0]));
    dimensionedScalar rho1 = sD1.lookup("rho");	
    dictionary sD2(tp.subDict(phases[1]));
    dimensionedScalar rho2 = sD2.lookup("rho");

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
    
    //Calculating velocity field
	scalarField eta = wave.getEta(x,t,phase);
	scalarField u = wave.getU(x,y,t,phase);
	scalarField v = wave.getV(x,y,t,phase);
    scalarField umax = wave.getU(x,eta,t,phase);
    scalarField vmax = wave.getV(x,eta,t,phase);

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
    vector uAir = vector(waveProps.lookupOrDefault<vector>("uAir",vector::zero));
    
    forAll(U,ci)
    {
        u[ci] = sign(u[ci])*min(mag(u[ci]),mag(umax[ci]));
        v[ci] = sign(v[ci])*min(mag(v[ci]),mag(vmax[ci]));
        U[ci] = (
                    alpha1[ci]*rho1.value()*(u[ci]*K/mag(K) + v[ci]*up)
                    + (1-alpha1[ci])*rho2.value()*uAir
                )/rho[ci];
    }
    U.write();

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
        linearInterpolate(U) & mesh.Sf()
    );
    
    //Calculating p and p_rgh

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
        rho*(g & mesh.C())
    );

    scalarField pByRho = wave.getPbyRho(x, y, t, phase);
    forAll(p, ci)
    {
        p[ci] = rho[ci]*pByRho[ci];
    }    

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
    p_rgh = p - rho*(g & mesh.C());
    
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
        p_rgh = p - rho*(g & mesh.C());
    }
    p.write();
    p_rgh.write();

    fv::IOoptionList fvOptions(mesh);