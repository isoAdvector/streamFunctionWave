/*---------------------------------------------------------------------------*\
Global
    setWaveFields

Description
    Sets alpha1, U, phi, p and p_rgh to the stream function wave solution at
    the current time.

\*---------------------------------------------------------------------------*/

{
    const scalar t = runTime.value();

    // Signed distance from the free surface to the mesh points, positive in
    // the water
    scalarField fp(wave.getEta(xp, t, phase) - yp);

    // Volume fraction from the geometric intersection of the cells with the
    // free surface
    cutCellIso cutCell(mesh, fp);
    scalarField& alphaIn = alpha1.primitiveFieldRef();
    forAll(alphaIn, celli)
    {
        cutCell.calcSubCell(celli, 0.0);
        alphaIn[celli] = clamp(cutCell.VolumeOfFluid(), zero_one{});
    }
    alpha1.correctBoundaryConditions();

    alpha2 = scalar(1) - alpha1;
    rho = alpha1*rho1 + alpha2*rho2;

    // Wave kinematics in the cell centres. Above the free surface the Fourier
    // series diverges, so the velocity there is limited by its value at the
    // surface.
    const scalarField etac(wave.getEta(xc, t, phase));
    const scalarField uc(wave.getU(xc, yc, t, phase));
    const scalarField vc(wave.getV(xc, yc, t, phase));
    const scalarField uSurf(wave.getU(xc, etac, t, phase));
    const scalarField vSurf(wave.getV(xc, etac, t, phase));

    const scalarField& rhoIn = rho.primitiveField();
    vectorField& UIn = U.primitiveFieldRef();
    forAll(UIn, celli)
    {
        const scalar u = sign(uc[celli])*min(mag(uc[celli]), mag(uSurf[celli]));
        const scalar v = sign(vc[celli])*min(mag(vc[celli]), mag(vSurf[celli]));

        // Momentum weighted mean of the water and the air velocity
        UIn[celli] =
        (
            alphaIn[celli]*rho1.value()*((u - uFrame)*kDir + v*up)
          + (scalar(1) - alphaIn[celli])*rho2.value()*uAir
        )/rhoIn[celli];
    }
    U.correctBoundaryConditions();

    phi = fvc::flux(U);

    // Wave pressure, Fenton (1999), eq. (3.28), scaled with the local density
    p.primitiveFieldRef() = rhoIn*wave.getPbyRho(xc, yc, t, phase);

    volScalarField::Boundary& pBf = p.boundaryFieldRef();
    forAll(pBf, patchi)
    {
        pBf[patchi] =
            rho.boundaryField()[patchi]
           *wave.getPbyRho(xb[patchi], yb[patchi], t, phase);
    }

    p_rgh = p - rho*gh;
}


// ************************************************************************* //
