/*---------------------------------------------------------------------------*\
Global
    alphaCourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers of the cells
    containing the interface.

    Same as the alphaCourantNo.H of the VoF solvers, except that the
    interface cells are found from alpha1 directly, avoiding the dependency
    on interfaceProperties.

\*---------------------------------------------------------------------------*/

scalar maxAlphaCo
(
    runTime.controlDict().getOrDefault<scalar>("maxAlphaCo", maxCo)
);

scalar alphaCoNum = 0.0;
scalar meanAlphaCoNum = 0.0;

if (mesh.nInternalFaces())
{
    const scalarField sumPhi
    (
        pos0(alpha1.primitiveField() - 0.01)
       *pos0(0.99 - alpha1.primitiveField())
       *fvc::surfaceSum(mag(phi))().primitiveField()
    );

    alphaCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    meanAlphaCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Interface Courant Number mean: " << meanAlphaCoNum
    << " max: " << alphaCoNum << endl;

// ************************************************************************* //
