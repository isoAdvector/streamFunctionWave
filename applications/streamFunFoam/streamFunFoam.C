/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is an extenstion of OpenFOAM.

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

Application
    streamFunFoam calculates the theoretical alpha1 and U for a stream function
    wave. This application is used to validate stream function waves generated
    with CFD.

Author
    Johan Roenby, DHI

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "streamFunctionWave.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "cutCellIso.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createTimeControls.H"

    #include "createFields.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        t = runTime.time().value();
        etap = wave.getEta(xp, t, phase) - yp;
        {
            alpha1 = 0.0;
            scalarField& alphaIn = alpha1.primitiveFieldRef();
            cutCellIso cci(mesh, etap);
            forAll(alphaIn, cellI)
            {
                cci.calcSubCell(cellI, 0.0);
                alphaIn[cellI] = cci.VolumeOfFluid();
            }
        }

        eta = wave.getEta(x,t,phase);
        u = wave.getU(x,y,t,phase);
        v = wave.getV(x,y,t,phase);
        umax = wave.getU(x,eta,t,phase);
        vmax = wave.getV(x,eta,t,phase);
        rho = alpha1*rho1 + (scalar(1) - alpha1)*rho2;
        pByRho = wave.getPbyRho(x, y, t, phase);
        forAll(p, ci)
        {
            p[ci] = rho[ci]*pByRho[ci];
        }
        p_rgh = p - rho*(g & mesh.C());

        forAll(U,ci)
        {
            u[ci] = sign(u[ci])*min(mag(u[ci]),mag(umax[ci]));
            v[ci] = sign(v[ci])*min(mag(v[ci]),mag(vmax[ci]));
            U[ci] = (
                    alpha1[ci]*rho1.value()*((u[ci] - comoving*c)*K/mag(K)
                        + v[ci]*up) + (1-alpha1[ci])*rho2.value()*uAir
                    )/rho[ci];
        }
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
