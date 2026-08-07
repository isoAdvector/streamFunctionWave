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

Application
    streamFunFoam

Description
    Sets the alpha, U, phi, p and p_rgh fields of an interFoam type case to
    the analytical solution of a stream function theory wave, at every time
    step, and writes them at every write time.

    The wave is calculated with the streamFunctionWave class following:

    \verbatim
        Fenton, J. D. (1999).
        Numerical methods for nonlinear waves.
        Advances in Coastal and Ocean Engineering, 5, 241-324.
    \endverbatim

    No equations are solved - the fields are explicit functions of position
    and time - so the utility serves two purposes:

      -# generating consistent initial conditions for a wave simulation, and
      -# generating a reference solution to compare simulation results with.

    The fields are updated at every time step, with the time step controlled
    by maxCo, maxAlphaCo and maxDeltaT exactly as in the VoF solvers, so that
    function objects sample the wave at the same rate as in the simulation
    being compared with.

Usage
    Requires the following case files:

      - constant/waveProperties: wave input, see the waveProperties file of
        the calcStreamFunctionWave application. In addition this utility
        reads:
          - direction:     wave propagation direction (mandatory)
          - waveMakerPos:  origin of the horizontal coordinate
                           (default: lowest x of the mesh in the wave
                           direction)
          - seabedHeight:  origin of the vertical coordinate
                           (default: lowest point of the mesh)
          - phase:         wave phase in wave lengths
                           (default: eta = d at x = 0, t = 0)
          - uAir:          uniform velocity assigned to the air phase
                           (default: (0 0 0))
          - comoving:      subtract the wave celerity from the water velocity
                           (default: false)
          - c:             wave celerity used when comoving is true
                           (default: 0)

      - constant/transportProperties: phase names and densities
      - constant/g
      - 0/alpha.<phase1>, 0/U and 0/p_rgh

Author
    Johan Roenby, STROMNING

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressibleTwoPhaseMixture.H"
#include "cutCellIso.H"
#include "streamFunctionWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Set the two-phase fields of an interFoam type case to the"
        " analytical solution of a stream function theory wave"
        " and write them at every write time."
    );

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createTimeControls.H"
    #include "createFields.H"

    // Fields at the start time
    #include "setWaveFields.H"

    Info<< "Writing fields for time = " << runTime.timeName() << nl << endl;
    alpha1.write();
    U.write();
    phi.write();
    p.write();
    p_rgh.write();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        // Follow the time step the VoF solvers would take for this wave
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "setWaveFields.H"

        if (runTime.writeTime())
        {
            runTime.write();
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
