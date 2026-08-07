/*---------------------------------------------------------------------------*\
Global
    setDeltaT

Description
    Reset the timestep to maintain a constant maximum Courant number.
    Reduction of the time step is immediate, but increase is damped to avoid
    unstable oscillations.

    Same as the setDeltaT.H of the VoF solvers, kept locally since the
    setDeltaT.H of finiteVolume does not account for maxAlphaCo.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    const scalar maxDeltaTFact =
        min(maxCo/(CoNum + SMALL), maxAlphaCo/(alphaCoNum + SMALL));

    const scalar deltaTFact =
        min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}

// ************************************************************************* //
