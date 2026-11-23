/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Johan Roenby
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

Description
    The stream function wave table compiled into the library, written by
    generateStreamFunctionWaveTable. A plain array of constant numbers, so
    it is placed in read-only memory: shared between the processes on a
    node and loaded only where used.

\*---------------------------------------------------------------------------*/

#include "streamFunctionWaveTable.H"

// * * * * * * * * * * * * * * * Local Data  * * * * * * * * * * * * * * * * //

namespace
{

const double tableData[] =
{
    #include "table/streamFunctionWaveTable.dat"
};

} // End anonymous namespace


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::streamFunctionWaveTable& Foam::streamFunctionWaveTable::builtIn()
{
    static const streamFunctionWaveTable table
    (
        tableData,
        sizeof(tableData)/sizeof(double)
    );
    return table;
}


// ************************************************************************* //
