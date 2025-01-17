/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Modified code Copyright (C) 2022 henning Scheufler
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "standard.H"

namespace Foam
{
    defineTypeNameAndDebug(standard, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::standard::standard
(
    const fvMesh& mesh,
    const dictionary& dict
)
{

}


void Foam::standard::update(algebraicVoF& newSurf,algebraicVoF& oldSurf,Foam::timeState state)
{
    Info << "update surface"  << endl;
}
// ************************************************************************* //
