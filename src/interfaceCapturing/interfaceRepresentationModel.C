/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 DLR
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

#include "interfaceRepresentationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Surface, class SurfaceStrategy>
Foam::interfaceRepresentationModel<Surface,SurfaceStrategy>::interfaceRepresentationModel
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    interfaceRepresentation(alpha1,phi,U,dict),
    interfaceRepresentationModelCoeffs_(dict),
    oldSurface_(alpha1,phi,U),
    newSurface_(alpha1,phi,U),
    surfaceStrat_(alpha1.mesh())
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- update the interface for the given time
template<class Surface, class SurfaceStrategy>
const Surface& Foam::interfaceRepresentationModel<Surface,SurfaceStrategy>::surface(timeState state)
{
    surfaceStrat_.update(newSurface_,oldSurface_,state);
    if (state == timeState::newState)
    {
        return newSurface_;
    }
    return oldSurface_;
};


// ************************************************************************* //
