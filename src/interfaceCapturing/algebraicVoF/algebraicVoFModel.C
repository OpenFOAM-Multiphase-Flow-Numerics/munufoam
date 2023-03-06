/*---------------------------------------------------------------------------*\
    Modified work | Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
    unofficial extension to OpenFOAM.

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

#include "algebraicVoFModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class Surface, class SurfaceStrategy, class AdvectionStrategy>
Foam::algebraicVoFModel<Surface,SurfaceStrategy,AdvectionStrategy>::algebraicVoFModel
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    algebraicVoFMethod(alpha1,phi,U),
    oldSurface_(alpha1.oldTime(),phi,U),
    newSurface_(alpha1,phi,U),
    surfaceStrat_(alpha1.mesh(),alpha1.mesh().solverDict(alpha1.name())),
    advectionStrat_(alpha1,alphaPhi_,phi,U)
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


//- advection of the interface
template<class Surface, class SurfaceStrategy, class AdvectionStrategy>
void Foam::algebraicVoFModel<Surface,SurfaceStrategy,AdvectionStrategy>::advect
(
    const volScalarField::Internal& Sp,
    const volScalarField::Internal& Su
)
{
    advectionStrat_.advect(newSurface_,oldSurface_,Sp,Su);
};

template<class Surface, class SurfaceStrategy, class AdvectionStrategy>
const Surface&  Foam::algebraicVoFModel<Surface,SurfaceStrategy,AdvectionStrategy>::surface
(
    timeState state
)
{
    surfaceStrat_.update(newSurface_,oldSurface_,state);
    if (state == timeState::newState)
    {
        return newSurface_;
    }
    return oldSurface_;
};

// ************************************************************************* //
