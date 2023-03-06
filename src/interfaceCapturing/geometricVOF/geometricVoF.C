/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
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

#include "geometricVoF.H"

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"

namespace Foam
{
    defineTypeNameAndDebug(geometricVoF, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::geometricVoF::geometricVoF
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    timeState state
)
:
    surfaceBase(alpha1,phi,U),
    state(state),
    alpha_(alpha1),
    normal
    (
        IOobject
        (
            IOobject::groupName("interfaceNormal", state == timeState::oldState ? word(alpha1.group() + "_old") : word(alpha1.group() + "_new")),
            alpha1.mesh().time().timeName(),
            alpha1.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1.mesh(),
        dimensionedVector(dimArea, Zero)
    ),
    centre
    (
        IOobject
        (
            IOobject::groupName("interfaceCentre", state == timeState::oldState ? word(alpha1.group() + "_old") : word(alpha1.group() + "_new")),
            alpha1.mesh().time().timeName(),
            alpha1.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1.mesh(),
        dimensionedVector(dimLength, Zero)
    )
{
    
}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //



// ************************************************************************* //
